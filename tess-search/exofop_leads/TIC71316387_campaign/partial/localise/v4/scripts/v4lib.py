"""v4 (independent re-derivation) TESS pixel-level library for TIC 71316387.

Written from scratch for v4; only the downloaded TESScut cubes (../data/tesscut) and PRF files (../prf_files) are reused.

Method differences w.r.t. v3 (deliberate, for independence):
  * per-window, per-pixel MATCHED-FILTER regression  f(t) = b0 + b1*(t-tc) - A * S(t-tc)  over |t-tc| <= W*T14,
    S = exposure-integrated eclipse template (unit peak) -> 'amplitude image' A (e/s at eclipse peak) and its error;
    all in-eclipse cadences (incl. ingress/egress) are used with template weighting;
  * PRF positions from a vectorised grid search over precomputed PRF images (0.02 px) + 2-D quadratic refinement;
  * own pointing-excursion cut (moment centroid, 0.25-d running median, 5 MAD) and own Gaia query (v4/gaia).
Conventions: 0-based pixel coordinates, integer = pixel centre (PRF.locate convention); images [row, col] = [y, x];
sky offsets (dE, dN) in arcsec with dE = dRA cos(dec); BTJD = BJD_TDB - 2457000.
"""
import os, glob, json
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table

BASE = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise'
V4 = BASE + '/v4'
SECTORS = [1, 28, 68, 95, 105, 106]
P = 1.3036712
T0 = 3165.6697
T14 = float(os.environ.get('V4_T14H', 1.54)) / 24.0
TARGET = 6571343840265900544
RA0, DE0 = 330.887948, -42.835511
WHALF = float(os.environ.get('V4_W', 2.6))        # window half-width in T14
INMIN = 0.5                                       # in-eclipse coverage requirement (fraction of |dt|<T14/2 cadences)
SIDEMIN = 0.4                                     # each OOT side coverage requirement
BOX = int(os.environ.get('V4_BOX', 3))            # fit box half size -> (2*BOX+1)^2
BASEDEG = int(os.environ.get('V4_BASEDEG', 1))    # baseline polynomial degree
SHAPE = os.environ.get('V4_SHAPE', 'V')           # 'V' triangle | 'U' trapezoid (T12 = 0.3*T14)


def to_json(obj, path):
    def conv(o):
        if isinstance(o, np.floating): return float(o)
        if isinstance(o, np.integer): return int(o)
        if isinstance(o, np.ndarray): return o.tolist()
        if isinstance(o, np.bool_): return bool(o)
        raise TypeError(type(o))
    with open(path, 'w') as f:
        json.dump(obj, f, indent=1, default=conv)


def tmag(g, c):
    """Stassun et al. (2019) TIC-8 T(G, BP-RP)."""
    g = np.asarray(g, float); c = np.asarray(c, float)
    t = g - 0.00522555 * c ** 3 + 0.0891337 * c ** 2 - 0.633923 * c + 0.0324473
    return np.where(np.isfinite(c), t, g - 0.430)


def jyear(btjd):
    return 2000.0 + (np.asarray(btjd) + 2457000.0 - 2451545.0) / 365.25


_GAIA = None


def gaia():
    global _GAIA
    if _GAIA is None:
        t = Table.read(V4 + '/gaia/field_8arcmin.ecsv')
        d = {c: np.array(t[c].filled(np.nan) if hasattr(t[c], 'filled') else t[c], dtype=float) for c in
             ('ra', 'dec', 'pmra', 'pmdec', 'phot_g_mean_mag', 'bp_rp', 'parallax', 'ruwe')}
        d['source_id'] = np.array(t['source_id'], dtype=np.int64)
        d['T'] = tmag(d['phot_g_mean_mag'], d['bp_rp'])
        _GAIA = d
    return _GAIA


def template(dt, texp, T=None, nsub=21):
    """Unit-peak eclipse template integrated over the exposure; V = triangle of total width T14, U = trapezoid."""
    T = T14 if T is None else T
    s = (np.arange(nsub) + 0.5) / nsub - 0.5
    a = np.abs(np.asarray(dt)[:, None] + s[None, :] * texp)
    if SHAPE == 'U':
        t12 = 0.3 * T
        v = np.clip((0.5 * T - a) / t12, 0, 1)
    else:
        v = np.clip(1 - a / (0.5 * T), 0, None)
    return v.mean(1)


class Sector:
    def __init__(self, sec, cut=True):
        self.sec = sec
        fn = glob.glob(f'{BASE}/data/tesscut/s{sec}/*astrocut.fits')[0]
        with fits.open(fn) as h:
            self.cam = int(h[0].header['CAMERA']); self.ccd = int(h[0].header['CCD'])
            self.col0 = int(h[1].header['1CRV4P']); self.row0 = int(h[1].header['2CRV4P'])
            d = h[1].data
            self.t = np.array(d['TIME'], float)
            self.f = np.array(d['FLUX'], np.float64)
            self.q = np.array(d['QUALITY'], int)
            self.wcs = WCS(h[2].header)
            self.ffi = np.array(d['FFI_FILE'])
        self.ny, self.nx = self.f.shape[1:]
        fin = np.isfinite(self.t) & np.isfinite(self.f).reshape(len(self.t), -1).all(1)
        self.f[~fin] = 0.0
        self.good = fin & (self.q == 0)
        self.texp = float(np.nanmedian(np.diff(self.t[self.good])))
        self.epoch = float(jyear(np.median(self.t[self.good])))
        g = gaia()
        dty = self.epoch - 2016.0
        dec = g['dec'] + np.nan_to_num(g['pmdec']) * dty / 3.6e6
        ra = g['ra'] + np.nan_to_num(g['pmra']) * dty / 3.6e6 / np.cos(np.radians(g['dec']))
        x, y = self.wcs.all_world2pix(ra, dec, 0)
        self.gx, self.gy = x, y
        self.it = int(np.where(g['source_id'] == TARGET)[0][0])
        self.xt, self.yt = float(x[self.it]), float(y[self.it])
        self.ix, self.iy = int(round(self.xt)), int(round(self.yt))
        self.J = self.jacobian(self.xt, self.yt)
        self.ncut = 0
        if cut:
            self.pointing_cut()

    def jacobian(self, x, y, h=0.05):
        """[dE, dN] (arcsec) = J @ [dx, dy] (px), central differences of the cutout WCS."""
        w = self.wcs
        def sky(px, py):
            r = np.array(w.all_pix2world([[px, py]], 0))[0]
            return r
        r0 = sky(x, y)
        c = np.cos(np.radians(r0[1]))
        cols = []
        for dx, dy in ((h, 0), (0, h)):
            a = sky(x + dx, y + dy); b = sky(x - dx, y - dy)
            da = (a[0] - b[0] + 180) % 360 - 180
            cols.append(np.array([da * c * 3600, (a[1] - b[1]) * 3600]) / (2 * h))
        return np.column_stack(cols)

    def pointing_cut(self, nmad=5.0, win=0.25):
        sl = (slice(self.iy - 2, self.iy + 3), slice(self.ix - 2, self.ix + 3))
        st = self.f[:, sl[0], sl[1]]
        st = st - np.percentile(self.f.reshape(len(st), -1), 20, axis=1)[:, None, None]   # per-cadence sky level
        yy, xx = np.mgrid[sl]
        s = st.sum((1, 2))
        with np.errstate(all='ignore'):
            cx = (st * xx).sum((1, 2)) / s; cy = (st * yy).sum((1, 2)) / s
        g = np.where(self.good)[0]
        tg = self.t[g]
        lo = np.searchsorted(tg, tg - win); hi = np.searchsorted(tg, tg + win)
        bad = np.zeros(len(self.t), bool)
        for c in (cx, cy):
            cg = c[g]
            r = np.array([cg[k] - np.median(cg[lo[k]:hi[k]]) for k in range(len(g))])
            mad = 1.4826 * np.median(np.abs(r))
            bad[g[np.abs(r) > nmad * mad]] = True
        self.ncut = int(bad.sum())
        self.good &= ~bad

    # ---------------- eclipse windows -----------------
    def event_times(self, phase):
        tg = self.t[self.good]
        n = np.arange(np.floor((tg.min() - T0) / P) - 1, np.ceil((tg.max() - T0) / P) + 2)
        return T0 + (n + phase) * P, n.astype(int)

    def window_fit(self, tc):
        """Matched-filter regression of every pixel in the window around tc.
        Returns dict(A, sA, B, ok, n_in, ...) or None if coverage is insufficient."""
        dt = self.t - tc
        m = self.good & (np.abs(dt) <= WHALF * T14)
        ii = np.where(m)[0]
        dd = dt[ii]
        nin_exp = T14 / self.texp
        n_in = int((np.abs(dd) < 0.5 * T14).sum())
        side_exp = (WHALF - 0.6) * T14 / self.texp
        n_l = int((dd < -0.6 * T14).sum()); n_r = int((dd > 0.6 * T14).sum())
        if n_in < max(1, INMIN * nin_exp) or n_l < max(3, SIDEMIN * side_exp) or n_r < max(3, SIDEMIN * side_exp):
            return None
        S = template(dd, self.texp)
        X = np.column_stack([dd ** k for k in range(BASEDEG + 1)] + [-S])
        Y = self.f[ii].reshape(len(ii), -1)
        XtX = X.T @ X
        Ci = np.linalg.inv(XtX)
        coef = Ci @ (X.T @ Y)
        res = Y - X @ coef
        dof = len(ii) - X.shape[1]
        rms = np.sqrt((res ** 2).sum(0) / dof)
        A = coef[-1]; sA = rms * np.sqrt(Ci[-1, -1])
        shp = (self.ny, self.nx)
        return dict(tc=float(tc), A=A.reshape(shp), sA=sA.reshape(shp), B=coef[0].reshape(shp), n=len(ii), n_in=n_in,
                    idx=ii, X=X, Ci=Ci)

    def all_windows(self, fake_phases=None):
        """Real deep (0) and shallow (0.5) windows, plus fake windows at off-eclipse phases (never overlapping eclipses)."""
        out = dict(A=[], sA=[], B=[], kind=[], tc=[], ep=[], n_in=[])
        reals = []
        for kind, ph in (('deep', 0.0), ('shallow', 0.5)):
            tcs, ns = self.event_times(ph)
            reals += list(tcs)
            for tc, n in zip(tcs, ns):
                w = self.window_fit(tc)
                if w is None:
                    continue
                for k, v in (('A', w['A']), ('sA', w['sA']), ('B', w['B']), ('kind', kind), ('tc', tc), ('ep', n), ('n_in', w['n_in'])):
                    out[k].append(v)
        reals = np.array(reals)
        if fake_phases is None:
            fake_phases = np.r_[np.arange(0.15, 0.3501, 0.025), np.arange(0.65, 0.8501, 0.025)]
        tcs, ns = self.event_times(0.0)
        for tc0, n in zip(tcs, ns):
            for ph in fake_phases:
                tc = tc0 + ph * P
                if np.min(np.abs(tc - reals)) < (WHALF + 0.55) * T14:
                    continue
                w = self.window_fit(tc)
                if w is None:
                    continue
                for k, v in (('A', w['A']), ('sA', w['sA']), ('B', w['B']), ('kind', 'fake'), ('tc', tc), ('ep', n), ('n_in', w['n_in'])):
                    out[k].append(v)
        return {k: np.array(v) for k, v in out.items()}


def stack(W, idx, rel=1.0):
    """Weighted stack of amplitude images (scalar weight per window, preserves PSF shape).
    rel: relative depth used to normalise (1 for deep; ratio for shallow when merging)."""
    idx = np.asarray(idx)
    rel = np.broadcast_to(np.asarray(rel, float), idx.shape)
    v = np.array([np.median(s[s > 0] ** 2) for s in W['sA'][idx]]) / rel ** 2
    w = 1.0 / v
    An = W['A'][idx] / rel[:, None, None]
    sn = W['sA'][idx] / rel[:, None, None]
    As = (w[:, None, None] * An).sum(0) / w.sum()
    ss = np.sqrt((w[:, None, None] ** 2 * sn ** 2).sum(0)) / w.sum()
    wb = 1.0 / np.array([np.median(s[s > 0] ** 2) for s in W['sA'][idx]])
    Bs = (wb[:, None, None] * W['B'][idx]).sum(0) / wb.sum()
    return As, ss, Bs


class PRFGrid:
    """Precomputed TESS PRF (PRF package) images on a sub-pixel grid around the target, restricted to the fit box."""
    def __init__(self, S, half=1.3, step=0.02, box=BOX):
        import PRF
        sub = 's0001' if S.sec < 4 else 's0004'
        self.S = S
        self.prf = PRF.TESS_PRF(S.cam, S.ccd, S.sec, S.col0 + S.xt, S.row0 + S.yt, localdatadir=f'{BASE}/prf_files/{sub}/')
        self.shape = (S.ny, S.nx)
        self.box = box
        self.sl = (slice(S.iy - box, S.iy + box + 1), slice(S.ix - box, S.ix + box + 1))
        self.gxv = S.xt + np.arange(-half, half + 1e-9, step)
        self.gyv = S.yt + np.arange(-half, half + 1e-9, step)
        self.step = step
        GX, GY = np.meshgrid(self.gxv, self.gyv)
        self.GX, self.GY = GX.ravel(), GY.ravel()
        self.M = np.array([self.img(x, y)[self.sl].ravel() for x, y in zip(self.GX, self.GY)])     # (K, npix)
        # tied Gaia neighbours (T < 19) that put light into the box: relative flux and offset
        g = gaia()
        d = np.hypot(S.gx - S.xt, S.gy - S.yt)
        sel = np.where((g['T'] < 19.0) & (d < box + 6) & (np.arange(len(d)) != S.it))[0]
        self.nb = [(float(S.gx[i] - S.xt), float(S.gy[i] - S.yt), float(10 ** (-0.4 * (g['T'][i] - g['T'][S.it]))), int(g['source_id'][i])) for i in sel]
        self.N = np.array([self.nbimg(x, y)[self.sl].ravel() for x, y in zip(self.GX, self.GY)])

    def img(self, x, y):
        return self.prf.locate(x, y, self.shape)

    def nbimg(self, x, y):
        m = np.zeros(self.shape)
        for dx, dy, r, _ in self.nb:
            m += r * self.img(x + dx, y + dy)
        return m

    def fit(self, im, sig, direct=False, plane=False):
        """Weighted LSQ  im = a*PRF(x,y) [+ tied neighbours] + c [+ plane] on the box, grid search + quadratic refinement.
        Returns x, y, a, c, chi2, and the chi2 grid minimum index."""
        b = im[self.sl].ravel(); w = 1.0 / sig[self.sl].ravel()
        M = self.M + (self.N if direct else 0.0)
        yy, xx = np.mgrid[self.sl]
        extra = [np.ones_like(b)]
        if plane:
            extra += [(xx.ravel() - self.S.xt), (yy.ravel() - self.S.yt)]
        E = np.column_stack(extra) * w[:, None]           # (npix, q)
        bw = b * w
        # project out the nuisance columns: chi2(k) = |P bw|^2 - (m_k^T P bw)^2 / (m_k^T P m_k)
        Q, _ = np.linalg.qr(E)
        Pb = bw - Q @ (Q.T @ bw)
        Mw = M * w[None, :]
        PM = Mw - (Mw @ Q) @ Q.T
        num = PM @ Pb
        den = (PM * PM).sum(1)
        chi2 = Pb @ Pb - num ** 2 / den
        k = int(np.argmin(chi2))
        ny, nx = len(self.gyv), len(self.gxv)
        j, i = divmod(k, nx)
        x, y = self.GX[k], self.GY[k]
        if 0 < i < nx - 1 and 0 < j < ny - 1:
            C = chi2.reshape(ny, nx)[j - 1:j + 2, i - 1:i + 2]
            # 2-D quadratic through the 3x3 neighbourhood
            uu, vv = np.meshgrid([-1, 0, 1], [-1, 0, 1])
            A = np.column_stack([np.ones(9), uu.ravel(), vv.ravel(), uu.ravel() ** 2, uu.ravel() * vv.ravel(), vv.ravel() ** 2])
            cq, *_ = np.linalg.lstsq(A, C.ravel(), rcond=None)
            H = np.array([[2 * cq[3], cq[4]], [cq[4], 2 * cq[5]]])
            if np.all(np.linalg.eigvalsh(H) > 0):
                du = -np.linalg.solve(H, cq[1:3])
                if np.all(np.abs(du) <= 1.0):
                    x += du[0] * self.step; y += du[1] * self.step
        a = num[k] / den[k]
        edge = (i == 0) or (i == nx - 1) or (j == 0) or (j == ny - 1)
        return dict(x=float(x), y=float(y), a=float(a), chi2=float(chi2[k]), npix=len(b), edge=bool(edge), k=k)

    def fit_cont(self, im, sig, x0, y0, direct=False):
        """Continuous refinement (Nelder-Mead on the exact PRF interpolation) for checks."""
        from scipy.optimize import minimize
        b = im[self.sl].ravel(); w = 1.0 / sig[self.sl].ravel()
        def chi(p):
            m = self.img(p[0], p[1])
            if direct:
                m = m + self.nbimg(p[0], p[1])
            A = np.column_stack([m[self.sl].ravel(), np.ones_like(b)]) * w[:, None]
            c, *_ = np.linalg.lstsq(A, b * w, rcond=None)
            r = A @ c - b * w
            return float(r @ r)
        r = minimize(chi, [x0, y0], method='Nelder-Mead', options=dict(xatol=1e-4, fatol=1e-6, maxiter=400))
        return dict(x=float(r.x[0]), y=float(r.x[1]), chi2=float(r.fun))

    def empirical(self, B, sig):
        """Unit-flux empirical target image from a direct image: (B - c - a*neighbours)/a."""
        r = self.fit(B, sig, direct=True)
        b = B[self.sl].ravel(); w = 1.0 / sig[self.sl].ravel()
        m = self.img(r['x'], r['y']); n = self.nbimg(r['x'], r['y'])
        A = np.column_stack([(m + n)[self.sl].ravel(), np.ones_like(b)]) * w[:, None]
        c, *_ = np.linalg.lstsq(A, b * w, rcond=None)
        E = (B - c[1] - c[0] * n) / c[0]
        return E, r, c
