"""v3 (independent re-derivation) library: TESScut cubes, eclipse windows, difference images, TESS-PRF fits.

Conventions
  pixel coords 0-based, integer = pixel centre (astropy origin 0; PRF.locate uses the same convention)
  images [row, col] = [y, x]
  difference image D = baseline(predicted at in-eclipse cadences) - in-eclipse  -> fading source = positive PSF
  sky offsets (dE, dN) arcsec, dE = dRA*cos(dec)
  BTJD = BJD_TDB - 2457000
"""
import os, glob, json
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares

BASE = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise'
V3 = BASE + '/v3'
SECTORS = [1, 28, 68, 95, 105, 106]
P = 1.3036712
T0 = 3165.6697
T14 = 1.55 / 24.0
TARGET = 6571343840265900544
RA0, DE0 = 330.887948, -42.835511

# window geometry (fractions of T14), overridable for variants
IN_HALF = float(os.environ.get('V3_IN_HALF', 0.30)) * T14
OOT_A = float(os.environ.get('V3_OOT_A', 0.9)) * T14
OOT_B = float(os.environ.get('V3_OOT_B', 2.6)) * T14
BOX = int(os.environ.get('V3_BOX', 3))          # fit box half-size (px): (2*BOX+1)^2 pixels


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
    """Stassun et al. (2019) TIC-8 relation T(G, BP-RP)."""
    g = np.asarray(g, float); c = np.asarray(c, float)
    t = g - 0.00522555 * c**3 + 0.0891337 * c**2 - 0.633923 * c + 0.0324473
    return np.where(np.isfinite(c), t, g - 0.430)


def jyear(btjd):
    return 2000.0 + (np.asarray(btjd) + 2457000.0 - 2451545.0) / 365.25


def gaia_field():
    d = pd.read_csv(V3 + '/gaia/field_6arcmin.csv')
    d['T'] = tmag(d.phot_g_mean_mag, d.bp_rp)
    return d


class Cube:
    def __init__(self, sector, centroid_cut=True):
        self.sector = sector
        fn = glob.glob(f'{BASE}/data/tesscut/s{sector}/*astrocut.fits')[0]
        with fits.open(fn) as h:
            self.h0 = h[0].header; self.h1 = h[1].header; self.h2 = h[2].header
            d = h[1].data
            self.t = np.array(d['TIME'], float)
            self.f = np.array(d['FLUX'], float)
            self.fe = np.array(d['FLUX_ERR'], float)
            self.q = np.array(d['QUALITY'], int)
        self.cam = int(self.h0['CAMERA']); self.ccd = int(self.h0['CCD'])
        self.col0 = int(self.h1['1CRV4P']); self.row0 = int(self.h1['2CRV4P'])
        self.ny, self.nx = self.f.shape[1:]
        fin = np.isfinite(self.t) & np.all(np.isfinite(self.f.reshape(len(self.t), -1)), axis=1)
        self.good = (self.q == 0) & fin
        self.f[~fin] = 0.0
        self.dt = float(np.nanmedian(np.diff(self.t)))
        self.wcs = WCS(self.h2)
        self.epoch = float(jyear(np.nanmedian(self.t[self.good])))
        g = gaia_field()
        dtyr = self.epoch - g.ref_epoch.values
        dec = g.dec.values + np.nan_to_num(g.pmdec.values) * dtyr / 3.6e6
        ra = g.ra.values + np.nan_to_num(g.pmra.values) * dtyr / 3.6e6 / np.cos(np.radians(g.dec.values))
        x, y = self.wcs.all_world2pix(ra, dec, 0)
        g['x'] = x; g['y'] = y; g['ra_e'] = ra; g['dec_e'] = dec
        self.it = int(np.where(g.source_id.values == TARGET)[0][0])
        self.g = g
        self.xt, self.yt = float(x[self.it]), float(y[self.it])
        self.jac = self._jac(self.xt, self.yt)
        self.ix, self.iy = int(round(self.xt)), int(round(self.yt))
        self.ncut = 0
        if centroid_cut:
            self._centroid_cut()

    def _jac(self, x, y, h=0.01):
        w = self.wcs
        r0 = np.array(w.all_pix2world([[x, y]], 0))[0]
        rx = np.array(w.all_pix2world([[x + h, y]], 0))[0]
        ry = np.array(w.all_pix2world([[x, y + h]], 0))[0]
        c = np.cos(np.radians(r0[1]))
        f = lambda r: np.array([(r[0] - r0[0]) * c * 3600, (r[1] - r0[1]) * 3600])
        return np.column_stack([f(rx) / h, f(ry) / h])    # [dE, dN] = J @ [dx, dy]  (arcsec / px)

    def _centroid_cut(self, nmad=6.0, win=0.5):
        """Remove pointing excursions: 3x3-px flux-weighted centroid deviating > nmad*MAD from a 0.5-d running median.
        (A 0.6% eclipse of a source <= 5" away moves this centroid by < 2 mpx, far below the threshold.)"""
        sl = (slice(self.iy - 1, self.iy + 2), slice(self.ix - 1, self.ix + 2))
        yy, xx = np.mgrid[sl]
        st = self.f[:, sl[0], sl[1]]
        s = st.sum((1, 2))
        with np.errstate(all='ignore'):
            cx = (st * xx).sum((1, 2)) / s; cy = (st * yy).sum((1, 2)) / s
        g = np.where(self.good)[0]
        bad = np.zeros(len(self.t), bool)
        for c in (cx, cy):
            r = np.full(len(self.t), np.nan)
            tg = self.t[g]; cg = c[g]
            # running median via sorted window indices
            lo = np.searchsorted(tg, tg - win); hi = np.searchsorted(tg, tg + win)
            for k, i in enumerate(g):
                r[i] = c[i] - np.median(cg[lo[k]:hi[k]])
            mad = 1.4826 * np.nanmedian(np.abs(r[g]))
            bad |= np.abs(r) > nmad * mad
        self.ncut = int((bad & self.good).sum())
        self.good &= ~bad
        self.cx, self.cy = cx, cy

    def epochs(self, phase):
        tg = self.t[self.good]
        n = np.arange(np.floor((tg.min() - T0) / P) - 1, np.ceil((tg.max() - T0) / P) + 2)
        return T0 + (n + phase) * P, n.astype(int)

    def window(self, tc):
        dt = self.t - tc
        g = self.good
        ih = max(IN_HALF, 0.55 * self.dt)
        inn = np.where(g & (np.abs(dt) <= ih))[0]
        lft = np.where(g & (dt <= -OOT_A) & (dt >= -OOT_B))[0]
        rgt = np.where(g & (dt >= OOT_A) & (dt <= OOT_B))[0]
        nside = (OOT_B - OOT_A) / self.dt
        nin = 2 * ih / self.dt
        ok = len(inn) >= max(1, 0.5 * nin) and len(lft) >= max(2, 0.4 * nside) and len(rgt) >= max(2, 0.4 * nside)
        return dict(tc=float(tc), inn=inn, lft=lft, rgt=rgt, ok=bool(ok))


def vshape(t, tc, dt, nsub=15):
    """Exposure-integrated triangular (V) eclipse profile with unit peak and total duration T14."""
    s = (np.arange(nsub) + 0.5) / nsub - 0.5
    a = np.abs(t[:, None] + s[None, :] * dt - tc)
    return np.clip(1 - a / (0.5 * T14), 0, None).mean(1)


def diff_image(C, w):
    """Per-pixel straight-line baseline fitted to both OOT sides, evaluated at the in-eclipse cadences.
    Returns D (baseline - in), sigma_D (from OOT residual scatter, propagated), B (baseline image), sbar."""
    oot = np.r_[w['lft'], w['rgt']]; inn = w['inn']; tc = w['tc']
    A = np.column_stack([np.ones(len(oot)), C.t[oot] - tc])
    Ai = np.column_stack([np.ones(len(inn)), C.t[inn] - tc])
    Y = C.f[oot].reshape(len(oot), -1)
    coef, *_ = np.linalg.lstsq(A, Y, rcond=None)
    res = Y - A @ coef
    rms = np.sqrt((res ** 2).sum(0) / max(len(oot) - 2, 1))
    pred = Ai @ coef
    D = (pred - C.f[inn].reshape(len(inn), -1)).mean(0)
    cinv = np.linalg.inv(A.T @ A)
    ab = Ai.mean(0)
    sig = rms * np.sqrt(1.0 / len(inn) + ab @ cinv @ ab)
    B = pred.mean(0)
    sb = vshape(C.t[inn], tc, C.dt).mean()
    shp = (C.ny, C.nx)
    return D.reshape(shp), sig.reshape(shp), B.reshape(shp), float(sb)


def build_windows(C, fake_phases=None):
    """All real deep (phase 0) / shallow (0.5) windows and fake windows at off-eclipse phases."""
    out = dict(D=[], S=[], B=[], sb=[], kind=[], tc=[], n=[])
    reals = []
    for kind, ph in (('deep', 0.0), ('shallow', 0.5)):
        tcs, ns = C.epochs(ph)
        reals += list(tcs)
        for tc, n in zip(tcs, ns):
            w = C.window(tc)
            if not w['ok']:
                continue
            D, s, B, sb = diff_image(C, w)
            for k, v in zip(('D', 'S', 'B', 'sb', 'kind', 'tc', 'n'), (D, s, B, sb, kind, tc, n)):
                out[k].append(v)
    reals = np.array(reals)
    if fake_phases is None:
        fake_phases = np.r_[np.arange(0.15, 0.351, 0.025), np.arange(0.65, 0.851, 0.025)]
    tcs, ns = C.epochs(0.0)
    for tc0, n in zip(tcs, ns):
        for ph in fake_phases:
            tc = tc0 + ph * P
            if np.min(np.abs(tc - reals)) < OOT_B + 0.5 * T14 + 0.003:
                continue
            w = C.window(tc)
            if not w['ok']:
                continue
            D, s, B, sb = diff_image(C, w)
            for k, v in zip(('D', 'S', 'B', 'sb', 'kind', 'tc', 'n'), (D, s, B, sb, 'fake', tc, n)):
                out[k].append(v)
    return {k: np.array(v) for k, v in out.items()}


def stack(W, idx, rel):
    """Inverse-variance stack of per-window D normalised by (rel * sbar) -> units of the peak amplitude of the
    reference (deep) eclipse.  rel = relative depth (1 deep, ~0.47 shallow).  B stacked with plain 1/var weights."""
    idx = np.asarray(idx)
    a = np.asarray(rel, float) * W['sb'][idx]
    v = np.array([np.median(s ** 2) for s in W['S'][idx]])
    wk = a ** 2 / v
    Dn = W['D'][idx] / a[:, None, None]
    Sn = W['S'][idx] / a[:, None, None]
    Ds = (wk[:, None, None] * Dn).sum(0) / wk.sum()
    Ss = np.sqrt((wk[:, None, None] ** 2 * Sn ** 2).sum(0)) / wk.sum()
    wb = 1 / v
    Bs = (wb[:, None, None] * W['B'][idx]).sum(0) / wb.sum()
    return Ds, Ss, Bs


class Model:
    """TESS PRF (PRF package) for the cutout; fits restricted to a (2*BOX+1)^2 box around the target pixel."""
    def __init__(self, C, kind='prf', gsig=None):
        import PRF
        sub = 's0001' if C.sector < 4 else 's0004'
        self.C = C
        self.kind = kind
        if kind == 'prf':
            self.p = PRF.TESS_PRF(C.cam, C.ccd, C.sector, C.col0 + C.xt, C.row0 + C.yt,
                                  localdatadir=f'{BASE}/prf_files/{sub}/')
        self.gsig = gsig
        self.shape = (C.ny, C.nx)
        b = BOX
        self.sl = (slice(C.iy - b, C.iy + b + 1), slice(C.ix - b, C.ix + b + 1))
        yy, xx = np.mgrid[self.sl]
        self.xx = xx.ravel(); self.yy = yy.ravel()
        # Gaia neighbours (T < 19) that can put light in the box, tied to the target with Gaia flux ratios
        g = C.g
        dpx = np.hypot(g.x - C.xt, g.y - C.yt)
        sel = np.where((g['T'] < 19.0) & (dpx < b + 5) & (np.arange(len(g)) != C.it))[0]
        self.nb = [(float(g.x.values[i] - C.xt), float(g.y.values[i] - C.yt),
                    float(10 ** (-0.4 * (g['T'].values[i] - g['T'].values[C.it])))) for i in sel]

    def img(self, x, y):
        if self.kind == 'prf':
            return self.p.locate(x, y, self.shape)
        from scipy.special import erf
        sx, sy = self.gsig
        xe = np.arange(self.shape[1] + 1) - 0.5; ye = np.arange(self.shape[0] + 1) - 0.5
        cx = 0.5 * np.diff(erf((xe - x) / (np.sqrt(2) * sx))); cy = 0.5 * np.diff(erf((ye - y) / (np.sqrt(2) * sy)))
        return np.outer(cy, cx)

    def box(self, im):
        return im[self.sl].ravel()

    def nbimg(self, x, y):
        m = np.zeros(self.shape)
        for dx, dy, r in self.nb:
            m += r * self.img(x + dx, y + dy)
        return m

    def fit(self, im, sig, direct=False, x0=None, y0=None, plane=False, maxoff=2.5):
        """Variable-projection PRF fit: im = a * [PRF(x,y) (+ tied neighbours if direct)] + c (+ plane)."""
        b = self.box(im); w = 1.0 / self.box(sig)
        x0 = self.C.xt if x0 is None else x0; y0 = self.C.yt if y0 is None else y0
        ones = np.ones_like(b)
        extra = [ones] + ([self.xx - self.C.xt, self.yy - self.C.yt] if plane else [])

        def design(p):
            m = self.img(p[0], p[1])
            if direct:
                m = m + self.nbimg(p[0], p[1])
            return np.column_stack([self.box(m)] + extra) * w[:, None]

        def res(p):
            A = design(p)
            c, *_ = np.linalg.lstsq(A, b * w, rcond=None)
            return A @ c - b * w

        best = None
        for gx in (-0.5, 0, 0.5):
            for gy in (-0.5, 0, 0.5):
                p = np.array([x0 + gx, y0 + gy]) if not direct else np.array([x0 + 0.4 * gx, y0 + 0.4 * gy])
                c2 = (res(p) ** 2).sum()
                if best is None or c2 < best[0]:
                    best = (c2, p)
        lo = [x0 - maxoff, y0 - maxoff]; hi = [x0 + maxoff, y0 + maxoff]
        r = least_squares(res, best[1], bounds=(lo, hi), diff_step=2e-4, x_scale=0.05, xtol=1e-10, ftol=1e-12, gtol=1e-12)
        A = design(r.x)
        c, *_ = np.linalg.lstsq(A, b * w, rcond=None)
        try:
            cov = np.linalg.inv(r.jac.T @ r.jac)
        except np.linalg.LinAlgError:
            cov = np.full((2, 2), np.nan)
        return dict(x=float(r.x[0]), y=float(r.x[1]), amp=float(c[0]), const=float(c[1]), chi2=float((r.fun ** 2).sum()),
                    npix=len(b), cov=cov, bound=bool(np.any(r.active_mask != 0)))

    def chi2_fixed(self, im, sig, x, y):
        b = self.box(im); w = 1.0 / self.box(sig)
        A = np.column_stack([self.box(self.img(x, y)), np.ones_like(b)]) * w[:, None]
        c, *_ = np.linalg.lstsq(A, b * w, rcond=None)
        return float(((A @ c - b * w) ** 2).sum()), float(c[0])

    def empirical(self, B, sig):
        """Unit-flux empirical target image from a direct image (background and tied neighbours removed)."""
        r = self.fit(B, sig, direct=True)
        E = (B - r['const'] - r['amp'] * self.nbimg(r['x'], r['y'])) / r['amp']
        return E, r
