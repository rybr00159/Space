"""Pixel-level localisation of the P = 1.3037 d eclipsing binary blended with TIC 71316387.

Adapted from /home/user/Space/tess-search/tic87055304/analysis/loc_independent/{common,pipeline}.py
(read-only there; this is a modified copy).

Conventions
-----------
* Pixel coords 0-based (x = column, y = row); integers = pixel centres (astropy WCS origin 0, PRF.locate).
* Images indexed [row, col] = [y, x].
* Difference image D = (OOT baseline predicted at in-eclipse times) - (in-eclipse flux) -> dimming star = positive PRF.
* Sky offsets (dE, dN) in arcsec, E = +RA*cos(dec).
* Time BTJD = BJD_TDB - 2457000.
"""
import os, glob, json
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares, lsq_linear

L = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise'
R = L + '/v2'
SECTORS = [1, 28, 68, 95, 105, 106]
P = 1.3036712
T0 = 3165.6697            # deep eclipse (BTJD)
T14 = 1.60 / 24.0         # eclipse duration (d); triangle fits to the TESScut aperture LC give 1.55-1.66 h
IN_HALF = float(os.environ.get('LOC_IN_HALF', 0.35)) * T14
OOT_IN = float(os.environ.get('LOC_OOT_IN', 1.0)) * T14
OOT_OUT = float(os.environ.get('LOC_OOT_OUT', 2.5)) * T14
BASELINE = os.environ.get('LOC_BASELINE', 'linear')
RFIT = float(os.environ.get('LOC_RFIT', 3.0))
TARGET_ID = 6571343840265900544
NB18_ID = 6571343844562758144
DEC0 = -42.835511
BAD = {1: [(1347.0, 1349.7)]}   # S1 pointing-jitter episode (see plots/s00_s1_raw.png)


def cutout_file(sector):
    return glob.glob(f'{L}/data/tesscut/s{sector}/*astrocut.fits')[0]


def tmag_from_gaia(g, bprp):
    g = np.asarray(g, float); c = np.asarray(bprp, float)
    t = g - 0.00522555 * c**3 + 0.0891337 * c**2 - 0.633923 * c + 0.0324473
    return np.where(np.isfinite(c), t, g - 0.430)


def btjd_to_jyear(t):
    return 2000.0 + (np.asarray(t) + 2457000.0 - 2451545.0) / 365.25


def load_gaia():
    d = pd.read_csv(f'{L}/gaia/field_378arcsec_G19.5.csv')
    d['tmag'] = tmag_from_gaia(d.phot_g_mean_mag, d.bp_rp)
    return d.reset_index(drop=True)


def propagate(d, epoch_year):
    dt = epoch_year - d.ref_epoch.values
    pmra = np.nan_to_num(d.pmra.values); pmdec = np.nan_to_num(d.pmdec.values)
    dec = d.dec.values + pmdec * dt / 3.6e6
    ra = d.ra.values + pmra * dt / 3.6e6 / np.cos(np.deg2rad(d.dec.values))
    return ra, dec


class Sector:
    def __init__(self, sector):
        self.sector = sector
        with fits.open(cutout_file(sector)) as h:
            self.hdr0 = h[0].header; self.hdr1 = h[1].header; self.hdr2 = h[2].header
            d = h[1].data
            self.t = np.array(d['TIME'], float)
            self.flux = np.array(d['FLUX'], float)
            self.q = np.array(d['QUALITY'], int)
        self.camera = int(self.hdr0['CAMERA']); self.ccd = int(self.hdr0['CCD'])
        self.col0 = int(self.hdr1['1CRV4P']); self.row0 = int(self.hdr1['2CRV4P'])
        self.wcs = WCS(self.hdr2)
        self.ny, self.nx = self.flux.shape[1:]
        fin = np.all(np.isfinite(self.flux.reshape(len(self.t), -1)), axis=1) & np.isfinite(self.t)
        self.good = (self.q == 0) & fin
        for a, b in BAD.get(sector, []):
            self.good &= ~((self.t > a) & (self.t < b))
        self.flux[~fin] = 0.0
        self.nbad_centroid = 0
        self.exptime = float(np.nanmedian(np.diff(self.t)))
        self.epoch = btjd_to_jyear(0.5 * (np.nanmin(self.t) + np.nanmax(self.t)))
        g = load_gaia()
        ra, dec = propagate(g, self.epoch)
        g['ra_ep'] = ra; g['dec_ep'] = dec
        x, y = self.wcs.all_world2pix(ra, dec, 0)
        g['x_wcs'] = x; g['y_wcs'] = y
        self.it = int(np.where(g.source_id == TARGET_ID)[0][0])
        self.ra_t, self.dec_t = ra[self.it], dec[self.it]
        g['de'] = (ra - self.ra_t) * np.cos(np.deg2rad(self.dec_t)) * 3600
        g['dn'] = (dec - self.dec_t) * 3600
        g['sep_t'] = np.hypot(g.de, g.dn)
        self.gaia = g
        self.jac = self.pix2sky_jacobian(x[self.it], y[self.it])
        self.jinv = np.linalg.inv(self.jac)
        self.xt_wcs, self.yt_wcs = float(x[self.it]), float(y[self.it])
        if os.environ.get('LOC_CENTCUT', '1') == '1':
            self._centroid_cut()

    def _centroid_cut(self, nmad=6.0, win=0.5):
        """Flag cadences whose target flux-weighted centroid deviates by > nmad*MAD from a running median
        (pointing-jitter excursions). An eclipse on a source ~1-5\" away moves the centroid by <~1 mpx,
        far below this threshold, so the cut cannot select on the eclipse itself."""
        m = circle_mask((self.ny, self.nx), self.xt_wcs, self.yt_wcs, 2.0)
        yy, xx = np.mgrid[0:self.ny, 0:self.nx]
        f = self.flux[:, m]
        fs = f.sum(1)
        with np.errstate(invalid='ignore', divide='ignore'):
            cx = (f * xx[m]).sum(1) / fs; cy = (f * yy[m]).sum(1) / fs
        bad = np.zeros(len(self.t), bool)
        g = np.where(self.good)[0]
        for c in (cx, cy):
            rm = np.full(len(self.t), np.nan)
            tg = self.t[g]; cg = c[g]
            for j, i in enumerate(g):
                s_ = np.abs(tg - self.t[i]) < win
                rm[i] = np.median(cg[s_])
            r = c - rm
            mad = 1.4826 * np.nanmedian(np.abs(r[g]))
            bad |= np.abs(r) > nmad * mad
        self.nbad_centroid = int((bad & self.good).sum())
        self.good &= ~bad
        self.cx, self.cy = cx, cy

    def pix2sky_jacobian(self, x, y, h=0.01):
        r0 = np.array(self.wcs.all_pix2world([[x, y]], 0))[0]
        rx = np.array(self.wcs.all_pix2world([[x + h, y]], 0))[0]
        ry = np.array(self.wcs.all_pix2world([[x, y + h]], 0))[0]
        def off(r):
            return np.array([(r[0] - r0[0]) * np.cos(np.deg2rad(r0[1])) * 3600, (r[1] - r0[1]) * 3600])
        return np.column_stack([off(rx) / h, off(ry) / h])   # [dE,dN] = J @ [dx,dy]

    def eclipse_epochs(self, phase=0.0, period=P, t0=T0):
        t = self.t[self.good]
        n = np.arange(np.ceil((t.min() - t0) / period) - 1, np.floor((t.max() - t0) / period) + 2)
        return t0 + (n + phase) * period, n.astype(int)

    def window(self, tc, min_in_frac=0.5, min_side_frac=0.35):
        dt = self.t - tc
        g = self.good
        inh = max(IN_HALF, 0.75 * self.exptime)      # 30-min S1: keep >=1-2 cadences in eclipse
        inn = np.where(g & (np.abs(dt) < inh))[0]
        lft = np.where(g & (dt <= -OOT_IN) & (dt >= -OOT_OUT))[0]
        rgt = np.where(g & (dt >= OOT_IN) & (dt <= OOT_OUT))[0]
        n_in_nom = 2 * inh / self.exptime
        n_side_nom = (OOT_OUT - OOT_IN) / self.exptime
        ok = (len(inn) >= max(1, min_in_frac * n_in_nom) and len(lft) >= max(2, min_side_frac * n_side_nom)
              and len(rgt) >= max(2, min_side_frac * n_side_nom))
        return dict(tc=float(tc), inn=inn, lft=lft, rgt=rgt, ok=bool(ok), n_in=len(inn), n_l=len(lft), n_r=len(rgt))


def eclipse_shape(t, tc, exptime, nsub=11):
    """Normalised V-shaped (triangular) eclipse profile, integrated over the exposure. Peak = 1."""
    sub = np.linspace(-0.5, 0.5, nsub) * exptime
    a = np.abs(t[:, None] + sub[None, :] - tc)
    return np.clip(1.0 - a / (T14 / 2), 0, 1).mean(1)


def difference_image(t, cube, win, regressors=None):
    """Per-pixel linear baseline on both OOT sides, evaluated at in-eclipse times.
    Returns D, sigma(D), baseline image at the in-eclipse times (direct image), OOT rms."""
    oot = np.concatenate([win['lft'], win['rgt']]); inn = win['inn']; tc = win['tc']
    A_oot = [np.ones(len(oot))] + ([t[oot] - tc] if BASELINE == 'linear' else [])
    A_in = [np.ones(len(inn))] + ([t[inn] - tc] if BASELINE == 'linear' else [])
    if regressors is not None:
        for k in range(regressors.shape[1]):
            A_oot.append(regressors[oot, k]); A_in.append(regressors[inn, k])
    A_oot = np.column_stack(A_oot); A_in = np.column_stack(A_in)
    Y = cube[oot].reshape(len(oot), -1)
    coef, *_ = np.linalg.lstsq(A_oot, Y, rcond=None)
    res = Y - A_oot @ coef
    dof = max(len(oot) - A_oot.shape[1], 1)
    rms = np.sqrt((res**2).sum(0) / dof)
    pred_in = (A_in @ coef).mean(0)
    obs_in = cube[inn].reshape(len(inn), -1).mean(0)
    D = pred_in - obs_in
    ainv = np.linalg.pinv(A_oot.T @ A_oot)
    abar = A_in.mean(0)
    var_pred = abar @ ainv @ abar
    sig = rms * np.sqrt(1.0 / len(inn) + var_pred)
    shp = cube.shape[1:]
    return D.reshape(shp), sig.reshape(shp), pred_in.reshape(shp), rms.reshape(shp)


class TessPRF:
    def __init__(self, S, xc=None, yc=None):
        import PRF
        sub = 's0001' if S.sector < 4 else 's0004'
        xc = S.xt_wcs if xc is None else xc; yc = S.yt_wcs if yc is None else yc
        self.p = PRF.TESS_PRF(S.camera, S.ccd, S.sector, S.col0 + xc, S.row0 + yc,
                              localdatadir=f'{L}/prf_files/{sub}/')
        self.shape = (S.ny, S.nx)
        self.name = 'TESS_PRF'

    def __call__(self, x, y):
        if x < -6.4 or y < -6.4 or x > self.shape[1] + 5.4 or y > self.shape[0] + 5.4:
            return np.zeros(self.shape)
        return self.p.locate(x, y, self.shape)


class GaussPSF:
    def __init__(self, S, sx=0.75, sy=None):
        from scipy.special import erf
        self.erf = erf
        self.shape = (S.ny, S.nx); self.sx = sx; self.sy = sx if sy is None else sy
        self.xe = np.arange(self.shape[1] + 1) - 0.5; self.ye = np.arange(self.shape[0] + 1) - 0.5
        self.name = f'Gauss({self.sx:.2f},{self.sy:.2f})'

    def __call__(self, x, y):
        cx = 0.5 * np.diff(self.erf((self.xe - x) / (np.sqrt(2) * self.sx)))
        cy = 0.5 * np.diff(self.erf((self.ye - y) / (np.sqrt(2) * self.sy)))
        return np.outer(cy, cx)


def circle_mask(shape, xc, yc, rad):
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    return np.hypot(xx - xc, yy - yc) <= rad


def fit_free(img, sig, psf, x0, y0, mask, fixed_img=None, fixed_free_cols=None, box=1.5, ngrid=5, grid=0.6):
    """One PSF at a free position + amplitude + constant [+ fixed_img*free scale? no: fixed_img is
    subtracted as a known component scaled by the same amplitude as the main PSF (tied flux ratio)]
    [+ fixed_free_cols: extra PSF images with free non-negative-agnostic linear amplitudes].
    Variable projection over (x, y)."""
    w = 1.0 / sig[mask]
    b = img[mask] * w
    extra = [] if fixed_free_cols is None else [c[mask] for c in fixed_free_cols]

    def lin(p):
        main = psf(p[0], p[1])
        if fixed_img is not None:
            main = main + fixed_img(p[0], p[1])
        cols = [main[mask]] + extra + [np.ones(mask.sum())]
        A = np.column_stack(cols) * w[:, None]
        coef, *_ = np.linalg.lstsq(A, b, rcond=None)
        return A, coef

    def resid(p):
        A, coef = lin(p)
        return A @ coef - b

    best = None
    for gx in np.linspace(-grid, grid, ngrid):
        for gy in np.linspace(-grid, grid, ngrid):
            p = np.array([x0 + gx, y0 + gy])
            c2 = (resid(p) ** 2).sum()
            if best is None or c2 < best[0]:
                best = (c2, p)
    lo = [x0 - box, y0 - box]; hi = [x0 + box, y0 + box]
    r = least_squares(resid, best[1], bounds=(lo, hi), diff_step=1e-4, x_scale=0.05, xtol=1e-10, ftol=1e-12)
    A, coef = lin(r.x)
    J = r.jac
    try:
        cov = np.linalg.inv(J.T @ J)
    except np.linalg.LinAlgError:
        cov = np.full((2, 2), np.nan)
    return dict(x=float(r.x[0]), y=float(r.x[1]), cov=cov, amp=float(coef[0]), coef=coef,
                chi2=float((r.fun ** 2).sum()), npix=int(mask.sum()), at_bound=bool(np.any(r.active_mask != 0)))


def to_json(obj, path):
    def conv(o):
        if isinstance(o, (np.floating,)): return float(o)
        if isinstance(o, (np.integer,)): return int(o)
        if isinstance(o, np.ndarray): return o.tolist()
        if isinstance(o, (np.bool_,)): return bool(o)
        raise TypeError(type(o))
    with open(path, 'w') as f:
        json.dump(obj, f, indent=1, default=conv)
