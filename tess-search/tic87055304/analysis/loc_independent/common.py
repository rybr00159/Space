"""Shared helpers for the independent pixel-level localisation of the TIC 87055304 transit.

Conventions
-----------
* Pixel coordinates are 0-based (x = column, y = row); integer values are pixel centres.
  This matches astropy.wcs (origin 0) and PRF.TESS_PRF.locate.
* Image arrays are indexed [row, col] = [y, x].
* Difference image D = (out-of-transit baseline predicted at in-transit times) - (in-transit flux),
  so a star that dims gives a POSITIVE PRF-shaped bump.
* Time is BTJD = BJD_TDB - 2457000.
"""
import os
import glob
import json
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.special import erf
from scipy.optimize import least_squares, lsq_linear

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = '/home/user/work/tic87/data'
TESSCUT = sorted(glob.glob(f'{DATA}/tesscut/tess-s*_21x21_astrocut.fits'))
SECTORS = [13, 39, 66, 93, 104]

P = 4.136195
T0 = 1658.41325
T14 = 0.1875
IN_HALF = float(os.environ.get('LOC_IN_HALF', 0.35)) * T14   # central 70 % of T14
OOT_IN = float(os.environ.get('LOC_OOT_IN', 1.0)) * T14      # T14/2 + 0.5*T14 gap from the transit edge
OOT_OUT = float(os.environ.get('LOC_OOT_OUT', 2.5)) * T14    # baseline extends 1.5*T14 beyond the gap on each side
BASELINE = os.environ.get('LOC_BASELINE', 'linear')          # 'linear' or 'const' per-pixel OOT baseline
BAD_EPISODE = (4178.8, 4179.2)  # S104 scattered-light episode

TARGET_ID = 6720986204578421888
NEIGH_ID = 6720986208895208960
C3_ID = 6720986208880830336   # G=17.96 (T~17.35) star 11.3" from target, 3.4" from the neighbour
DEC0 = -43.995371
COSD = np.cos(np.deg2rad(DEC0))


def cutout_file(sector):
    return [f for f in TESSCUT if f'tess-s{sector:04d}-' in f][0]


def tmag_from_gaia(g, bprp):
    """Stassun et al. (2019) TIC-8 relation; G-0.430 when no colour."""
    g = np.asarray(g, float)
    c = np.asarray(bprp, float)
    t = g - 0.00522555 * c**3 + 0.0891337 * c**2 - 0.633923 * c + 0.0324473
    return np.where(np.isfinite(c), t, g - 0.430)


def btjd_to_jyear(t):
    return 2000.0 + (np.asarray(t) + 2457000.0 - 2451545.0) / 365.25


def load_gaia(include_wide=True):
    """Gaia DR3 list: full 90" CSV, plus G<17.5 stars from 90"-330" (wide query) if present."""
    d = pd.read_csv(f'{DATA}/gaia_dr3_90arcsec.csv')
    cols = ['source_id', 'ra', 'dec', 'pmra', 'pmdec', 'parallax', 'ref_epoch', 'phot_g_mean_mag',
            'phot_bp_mean_mag', 'phot_rp_mean_mag', 'bp_rp', 'ruwe']
    d = d[cols].copy()
    d['src'] = 'csv90'
    wfile = f'{HERE}/gaia_wide/gaia_dr3_330arcsec_G17p5_vizier.csv'
    if include_wide and os.path.exists(wfile):
        w = pd.read_csv(wfile)
        w = w[~w.source_id.isin(d.source_id)][cols].copy()
        w['src'] = 'wide330'
        d = pd.concat([d, w], ignore_index=True)
    d['tmag'] = tmag_from_gaia(d.phot_g_mean_mag, d.bp_rp)
    return d.reset_index(drop=True)


def propagate(d, epoch_year):
    """Propagate Gaia positions from ref_epoch to epoch_year (linear, pm only)."""
    dt = epoch_year - d.ref_epoch.values
    pmra = np.nan_to_num(d.pmra.values)
    pmdec = np.nan_to_num(d.pmdec.values)
    dec = d.dec.values + pmdec * dt / 3.6e6
    ra = d.ra.values + pmra * dt / 3.6e6 / np.cos(np.deg2rad(d.dec.values))
    return ra, dec


def sky_offset_arcsec(ra, dec, ra0, dec0):
    """Tangent-plane offsets (east, north) in arcsec of (ra,dec) from (ra0,dec0)."""
    de = (np.asarray(ra) - ra0) * np.cos(np.deg2rad(dec0)) * 3600.0
    dn = (np.asarray(dec) - dec0) * 3600.0
    return de, dn


class Sector:
    """TESScut cube, quality mask, WCS and Gaia star pixel positions for one sector."""

    def __init__(self, sector, include_wide=True):
        self.sector = sector
        f = cutout_file(sector)
        with fits.open(f) as h:
            self.hdr0 = h[0].header
            self.hdr1 = h[1].header
            self.hdr2 = h[2].header
            d = h[1].data
            self.t = np.array(d['TIME'], float)
            self.flux = np.array(d['FLUX'], float)
            self.ferr = np.array(d['FLUX_ERR'], float)
            self.bkg = np.array(d['FLUX_BKG'], float)
            self.q = np.array(d['QUALITY'], int)
        self.camera = int(self.hdr0['CAMERA'])
        self.ccd = int(self.hdr0['CCD'])
        self.col0 = int(self.hdr1['1CRV4P'])   # CCD column of stamp pixel x=0
        self.row0 = int(self.hdr1['2CRV4P'])   # CCD row of stamp pixel y=0
        self.wcs = WCS(self.hdr2)
        self.ny, self.nx = self.flux.shape[1:]
        fin = np.all(np.isfinite(self.flux.reshape(len(self.t), -1)), axis=1) & np.isfinite(self.t)
        self.good = (self.q == 0) & fin
        self.good &= ~((self.t > BAD_EPISODE[0]) & (self.t < BAD_EPISODE[1]))
        self.exptime = float(np.nanmedian(np.diff(self.t)))
        self.epoch = btjd_to_jyear(0.5 * (np.nanmin(self.t) + np.nanmax(self.t)))
        self.gaia = load_gaia(include_wide)
        ra, dec = propagate(self.gaia, self.epoch)
        self.gaia['ra_ep'] = ra
        self.gaia['dec_ep'] = dec
        x, y = self.wcs.all_world2pix(ra, dec, 0)
        self.gaia['x_wcs'] = x
        self.gaia['y_wcs'] = y
        it = np.where(self.gaia.source_id == TARGET_ID)[0][0]
        ineigh = np.where(self.gaia.source_id == NEIGH_ID)[0][0]
        self.itarget, ineigh_ = it, ineigh
        self.ineigh = ineigh
        self.ra_t, self.dec_t = ra[it], dec[it]
        de, dn = sky_offset_arcsec(ra, dec, self.ra_t, self.dec_t)
        self.gaia['de'] = de
        self.gaia['dn'] = dn
        self.gaia['sep_t'] = np.hypot(de, dn)
        # local Jacobian: pixel -> sky (arcsec east, north) at the target
        self.jac = self.pix2sky_jacobian(x[it], y[it])

    def pix2sky_jacobian(self, x, y, h=0.01):
        r0 = np.array(self.wcs.all_pix2world([[x, y]], 0))[0]
        rx = np.array(self.wcs.all_pix2world([[x + h, y]], 0))[0]
        ry = np.array(self.wcs.all_pix2world([[x, y + h]], 0))[0]
        def off(r):
            return np.array([(r[0] - r0[0]) * np.cos(np.deg2rad(r0[1])) * 3600, (r[1] - r0[1]) * 3600])
        return np.column_stack([off(rx) / h, off(ry) / h])   # 2x2: [dE,dN] = J @ [dx,dy]

    def transit_epochs(self, period=P, t0=T0, shift=0.0):
        t = self.t[self.good]
        n = np.arange(np.ceil((t.min() - t0) / period) - 1, np.floor((t.max() - t0) / period) + 2)
        return t0 + (n + shift) * period, n.astype(int)

    def window(self, tc, min_in_frac=0.5, min_side_frac=0.25):
        """Index arrays of in-transit and out-of-transit (both sides) cadences, or None."""
        dt = self.t - tc
        g = self.good
        inn = np.where(g & (np.abs(dt) < IN_HALF))[0]
        lft = np.where(g & (dt <= -OOT_IN) & (dt >= -OOT_OUT))[0]
        rgt = np.where(g & (dt >= OOT_IN) & (dt <= OOT_OUT))[0]
        n_in_nom = 2 * IN_HALF / self.exptime
        n_side_nom = (OOT_OUT - OOT_IN) / self.exptime
        ok = (len(inn) >= min_in_frac * n_in_nom and len(lft) >= min_side_frac * n_side_nom
              and len(rgt) >= min_side_frac * n_side_nom)
        return dict(tc=tc, inn=inn, lft=lft, rgt=rgt, ok=bool(ok),
                    n_in=len(inn), n_l=len(lft), n_r=len(rgt))


def difference_image(t, cube, win, extra_cube=None, regressors=None):
    """Per-pixel linear-in-time baseline fitted to both OOT sides, evaluated at the
    in-transit cadences; D = mean(baseline_in) - mean(flux_in).
    Returns D, sigma(D) (per pixel, from OOT residual rms) and the OOT residual rms.
    regressors: optional (ncad, k) array of extra baseline regressors (e.g. pointing)."""
    oot = np.concatenate([win['lft'], win['rgt']])
    inn = win['inn']
    tc = win['tc']
    c = cube if extra_cube is None else cube + extra_cube
    A_oot = [np.ones(len(oot))] + ([t[oot] - tc] if BASELINE == 'linear' else [])
    A_in = [np.ones(len(inn))] + ([t[inn] - tc] if BASELINE == 'linear' else [])
    if regressors is not None:
        for k in range(regressors.shape[1]):
            r = regressors[:, k]
            A_oot.append(r[oot]); A_in.append(r[inn])
    A_oot = np.column_stack(A_oot)
    A_in = np.column_stack(A_in)
    Y = c[oot].reshape(len(oot), -1)
    coef, *_ = np.linalg.lstsq(A_oot, Y, rcond=None)
    res = Y - A_oot @ coef
    dof = max(len(oot) - A_oot.shape[1], 1)
    rms = np.sqrt((res**2).sum(0) / dof)
    pred_in = (A_in @ coef).mean(0)
    obs_in = c[inn].reshape(len(inn), -1).mean(0)
    D = pred_in - obs_in
    # variance: in-transit mean noise + baseline prediction noise
    ainv = np.linalg.pinv(A_oot.T @ A_oot)
    abar = A_in.mean(0)
    var_pred = abar @ ainv @ abar
    sig = rms * np.sqrt(1.0 / len(inn) + var_pred)
    shp = cube.shape[1:]
    return D.reshape(shp), sig.reshape(shp), rms.reshape(shp)


# ---------------------------------------------------------------- PSF models
class TessPRF:
    """Wrapper around PRF.TESS_PRF (Bell & Higgins) using locally downloaded MAST PRF files."""

    def __init__(self, sec):
        import PRF
        self.p = PRF.TESS_PRF(sec.camera, sec.ccd, sec.sector, sec.col0 + sec.nx / 2, sec.row0 + sec.ny / 2,
                              localdatadir=f'{HERE}/prf_files/')
        self.shape = (sec.ny, sec.nx)
        self.name = 'TESS_PRF'

    def __call__(self, x, y):
        # locate() fails if the source is far outside the stamp; clip to a safe range
        if x < -6.4 or y < -6.4 or x > self.shape[1] + 5.4 or y > self.shape[0] + 5.4:
            return np.zeros(self.shape)
        return self.p.locate(x, y, self.shape)


class GaussPSF:
    """Pixel-integrated circular (or elliptical) Gaussian PSF; integer coordinates = pixel centres."""

    def __init__(self, sec, sx=0.75, sy=None, rho=0.0):
        self.shape = (sec.ny, sec.nx)
        self.sx = sx
        self.sy = sx if sy is None else sy
        self.name = f'Gauss(sx={self.sx:.3f},sy={self.sy:.3f})'
        self.xe = np.arange(self.shape[1] + 1) - 0.5
        self.ye = np.arange(self.shape[0] + 1) - 0.5

    def __call__(self, x, y):
        cx = 0.5 * np.diff(erf((self.xe - x) / (np.sqrt(2) * self.sx)))
        cy = 0.5 * np.diff(erf((self.ye - y) / (np.sqrt(2) * self.sy)))
        return np.outer(cy, cx)


# ---------------------------------------------------------------- fitting
def fit_fixed_positions(img, sig, psf, xs, ys, mask, bkg=True, nonneg=None):
    """Linear least squares for amplitudes at fixed positions (+ constant).
    Returns amplitudes, covariance, chi2, model."""
    cols = [psf(x, y)[mask] for x, y in zip(xs, ys)]
    if bkg:
        cols.append(np.ones(mask.sum()))
    A = np.column_stack(cols)
    w = 1.0 / sig[mask]
    Aw = A * w[:, None]
    bw = img[mask] * w
    if nonneg is not None:
        lb = np.full(A.shape[1], -np.inf)
        lb[nonneg] = 0.0
        r = lsq_linear(Aw, bw, bounds=(lb, np.inf), method='bvls')
        coef = r.x
    else:
        coef, *_ = np.linalg.lstsq(Aw, bw, rcond=None)
    cov = np.linalg.pinv(Aw.T @ Aw)
    model_m = A @ coef
    chi2 = float((((img[mask] - model_m) * w) ** 2).sum())
    full = np.zeros_like(img)
    for k, (x, y) in enumerate(zip(xs, ys)):
        full += coef[k] * psf(x, y)
    if bkg:
        full += coef[-1]
    return coef, cov, chi2, full


def fit_free_position(img, sig, psf, x0, y0, mask, bkg=True, fixed=None, box=2.5, plane=False):
    """One PSF at a free position (+ amplitude, constant, and optional fixed-position
    nuisance sources with free amplitudes). Variable projection: linear amplitudes are
    solved for each trial position.  Returns dict with x, y, cov (x,y), amp, chi2."""
    w = 1.0 / sig[mask]
    b = img[mask] * w
    fixed_cols = []
    if plane:
        yy, xx = np.mgrid[0:img.shape[0], 0:img.shape[1]]
        fixed_cols += [(xx - x0)[mask], (yy - y0)[mask]]
    if fixed is not None:
        for (xf, yf) in fixed:
            fixed_cols.append(psf(xf, yf)[mask])

    def lin(p):
        cols = [psf(p[0], p[1])[mask]] + fixed_cols
        if bkg:
            cols.append(np.ones(mask.sum()))
        A = np.column_stack(cols) * w[:, None]
        coef, *_ = np.linalg.lstsq(A, b, rcond=None)
        return A, coef

    def resid(p):
        A, coef = lin(p)
        return A @ coef - b

    # coarse grid start to avoid local minima
    best = None
    gr = min(1.5, 0.9 * box)
    for gx in np.linspace(-gr, gr, 7):
        for gy in np.linspace(-gr, gr, 7):
            p = np.array([x0 + gx, y0 + gy])
            c2 = (resid(p) ** 2).sum()
            if best is None or c2 < best[0]:
                best = (c2, p)
    lo = [x0 - box, y0 - box]
    hi = [x0 + box, y0 + box]
    r = least_squares(resid, best[1], bounds=(lo, hi), diff_step=1e-3, x_scale=0.1)
    A, coef = lin(r.x)
    J = r.jac
    try:
        cov = np.linalg.inv(J.T @ J)
    except np.linalg.LinAlgError:
        cov = np.full((2, 2), np.nan)
    return dict(x=float(r.x[0]), y=float(r.x[1]), cov=cov, amp=float(coef[0]), coef=coef,
                chi2=float((r.fun ** 2).sum()), npix=int(mask.sum()), at_bound=bool(np.any(r.active_mask != 0)))


def circle_mask(shape, xc, yc, rad):
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    return np.hypot(xx - xc, yy - yc) <= rad


def to_json(obj, path):
    def conv(o):
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.bool_,)):
            return bool(o)
        raise TypeError(type(o))
    with open(path, 'w') as f:
        json.dump(obj, f, indent=1, default=conv)
