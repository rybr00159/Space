"""Referee helper library (independent of the original vetting scripts).

Loaders for cached light curves and TESScut cutouts, a PRF model wrapper,
Gaia neighbour helpers and a few small fitting utilities.
"""
import os, glob, json
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS

BASE = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/five_referee'
LC_DIR = '/home/user/work/lc_cache'
CUT_DIR = '/home/user/work/tesscut5'
PRF_DIR = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/prf_files/s0004/'

GAIA_ID = {90086266: 6722830841491505536, 61332759: 6711181721954396928,
           214311906: 5925373268528684288, 55829034: 4038825960551504768,
           61141412: 6712312986277027968}


def save_json(obj, fn):
    def conv(o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.bool_,)):
            return bool(o)
        raise TypeError(type(o))
    with open(fn, 'w') as f:
        json.dump(obj, f, indent=1, default=conv)


# ----------------------------------------------------------------- light curves
def lc_path(tic, sector, prod):
    if prod == 'QLP':
        p = f'{LC_DIR}/hlsp_qlp_tess_ffi_s{sector:04d}-{tic:016d}_tess_v01_llc.fits'
    elif prod == 'TESS-SPOC':
        p = f'{LC_DIR}/hlsp_tess-spoc_tess_phot_{tic:016d}-s{sector:04d}_tess_v1_lc.fits'
    elif prod == 'SPOC':
        g = glob.glob(f'{LC_DIR}/tess*-s{sector:04d}-{tic:016d}-*-s_lc.fits')
        p = g[0] if g else ''
    else:
        raise ValueError(prod)
    return p if os.path.exists(p) else None


def load_lc(tic, sector, prod, col):
    """Return t, f (normalised by median), e, extra dict. QUALITY==0 and finite only."""
    p = lc_path(tic, sector, prod)
    if p is None:
        return None
    with fits.open(p) as h:
        d = h[1].data
        names = d.columns.names
        if col not in names:
            # QLP old sectors use KSPSAP instead of DET
            alt = {'DET_FLUX': 'KSPSAP_FLUX', 'KSPSAP_FLUX': 'DET_FLUX'}.get(col)
            if alt in names:
                col = alt
            else:
                return None
        t = np.array(d['TIME'], float)
        f = np.array(d[col], float)
        q = np.array(d['QUALITY'])
        ecol = col + '_ERR'
        e = np.array(d[ecol], float) if ecol in names else np.full_like(f, np.nan)
        extra = {}
        for c in ['SAP_BKG', 'SAP_X', 'SAP_Y', 'MOM_CENTR1', 'MOM_CENTR2', 'POS_CORR1', 'POS_CORR2', 'ORBITID']:
            if c in names:
                extra[c] = np.array(d[c], float)
    m = (q == 0) & np.isfinite(t) & np.isfinite(f)
    med = np.nanmedian(f[m])
    t, f, e = t[m], f[m] / med, e[m] / med
    if not np.all(np.isfinite(e)):
        # robust point-to-point scatter as error
        s = 1.4826 * np.nanmedian(np.abs(np.diff(f))) / np.sqrt(2)
        e = np.full_like(f, s)
    for k in extra:
        extra[k] = extra[k][m]
    extra['col'] = col
    extra['median'] = med
    return t, f, e, extra


# ----------------------------------------------------------------- cutouts
def cut_path(tic, sector):
    g = glob.glob(f'{CUT_DIR}/{tic}/tess-s{sector:04d}-*_astrocut.fits')
    return g[0] if g else None


def load_cut(tic, sector, qual0=True):
    p = cut_path(tic, sector)
    with fits.open(p) as h:
        d = h[1].data
        hd0 = h[0].header
        hd1 = h[1].header
        wcs = WCS(h[2].header)
        out = dict(time=np.array(d['TIME'], float), flux=np.array(d['FLUX'], float),
                   ferr=np.array(d['FLUX_ERR'], float), bkg=np.array(d['FLUX_BKG'], float),
                   quality=np.array(d['QUALITY']), pc1=np.array(d['POS_CORR1'], float),
                   pc2=np.array(d['POS_CORR2'], float), cadno=np.array(d['CADENCENO']))
        out['cam'] = hd0['CAMERA']
        out['ccd'] = hd0['CCD']
        out['sector'] = hd0['SECTOR']
        out['col0'] = hd1['1CRV4P']  # CCD column of 0-based pixel x=0
        out['row0'] = hd1['2CRV4P']
        out['wcs'] = wcs
    good = np.isfinite(out['time']) & np.all(np.isfinite(out['flux'].reshape(len(out['time']), -1)), axis=1)
    if qual0:
        good &= out['quality'] == 0
    for k in ['time', 'flux', 'ferr', 'bkg', 'quality', 'pc1', 'pc2', 'cadno']:
        out[k] = out[k][good]
    return out


# ----------------------------------------------------------------- PRF
class PRFModel:
    def __init__(self, cut):
        import PRF
        ny, nx = cut['flux'].shape[1:]
        self.shape = (ny, nx)
        cc = cut['col0'] + nx / 2.
        rr = cut['row0'] + ny / 2.
        self.prf = PRF.TESS_PRF(cut['cam'], cut['ccd'], 4, cc, rr, localdatadir=PRF_DIR)

    def __call__(self, x, y):
        """x, y: 0-based cutout pixel coords (integer = pixel centre)."""
        return self.prf.locate(x, y, self.shape)


# ----------------------------------------------------------------- Gaia helpers
def tmag_from_gaia(g, bprp):
    g = np.asarray(g, float)
    c = np.asarray(bprp, float)
    t = g - 0.00522555 * c ** 3 + 0.0891337 * c ** 2 - 0.633923 * c + 0.0324473
    t = np.where(np.isfinite(c), t, g - 0.430)
    return t


def gaia_table(tic, epoch_jyear):
    """Gaia neighbours with positions propagated to epoch (Julian year)."""
    d = pd.read_csv(f'{BASE}/data/gaia_{tic}.csv')
    dt = epoch_jyear - d['ref_epoch'].values
    pmra = np.nan_to_num(d['pmra'].values)
    pmdec = np.nan_to_num(d['pmdec'].values)
    dec = d['dec'].values + pmdec * dt / 3.6e6
    ra = d['ra'].values + pmra * dt / 3.6e6 / np.cos(np.radians(d['dec'].values))
    d['ra_ep'] = ra
    d['dec_ep'] = dec
    d['tmag'] = tmag_from_gaia(d['phot_g_mean_mag'], d['bp_rp'])
    tgt = d[d.source_id == GAIA_ID[tic]].iloc[0]
    d['dE'] = (d['ra_ep'] - tgt['ra_ep']) * np.cos(np.radians(tgt['dec_ep'])) * 3600.
    d['dN'] = (d['dec_ep'] - tgt['dec_ep']) * 3600.
    d['sep_t'] = np.hypot(d['dE'], d['dN'])
    d = d.sort_values('sep_t').reset_index(drop=True)
    return d


def btjd_to_jyear(t):
    return 2000.0 + (t + 2457000.0 - 2451545.0) / 365.25


# ----------------------------------------------------------------- misc fitting
def robust_std(x):
    x = np.asarray(x)
    x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def trapezoid(t, t0, depth, T14, T12):
    """Symmetric trapezoid; T12 = ingress duration (T12 = T14/2 gives a V)."""
    T12 = np.clip(T12, 1e-5, T14 / 2.)
    x = np.abs(t - t0)
    m = np.zeros_like(t)
    full = x <= T14 / 2. - T12
    ing = (x > T14 / 2. - T12) & (x < T14 / 2.)
    m[full] = depth
    m[ing] = depth * (T14 / 2. - x[ing]) / T12
    return 1. - m


def supersample(func, t, texp, n=7, **kw):
    if texp is None or texp < 1e-4:
        return func(t, **kw)
    off = (np.arange(n) - (n - 1) / 2.) / n * texp
    tt = (t[:, None] + off[None, :]).ravel()
    return func(tt, **kw).reshape(len(t), n).mean(axis=1)


def aperture_lc(tic, sector, r=1.5, center=None, sky_r=(6, 9), qual0=True):
    """Simple circular-aperture photometry on the TESScut cutout around the target
    (Gaia position propagated to the sector epoch) or a given (ra, dec) centre.
    Returns t, normalised flux, error, extra dict (bkg, centroids, pos_corr, raw flux)."""
    cut = load_cut(tic, sector, qual0=qual0)
    if center is None:
        g = gaia_table(tic, btjd_to_jyear(np.median(cut['time'])))
        tgt = g.iloc[0]
        center = (tgt['ra_ep'], tgt['dec_ep'])
    x, y = [float(v) for v in cut['wcs'].world_to_pixel_values(center[0], center[1])]
    ny, nx = cut['flux'].shape[1:]
    yy, xx = np.mgrid[:ny, :nx]
    rr = np.hypot(xx - x, yy - y)
    ap = rr <= r
    sky = (rr > sky_r[0]) & (rr < sky_r[1])
    fl = cut['flux']
    bk = np.nanmedian(fl[:, sky], axis=1)
    f = np.nansum(fl[:, ap], axis=1) - bk * ap.sum()
    ap2 = rr <= 2.5
    w = fl[:, ap2] - bk[:, None]
    cx = np.sum(w * xx[ap2], axis=1) / np.sum(w, axis=1)
    cy = np.sum(w * yy[ap2], axis=1) / np.sum(w, axis=1)
    med = np.median(f)
    extra = dict(SAP_BKG=bk + np.nanmedian(cut['bkg'][:, ap], axis=1), MOM_CENTR1=cx, MOM_CENTR2=cy,
                 POS_CORR1=cut['pc1'], POS_CORR2=cut['pc2'], col='CUT_r%.1f' % r, median=med, raw=f,
                 quality=cut['quality'], npix=int(ap.sum()), xy=(x, y))
    e = np.full_like(f, robust_std(np.diff(f / med)) / np.sqrt(2))
    return cut['time'], f / med, e, extra
