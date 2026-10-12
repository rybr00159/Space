"""Shared helpers for the TIC 87055304 PSF-photometry investigation (psf_photometry/).

- load_sector: TESScut cube, quality mask, WCS, CCD geometry.
- gaia_stars: Gaia DR3 stars (VizieR, epoch 2016.0) propagated to the sector epoch with proper motions and converted to
  cutout pixel coordinates (0-based, integer = pixel centre; same convention as PRF.TESS_PRF.locate).
- PRFModel: TESS_PRF (Bell & Higgins) interpolated to the target's CCD position, optionally smeared by a Gaussian
  (pointing jitter / focus) of width sig_px applied to the 9x-supersampled PRF, i.e. a convolution in source position.
- template: unit-depth transit template (batman, shape from the SAP fit in CONTEXT.md).
"""
import os, glob
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter
import PRF

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = "/home/user/work/tic87/psf_photometry"  # read-only: prf_files/, gaia_330as_G19.csv
DATA = "/home/user/work/tic87/data"
SECTORS = [13, 39, 66, 93, 104]
P, T0, DUR = 4.136195, 1658.41325, 0.1875
TGT = 6720986204578421888
NB = 6720986208895208960
PRFDIR = os.path.join(ROOT, "prf_files")


def tmag_from_gaia(G, bp, rp):
    """Stassun et al. (2019) TIC-8 relation; G - 0.43 when no colour."""
    c = bp - rp
    T = G - 0.00522555 * c ** 3 + 0.0891337 * c ** 2 - 0.633923 * c + 0.0324473
    return np.where(np.isfinite(c) & (c > -1) & (c < 6), T, G - 0.43)


def load_sector(s):
    f = glob.glob(f"{DATA}/tesscut/tess-s{s:04d}-*_21x21_astrocut.fits")[0]
    with fits.open(f) as h:
        d = h[1].data
        t = np.array(d["TIME"], float)
        F = np.array(d["FLUX"], float)
        E = np.array(d["FLUX_ERR"], float)
        q = np.array(d["QUALITY"])
        cad = np.array(d["CADENCENO"])
        h1 = h[1].header
        cam, ccd = h[0].header["CAMERA"], h[0].header["CCD"]
        col0, row0 = h1["1CRV4P"], h1["2CRV4P"]
        w = WCS(h[2].header)
    good = (q == 0) & np.isfinite(t) & np.all(np.isfinite(F.reshape(len(t), -1)), axis=1)
    exptime = np.nanmedian(np.diff(t[good])) * 86400
    return dict(sector=s, t=t[good], F=F[good], E=E[good], q=q[good], cad=cad[good], cam=cam, ccd=ccd,
                col0=col0, row0=row0, wcs=w, exptime=exptime, file=f)


def epoch_year(t):
    return 2000.0 + (np.nanmean(t) + 2457000.0 - 2451545.0) / 365.25


def gaia_stars(sec, gmax=19.0):
    g = pd.read_csv(os.path.join(ROOT, "gaia_330as_G19.csv"))
    g = g[g.Gmag < gmax].copy()
    dt = epoch_year(sec["t"]) - 2016.0
    pmra = np.nan_to_num(g.pmRA.values); pmde = np.nan_to_num(g.pmDE.values)
    ra = g.RA_ICRS.values + pmra * dt / 3.6e6 / np.cos(np.radians(g.DE_ICRS.values))
    de = g.DE_ICRS.values + pmde * dt / 3.6e6
    x, y = sec["wcs"].world_to_pixel_values(ra, de)
    g["ra_ep"], g["de_ep"], g["x"], g["y"] = ra, de, x, y
    g["T"] = tmag_from_gaia(g.Gmag.values, g.BPmag.values, g.RPmag.values)
    it = np.where(g.Source.values == TGT)[0][0]
    g["dx_t"] = g.x - g.x.values[it]; g["dy_t"] = g.y - g.y.values[it]
    g["r_t"] = np.hypot(g.dx_t, g.dy_t)
    # keep stars whose PRF can touch the 21x21 stamp
    g = g[(g.x > -6) & (g.x < 26) & (g.y > -6) & (g.y < 26)].reset_index(drop=True)
    return g


class PRFModel:
    """TESS_PRF at the target's CCD position; optional Gaussian smear sig_px (pixels) of the source position."""

    def __init__(self, sec, sig_px=0.0, ccdpos=None, xt=10.0, yt=10.0):
        col = sec["col0"] + xt if ccdpos is None else ccdpos[0]
        row = sec["row0"] + yt if ccdpos is None else ccdpos[1]
        # TESS_PRF drops any file whose path contains 'phot' (our folder name does), so pass a relative path
        cwd = os.getcwd(); os.chdir(ROOT)
        try:
            self.base = PRF.TESS_PRF(sec["cam"], sec["ccd"], sec["sector"], col, row, localdatadir="prf_files")
        finally:
            os.chdir(cwd)
        self.raw = self.base.prf.copy()
        self.set_sigma(sig_px)

    def set_sigma(self, sig_px):
        self.sig = float(sig_px)
        prf = gaussian_filter(self.raw, self.sig * 9.0, mode="constant") if self.sig > 0 else self.raw
        self.base.prf = prf
        reshaped = np.zeros((11, 11, 13, 13))
        for i in range(9):
            for j in range(9):
                reshaped[i + 1, j + 1, :, :] = prf[8 - i::9, 8 - j::9]
        for j in range(1, 10):
            reshaped[j, 0, :, :] = np.append(reshaped[j, -2, :, 1:], np.zeros((13, 1)), axis=1)
            reshaped[j, -1, :, :] = np.append(np.zeros((13, 1)), reshaped[j, 1, :, :-1], axis=1)
        for i in range(0, 11):
            reshaped[0, i, :, :] = np.append(reshaped[-2, i, 1:, :], np.zeros((1, 13)), axis=0)
            reshaped[-1, i, :, :] = np.append(np.zeros((1, 13)), reshaped[1, i, :-1, :], axis=0)
        self.base.reshaped = reshaped

    def img(self, x, y, shape=(21, 21)):
        """Unit-flux PRF image of a star at cutout position (x, y). Stars up to 6 px off the stamp are handled."""
        return self.base.locate(float(x), float(y), shape)

    def imgs(self, xs, ys, shape=(21, 21)):
        return np.array([self.img(x, y, shape) for x, y in zip(xs, ys)])


class GaussPSF:
    """Pixel-integrated circular Gaussian (alternative PSF model for robustness tests)."""

    def __init__(self, sig_px=0.9):
        self.sig = sig_px

    def img(self, x, y, shape=(21, 21)):
        from scipy.special import erf
        ny, nx = shape
        xe = np.arange(nx + 1) - 0.5; ye = np.arange(ny + 1) - 0.5
        s = self.sig * np.sqrt(2)
        cx = np.diff(0.5 * erf((xe - x) / s)); cy = np.diff(0.5 * erf((ye - y) / s))
        return np.outer(cy, cx)

    def imgs(self, xs, ys, shape=(21, 21)):
        return np.array([self.img(x, y, shape) for x, y in zip(xs, ys)])


def template(t, exptime, P_=P, t0=T0):
    """Unit-depth transit template: batman shape (rp 0.0668, a/R* 4.30, b 0.85) normalised to depth 1 at mid-transit.
    Returns (1 - lc)/(1 - lc_min) so that flux = baseline * (1 - A * template)."""
    import batman
    p = batman.TransitParams()
    a, b, rp = 4.30, 0.85, 0.0668
    p.t0, p.per, p.rp, p.a, p.inc, p.ecc, p.w = t0, P_, rp, a, np.degrees(np.arccos(b / a)), 0.0, 90.0
    p.u, p.limb_dark = [0.45, 0.2], "quadratic"
    ss = max(1, int(round(exptime / 120)))
    lc = batman.TransitModel(p, np.asarray(t, float), supersample_factor=ss, exp_time=exptime / 86400).light_curve(p)
    tc = np.array([t0])
    lc0 = batman.TransitModel(p, tc).light_curve(p)[0]
    return (1 - lc) / (1 - lc0)


def phase_dt(t, P_=P, t0=T0):
    ep = np.round((t - t0) / P_)
    return t - t0 - ep * P_, ep
