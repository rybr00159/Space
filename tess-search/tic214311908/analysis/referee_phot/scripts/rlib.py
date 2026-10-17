"""Referee library (independent of the phot investigator's photlib): cutout loading, Gaia propagation, SPOC PRF via
PRF.TESS_PRF.locate (package routine, not a re-implementation), simple helpers."""
import os, glob, json, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS

R = "tic214311908/analysis/referee_phot"
CUT = "/home/user/work/tesscut5/214311906"
PRFDIR = "<local>/tic713/localise/prf_files"
SECTORS = (12, 39, 66, 93, 103, 104)
GID = dict(Mira=5925373268528684288, M3V=5925373272850784768, pA=5925373272850782208, pB=5925373272827253632,
           bg3=5925373272827260544, ref=5925373204131295872)
P0, T00 = 5.6914968, 3840.9153          # task ephemeris (QLP)


def tflux(T):
    return 15000.0 * 10 ** (-0.4 * (np.asarray(T, float) - 10.0))


def rstd(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def gaia_T(g):
    """Stassun+19 Gaia->TESS relation (colour clipped to [-0.2, 3.0]); G-0.43 without colour."""
    c = np.clip(g.bp_rp.values, -0.2, 3.0)
    T = g.phot_g_mean_mag - 0.00522555 * c ** 3 + 0.0891337 * c ** 2 - 0.633923 * c + 0.0324473
    return np.where(np.isfinite(g.bp_rp), T, g.phot_g_mean_mag - 0.43)


class Cut:
    def __init__(self, s):
        f = glob.glob(f"{CUT}/tess-s{s:04d}-*_astrocut.fits")[0]
        with fits.open(f) as h:
            self.h0 = h[0].header; h1 = h[1].header; self.wcs = WCS(h[2].header)
            d = h[1].data
            self.t = np.array(d["TIME"], float); self.tcorr = np.array(d["TIMECORR"], float)
            self.f = np.array(d["FLUX"], float); self.q = np.array(d["QUALITY"], int)
            self.bkg_spoc = np.array(d["FLUX_BKG"], float) if "FLUX_BKG" in d.names else None
        self.s = s; self.cam = int(self.h0["CAMERA"]); self.ccd = int(self.h0["CCD"])
        self.col0 = int(h1["1CRV4P"]); self.row0 = int(h1["2CRV4P"])
        self.n, self.ny, self.nx = self.f.shape
        self.fin = np.isfinite(self.t) & np.all(np.isfinite(self.f.reshape(self.n, -1)), 1)
        self.dt = float(np.nanmedian(np.diff(self.t[self.fin])))
        self.epoch = 2000.0 + (np.nanmedian(self.t[self.fin]) + 2457000 - 2451545.0) / 365.25
        self._prf = None

    def stars(self, Tmax=19.0):
        g = pd.read_csv(f"{R}/data/gaia330.csv")
        g["T"] = gaia_T(g); g = g[g["T"] < Tmax].copy()
        dt = self.epoch - 2016.0
        ra = g.ra + np.nan_to_num(g.pmra) * dt / 3.6e6 / np.cos(np.radians(g.dec))
        de = g.dec + np.nan_to_num(g.pmdec) * dt / 3.6e6
        x, y = self.wcs.all_world2pix(ra.values, de.values, 0)
        g["x"], g["y"] = x, y
        return g.reset_index(drop=True)

    def prf(self):
        if self._prf is None:
            import PRF
            sub = "s0001" if self.s < 4 else "s0004"
            self._prf = PRF.TESS_PRF(self.cam, self.ccd, self.s, self.col0 + 10, self.row0 + 10, localdatadir=f"{PRFDIR}/{sub}/")
        return self._prf

    def psf(self, x, y):
        """Unit-sum (over 13x13) SPOC PRF placed at 0-based (x=col, y=row) in the cutout frame."""
        return self.prf().locate(x, y, stampsize=(self.ny, self.nx))

    def en(self, x, y, dx, dy):
        """(dE, dN) arcsec of a pixel offset (dx, dy) at (x, y)."""
        r0 = np.array(self.wcs.all_pix2world([[x, y]], 0))[0]; r1 = np.array(self.wcs.all_pix2world([[x + dx, y + dy]], 0))[0]
        return np.array([(r1[0] - r0[0]) * np.cos(np.radians(r0[1])) * 3600, (r1[1] - r0[1]) * 3600])


def save_json(o, fn):
    def conv(x):
        if isinstance(x, (np.integer,)): return int(x)
        if isinstance(x, (np.floating,)): return float(x)
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, np.bool_): return bool(x)
        raise TypeError(type(x))
    json.dump(o, open(fn, "w"), indent=1, default=conv)
