"""Shared SED machinery: SVO filter curves + zero points, Castelli & Kurucz (2004) ATLAS9 surface fluxes,
Fitzpatrick (1999) R_V=3.1 extinction applied to the spectrum before band integration (photon-counting
band averages), and a 4-D interpolation grid (Teff, logg, [M/H], A_V) of band-averaged surface fluxes.
"""
import os, re, glob
import numpy as np
from astropy.io import fits
import extinction
from scipy.interpolate import RegularGridInterpolator

ROOT = "/home/user/work/tic87/stellar"  # read-only (filters, ck04)
GRID = "/home/user/work/tic87/neighbour/data/sed_grid.npz"
FDIR = os.path.join(ROOT, "cache/filters"); CKDIR = os.path.join(ROOT, "cache/ck04")
C_AA = 2.99792458e18  # Angstrom/s

FILTERS = {  # name: (svo id, system)
    "G": "GAIA/GAIA3.G", "BP": "GAIA/GAIA3.Gbp", "RP": "GAIA/GAIA3.Grp",
    "J": "2MASS/2MASS.J", "H": "2MASS/2MASS.H", "Ks": "2MASS/2MASS.Ks",
    "W1": "WISE/WISE.W1", "W2": "WISE/WISE.W2", "W3": "WISE/WISE.W3",
    "SM_u": "SkyMapper/SkyMapper.u", "SM_v": "SkyMapper/SkyMapper.v", "SM_g": "SkyMapper/SkyMapper.g",
    "SM_r": "SkyMapper/SkyMapper.r", "SM_i": "SkyMapper/SkyMapper.i", "SM_z": "SkyMapper/SkyMapper.z",
    "BT": "TYCHO/TYCHO.B", "VT": "TYCHO/TYCHO.V", "B": "Generic/Bessell.B", "V": "Generic/Bessell.V",
    "gp": "SLOAN/SDSS.g", "rp": "SLOAN/SDSS.r", "ip": "SLOAN/SDSS.i", "TESS": "TESS/TESS.Red",
    "U": "Generic/Bessell.U", "R": "Generic/Bessell.R", "I": "Generic/Bessell.I",
}


def load_filter(name):
    fid = FILTERS[name]
    txt = open(os.path.join(FDIR, fid.replace("/", "_") + ".xml")).read()
    zp = float(re.search(r'name="ZeroPoint"[^>]*value="([^"]+)"', txt).group(1))
    lp = float(re.search(r'name="WavelengthPivot"[^>]*value="([^"]+)"', txt).group(1))
    rows = re.findall(r"<TR>\s*<TD>([^<]+)</TD>\s*<TD>([^<]+)</TD>\s*</TR>", txt)
    w = np.array([float(a) for a, b in rows]); t = np.array([float(b) for a, b in rows])
    o = np.argsort(w)
    zp_flam = zp * 1e-23 * C_AA / lp ** 2   # erg/s/cm2/A for mag = 0
    return dict(w=w[o], t=np.clip(t[o], 0, None), zp_jy=zp, lpiv=lp, zp_flam=zp_flam)


FILT = {k: load_filter(k) for k in FILTERS}


def band_mean(wl, fl, f):
    """photon-counting band average <F_lambda> = int F T lambda dlam / int T lambda dlam, evaluated on the union of
    the model and filter wavelength grids (filter linear, model log-log interpolated: the CK04 grid is sparse in the
    mid-IR, e.g. 100200 -> 200000 A)."""
    lo, hi = f["w"][0], f["w"][-1]
    w = np.unique(np.concatenate([f["w"], wl[(wl > lo) & (wl < hi)]]))
    t = np.interp(w, f["w"], f["t"], left=0, right=0)
    pos = fl > 0
    F = 10 ** np.interp(np.log10(w), np.log10(wl[pos]), np.log10(fl[pos]))
    return np.trapezoid(F * t * w, w) / np.trapezoid(t * w, w)


def ck04_files():
    out = {}
    for fn in glob.glob(os.path.join(CKDIR, "ck*_*.fits")):
        b = os.path.basename(fn)[:-5]
        z, teff = b.split("_")
        mh = (-1 if z[2] == "m" else 1) * int(z[3:]) / 10.
        out[(int(teff), mh)] = fn
    return out


def load_spec(fn, logg):
    d = fits.open(fn)[1].data
    col = "g%02d" % int(round(logg * 10))
    return np.array(d["WAVELENGTH"], float), np.array(d[col], float)


TEFFS = np.arange(3500, 8001, 250)
LOGGS = np.array([3.0, 3.5, 4.0, 4.5, 5.0])
MHS = np.array([-1.5, -1.0, -0.5, 0.0, 0.2, 0.5])
AVS = np.round(np.arange(0, 1.51, 0.05), 3)


def build_grid(bands, path=GRID):
    if os.path.exists(path):
        z = np.load(path, allow_pickle=True)
        if list(z["bands"]) == list(bands):
            return z["grid"], z["fbol"]
    files = ck04_files()
    grid = np.full((len(TEFFS), len(LOGGS), len(MHS), len(AVS), len(bands)), np.nan)
    fbol = np.full((len(TEFFS), len(LOGGS), len(MHS)), np.nan)
    for i, T in enumerate(TEFFS):
        for k, mh in enumerate(MHS):
            fn = files.get((int(T), mh))
            if fn is None: continue
            for j, lg in enumerate(LOGGS):
                wl, fl = load_spec(fn, lg)
                if not np.any(fl > 0): continue
                fbol[i, j, k] = np.trapezoid(fl, wl)
                m = (wl > 900) & (wl < 300000)
                wl_, fl_ = wl[m], fl[m]
                ext1 = extinction.fitzpatrick99(wl_, 1.0, 3.1)   # A_lambda for A_V = 1
                for l, av in enumerate(AVS):
                    fr = fl_ * 10 ** (-0.4 * ext1 * av)
                    for b, bn in enumerate(bands):
                        grid[i, j, k, l, b] = band_mean(wl_, fr, FILT[bn])
    np.savez(path, grid=grid, fbol=fbol, bands=np.array(bands))
    return grid, fbol


# Offsets added to synthetic magnitudes so that CALSPEC Vega has the magnitude defined by each system
# (validate_synphot.py gives the synthetic Vega mags with SVO zero points):
#   2MASS (Cohen et al. 2003): J=-0.001, H=+0.019, Ks=-0.017 ; synthetic 0.020, 0.017, 0.023
#   WISE (Wright et al. 2010): Vega = 0 in W1, W2 ; synthetic 0.049, 0.042 (W3 0.207 -> W3 not used in fits)
#   Johnson/Bessell B, V (Vega V = B = 0.03): synthetic -0.026, 0.015
VEGA_OFFSET = {"J": -0.021, "H": 0.002, "Ks": -0.040, "W1": -0.049, "W2": -0.042, "B": 0.056, "V": 0.015}


class SEDModel:
    def __init__(self, bands):
        self.bands = list(bands)
        grid, fbol = build_grid(self.bands)
        self.ok = np.isfinite(grid).all(axis=(3, 4))
        lg = np.log10(np.where(grid > 0, grid, np.nan))
        # fill missing models (e.g. logg/Teff corners) by nearest valid logg to keep interpolator defined
        for i in range(lg.shape[0]):
            for k in range(lg.shape[2]):
                for j in range(lg.shape[1]):
                    if not np.isfinite(lg[i, j, k]).all():
                        good = [jj for jj in range(lg.shape[1]) if np.isfinite(lg[i, jj, k]).all()]
                        if good:
                            jj = min(good, key=lambda x: abs(x - j)); lg[i, j, k] = lg[i, jj, k]
        self.interp = RegularGridInterpolator((TEFFS, LOGGS, MHS, AVS), lg, bounds_error=False, fill_value=np.nan)
        self.fbol_interp = RegularGridInterpolator((TEFFS, LOGGS, MHS), np.log10(fbol), bounds_error=False, fill_value=np.nan)
        self.zp = np.array([FILT[b]["zp_flam"] for b in self.bands])
        self.off = np.array([VEGA_OFFSET.get(b, 0.0) for b in self.bands])

    def mags(self, teff, logg, mh, av, theta_mas):
        """apparent magnitudes for angular diameter theta (mas)"""
        th = theta_mas / 206264806.247  # radians
        lf = self.interp([[teff, logg, mh, av]])[0]
        return -2.5 * (lf + 2 * np.log10(th / 2)) + 2.5 * np.log10(self.zp) + self.off

    def fbol(self, teff, logg, mh, theta_mas):
        th = theta_mas / 206264806.247
        return 10 ** self.fbol_interp([[teff, logg, mh]])[0] * (th / 2) ** 2
