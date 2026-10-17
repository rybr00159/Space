"""Shared helpers for the deblended PSF photometry of the TIC 214311906/08/10 field.

Conventions: pixel coordinates 0-based, images indexed [row, col] = [y, x]; integer coordinates are pixel centres.
Star positions come from Gaia DR3 propagated to the sector epoch and the TESScut WCS, plus a per-sector WCS offset
calibrated on normal-colour field stars (never on the Mira).
"""
import os, glob, json, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS

ROOT = "tic214311908/analysis/phot"
CUTDIR = "/home/user/work/tesscut5/214311906"
PRFDIR = "<local>/tic713/localise/prf_files"
LCDIR = "/home/user/work/lc_cache"
GAIA = f"{ROOT}/data/gaia/gaia_360arcsec.csv"
SECTORS = (12, 39, 66, 93, 103, 104)
IDS = dict(Mira=5925373268528684288, M3V=5925373272850784768, pairA=5925373272850782208, pairB=5925373272827253632,
           ref=5925373204131295872)
TICS = dict(Mira=214311906, M3V=214311908, pair=214311910)
P_EPH, T0_EPH = 5.6914968266065875, 3840.915307691422      # ephemeris_v2.json (QLP, corrected S103/S104 times)
T14_EPH = 0.067
RA_M3V, DEC_M3V = 260.970485, -51.324202                     # TIC 214311908 (J2000, epoch 2000 position irrelevant at ms level)


def tmag_flux(T):
    """TESS e-/s for a TESS magnitude (zero point 15000 e/s at T=10, Vanderspek et al.)."""
    return 15000.0 * 10 ** (-0.4 * (np.asarray(T) - 10.0))


class Cutout:
    def __init__(self, sector):
        f = glob.glob(f"{CUTDIR}/tess-s{sector:04d}-*_astrocut.fits")[0]
        with fits.open(f) as h:
            self.hdr0 = h[0].header; h1 = h[1].header; self.wcs = WCS(h[2].header)
            d = h[1].data
            self.t = np.array(d["TIME"], float); self.timecorr = np.array(d["TIMECORR"], float)
            self.flux = np.array(d["FLUX"], float); self.ferr = np.array(d["FLUX_ERR"], float)
            self.q = np.array(d["QUALITY"], int); self.cadn = np.array(d["CADENCENO"], int)
        self.file = f; self.sector = sector
        self.camera = int(self.hdr0["CAMERA"]); self.ccd = int(self.hdr0["CCD"])
        self.col0 = int(h1["1CRV4P"]); self.row0 = int(h1["2CRV4P"])
        self.ny, self.nx = self.flux.shape[1:]
        self.fin = np.all(np.isfinite(self.flux.reshape(len(self.t), -1)), axis=1) & np.isfinite(self.t)
        self.exptime = float(np.nanmedian(np.diff(self.t[self.fin])))
        self.epoch = 2000.0 + (np.nanmedian(self.t[self.fin]) + 2457000.0 - 2451545.0) / 365.25

    def gaia(self, Tmax=19.0):
        g = pd.read_csv(GAIA)
        g = g[g["T"] < Tmax].copy()
        dt = self.epoch - 2016.0
        g["ra_ep"] = g.ra + np.nan_to_num(g.pmra) * dt / 3.6e6 / np.cos(np.radians(g.dec))
        g["dec_ep"] = g.dec + np.nan_to_num(g.pmdec) * dt / 3.6e6
        x, y = self.wcs.all_world2pix(g.ra_ep.values, g.dec_ep.values, 0)
        g["x"], g["y"] = x, y
        return g.reset_index(drop=True)

    def jacobian(self, x, y):
        """d(dE, dN)/d(x, y) in arcsec per pixel at (x, y)."""
        r0 = np.array(self.wcs.all_pix2world([[x, y]], 0))[0]
        def off(r): return np.array([(r[0] - r0[0]) * np.cos(np.radians(r0[1])) * 3600, (r[1] - r0[1]) * 3600])
        rx = np.array(self.wcs.all_pix2world([[x + 0.01, y]], 0))[0]
        ry = np.array(self.wcs.all_pix2world([[x, y + 0.01]], 0))[0]
        return np.column_stack([off(rx) / 0.01, off(ry) / 0.01])


class PRFModel:
    """TESS SPOC PRF (via the PRF package's model files) at one detector location, with fast evaluation of
    a star image at sub-pixel position and optional Gaussian blur (extra jitter / focus) in pixels."""
    def __init__(self, C, xc=10.0, yc=10.0):
        import PRF
        sub = "s0001" if C.sector < 4 else "s0004"
        p = PRF.TESS_PRF(C.camera, C.ccd, C.sector, C.col0 + xc, C.row0 + yc, localdatadir=f"{PRFDIR}/{sub}/")
        self.reshaped = p.reshaped          # [subrow, subcol, 13, 13]
        self.shape = (C.ny, C.nx)
        self.samples = np.arange(-1 / 18, 19.1 / 18, 1 / 9)
        self.blur = 0.0

    def stamp(self, x, y):
        colint = np.floor(x + 0.5); cf = (x + 0.5) % 1
        rowint = np.floor(y + 0.5); rf = (y + 0.5) % 1
        s = self.samples
        cb = np.max(np.where(s < cf)[0]); ca = cb + 1
        rb = np.max(np.where(s < rf)[0]); ra = rb + 1
        wx = (cf - s[cb]) / (s[ca] - s[cb]); wy = (rf - s[rb]) / (s[ra] - s[rb])
        R = self.reshaped
        st = (R[rb, cb] * (1 - wx) * (1 - wy) + R[rb, ca] * wx * (1 - wy) + R[ra, cb] * (1 - wx) * wy + R[ra, ca] * wx * wy)
        st = st / st.sum()
        return st, int(colint), int(rowint)

    def image(self, x, y, shape=None, blur=None):
        shape = shape or self.shape
        blur = self.blur if blur is None else blur
        if blur > 0:   # average of shifted PRFs over a Gaussian kernel (Gauss-Hermite style 5x5 grid)
            n = np.array([-2, -1, 0, 1, 2]) * blur
            w1 = np.exp(-0.5 * (n / blur) ** 2); w1 /= w1.sum()
            out = np.zeros(shape)
            for i, dx in enumerate(n):
                for j, dy in enumerate(n):
                    out += w1[i] * w1[j] * self.image(x + dx, y + dy, shape, blur=0.0)
            return out
        out = np.zeros(shape)
        if x < -7 or y < -7 or x > shape[1] + 6 or y > shape[0] + 6: return out
        st, ci, ri = self.stamp(x, y)
        m = 6
        r0, r1 = max(0, ri - m), min(shape[0], ri + m + 1)
        c0, c1 = max(0, ci - m), min(shape[1], ci + m + 1)
        if r1 <= r0 or c1 <= c0: return out
        out[r0:r1, c0:c1] = st[r0 - (ri - m):r1 - (ri - m), c0 - (ci - m):c1 - (ci - m)]
        return out

    def grads(self, x, y, h=0.03, blur=None):
        gx = (self.image(x + h, y, blur=blur) - self.image(x - h, y, blur=blur)) / (2 * h)
        gy = (self.image(x, y + h, blur=blur) - self.image(x, y - h, blur=blur)) / (2 * h)
        return gx, gy


def roemer_target(sector, t_btjd, timecorr, ra=RA_M3V, dec=DEC_M3V, cache=f"{ROOT}/data/horizons"):
    """Seconds to ADD to TESScut TIME to barycentre to the target direction (TESScut TIME is barycentred for the FFI
    reference direction). Uses JPL Horizons barycentric positions of TESS (id -95). Cached per sector."""
    os.makedirs(cache, exist_ok=True)
    fn = f"{cache}/horizons_s{sector}.npz"
    tsc = t_btjd - timecorr
    if os.path.exists(fn):
        z = np.load(fn); jd, xyz = z["jd"], z["xyz"]
    else:
        from astroquery.jplhorizons import Horizons
        t0, t1 = np.nanmin(tsc), np.nanmax(tsc)
        o = Horizons(id="-95", location="@0", epochs={"start": f"JD{t0 + 2457000 - 0.5:.4f}", "stop": f"JD{t1 + 2457000 + 0.5:.4f}", "step": "30m"})
        v = o.vectors(refplane="earth")
        jd = np.array(v["datetime_jd"], float); xyz = np.vstack([np.asarray(v[c], float) for c in ("x", "y", "z")]).T
        np.savez(fn, jd=jd, xyz=xyz)
    RA, DE = np.radians(ra), np.radians(dec)
    N = np.array([np.cos(DE) * np.cos(RA), np.cos(DE) * np.sin(RA), np.sin(DE)])
    x = jd - 2457000
    r = np.vstack([np.interp(tsc, x, xyz[:, i]) for i in range(3)]).T
    AU_S = 149597870.700e3 / 299792458.0
    return r @ N * AU_S - timecorr * 86400.0


def robust_std(x):
    x = np.asarray(x); x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def trap(t, tc, depth, T14, T23):
    a = np.abs(t - tc)
    y = np.where(a <= T23 / 2, 1.0, np.where(a >= T14 / 2, 0.0, (T14 / 2 - a) / max((T14 - T23) / 2, 1e-9)))
    return 1 - depth * y


def save_json(obj, fn):
    def conv(o):
        if isinstance(o, np.integer): return int(o)
        if isinstance(o, np.floating): return float(o)
        if isinstance(o, np.ndarray): return o.tolist()
        if isinstance(o, np.bool_): return bool(o)
        raise TypeError(type(o))
    with open(fn, "w") as f: json.dump(obj, f, indent=1, default=conv)
