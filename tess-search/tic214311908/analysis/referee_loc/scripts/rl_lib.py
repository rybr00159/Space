"""Referee library (independent of the 'loc' agent's l0lib/l1img/l2fit).

Conventions: cutout arrays [row, col]; positions (x=col, y=row), 0-based, integer = pixel centre (PRF package
convention).  Sky offsets (dE, dN) arcsec.  Times: BJD_TDB - 2457000 for the TARGET direction (TESScut TIME minus its
TIMECORR plus my own Roemer delay toward the M3V; geocentre approximation, <~1.3 s).
"""
import os, glob, json, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS

ROOT = 'tic214311908/analysis/referee_loc'
CUT = '/home/user/work/tesscut5/214311906'
PRFDIR = '<local>/tic713/localise/prf_files'
RA0, DE0 = 260.977650, -51.323690          # cutout centre (Mira)
P, T0 = 5.6914968, 3840.9153               # BJD_TDB - 2457000
T14 = 1.63 / 24
FLAT = 0.70
SECTORS = (12, 39, 66, 93, 103, 104)
ID = dict(M3V=5925373272850784768, pairA=5925373272850782208, pairB=5925373272827253632,
          Mira=5925373268528684288, ref=5925373204131295872, G175=5925373272827260544,
          G177=5925373272827262208, G170=5925373272827256576)
NAME = {v: k for k, v in ID.items()}


def gaia():
    g = pd.read_csv(f'{ROOT}/data/gaia_400.csv')
    c = np.cos(np.radians(DE0))
    g['dE'] = (g.ra - RA0) * c * 3600
    g['dN'] = (g.dec - DE0) * 3600
    g['C'] = g.bp_rp
    g.loc[g.source_id == ID['pairB'], 'C'] = 3.3        # no Gaia colour; assume ~M4 (variants test this)
    g['C'] = g['C'].fillna(1.2)
    # TESS magnitude: Stassun et al. 2019 (TIC-8) Gaia relation
    C = g.C.values
    g['T'] = g.phot_g_mean_mag - 0.00522555 * C ** 3 + 0.0891337 * C ** 2 - 0.633923 * C + 0.0324473
    g['F'] = 15000.0 * 10 ** (-0.4 * (g['T'] - 10.0))
    g['name'] = [NAME.get(s, str(s)) for s in g.source_id]
    return g


def bary_target(tsc):
    """Roemer delay (s) toward the M3V for spacecraft TDB times tsc (BTJD); geocentric approximation."""
    from astropy.time import Time
    from astropy.coordinates import get_body_barycentric
    import astropy.units as u
    ra, de = np.radians(260.97, ), np.radians(-51.3238)
    n = np.array([np.cos(de) * np.cos(ra), np.cos(de) * np.sin(ra), np.sin(de)])
    pos = get_body_barycentric('earth', Time(tsc + 2457000.0, format='jd', scale='tdb'))
    r = np.vstack([pos.x.to(u.m).value, pos.y.to(u.m).value, pos.z.to(u.m).value]).T
    return r @ n / 299792458.0


# lateral-colour models (radial shift in arcsec, NEGATIVE = away from boresight), relative to C=0.786
CURVE_C = np.array([0.0, 0.79, 0.90, 1.09, 1.27, 1.54, 1.83, 2.34, 2.48, 2.94])
CURVE_R = np.array([1.20, 0.0, -0.15, -0.47, -0.74, -1.20, -1.79, -2.80, -3.42, -3.95])


def colour_shift(C, model='curve'):
    C = np.asarray(C, float)
    if model == 'none':
        return np.zeros_like(C)
    if model == 'linear':                 # simple 1.5"/mag everywhere (task text)
        return -1.5 * (C - 0.786)
    r = np.interp(C, CURVE_C, CURVE_R)    # campaign curve, extrapolated at 1.5"/mag beyond 2.94
    return np.where(C > CURVE_C[-1], CURVE_R[-1] - 1.5 * (C - CURVE_C[-1]), r)


class Sec:
    def __init__(self, sec, cube=True):
        f = glob.glob(f'{CUT}/tess-s{sec:04d}-*_astrocut.fits')[0]
        with fits.open(f) as h:
            H0 = h[0].header; H1 = h[1].header
            self.wcs = WCS(h[2].header)
            d = h[1].data
            t = np.array(d['TIME'], float); tc = np.array(d['TIMECORR'], float)
            self.q = np.array(d['QUALITY'], int)
            if cube:
                self.flux = np.array(d['FLUX'], np.float64)
                self.ferr = np.array(d['FLUX_ERR'], np.float64)
                self.bkg = np.array(d['FLUX_BKG'], np.float32)
        self.sec, self.cam, self.ccd = sec, int(H0['CAMERA']), int(H0['CCD'])
        self.col0, self.row0 = int(H1['1CRV4P']), int(H1['2CRV4P'])
        ok = np.isfinite(t) & (t > 0)
        self.t_raw = t.copy()
        tt = np.full_like(t, np.nan)
        tt[ok] = t[ok] - tc[ok] + bary_target(t[ok] - tc[ok]) / 86400.0
        self.t = tt
        self.dcorr_s = np.full_like(t, np.nan); self.dcorr_s[ok] = (tt[ok] - t[ok]) * 86400
        self.good = ok & (self.q == 0)
        if cube:
            self.good &= np.isfinite(self.flux).all(axis=(1, 2))
        self.ny, self.nx = (self.flux.shape[1:] if cube else (21, 21))
        self.dt = float(np.nanmedian(np.diff(self.t[self.good])))
        self.epoch = 2016.0 + (np.nanmedian(self.t[self.good]) + 2457000 - 2457389.0) / 365.25
        # local Jacobian pixel->(dE,dN) at the cutout centre
        xc, yc = self.wcs.all_world2pix([[RA0, DE0]], 0)[0]
        def sky(x, y):
            r, d = self.wcs.all_pix2world([[x, y]], 0)[0]
            return np.array([(r - RA0) * np.cos(np.radians(DE0)) * 3600, (d - DE0) * 3600])
        b = sky(xc, yc)
        self.J = np.column_stack([(sky(xc + 0.01, yc) - b) / 0.01, (sky(xc, yc + 0.01) - b) / 0.01])
        self.Ji = np.linalg.inv(self.J)
        self.xc, self.yc = xc, yc
        # boresight via tess-point
        import tess_stars2px as t2p
        res = t2p.tess_stars2px_function_entry(1, RA0, DE0, trySector=sec)
        sc = res[-1]
        bra, bde = float(sc.camRa[self.cam - 1][0]), float(sc.camDec[self.cam - 1][0])
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        c0 = SkyCoord(RA0 * u.deg, DE0 * u.deg); cb = SkyCoord(bra * u.deg, bde * u.deg)
        pa = np.radians(float(c0.position_angle(cb).deg))
        self.u = np.array([np.sin(pa), np.cos(pa)])     # (E,N) unit vector toward the boresight
        self.theta = float(c0.separation(cb).deg)
        self.off = np.zeros(2)                          # WCS offset (px) - set after calibration
        self._prf = None

    def xy(self, g, cmodel='curve', cscale=1.0, use_colour=True):
        """Apparent pixel positions of catalogue rows."""
        dt = self.epoch - 2016.0
        ra = g.ra.values + np.nan_to_num(g.pmra.values) * dt / 3.6e6 / np.cos(np.radians(g.dec.values))
        de = g.dec.values + np.nan_to_num(g.pmdec.values) * dt / 3.6e6
        x, y = self.wcs.all_world2pix(ra, de, 0)
        p = np.column_stack([x, y]) + self.off
        if use_colour:
            s = cscale * colour_shift(g.C.values, cmodel)
            p = p + (self.Ji @ (s[:, None] * self.u[None, :]).T).T
        return p

    def prf(self):
        if self._prf is None:
            import PRF
            p = PRF.TESS_PRF(self.cam, self.ccd, self.sec, self.col0 + 10, self.row0 + 10,
                             localdatadir=f'{PRFDIR}/s0004/')
            shp = (self.ny, self.nx)
            def f(x, y):
                if x < -6 or y < -6 or x > shp[1] + 5 or y > shp[0] + 5:
                    return np.zeros(shp)
                return p.locate(x, y, shp)
            self._prf = f
        return self._prf


def trapezoid(t, tc, dt, t14=T14, flat=FLAT, nsub=15):
    """Exposure-averaged unit-depth trapezoid."""
    o = (np.arange(nsub) + 0.5) / nsub - 0.5
    tt = np.abs(t[:, None] + o[None, :] * dt - tc)
    t23 = flat * t14
    v = np.clip((t14 / 2 - tt) / ((t14 - t23) / 2), 0, 1)
    return v.mean(1)


def jdump(o, fn):
    def cv(x):
        if isinstance(x, np.integer): return int(x)
        if isinstance(x, np.floating): return float(x)
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, np.bool_): return bool(x)
        raise TypeError(type(x))
    json.dump(o, open(fn, 'w'), indent=1, default=cv)
