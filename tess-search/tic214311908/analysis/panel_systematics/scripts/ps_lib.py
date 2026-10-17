"""Systematics-panel library (own code). Reads the TESScut cubes read-only.

Inputs shared with earlier steps (calibration only): the Gaia DR3 catalogue saved by referee_loc
(data/gaia_400.csv), the per-sector WCS offsets measured by referee_loc (out/wcs_mira.json, offset_px), the SPOC PRF
model files and the lateral-colour curve (radial shift vs BP-RP).
Times: TESScut TIME - TIMECORR + Roemer delay toward the M3V (BJD_TDB - 2457000).
Ephemeris: adopted T0 = 3840.9150, P = 5.691506 d.
"""
import glob, json, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS

TIC = 'tic214311908/analysis'
ROOT = f'{TIC}/panel_systematics'
CUT = '/home/user/work/tesscut5/214311906'
PRFDIR = '<local>/tic713/localise/prf_files'
RA0, DE0 = 260.977650, -51.323690
T0, P = 3840.9150, 5.691506
T14, FLAT = 1.70 / 24, 0.72
SECTORS = (12, 39, 66, 93, 103, 104)
ID = dict(M3V=5925373272850784768, pairA=5925373272850782208, pairB=5925373272827253632,
          Mira=5925373268528684288, G175=5925373272827260544, G177=5925373272827262208, G170=5925373272827256576)
FM3V = 979.0
CURVE_C = np.array([0.0, 0.79, 0.90, 1.09, 1.27, 1.54, 1.83, 2.34, 2.48, 2.94])
CURVE_R = np.array([1.20, 0.0, -0.15, -0.47, -0.74, -1.20, -1.79, -2.80, -3.42, -3.95])


def colour_shift(C):
    C = np.asarray(C, float)
    r = np.interp(C, CURVE_C, CURVE_R)
    return np.where(C > CURVE_C[-1], CURVE_R[-1] - 1.5 * (C - CURVE_C[-1]), r)


def roemer(t):
    from astropy.time import Time
    from astropy.coordinates import get_body_barycentric
    import astropy.units as u
    ra, de = np.radians(260.9697), np.radians(-51.3242)
    n = np.array([np.cos(de) * np.cos(ra), np.cos(de) * np.sin(ra), np.sin(de)])
    pos = get_body_barycentric('earth', Time(t + 2457000.0, format='jd', scale='tdb'))
    r = np.vstack([pos.x.to(u.m).value, pos.y.to(u.m).value, pos.z.to(u.m).value]).T
    return r @ n / 299792458.0 / 86400.0


def gaia():
    g = pd.read_csv(f'{TIC}/referee_loc/data/gaia_400.csv')
    g['C'] = g.bp_rp
    g.loc[g.source_id == ID['pairB'], 'C'] = 3.3
    g['C'] = g['C'].fillna(1.2)
    C = g.C.values
    g['T'] = g.phot_g_mean_mag - 0.00522555 * C ** 3 + 0.0891337 * C ** 2 - 0.633923 * C + 0.0324473
    g['F'] = 15000.0 * 10 ** (-0.4 * (g['T'] - 10.0))
    return g


def trap(t, tc, dt, t14=T14, flat=FLAT, nsub=11):
    o = (np.arange(nsub) + 0.5) / nsub - 0.5
    tc = np.asarray(tc, float)
    tt = np.abs(t[:, None] + o[None, :] * dt - (tc[:, None] if tc.ndim else tc))
    t23 = flat * t14
    v = np.clip((t14 / 2 - tt) / ((t14 - t23) / 2), 0, 1)
    return v.mean(1)


class Cube:
    def __init__(self, sec):
        f = glob.glob(f'{CUT}/tess-s{sec:04d}-*_astrocut.fits')[0]
        with fits.open(f) as h:
            self.cam, self.ccd = int(h[0].header['CAMERA']), int(h[0].header['CCD'])
            H1 = h[1].header
            self.col0, self.row0 = int(H1['1CRV4P']), int(H1['2CRV4P'])
            self.wcs = WCS(h[2].header)
            d = h[1].data
            t = np.array(d['TIME'], float); tc = np.array(d['TIMECORR'], float)
            self.q = np.array(d['QUALITY'], int)
            self.flux = np.array(d['FLUX'], np.float64)
            self.bkg = np.array(d['FLUX_BKG'], np.float64)
        ok = np.isfinite(t) & (t > 0)
        tt = np.full_like(t, np.nan)
        tt[ok] = t[ok] - tc[ok] + roemer(t[ok] - tc[ok])
        self.t = tt
        self.sec = sec
        self.good = ok & (self.q == 0) & np.isfinite(self.flux).all(axis=(1, 2)) & np.isfinite(self.bkg).all(axis=(1, 2))
        self.dt = float(np.nanmedian(np.diff(self.t[self.good])))
        self.epoch = 2016.0 + (np.nanmedian(self.t[self.good]) + 2457000 - 2457389.0) / 365.25
        self.ny, self.nx = self.flux.shape[1:]
        off = json.load(open(f'{TIC}/referee_loc/out/wcs_mira.json'))[str(sec)]['offset_px']
        self.off = np.array(off)
        # boresight direction for the lateral colour shift
        import tess_stars2px as t2p
        res = t2p.tess_stars2px_function_entry(1, RA0, DE0, trySector=sec)
        sc = res[-1]
        bra, bde = float(sc.camRa[self.cam - 1][0]), float(sc.camDec[self.cam - 1][0])
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        pa = np.radians(float(SkyCoord(RA0 * u.deg, DE0 * u.deg).position_angle(SkyCoord(bra * u.deg, bde * u.deg)).deg))
        self.u = np.array([np.sin(pa), np.cos(pa)])
        xc, yc = self.wcs.all_world2pix([[RA0, DE0]], 0)[0]

        def sky(x, y):
            r, dd = self.wcs.all_pix2world([[x, y]], 0)[0]
            return np.array([(r - RA0) * np.cos(np.radians(DE0)) * 3600, (dd - DE0) * 3600])
        b = sky(xc, yc)
        self.J = np.column_stack([(sky(xc + 0.01, yc) - b) / 0.01, (sky(xc, yc + 0.01) - b) / 0.01])
        self.Ji = np.linalg.inv(self.J)
        self._prf = None

    def xy(self, g, colour=True):
        dt = self.epoch - 2016.0
        ra = g.ra.values + np.nan_to_num(g.pmra.values) * dt / 3.6e6 / np.cos(np.radians(g.dec.values))
        de = g.dec.values + np.nan_to_num(g.pmdec.values) * dt / 3.6e6
        x, y = self.wcs.all_world2pix(ra, de, 0)
        p = np.column_stack([x, y]) + self.off
        if colour:
            s = colour_shift(g.C.values)
            p = p + (self.Ji @ (s[:, None] * self.u[None, :]).T).T
        return p

    def prf(self, x, y):
        if self._prf is None:
            import PRF
            self._prf = PRF.TESS_PRF(self.cam, self.ccd, self.sec, self.col0 + 10, self.row0 + 10,
                                     localdatadir=f'{PRFDIR}/s0004/')
        return self._prf.locate(x, y, (self.ny, self.nx))

    def momentum_dumps(self):
        # QUALITY bit 32 (reaction-wheel desaturation) and bit 4 (coarse point)
        m = np.isfinite(self.t) & ((self.q & 32) > 0)
        return self.t[m]


def jdump(o, fn):
    def cv(x):
        if isinstance(x, np.integer): return int(x)
        if isinstance(x, np.floating): return float(x)
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, np.bool_): return bool(x)
        raise TypeError(type(x))
    json.dump(o, open(fn, 'w'), indent=1, default=cv)
