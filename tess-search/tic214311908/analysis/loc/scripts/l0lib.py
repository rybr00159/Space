"""Shared library for the independent pixel-level attribution of the 5.69-d transit near TIC 214311906/08/10.

Conventions: cutout pixel coords 0-based, image arrays [row(y), col(x)]; sky offsets (dE, dN) in arcsec relative to
the Mira's Gaia position at the sector epoch.  Flux in e-/s.  Times BTJD (TDB), with the TESScut FFI-direction
barycentric correction replaced by the Roemer delay towards the target (Earth position; error < 1.3 s).
"""
import os, glob, json, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS

BASE = 'tic214311908/analysis/loc'
CUT = '/home/user/work/tesscut5/214311906'
PRFDIR = '<local>/tic713/localise/prf_files'
RA0, DE0 = 260.977650, -51.323690
# ephemeris (13 QLP transits, S103/104 corrected): eph214v2 / ephemeris_v2.json
P, T0 = 5.6914968266065875, 3840.915307691422
T14 = 1.63 / 24.0          # total duration (d)
FLAT = 0.70                # T23/T14
SECTORS = (12, 39, 66, 93, 103, 104)
IDS = dict(M3V=5925373272850784768, pairA=5925373272850782208, pairB=5925373272827253632, Mira=5925373268528684288,
           ref=5925373204131295872)
C0 = 0.786
# lateral-colour curve (radial shift in arcsec, negative = away from boresight, relative to colour C0), from the
# TIC 71316387 campaign field-star fits (v4/f03 whole-CCD curve at field angle ~10.5 deg; theta^0.35 scaling)
COLCURVE_C = np.array([0.0, 0.79, 0.90, 1.09, 1.27, 1.54, 1.83, 2.34, 2.48, 2.94])
COLCURVE_R = np.array([1.20, 0.0, -0.15, -0.47, -0.74, -1.20, -1.79, -2.80, -3.42, -3.95])
THETA_REF, ALPHA_TH = 10.5, 0.35


def colour_shift(C, scale=1.0):
    """Radial colour shift (arcsec; negative = away from boresight) at the reference field angle."""
    C = np.asarray(C, float)
    r = np.interp(C, COLCURVE_C, COLCURVE_R)
    hi = C > COLCURVE_C[-1]
    r = np.where(hi, COLCURVE_R[-1] - 1.5 * (C - COLCURVE_C[-1]), r)
    return scale * r


def tmag(G, C):
    """Stassun+2019 Gaia->TESS relation (dwarf/giant common polynomial), valid to BP-RP ~ 5 roughly."""
    C = np.asarray(C, float)
    return G - 0.00522555 * C ** 3 + 0.0891337 * C ** 2 - 0.633923 * C + 0.0324473


def flux_from_T(T):
    return 15000.0 * 10 ** (-0.4 * (np.asarray(T, float) - 10.0))


def stars(maxsep=330.0):
    g = pd.read_csv(f'{BASE}/data/gaia_330.csv')
    g = g[g.sep < maxsep].copy()
    g['C'] = g.bp_rp
    g['C_assumed'] = ~np.isfinite(g.C)
    g.loc[g.source_id == IDS['pairB'], 'C'] = 3.3      # assumed: fainter M4.5 companion of pair A (no Gaia colour)
    g['C'] = g.C.fillna(1.2)
    g['G'] = g.phot_g_mean_mag
    g['T'] = tmag(g.G, g.C)
    g['F'] = flux_from_T(g['T'])
    return g.reset_index(drop=True)


def roemer_add(t, timecorr):
    """Seconds to ADD to TESScut TIME to obtain the target-specific BJD_TDB."""
    from astropy.time import Time
    from astropy.coordinates import get_body_barycentric
    import astropy.units as u
    ra, de = np.radians(RA0), np.radians(DE0)
    n = np.array([np.cos(de) * np.cos(ra), np.cos(de) * np.sin(ra), np.sin(de)])
    tsc = t - timecorr
    pos = get_body_barycentric('earth', Time(tsc + 2457000, format='jd', scale='tdb'))
    r = np.vstack([pos.x.to(u.au).value, pos.y.to(u.au).value, pos.z.to(u.au).value]).T
    return r @ n * 499.004784 - timecorr * 86400.0


class Sector:
    def __init__(self, sec, load_cube=True):
        f = glob.glob(f'{CUT}/tess-s{sec:04d}-*_astrocut.fits')[0]
        self.file = f
        with fits.open(f) as h:
            h0 = h[0].header; h1 = h[1].header
            self.wcs = WCS(h[2].header)
            d = h[1].data
            t = np.array(d['TIME'], float); tc = np.array(d['TIMECORR'], float); q = np.array(d['QUALITY'], int)
            if load_cube:
                self.flux = np.array(d['FLUX'], float)
                self.ferr = np.array(d['FLUX_ERR'], float)
        self.sec = sec; self.cam = int(h0['CAMERA']); self.ccd = int(h0['CCD'])
        self.col0 = int(h1['1CRV4P']); self.row0 = int(h1['2CRV4P'])
        ok = np.isfinite(t) & (t > 0)
        add = np.zeros_like(t); add[ok] = roemer_add(t[ok], tc[ok])
        self.t = t + add / 86400.0
        self.q = q
        self.good = ok & (q == 0)
        if load_cube:
            self.good &= np.all(np.isfinite(self.flux.reshape(len(t), -1)), axis=1)
            self.ny, self.nx = self.flux.shape[1:]
        else:
            self.ny = self.nx = 21
        self.exptime = float(np.nanmedian(np.diff(self.t[self.good])))
        self.epoch = 2000.0 + (np.nanmedian(self.t[self.good]) + 2457000.0 - 2451545.0) / 365.25
        # Jacobian pixel -> (dE, dN) arcsec at the Mira
        x0, y0 = self.wcs.all_world2pix([[RA0, DE0]], 0)[0]
        self.xm0, self.ym0 = float(x0), float(y0)
        def off(px):
            r = self.wcs.all_pix2world([px], 0)[0]
            return np.array([(r[0] - RA0) * np.cos(np.radians(DE0)) * 3600, (r[1] - DE0) * 3600])
        b = off([x0, y0])
        self.J = np.column_stack([(off([x0 + 0.01, y0]) - b) / 0.01, (off([x0, y0 + 0.01]) - b) / 0.01])
        self.Ji = np.linalg.inv(self.J)
        # boresight direction and field angle
        import tess_stars2px as t2p
        res = t2p.tess_stars2px_function_entry(1, RA0, DE0, trySector=sec)
        sc = res[-1]
        bra, bde = float(sc.camRa[self.cam - 1][0]), float(sc.camDec[self.cam - 1][0])
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        c0 = SkyCoord(RA0 * u.deg, DE0 * u.deg); cb = SkyCoord(bra * u.deg, bde * u.deg)
        self.pa_bs = float(c0.position_angle(cb).deg)      # PA (E of N) from the target towards the boresight
        self.theta = float(c0.separation(cb).deg)
        pa = np.radians(self.pa_bs)
        self.u_rad = np.array([np.sin(pa), np.cos(pa)])     # (E, N) unit vector towards the boresight
        self.colscale = (self.theta / THETA_REF) ** ALPHA_TH
        self.wcs_off = np.zeros(2)                          # set by calibration (pixels)
        cal = f'{BASE}/out/wcs_cal.json'
        if os.path.exists(cal):
            cc = json.load(open(cal)).get(str(sec))
            if cc: self.wcs_off = np.array(cc['offset_px'])
        self._prf = None

    # ---------------- positions
    def star_xy(self, g, colour=True, colscale=1.0, extra_shift=None):
        """Apparent pixel positions of Gaia stars: Gaia (propagated) -> WCS -> + calibrated WCS offset
        + lateral-colour shift for each star's colour."""
        dt = self.epoch - 2016.0
        ra = g.ra.values + np.nan_to_num(g.pmra.values) * dt / 3.6e6 / np.cos(np.radians(g.dec.values))
        de = g.dec.values + np.nan_to_num(g.pmdec.values) * dt / 3.6e6
        x, y = self.wcs.all_world2pix(ra, de, 0)
        xy = np.column_stack([x, y]) + self.wcs_off[None, :]
        if colour:
            sh = colour_shift(g.C.values, colscale * self.colscale)          # arcsec, radial (neg = away)
            dEN = sh[:, None] * self.u_rad[None, :]
            xy = xy + (self.Ji @ dEN.T).T
        if extra_shift is not None:
            xy = xy + extra_shift
        return xy

    def sky_to_xy(self, dEN, C=None, colscale=1.0):
        """Apparent pixel position of a source at sky offset dEN (arcsec, relative to the Mira's Gaia position at
        this epoch) with colour C (None -> no colour shift)."""
        dt = self.epoch - 2016.0
        g = pd.read_csv(f'{BASE}/data/gaia_330.csv', nrows=0) if False else None
        mira = MIRA_ROW
        ra = mira['ra'] + mira['pmra'] * dt / 3.6e6 / np.cos(np.radians(mira['dec']))
        de = mira['dec'] + mira['pmdec'] * dt / 3.6e6
        x, y = self.wcs.all_world2pix([[ra, de]], 0)[0]
        dEN = np.atleast_2d(dEN)
        xy = np.array([x, y])[None, :] + (self.Ji @ dEN.T).T + self.wcs_off[None, :]
        if C is not None:
            sh = colour_shift(C, colscale * self.colscale)
            xy = xy + (self.Ji @ (np.atleast_1d(sh)[:, None] * self.u_rad[None, :]).T).T
        return xy

    def xy_to_sky(self, xy, C=None, colscale=1.0):
        dt = self.epoch - 2016.0
        mira = MIRA_ROW
        ra = mira['ra'] + mira['pmra'] * dt / 3.6e6 / np.cos(np.radians(mira['dec']))
        de = mira['dec'] + mira['pmdec'] * dt / 3.6e6
        x, y = self.wcs.all_world2pix([[ra, de]], 0)[0]
        xy = np.atleast_2d(xy) - self.wcs_off[None, :] - np.array([x, y])[None, :]
        if C is not None:
            sh = colour_shift(C, colscale * self.colscale)
            xy = xy - (self.Ji @ (np.atleast_1d(sh)[:, None] * self.u_rad[None, :]).T).T
        return (self.J @ xy.T).T

    # ---------------- PRF
    def prf(self):
        if self._prf is None:
            import PRF
            sub = 's0001' if self.sec < 4 else 's0004'
            p = PRF.TESS_PRF(self.cam, self.ccd, self.sec, self.col0 + self.xm0, self.row0 + self.ym0,
                             localdatadir=f'{PRFDIR}/{sub}/')
            shape = (self.ny, self.nx)
            def f(x, y):
                if x < -6 or y < -6 or x > shape[1] + 5 or y > shape[0] + 5:
                    return np.zeros(shape)
                return p.locate(x, y, shape)
            self._prf = f
        return self._prf


def _mira_row():
    g = pd.read_csv(f'{BASE}/data/gaia_330.csv')
    r = g[g.source_id == IDS['Mira']].iloc[0]
    return dict(ra=float(r.ra), dec=float(r.dec), pmra=float(r.pmra), pmdec=float(r.pmdec))


MIRA_ROW = _mira_row()


def template(t, tc_list, exptime, t14=T14, flat=FLAT, nsub=11):
    """Exposure-integrated trapezoid (unit depth) for transits centred at tc_list."""
    t = np.asarray(t, float)
    off = (np.arange(nsub) - (nsub - 1) / 2) / nsub * exptime
    y = np.zeros_like(t)
    t23 = t14 * flat
    for tc in np.atleast_1d(tc_list):
        near = np.abs(t - tc) < t14
        if not near.any(): continue
        a = np.abs((t[near][:, None] + off[None, :]) - tc)
        v = np.where(a <= t23 / 2, 1.0, np.where(a >= t14 / 2, 0.0, (t14 / 2 - a) / ((t14 - t23) / 2)))
        y[near] += v.mean(1)
    return y


def real_tcs(t):
    n = np.arange(np.ceil((np.nanmin(t) - T0) / P) - 1, np.floor((np.nanmax(t) - T0) / P) + 2)
    return T0 + n * P


def to_json(obj, fn):
    def conv(o):
        if isinstance(o, (np.integer,)): return int(o)
        if isinstance(o, (np.floating,)): return float(o)
        if isinstance(o, np.ndarray): return o.tolist()
        if isinstance(o, (np.bool_,)): return bool(o)
        raise TypeError(type(o))
    with open(fn, 'w') as f:
        json.dump(obj, f, indent=1, default=conv)
