"""Referee-audit helpers (written from scratch for /home/user/work/tic87/referee_audit/).

Conventions (checked explicitly in r1_wcs_check.py):
  * FLUX cube is indexed [cadence, row, col] = [t, y, x].
  * Pixel coordinates are 0-based, integer = pixel centre (astropy origin 0).
  * WCS from the cutout APERTURE extension (ext 2). Cross-checked against the ext-1 column WCS and tess-point.
  * Time is BTJD = BJD_TDB - 2457000 (TESScut TIME, TIMEREF=SOLARSYSTEM, TIMESYS=TDB, BJDREFI=2457000).
  * Difference image D = baseline(in-transit times) - in-transit flux  => a dimming star gives D > 0.
"""
import glob
import json
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS

DATA = '/home/user/work/tic87/data'
OUT = '/home/user/work/tic87/referee_audit'
SECTORS = [13, 39, 66, 93, 104]
P, T0, T14 = 4.136195, 1658.41325, 0.1875
TID, NID = 6720986204578421888, 6720986208895208960
TIC_RA, TIC_DEC = 272.606557, -43.995371
S104_BAD = (4178.8, 4179.2)   # scattered-light episode at the start of S104 (CONTEXT.md)


def cutout_path(s):
    return glob.glob(f'{DATA}/tesscut/tess-s{s:04d}-*_21x21_astrocut.fits')[0]


def tmag(g, bprp):
    g = np.asarray(g, float); c = np.asarray(bprp, float)
    t = g - 0.00522555 * c**3 + 0.0891337 * c**2 - 0.633923 * c + 0.0324473
    return np.where(np.isfinite(c), t, g - 0.430)


def gaia():
    """Shared 90" Gaia DR3 CSV plus my own VizieR copy (G < 18 within 330", fetch_gaia.py) for the rest of the stamp."""
    d = pd.read_csv(f'{DATA}/gaia_dr3_90arcsec.csv')
    d = d[['source_id', 'ra', 'dec', 'pmra', 'pmdec', 'ref_epoch', 'phot_g_mean_mag', 'bp_rp']].copy()
    w = pd.read_csv(f'{OUT}/gaia/gaia_dr3_330as_G18_vizier.csv')
    w = w[~w.Source.isin(d.source_id)]
    w = pd.DataFrame(dict(source_id=w.Source, ra=w.RA_ICRS, dec=w.DE_ICRS, pmra=w.pmRA, pmdec=w.pmDE,
                          ref_epoch=2016.0, phot_g_mean_mag=w.Gmag, bp_rp=w['BP-RP']))
    d = pd.concat([d, w], ignore_index=True)
    d['tmag'] = tmag(d.phot_g_mean_mag, d.bp_rp)
    return d.copy()


def jyear(btjd):
    return 2000.0 + (np.asarray(btjd) + 2457000.0 - 2451545.0) / 365.25


def propagate(d, ep):
    dt = ep - d.ref_epoch.values
    pmra = np.nan_to_num(d.pmra.values); pmde = np.nan_to_num(d.pmdec.values)
    dec = d.dec.values + pmde * dt / 3.6e6
    ra = d.ra.values + pmra * dt / 3.6e6 / np.cos(np.radians(d.dec.values))
    return ra, dec


class Cube:
    def __init__(self, s, drop_bad=True):
        self.s = s
        with fits.open(cutout_path(s), memmap=False) as h:
            self.h0, self.h1, self.h2 = h[0].header, h[1].header, h[2].header
            d = h[1].data
            self.t = np.array(d['TIME'], float)
            self.f = np.array(d['FLUX'], float)
            self.fe = np.array(d['FLUX_ERR'], float)
            self.fb = np.array(d['FLUX_BKG'], float)
            self.q = np.array(d['QUALITY'], int)
            self.tcorr = np.array(d['TIMECORR'], float)
        self.cam, self.ccd = int(self.h0['CAMERA']), int(self.h0['CCD'])
        self.col0, self.row0 = int(self.h1['1CRV4P']), int(self.h1['2CRV4P'])
        self.w = WCS(self.h2)
        self.ny, self.nx = self.f.shape[1:]
        fin = np.isfinite(self.t) & np.all(np.isfinite(self.f.reshape(len(self.t), -1)), axis=1)
        self.good = (self.q == 0) & fin
        if drop_bad:
            self.good &= ~((self.t > S104_BAD[0]) & (self.t < S104_BAD[1]))
        self.dtc = float(np.nanmedian(np.diff(self.t)))
        self.epoch = float(jyear(np.nanmedian(self.t)))
        g = gaia()
        ra, de = propagate(g, self.epoch)
        x, y = self.w.all_world2pix(ra, de, 0)
        g = g.assign(ra_ep=ra, dec_ep=de, x=x, y=y)
        self.g = g
        it = int(np.where(g.source_id == TID)[0][0]); ineb = int(np.where(g.source_id == NID)[0][0])
        self.xt, self.yt = float(x[it]), float(y[it])
        self.xn, self.yn = float(x[ineb]), float(y[ineb])
        self.it, self.ineb = it, ineb
        v = np.array([self.xn - self.xt, self.yn - self.yt])
        self.sep_px = float(np.hypot(*v))
        self.u = v / self.sep_px                  # unit vector target -> neighbour (pixel frame)
        self.uperp = np.array([-self.u[1], self.u[0]])

    def epochs(self):
        tg = self.t[self.good]
        n = np.arange(np.ceil((tg.min() - T0) / P) - 1, np.floor((tg.max() - T0) / P) + 2).astype(int)
        return n

    def window(self, tc, in_half=0.375 * T14, oot_in=0.75 * T14, oot_out=2.0 * T14, min_in=0.5, min_side=0.25):
        dt = self.t - tc; g = self.good
        inn = np.where(g & (np.abs(dt) < in_half))[0]
        lft = np.where(g & (dt <= -oot_in) & (dt >= -oot_out))[0]
        rgt = np.where(g & (dt >= oot_in) & (dt <= oot_out))[0]
        n_in_nom = 2 * in_half / self.dtc; n_side_nom = (oot_out - oot_in) / self.dtc
        ok = len(inn) >= min_in * n_in_nom and len(lft) >= min_side * n_side_nom and len(rgt) >= min_side * n_side_nom
        return dict(tc=tc, inn=inn, lft=lft, rgt=rgt, ok=bool(ok))


def diff_image(t, cube, win, linear=True):
    """Per-pixel baseline fitted to both OOT sides; D = baseline(in) - mean flux(in). Returns D, sigma_D, O."""
    oot = np.concatenate([win['lft'], win['rgt']]); inn = win['inn']; tc = win['tc']
    Ao = np.column_stack([np.ones(len(oot))] + ([t[oot] - tc] if linear else []))
    Ai = np.column_stack([np.ones(len(inn))] + ([t[inn] - tc] if linear else []))
    Y = cube[oot].reshape(len(oot), -1)
    coef, *_ = np.linalg.lstsq(Ao, Y, rcond=None)
    res = Y - Ao @ coef
    rms = np.sqrt((res**2).sum(0) / max(len(oot) - Ao.shape[1], 1))
    D = (Ai @ coef).mean(0) - cube[inn].reshape(len(inn), -1).mean(0)
    abar = Ai.mean(0)
    var_pred = abar @ np.linalg.pinv(Ao.T @ Ao) @ abar
    sig = rms * np.sqrt(1.0 / len(inn) + var_pred)
    O = Y.mean(0)
    shp = cube.shape[1:]
    return D.reshape(shp), sig.reshape(shp), O.reshape(shp)


def dump(obj, path):
    def conv(o):
        if isinstance(o, np.floating): return float(o)
        if isinstance(o, np.integer): return int(o)
        if isinstance(o, np.ndarray): return o.tolist()
        if isinstance(o, np.bool_): return bool(o)
        raise TypeError(type(o))
    with open(path, 'w') as f:
        json.dump(obj, f, indent=1, default=conv)
