"""v4 SkyMapper Shallow Survey (unsaturated, FWHM 1.6-2.5"): template test of the VHS feature position.
Template = Gaia DR3 6571343707123804416 (G 13.64, 59" from the target; the brightest isolated star in the cutouts).
Joint fit: T = a*shift(P) + b + f*a*shift(P, + companion offset), primary shift re-optimised, at the VHS feature
position (1.1", PA 347) and at 7 control PAs at the same separation; injection calibration of f at control PAs."""
import sys, os, glob
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from scipy.ndimage import shift as ndshift
from scipy.optimize import minimize
import warnings; warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v4lib import V4, BASE, to_json

H = 10; PIX = 0.497; FLOOR = 0.01
g = Table.read(f'{V4}/gaia/wide_3p6deg_G17.fits').to_pandas()
IDS = {'T': 6571343840265900544, 'P': 6571343707123804416}


def pos(sid, ep):
    r = g[g.source_id == sid].iloc[0]
    return (r.ra + np.nan_to_num(r.pmra) * (ep - 2016) / 3.6e6 / np.cos(np.radians(r.dec)), r.dec + np.nan_to_num(r.pmdec) * (ep - 2016) / 3.6e6)


def stamp(d, w, sid, ep):
    x, y = w.all_world2pix([pos(sid, ep)], 0)[0]
    ix, iy = int(round(x)), int(round(y))
    big = d[iy - 30:iy + 31, ix - 30:ix + 31]; by, bx = np.mgrid[-30:31, -30:31]
    ann = big[np.hypot(bx, by) > 18]; sky = float(np.median(ann)); rms = float(1.4826 * np.median(np.abs(ann - sky)))
    raw = d[iy - H:iy + H + 1, ix - H:ix + H + 1].astype(float)
    return dict(s=raw - sky, rms=rms, fx=x - ix, fy=y - iy)


yy, xx = np.mgrid[-H:H + 1, -H:H + 1]


def solver(S, P, inj=None):
    s = S['s'] if inj is None else S['s'] + inj
    def solve(p, off=None):
        m = ndshift(P['s'], (p[1], p[0]), order=3, mode='constant', cval=0.0)
        cols = [m.ravel(), np.ones(m.size)]
        if off is not None:
            cols.append(ndshift(P['s'], (p[1] + off[1], p[0] + off[0]), order=3, mode='constant', cval=0.0).ravel())
        A = np.column_stack(cols)
        var = S['rms'] ** 2 + np.clip(s, 0, None).ravel() + (FLOOR * np.clip(m, 0, None).ravel() * (s.max() / max(P['s'].max(), 1))) ** 2
        wgt = 1 / np.sqrt(var)
        c, *_ = np.linalg.lstsq(A * wgt[:, None], s.ravel() * wgt, rcond=None)
        r = (A @ c - s.ravel()) * wgt
        return float(r @ r), c
    p0 = np.array([S['fx'] - P['fx'], S['fy'] - P['fy']])
    o = minimize(lambda p: solve(p)[0], p0, method='Nelder-Mead', options=dict(xatol=1e-3, fatol=1e-3))
    return o.x, solve


def fixed(S, P, off, inj=None):
    p0, solve = solver(S, P, inj)
    c20, _ = solve(p0)
    o = minimize(lambda q: solve(q, off)[0], p0, method='Nelder-Mead', options=dict(xatol=1e-3, fatol=1e-3))
    c2, c = solve(o.x, off)
    return float(c20 - c2), float(c[2] / c[0]), c20


out = []
for fn in sorted(glob.glob(f'{BASE}/v3/imaging/skymapper/*.fits')):
    d = fits.getdata(fn).astype(float); hd = fits.getheader(fn); w = WCS(hd)
    ep = 2014.66 if '2014' in os.path.basename(fn) else 2015.68
    S = stamp(d, w, IDS['T'], ep); P = stamp(d, w, IDS['P'], ep)
    x0, y0 = w.all_world2pix([pos(IDS['T'], ep)], 0)[0]
    r0 = w.all_pix2world([[x0, y0], [x0 + 1, y0], [x0, y0 + 1]], 0); cd = np.cos(np.radians(r0[0][1]))
    J = np.column_stack([[(r0[1][0] - r0[0][0]) * cd * 3600, (r0[1][1] - r0[0][1]) * 3600], [(r0[2][0] - r0[0][0]) * cd * 3600, (r0[2][1] - r0[0][1]) * 3600]])
    Ji = np.linalg.inv(J)
    res = {}
    for sep in (0.8, 1.1, 1.5):
        rows = []
        for pa in 347 + 45 * np.arange(8):
            off = Ji @ np.array([sep * np.sin(np.radians(pa)), sep * np.cos(np.radians(pa))])
            dch, f, c20 = fixed(S, P, off)
            rows.append(dict(pa=float(pa % 360), dchi2=dch, f=f))
        res[sep] = rows
    # injection calibration at 1.1" (control PAs), f_in = 0.025 (dm 4.0) and 0.06 (dm 3.05)
    p0, solve0 = solver(S, P)
    _, c0 = solve0(p0)
    cal = []
    for f_in in (0.025, 0.06):
        for pa in 347 + 45 * np.arange(1, 8):
            off = Ji @ np.array([1.1 * np.sin(np.radians(pa)), 1.1 * np.cos(np.radians(pa))])
            comp = c0[0] * f_in * ndshift(P['s'], (p0[1] + off[1], p0[0] + off[0]), order=3, mode='constant', cval=0.0)
            dch, f, _ = fixed(S, P, off, inj=comp)
            cal.append(dict(f_in=f_in, pa=float(pa % 360), dchi2=dch, f_out=f))
    out.append(dict(file=os.path.basename(fn), band=hd['FILTER'], chi2_0=c20, chi2r_0=c20 / (S['s'].size - 4), fixed=res, calib=cal))
    r11 = res[1.1]
    print(f"{os.path.basename(fn)} {hd['FILTER']}: chi2r0 {c20 / (S['s'].size - 4):.2f}; at 1.1\" feature PA 347: dchi2 {r11[0]['dchi2']:.1f} f {r11[0]['f']:+.4f}; "
          f"controls dchi2 {np.round([r['dchi2'] for r in r11[1:]], 1)} f {np.round([r['f'] for r in r11[1:]], 4)}", flush=True)
    for f_in in (0.025, 0.06):
        cc = [c for c in cal if c['f_in'] == f_in]
        print(f"     inj f={f_in}: recovered f {np.round([c['f_out'] for c in cc], 4)} dchi2 {np.round([c['dchi2'] for c in cc], 1)}")
to_json(out, f'{V4}/results/h05_skymapper.json')
