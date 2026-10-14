"""v4 VISTA/VHS closed-loop template test (independent of v3's code; uses the SODA cutouts v3 downloaded).

Stars on the same VIRCAM detector in every frame:
  T   = target (G 11.78, Ks 10.34)
  C1  = Gaia DR3 6571344458741189760 (G 11.83, BP-RP 0.90; 3.8' away)  -- equally saturated
  C2  = Gaia DR3 6571344737915960064 (G 13.39, BP-RP 0.94)             -- 1.6 mag fainter, (near-)linear core in Ks
For every ordered pair (star S, template P): S = a * shift(P, dx, dy) + b, pixels with raw counts > SAT in S or in the
shifted template masked (+1 px dilation); then a companion scan: + f * a * shift(P, dx + ex, dy + ey) on a 0.1" grid
(0.4-3.0"), f free (either sign) -> dchi2 map.  A real companion of T shows up as f > 0 at the same sky position in
(T|C1) and (T|C2), as f < 0 in (C1|T) and (C2|T), and NOT in (C1|C2) / (C2|C1).
Noise: var = sky_rms^2 + max(S,0)/GAIN + (FLOOR * model)^2.  dchi2 values are reported raw and divided by chi2_r.
Injection: fake companions (template-shaped, f = 10^(-0.4 dm)) added to T at 0.5-2.0", 8 PAs; recovered if the scan
peak lies within 0.25" of the injection and its dchi2 exceeds the largest dchi2 seen in the comparison-only pairs.
"""
import sys, os, glob, json
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from scipy.ndimage import shift as ndshift, binary_dilation
from scipy.optimize import minimize
import warnings; warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v4lib import V4, BASE, to_json

SAT = float(os.environ.get('H01_SAT', 22000.0)); GAIN = 4.2; FLOOR = float(os.environ.get('H01_FLOOR', 0.02))
H = 14; PIX = 0.339; EP = 2012.64
IDS = {'T': 6571343840265900544, 'C1': 6571344458741189760, 'C2': 6571344737915960064}
g = Table.read(f'{V4}/gaia/wide_3p6deg_G17.fits').to_pandas()
POS = {}
for k, sid in IDS.items():
    r = g[g.source_id == sid].iloc[0]
    POS[k] = (r.ra + np.nan_to_num(r.pmra) * (EP - 2016) / 3.6e6 / np.cos(np.radians(r.dec)), r.dec + np.nan_to_num(r.pmdec) * (EP - 2016) / 3.6e6)
FILES = sorted(glob.glob(f'{BASE}/v3/imaging/big/*.fz')) + sorted(glob.glob(f'{BASE}/v3/imaging/vhs_ADP.2014-12-22T14:5[67]*.fz'))
yy, xx = np.mgrid[-H:H + 1, -H:H + 1]


def stamp(d, w, k):
    x, y = w.all_world2pix([POS[k]], 0)[0]
    ix, iy = int(round(x)), int(round(y))
    big = d[iy - 40:iy + 41, ix - 40:ix + 41]
    by, bx = np.mgrid[-40:41, -40:41]
    ann = big[(np.hypot(bx, by) > 25)]
    sky = float(np.median(ann)); rms = float(1.4826 * np.median(np.abs(ann - sky)))
    raw = d[iy - H:iy + H + 1, ix - H:ix + H + 1].astype(float)
    return dict(raw=raw, s=raw - sky, sky=sky, rms=rms, fx=x - ix, fy=y - iy)


def fit_pair(S, P, ex=None, ey=None, f=None, inj=None):
    """Fit S with shifted template P; optional companion at (ex, ey) px offset from S's centre; returns chi2, params."""
    s = S['s'] if inj is None else S['s'] + inj
    rawS = S['raw'] if inj is None else S['raw'] + inj
    T = P['s']
    def model(p, comp=None):
        m = ndshift(T, (p[1], p[0]), order=3, mode='constant', cval=0.0)
        return m
    def solve(p, comp_off=None):
        m = model(p)
        satm = ndshift((P['raw'] > SAT).astype(float), (p[1], p[0]), order=1, mode='constant', cval=0.0) > 0.01
        mask = ~binary_dilation((rawS > SAT) | satm, iterations=1) & (np.hypot(xx, yy) <= H)
        cols = [m[mask], np.ones(mask.sum())]
        if comp_off is not None:
            mc = ndshift(T, (p[1] + comp_off[1], p[0] + comp_off[0]), order=3, mode='constant', cval=0.0)
            cols.append(mc[mask])
        A = np.column_stack(cols)
        var = S['rms'] ** 2 + np.clip(s[mask], 0, None) / GAIN
        # floor on the template model (iterate once with a provisional amplitude)
        c, *_ = np.linalg.lstsq(A / np.sqrt(var)[:, None], s[mask] / np.sqrt(var), rcond=None)
        var = var + (FLOOR * (A @ c)) ** 2
        w = 1 / np.sqrt(var)
        c, *_ = np.linalg.lstsq(A * w[:, None], s[mask] * w, rcond=None)
        r = (A @ c - s[mask]) * w
        return float(r @ r), c, int(mask.sum())
    # centre: S centre - P centre (in px), refined
    p0 = np.array([S['fx'] - P['fx'], S['fy'] - P['fy']])
    o = minimize(lambda p: solve(p)[0], p0, method='Nelder-Mead', options=dict(xatol=1e-3, fatol=1e-3, maxiter=200))
    c2, c, npx = solve(o.x)
    return o.x, c2, c, npx, solve


def scan(S, P, inj=None, rmin=0.4, rmax=3.0, step=0.1):
    p, c2, c, npx, solve = fit_pair(S, P, inj=inj)
    chi2r = c2 / (npx - 4)
    best = None; rows = []
    for ex in np.arange(-rmax, rmax + 1e-9, step):
        for ey in np.arange(-rmax, rmax + 1e-9, step):
            r = np.hypot(ex, ey)
            if r < rmin or r > rmax:
                continue
            c2c, cc, _ = solve(p, comp_off=(ex / PIX, ey / PIX))
            d = c2 - c2c
            f = cc[2] / cc[0]
            rows.append((ex, ey, d, f))
            if best is None or d > best[2]:
                best = (ex, ey, d, f)
    return dict(p=p.tolist(), chi2=c2, chi2r=chi2r, npix=npx, amp=float(c[0]), best=best, grid=np.array(rows))


def sky_of(ex, ey, w, xpix):
    """pixel offset (east/north in arcsec) conversion: use WCS local Jacobian to map pixel (dx, dy) -> (dE, dN)."""
    pass


def main():
    out = []
    for fn in FILES:
        hd0 = fits.getheader(fn, 0)
        for ext in range(1, len(fits.open(fn))):
            d = fits.getdata(fn, ext).astype(float); hd = fits.getheader(fn, ext)
            if min(d.shape) < 200:
                continue
            w = WCS(hd)
            band = hd0.get('FILTER') or hd.get('FILTER') or '?'
            st = {k: stamp(d, w, k) for k in IDS}
            # local pixel->sky: arcsec per px in E and N for +x, +y
            x0, y0 = w.all_world2pix([POS['T']], 0)[0]
            r0 = w.all_pix2world([[x0, y0], [x0 + 1, y0], [x0, y0 + 1]], 0)
            cd = np.cos(np.radians(r0[0][1]))
            Jx = np.array([(r0[1][0] - r0[0][0]) * cd * 3600, (r0[1][1] - r0[0][1]) * 3600])
            Jy = np.array([(r0[2][0] - r0[0][0]) * cd * 3600, (r0[2][1] - r0[0][1]) * 3600])
            J = np.column_stack([Jx, Jy]); Ji = np.linalg.inv(J)
            peaks = {k: float(v['raw'].max()) for k, v in st.items()}
            res = {}
            for a, b in [('T', 'C1'), ('C1', 'T'), ('T', 'C2'), ('C2', 'T'), ('C1', 'C2'), ('C2', 'C1')]:
                r = scan(st[a], st[b])
                ex, ey, dch, f = r['best']
                sky_off = J @ np.array([ex, ey]) / PIX
                # value at the candidate position of interest (1.12" PA 350) and at its mirror
                g_ = r['grid']
                def at(sep, pa):
                    e_, n_ = sep * np.sin(np.radians(pa)), sep * np.cos(np.radians(pa))
                    pxo = Ji @ np.array([e_, n_]) * PIX
                    k = np.argmin(np.hypot(g_[:, 0] - pxo[0], g_[:, 1] - pxo[1]))
                    return float(g_[k, 2]), float(g_[k, 3])
                d1, f1 = at(1.12, 350.0)
                d2, f2 = at(0.6, 335.0)
                res[f'{a}|{b}'] = dict(chi2r=r['chi2r'], npix=r['npix'], best_dchi2=dch, best_dchi2_over_chi2r=dch / r['chi2r'],
                                       best_sep=float(np.hypot(*sky_off)), best_pa=float(np.degrees(np.arctan2(*sky_off)) % 360), best_f=f,
                                       at112_350=dict(dchi2=d1, f=f1), at060_335=dict(dchi2=d2, f=f2))
                print(f"{os.path.basename(fn)[4:30]} e{ext} {band:2s} {a}|{b}: chi2r {r['chi2r']:6.2f} best dchi2 {dch:8.1f} ({dch / r['chi2r']:6.1f}) "
                      f"sep {np.hypot(*sky_off):.2f}\" PA {np.degrees(np.arctan2(*sky_off)) % 360:5.1f} f {f:+.4f} | @1.12\"/350: dchi2 {d1:7.1f} f {f1:+.4f} "
                      f"| @0.6\"/335: dchi2 {d2:7.1f} f {f2:+.4f}", flush=True)
            # injection into T (fitted with template C1): fake companion = C2 image (linear core) scaled to f * (T flux)
            # detected if dchi2 AT the injected position exceeds the largest dchi2 anywhere in the comparison-only scans
            thr = max(res['C1|C2']['best_dchi2'], res['C2|C1']['best_dchi2'])
            p2, _, c2c, _, _ = fit_pair(st['T'], st['C2'])
            pT, _, _, _, solveT = fit_pair(st['T'], st['C1'])
            inj = []
            for sep in (0.5, 0.75, 1.0, 1.5, 2.0):
                for dm in (1.5, 2.0, 2.5, 3.0, 3.5, 4.0):
                    ok = 0; n = 0
                    for pa in np.arange(10, 360, 45.0):
                        e_, n_ = sep * np.sin(np.radians(pa)), sep * np.cos(np.radians(pa))
                        pxo = Ji @ np.array([e_, n_])
                        comp = c2c[0] * 10 ** (-0.4 * dm) * ndshift(st['C2']['s'], (p2[1] + pxo[1], p2[0] + pxo[0]), order=3, mode='constant', cval=0.0)
                        p_, c0, c_, npx_, solve_ = fit_pair(st['T'], st['C1'], inj=comp)
                        c1, cc, _ = solve_(p_, comp_off=(pxo[0], pxo[1]))
                        n += 1
                        if (c0 - c1) > thr and cc[2] > 0:
                            ok += 1
                    inj.append(dict(sep=sep, dm=dm, frac=ok / n))
                print(f"   inj sep {sep:.2f}\": " + '  '.join(f"dm{r['dm']}:{r['frac']:.2f}" for r in inj if r['sep'] == sep) + f"  (thr {thr:.0f})", flush=True)
            out.append(dict(file=os.path.basename(fn), ext=ext, band=band, fwhm=hd.get('PSF_FWHM'), abmagsat=hd.get('ABMAGSAT'),
                            peaks_raw=peaks, sky={k: v['sky'] for k, v in st.items()}, pairs=res, threshold=thr, injection=inj))
    to_json(out, f'{V4}/results/h01_vhs_loop.json')


if __name__ == '__main__':
    main()
