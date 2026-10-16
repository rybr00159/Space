"""TIC 90086266: robustness of the difference-image offset (windows, box sizes,
odd/even subsets, orbits) plus model-independent centroid comparison."""
import sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/five_referee/code')
import numpy as np
from reflib import *
from diffimg import *

TIC = 90086266
T14 = 3.0 / 24
P = 2.339474
EV = json.load(open(f'{BASE}/results/t1_events.json'))
OUT = {}


def centroid(img, x0, y0, box=2):
    xc, yc = int(round(x0)), int(round(y0))
    sub = img[yc - box:yc + box + 1, xc - box:xc + box + 1]
    yy, xx = np.mgrid[yc - box:yc + box + 1, xc - box:xc + box + 1]
    return np.sum(sub * xx) / np.sum(sub), np.sum(sub * yy) / np.sum(sub)


for sec in (13, 104):
    cut = load_cut(TIC, sec)
    prf = PRFModel(cut)
    evd = EV[f'S{sec}_CUT:r2.0']['events']
    tcs = np.array([e['tc'] for e in evd])
    Eabs = np.array([e['E'] for e in evd])
    g = gaia_table(TIC, btjd_to_jyear(np.median(cut['time'])))
    tgt = g.iloc[0]
    xg, yg = [float(v) for v in cut['wcs'].world_to_pixel_values(tgt['ra_ep'], tgt['dec_ep'])]
    print(f'S{sec} WCS target pixel: x={xg:.3f} y={yg:.3f}; PRF built at cam{cut["cam"]} ccd{cut["ccd"]} col0={cut["col0"]} row0={cut["row0"]}')
    rows = []
    for fin, gap, wout in [(0.3, 0.6, 0.15), (0.2, 0.6, 0.25), (0.4, 0.7, 0.2), (0.3, 0.8, 0.3)]:
        evs, direct = event_diff_images(cut, tcs, T14, fin=fin, gap=gap, wout=wout)
        derr = np.sqrt((0.002 * np.abs(direct)) ** 2 + 5.0 ** 2)
        for box in (2, 3, 4):
            D, E = combine(evs)
            fd = fit_prf(direct, derr, prf, xg, yg, box=box)
            fdf = fit_prf(D, E, prf, xg, yg, box=box)
            off = pix_to_EN(cut, fdf['x'], fdf['y'], fd['x'], fd['y'])
            rows.append(('all', fin, gap, wout, box, off[0], off[1], fdf['chi2'] / fdf['dof']))
        # subsets with box 3
        for name, sel in [('even', Eabs % 2 == 0), ('odd', Eabs % 2 == 1), ('orbit1', tcs < np.median(cut['time'])), ('orbit2', tcs >= np.median(cut['time']))]:
            idx = [i for i in range(len(evs)) if sel[np.argmin(np.abs(tcs - evs[i][0]))]]
            if len(idx) < 2:
                continue
            D, E = combine(evs, idx)
            fd = fit_prf(direct, derr, prf, xg, yg, box=3)
            fdf = fit_prf(D, E, prf, xg, yg, box=3)
            off = pix_to_EN(cut, fdf['x'], fdf['y'], fd['x'], fd['y'])
            rows.append((name, fin, gap, wout, 3, off[0], off[1], fdf['chi2'] / fdf['dof']))
        # model-independent: flux-weighted centroids (5x5) of diff and direct
        D, E = combine(evs)
        cdx, cdy = centroid(D, xg, yg, 2)
        cix, ciy = centroid(direct, xg, yg, 2)
        offc = pix_to_EN(cut, cdx, cdy, cix, ciy)
        rows.append(('centroid5x5', fin, gap, wout, 2, offc[0], offc[1], np.nan))
    for r in rows:
        print('  %-11s fin=%.1f gap=%.1f wout=%.2f box=%d  offset dE=%+.2f" dN=%+.2f"  chi2/dof=%.2f' % r)
    OUT[f'S{sec}'] = [dict(subset=r[0], fin=r[1], gap=r[2], wout=r[3], box=r[4], dE=r[5], dN=r[6], redchi2=r[7]) for r in rows]
save_json(OUT, f'{BASE}/results/t1_loc_robust.json')
