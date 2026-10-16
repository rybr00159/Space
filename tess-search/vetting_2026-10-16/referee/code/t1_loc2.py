"""TIC 90086266: localisation calibrated with injections that use the real
direct image as the template (no PRF-model self-consistency), at several
orbital phases, for several window choices; plus PRF fits with a planar
background in the difference image."""
import sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/five_referee/code')
import numpy as np
from scipy.optimize import least_squares
from reflib import *
from diffimg import *

TIC = 90086266
T14 = 3.0 / 24
P = 2.339474
EV = json.load(open(f'{BASE}/results/t1_events.json'))
OUT = {}


def fit_prf_plane(img, err, prf, x0, y0, box=3):
    ny, nx = img.shape
    xc, yc = int(round(x0)), int(round(y0))
    ys, xs = slice(yc - box, yc + box + 1), slice(xc - box, xc + box + 1)
    sub, se = img[ys, xs], err[ys, xs]
    yy, xx = np.mgrid[ys, xs]

    def resid(p):
        x, y, A, c, gx, gy = p
        m = A * prf(x, y)[ys, xs] + c + gx * (xx - xc) + gy * (yy - yc)
        return ((sub - m) / se).ravel()
    A0 = np.nansum(sub)
    r = least_squares(resid, [x0, y0, A0, 0, 0, 0], x_scale=[0.05, 0.05, abs(A0) * 0.05, abs(A0) * 0.01, abs(A0) * 0.01, abs(A0) * 0.01],
                      diff_step=[1e-3, 1e-3, 1e-4, 1e-4, 1e-4, 1e-4])
    return r.x[0], r.x[1]


for sec in (13, 104):
    cut = load_cut(TIC, sec)
    prf = PRFModel(cut)
    tcs = np.array([e['tc'] for e in EV[f'S{sec}_CUT:r2.0']['events']])
    g = gaia_table(TIC, btjd_to_jyear(np.median(cut['time'])))
    tgt = g.iloc[0]
    xg, yg = [float(v) for v in cut['wcs'].world_to_pixel_values(tgt['ra_ep'], tgt['dec_ep'])]
    res = []
    for fin, gap, wout in [(0.3, 0.6, 0.15), (0.2, 0.6, 0.25), (0.4, 0.7, 0.2), (0.3, 0.8, 0.3)]:
        evs, direct = event_diff_images(cut, tcs, T14, fin=fin, gap=gap, wout=wout)
        derr = np.sqrt((0.002 * np.abs(direct)) ** 2 + 5.0 ** 2)
        D, E = combine(evs)
        fd = fit_prf(direct, derr, prf, xg, yg, box=3)
        fdf = fit_prf(D, E, prf, xg, yg, box=3)
        real = pix_to_EN(cut, fdf['x'], fdf['y'], fd['x'], fd['y'])
        xp, yp = fit_prf_plane(D, E, prf, xg, yg, box=3)
        realp = pix_to_EN(cut, xp, yp, fd['x'], fd['y'])
        depth = fdf['A'] / fd['A']
        # template = real direct image (background-subtracted already), V-shaped injection
        tmpl = direct - np.median(direct)
        injs, injsp = [], []
        for phs in (0.17, 0.25, 0.33, 0.5, 0.67, 0.75, 0.83):
            tci = tcs + phs * P
            tci = tci[(tci > cut['time'].min() + 0.4) & (tci < cut['time'].max() - 0.4)]
            # avoid real eclipses
            tci = np.array([x for x in tci if np.min(np.abs(x - tcs)) > 0.5])
            evi, di = event_diff_images(cut, tci, T14, fin=fin, gap=gap, wout=wout, inject=(tmpl, depth / (1 - fin / 1.0 * 0.5) ))
            if len(evi) < 4:
                continue
            Di, Ei = combine(evi)
            fdi = fit_prf(di, derr, prf, xg, yg, box=3)
            fi = fit_prf(Di, Ei, prf, xg, yg, box=3)
            injs.append(pix_to_EN(cut, fi['x'], fi['y'], fdi['x'], fdi['y']))
            xpi, ypi = fit_prf_plane(Di, Ei, prf, xg, yg, box=3)
            injsp.append(pix_to_EN(cut, xpi, ypi, fdi['x'], fdi['y']))
        injs, injsp = np.array(injs), np.array(injsp)
        rec = dict(fin=fin, gap=gap, wout=wout, real=list(real), real_plane=list(realp), inj_mean=injs.mean(0).tolist(), inj_std=injs.std(0).tolist(),
                   inj_plane_mean=injsp.mean(0).tolist(), inj_plane_std=injsp.std(0).tolist(), n_inj=len(injs))
        res.append(rec)
        print(f'S{sec} fin={fin} gap={gap} wout={wout}: real ({real[0]:+.2f},{real[1]:+.2f})"  plane ({realp[0]:+.2f},{realp[1]:+.2f})"  | '
              f'direct-template injections mean ({injs.mean(0)[0]:+.2f},{injs.mean(0)[1]:+.2f}) sd ({injs.std(0)[0]:.2f},{injs.std(0)[1]:.2f}) n={len(injs)}; '
              f'plane: mean ({injsp.mean(0)[0]:+.2f},{injsp.mean(0)[1]:+.2f}) sd ({injsp.std(0)[0]:.2f},{injsp.std(0)[1]:.2f})')
    OUT[f'S{sec}'] = res
save_json(OUT, f'{BASE}/results/t1_loc2.json')
