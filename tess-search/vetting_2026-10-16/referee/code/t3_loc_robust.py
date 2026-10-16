"""TIC 214311906: robustness of the ~19" W localisation to sky handling, windows,
fit box, and to which star is used as the reference (offset quoted relative to
the Gaia position of TIC 214311908 as well)."""
import sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/five_referee/code')
import numpy as np
from reflib import *
from diffimg import *

TIC = 214311906
P, T0, T14 = 5.69149, 3840.913, 1.55 / 24
A_ID = 5925373272850784768
OUT = {}
for variant in [dict(sky='p20', fin=0.35, gap=0.65, wout=0.12, box=2), dict(sky='none', fin=0.35, gap=0.65, wout=0.12, box=2),
                dict(sky='p20', fin=0.25, gap=0.8, wout=0.2, box=2), dict(sky='p20', fin=0.35, gap=0.65, wout=0.12, box=3),
                dict(sky='p05', fin=0.3, gap=0.7, wout=0.15, box=2)]:
    res = []
    for sec in (12, 39, 66, 93, 103, 104):
        cut = load_cut(TIC, sec)
        t = cut['time']
        if variant['sky'] != 'none':
            q = int(variant['sky'][1:])
            sky = np.percentile(cut['flux'].reshape(len(t), -1), q, axis=1)
            cut['flux'] = cut['flux'] - sky[:, None, None]
        n0, n1 = int(np.floor((t.min() - T0) / P)), int(np.ceil((t.max() - T0) / P))
        tcs = T0 + np.arange(n0, n1 + 1) * P
        tcs = tcs[(tcs > t.min() + 0.1) & (tcs < t.max() - 0.1)]
        evs, direct = event_diff_images(cut, tcs, T14, fin=variant['fin'], gap=variant['gap'], wout=variant['wout'])
        D, E = combine(evs)
        g = gaia_table(TIC, btjd_to_jyear(np.median(t)))
        tgt = g.iloc[0]; a = g[g.source_id == A_ID].iloc[0]
        xt, yt = [float(v) for v in cut['wcs'].world_to_pixel_values(tgt['ra_ep'], tgt['dec_ep'])]
        prf = PRFModel(cut)
        derr = np.sqrt((0.003 * np.abs(direct)) ** 2 + 10 ** 2)
        ft = fit_prf(direct, derr, prf, xt, yt, box=1)
        fd = fit_prf(D, E, prf, xt, yt, box=variant['box'])
        o = pix_to_EN(cut, fd['x'], fd['y'], ft['x'], ft['y'])
        res.append((sec, o[0], o[1], fd['A'], a['dE'], a['dN']))
    r = np.array(res)
    mE, mN = r[:, 1].mean(), r[:, 2].mean()
    sE, sN = r[:, 1].std(ddof=1) / np.sqrt(len(r)), r[:, 2].std(ddof=1) / np.sqrt(len(r))
    key = '_'.join(f'{k}{v}' for k, v in variant.items())
    OUT[key] = dict(per_sector=r.tolist(), mean=[mE, mN], sem=[sE, sN], dF_mean=float(r[:, 3].mean()))
    print(f"{key:45s} per-sector dE: {' '.join(f'{x:+5.1f}' for x in r[:,1])} | dN: {' '.join(f'{x:+5.1f}' for x in r[:,2])} | mean ({mE:+.1f}+-{sE:.1f}, {mN:+.1f}+-{sN:.1f})\" ; "
          f"TIC214311908 at ({r[0,4]:+.1f},{r[0,5]:+.1f})\"; dF {r[:,3].mean():.1f} e/s")
save_json(OUT, f'{BASE}/results/t3_loc_robust.json')
