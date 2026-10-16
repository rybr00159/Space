"""TIC 214311906 field: injection-recovery calibration of the difference-image
localisation. Inject the observed flux deficit (16.5 e-/s, flat 1.55-h shape)
on the Mira, on M dwarf B (TIC 214311910, ~7" W) and on M dwarf A (TIC 214311908,
~19" W) at non-eclipse phases, run the identical pipeline, and compare the
recovered offsets with the real one."""
import sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/five_referee/code')
import numpy as np
from reflib import *
from diffimg import *

TIC = 214311906
P, T0, T14 = 5.69149, 3840.913, 1.55 / 24
IDS = {'Mira': 5925373268528684288, 'B_7W': 5925373272850782208, 'A_19W': 5925373272850784768}
DF = 16.5
OUT = {}


def run(cut, tcs, prf, xt, yt, inject=None):
    c = dict(cut)
    if inject is not None:
        tmpl, pos_tcs = inject
        F = c['flux'].copy()
        t = c['time']
        for tc in pos_tcs:
            x = np.abs(t - tc)
            w = np.clip((T14 / 2 - x) / (0.2 * T14), 0, 1)  # flat-bottomed trapezoid, ingress 0.2*T14
            F -= w[:, None, None] * tmpl[None]
        c['flux'] = F
    evs, direct = event_diff_images(c, tcs, T14, fin=0.35, gap=0.65, wout=0.12)
    D, E = combine(evs)
    derr = np.sqrt((0.003 * np.abs(direct)) ** 2 + 10 ** 2)
    ft = fit_prf(direct, derr, prf, xt, yt, box=1)
    fd = fit_prf(D, E, prf, xt, yt, box=2)
    return pix_to_EN(cut, fd['x'], fd['y'], ft['x'], ft['y'])


for sec in (12, 39, 66, 93, 103, 104):
    cut = load_cut(TIC, sec)
    t = cut['time']
    sky = np.percentile(cut['flux'].reshape(len(t), -1), 20, axis=1)
    cut['flux'] = cut['flux'] - sky[:, None, None]
    g = gaia_table(TIC, btjd_to_jyear(np.median(t)))
    pix = {k: [float(v) for v in cut['wcs'].world_to_pixel_values(g[g.source_id == s].iloc[0]['ra_ep'], g[g.source_id == s].iloc[0]['dec_ep'])] for k, s in IDS.items()}
    gaia_off = {k: (float(g[g.source_id == s].iloc[0]['dE']), float(g[g.source_id == s].iloc[0]['dN'])) for k, s in IDS.items()}
    xt, yt = pix['Mira']
    prf = PRFModel(cut)
    n0, n1 = int(np.floor((t.min() - T0) / P)), int(np.ceil((t.max() - T0) / P))
    real_tcs = T0 + np.arange(n0, n1 + 1) * P
    real_tcs = real_tcs[(real_tcs > t.min() + 0.1) & (real_tcs < t.max() - 0.1)]
    real = run(cut, real_tcs, prf, xt, yt)
    rec = dict(real=list(real), gaia=gaia_off, inj={})
    for k in IDS:
        tmpl = prf(*pix[k]) * DF
        offs = []
        for phs in (0.25, 0.5, 0.75):
            tci = real_tcs + phs * P
            tci = tci[(tci > t.min() + 0.15) & (tci < t.max() - 0.15)]
            if len(tci) < 2:
                continue
            offs.append(run(cut, tci, prf, xt, yt, inject=(tmpl, tci)))
        offs = np.array(offs)
        rec['inj'][k] = dict(mean=offs.mean(0).tolist(), sd=offs.std(0).tolist(), n=len(offs))
    OUT[sec] = rec
    print(f"S{sec:3d}: real ({real[0]:+5.1f},{real[1]:+5.1f})\" | injected on Mira -> ({rec['inj']['Mira']['mean'][0]:+5.1f},{rec['inj']['Mira']['mean'][1]:+5.1f}); "
          f"on B {gaia_off['B_7W'][0]:+.1f},{gaia_off['B_7W'][1]:+.1f} -> ({rec['inj']['B_7W']['mean'][0]:+5.1f},{rec['inj']['B_7W']['mean'][1]:+5.1f}); "
          f"on A {gaia_off['A_19W'][0]:+.1f},{gaia_off['A_19W'][1]:+.1f} -> ({rec['inj']['A_19W']['mean'][0]:+5.1f},{rec['inj']['A_19W']['mean'][1]:+5.1f})  [sd of injections ~({rec['inj']['A_19W']['sd'][0]:.1f},{rec['inj']['A_19W']['sd'][1]:.1f})]")
# summary: mean over sectors of (real - injection) for each hypothesis
for k in IDS:
    d = np.array([np.array(OUT[s]['real']) - np.array(OUT[s]['inj'][k]['mean']) for s in OUT if isinstance(s, int)])
    m = d.mean(0); se = d.std(0, ddof=1) / np.sqrt(len(d))
    chi = np.sum((m / se) ** 2)
    OUT[f'summary_{k}'] = dict(mean_real_minus_inj=m.tolist(), sem=se.tolist(), chi2_2dof=float(chi))
    print(f'hypothesis {k:6s}: mean(real - injected) = ({m[0]:+.1f}+-{se[0]:.1f}, {m[1]:+.1f}+-{se[1]:.1f})"  chi2(2 dof)={chi:.1f}')
save_json(OUT, f'{BASE}/results/t3_inject.json')
