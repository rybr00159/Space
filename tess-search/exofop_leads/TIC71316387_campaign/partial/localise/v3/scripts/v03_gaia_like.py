"""Likelihood of the target's Gaia DR3 image diagnostics as a function of companion separation and contrast,
calibrated empirically on imaging-resolved pairs (v02 'pairs_gaia.csv') and imaging non-detections ('singles_gaia.csv').

Target diagnostics (D): ipd_frac_multi_peak = 41 (use 30-50), RUWE = 1.03 (< 1.4), no other Gaia source within 3"
(nearest 18.2"), BP blended in 47/47 transits (use >= 0.9).
For each (sep, dm) cell: L = k/n (pairs with G 9-14 meeting D), with a Jeffreys interval.  Single-star rate likewise.
Posterior over cells for flat priors in (log sep, dm) is reported for D and for D without the blending criterion,
for all distances and for parallax < 6 mas (to limit orbital-motion RUWE inflation in nearby calibrators).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from scipy.stats import beta
from v3lib import V3, to_json

P = pd.read_csv(V3 + '/speckle/pairs_gaia.csv'); S = pd.read_csv(V3 + '/speckle/singles_gaia.csv')
for d in (P, S):
    d['bpb'] = d.NBPblend / d.o_BPmag
SEPS = [0.1, 0.3, 0.5, 0.7, 0.9, 1.2, 1.6, 2.5]
DMS = [0, 2, 3, 4, 5, 7]


def jeff(k, n):
    return [float(beta.ppf(0.16, k + 0.5, n - k + 0.5)), float(beta.ppf(0.84, k + 0.5, n - k + 0.5))] if n else [None, None]


def crit(d, blend=True):
    m = (d.IPDfmp >= 30) & (d.IPDfmp <= 50) & (d.RUWE < 1.4)
    if 'n_gaia3' in d:
        m &= d.n_gaia3 == 0
    if blend:
        m &= d.bpb >= 0.9
    return m


out = {}
for dist_lab, dsel_p, dsel_s in [('all', np.ones(len(P), bool), np.ones(len(S), bool)), ('plx<6', P.Plx < 6, S.Plx < 6)]:
    for blend in (True, False):
        lab = f'{dist_lab}_{"blend" if blend else "noblend"}'
        gp = dsel_p & (P.Gmag > 9) & (P.Gmag < 14)
        gs = dsel_s & (S.Gmag > 9) & (S.Gmag < 14)
        ks = int((gs & crit(S, blend) & (S.RUWE < 1.4)).sum()); ns = int(gs.sum())
        cells = []
        for s0, s1 in zip(SEPS[:-1], SEPS[1:]):
            for d0, d1 in zip(DMS[:-1], DMS[1:]):
                m = gp & (P.sep >= s0) & (P.sep < s1) & (P.dm >= d0) & (P.dm < d1)
                n = int(m.sum()); k = int((m & crit(P, blend)).sum())
                Lm = (k + 0.5) / (n + 1) if n else np.nan      # Jeffreys posterior mean
                cells.append(dict(sep=[s0, s1], dm=[d0, d1], n=n, k=k, L=Lm, ci68=jeff(k, n)))
        # posterior for flat prior in log sep & dm: weight = L * dlogsep * ddm (cells with n >= 3 only)
        w = np.array([c['L'] * np.log(c['sep'][1] / c['sep'][0]) * (c['dm'][1] - c['dm'][0]) if c['n'] >= 3 else 0 for c in cells])
        w = w / w.sum()
        for c, ww in zip(cells, w):
            c['post'] = float(ww)
        # marginals
        msep = {f'{a}-{b}': float(sum(c['post'] for c in cells if c['sep'] == [a, b])) for a, b in zip(SEPS[:-1], SEPS[1:])}
        mdm = {f'{a}-{b}': float(sum(c['post'] for c in cells if c['dm'] == [a, b])) for a, b in zip(DMS[:-1], DMS[1:])}
        # conditional on dm 2-4 (TESS-implied contrast range)
        cond = [c for c in cells if c['dm'][0] >= 2 and c['dm'][1] <= 4 and c['n'] >= 3]
        wc = np.array([c['L'] * np.log(c['sep'][1] / c['sep'][0]) for c in cond]); wc = wc / wc.sum()
        msep_c = {}
        for c, ww in zip(cond, wc):
            key = f"{c['sep'][0]}-{c['sep'][1]}"; msep_c[key] = msep_c.get(key, 0) + float(ww)
        Lmax = max(c['L'] for c in cells if c['n'] >= 3)
        Ls = (ks + 0.5) / (ns + 1)
        out[lab] = dict(single=dict(k=ks, n=ns, L=Ls, ci68=jeff(ks, ns)), cells=cells, marg_sep=msep, marg_dm=mdm,
                        marg_sep_given_dm2to4=msep_c, Lmax_pair=Lmax, LR_bestpair_vs_single=Lmax / Ls)
        print(f'\n[{lab}] singles: {ks}/{ns} (L={Ls:.4f}); best pair-cell L={Lmax:.3f}; LR={Lmax / Ls:.0f}')
        print('  sep marginal:', {k: round(v, 3) for k, v in msep.items()})
        print('  dm  marginal:', {k: round(v, 3) for k, v in mdm.items()})
        print('  sep | dm 2-4:', {k: round(v, 3) for k, v in msep_c.items()})
        print('  cells k/n:')
        for s0, s1 in zip(SEPS[:-1], SEPS[1:]):
            print(f'   {s0:.1f}-{s1:.1f}: ' + '  '.join(f"{c['k']:2d}/{c['n']:<3d}" for c in cells if c['sep'] == [s0, s1]))
to_json(out, f'{V3}/results/v03_gaia_like.json')
