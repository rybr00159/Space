"""Collect referee fit results (out/fit_*.json, out/inj_*.npz) into out/summary.json and print tables.
For each run: real combined chi2(h) - chi2(M3V) for the main hypotheses, free-fit position, and from injections:
 - p_true(h) = fraction of injections ON h whose chi2(h)-chi2(M3V) >= real value (small = h excluded)
 - p_M3V(h)  = fraction of injections on the M3V whose chi2(h)-chi2(M3V) <= real value (two-sided consistency)
 - Gaussian-approx likelihood ratio L(M3V)/L(h) of the real statistic X_h = chi2(h)-chi2(M3V), using the mean/sd of
   X_h under each truth (an extrapolation where p_true = 0).
"""
import sys, glob; sys.path.insert(0, 'tic214311908/analysis/referee_loc/scripts')
from rl_lib import *
from scipy.stats import norm

S = {}
for f in sorted(glob.glob(f'{ROOT}/out/fit_*.json')):
    tag = f.split('/fit_')[-1][:-5]
    d = json.load(open(f)); MAIN = d['main']
    r = dict(real=d['real']['comb_dchi2'], free=d['real']['comb_free_pos'], free_dchi2=d['real']['comb_free_dchi2_vs_M3V'],
             A=d['real']['A_M3V_comb'], best=d['real']['comb_best'])
    npz = f'{ROOT}/out/inj_{tag}.npz'
    if os.path.exists(npz):
        Z = np.load(npz)
        pt, pm, lr = {}, {}, {}
        x3 = Z['M3V__comb_c'] - Z['M3V__comb_c'][:, [0]]
        for ii, h in enumerate(MAIN[1:], 1):
            xh = Z[f'{h}__comb_c'] - Z[f'{h}__comb_c'][:, [0]]
            xr = r['real'][h]
            pt[h] = float(np.mean(xh[:, ii] >= xr)); pm[h] = float(np.mean(x3[:, ii] <= xr))
            l3 = norm.logpdf(xr, x3[:, ii].mean(), x3[:, ii].std()); lh = norm.logpdf(xr, xh[:, ii].mean(), xh[:, ii].std())
            lr[h] = float((l3 - lh) / np.log(10))
        r.update(p_true=pt, p_M3V_le=pm, log10_LR_M3V_over_h=lr, n_inj=int(len(x3)),
                 free_M3V_inj_sd=d['inj']['M3V']['free_pos_sd'], free_M3V_inj_mean=d['inj']['M3V']['free_pos_mean'],
                 loc_thresholds=d['inj']['M3V']['true_minus_free_q'])
        if 'ring' in d: r['ring'] = d['ring']; r['ring_thresholds'] = d['ring_thresholds']
    if tag == 'base':
        r['per_sector'] = {s: dict(dchi2=v['dchi2'], best=v['best'], free=v['free_pos'], A=v['amp']['M3V'], A_err=v['amp_err_M3V'],
                                   chi2=v['chi2_M3V'], dof=v['dof']) for s, v in d['real']['per_sector'].items()}
    S[tag] = r
jdump(S, f'{ROOT}/out/summary.json')
H = ['pairA', 'pairB', 'Mira', 'G175', 'G177', 'G170']
print(f"{'run':16s} " + ' '.join(f'{h:>7s}' for h in H) + '   free(E,N)   dchi2free   A    | p_true(pairA,Mira,G175)  log10LR(G175)')
for tag, r in S.items():
    s = f"{tag:16s} " + ' '.join(f"{r['real'][h]:7.1f}" for h in H) + f"  ({r['free'][0]:+.1f},{r['free'][1]:+.1f})  {r['free_dchi2']:+6.2f}  {r['A']:5.2f}"
    if 'p_true' in r:
        s += f" | {r['p_true']['pairA']:.3f} {r['p_true']['Mira']:.3f} {r['p_true']['G175']:.3f}  {r['log10_LR_M3V_over_h']['G175']:.1f}  (n={r['n_inj']})"
    print(s)
