"""PSF-model mismatch tests: injections made with a PSF different from the fitting model ('alt' PRF of another CCD
location, 'blur' 0.15 px smear, 'emp' PRF + red field star residual pattern, 'miraemp' the Mira's empirical image) and
injections with each hypothesis's own best-fit amplitude ('ownamp'), compared with the baseline (same PRF).
For each true source: quantiles of combined chi2(pair)-chi2(M3V), chi2(Mira)-chi2(M3V), chi2(..260544)-chi2(M3V),
combined free-fit position, and the fraction of trials reaching the real-data values.
Usage: r05_mismatch.py METHOD(diff|amp)
Output: out/mismatch_{method}.json"""
import sys, json, os
sys.path.insert(0, 'tic214311908/analysis/loc/scripts')
import numpy as np
from l0lib import BASE, to_json
m = sys.argv[1]
base = np.load(f'{BASE}/out/sim_{m}_prf_all.npz', allow_pickle=True)
names = list(base['hyp_name']); idx = {n: i for i, n in enumerate(names)}
rc = base['real__chiH'].sum(0)
real = dict(pair=float(rc[idx['pairA']] - rc[idx['M3V']]), mira=float(rc[idx['Mira']] - rc[idx['M3V']]), f260544=float(rc[idx['5925373272827260544']] - rc[idx['M3V']]))
out = dict(real=real)
runs = [('baseline', f'{m}_prf_all')] + [(k, f'{m}_{k}_main') for k in ('alt', 'blur', 'emp')] + [('miraemp', f'{m}_miraemp')]
if m == 'diff': runs.append(('ownamp', 'diff_prf_main_ownamp'))
for lab, tag in runs:
    fn = f'{BASE}/out/sim_{tag}.npz'
    if not os.path.exists(fn): continue
    Z = np.load(fn, allow_pickle=True)
    for t in ('M3V', 'pairA', 'pairB', 'Mira', '5925373272827260544', '5925373272827262208'):
        k = f'inj_{t}__chiH'
        if k not in Z.files: continue
        c = Z[k].astype(float).sum(1)
        dp = c[:, idx['pairA']] - c[:, idx['M3V']]; dm = c[:, idx['Mira']] - c[:, idx['M3V']]; df = c[:, idx['5925373272827260544']] - c[:, idx['M3V']]
        own = c[:, idx[t]] - c[:, idx['M3V']]
        pos = Z[f'inj_{t}__gcomb_pos']
        r = dict(pair_minus_M3V_q=np.percentile(dp, [2.5, 50, 97.5]).tolist(), mira_minus_M3V_q=np.percentile(dm, [2.5, 50, 97.5]).tolist(),
                 f260544_minus_M3V_q=np.percentile(df, [2.5, 50, 97.5]).tolist(), true_minus_M3V_q=np.percentile(own, [2.5, 50, 97.5, 99.85]).tolist(),
                 free_pos_median=np.median(pos, 0).tolist(), M3V_best_frac=float((np.argmin(c, 1) == idx['M3V']).mean()),
                 true_best_frac=float((np.argmin(c, 1) == idx[t]).mean()))
        if t != 'M3V':
            rv = real['pair'] if t.startswith('pair') else real['mira'] if t == 'Mira' else real['f260544'] if t.endswith('260544') else float(rc[idx[t]] - rc[idx['M3V']])
            r['frac_reaching_real'] = float((own >= rv).mean()); r['real'] = rv
        out[f'{lab}:{t}'] = r
        print(f"{m} {lab:8s} true {t[-10:]:>10s}: chi2(true)-chi2(M3V) q2.5/50/97.5/99.85 {np.round(r['true_minus_M3V_q'], 1).tolist()}"
              + (f" (real {r['real']:.1f}; reached in {r['frac_reaching_real']*100:.1f}%)" if 'real' in r else '') +
              f"; pair-M3V median {r['pair_minus_M3V_q'][1]:.1f}, Mira-M3V {r['mira_minus_M3V_q'][1]:.1f}, 260544-M3V {r['f260544_minus_M3V_q'][1]:.1f}; "
              f"free fit {np.round(r['free_pos_median'], 1).tolist()}; M3V best {r['M3V_best_frac']*100:.0f}%, true best {r['true_best_frac']*100:.0f}%")
to_json(out, f'{BASE}/out/mismatch_{m}.json')
