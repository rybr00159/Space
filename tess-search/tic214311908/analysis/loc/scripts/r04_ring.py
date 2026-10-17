"""Uncatalogued-source test (item 5): injections of the observed signal at offsets r = 1-8" from the M3V (8 position
angles, M-dwarf colour).  For each r: fraction of injections in which the M3V position is excluded from the
combined localisation region at 95.4% / 99.73% (combined delta-chi2 at the M3V above the thresholds calibrated with
injections ON the M3V), and the coverage check of those thresholds (delta-chi2 at the true ring position).
Usage: r04_ring.py METHOD(diff|amp)
Output: out/ring_{method}.json"""
import sys, json
sys.path.insert(0, 'tic214311908/analysis/loc/scripts')
import numpy as np
from l0lib import BASE, to_json
m = sys.argv[1]
R = np.load(f'{BASE}/out/sim_{m}_ring.npz', allow_pickle=True)
P = json.load(open(f'{BASE}/out/post_{m}_prf_all.json'))
T = P['loc_thresholds_used']
out = dict(thresholds=T, per_r={})
allt = []
for r in (1.0, 2.0, 3.0, 4.0, 6.0, 8.0):
    dm, dt, fe = [], [], {}
    for pa in range(0, 360, 45):
        a = R[f'ring_{r}_{pa}__dchi_m3v']; b = R[f'ring_{r}_{pa}__dchi_true']
        dm.append(a); dt.append(b); fe[pa] = float((a > T[2]).mean())
    dm = np.concatenate(dm); dt = np.concatenate(dt); allt.append(dt)
    out['per_r'][str(r)] = dict(frac_M3V_excluded_95=float((dm > T[1]).mean()), frac_M3V_excluded_99p7=float((dm > T[2]).mean()),
                                median_dchi_M3V=float(np.median(dm)), coverage_true_in_95=float((dt <= T[1]).mean()),
                                coverage_true_in_99p7=float((dt <= T[2]).mean()), frac_excl_99p7_by_PA=fe)
    v = out['per_r'][str(r)]
    print(f"{m} r={r:3.0f}\": M3V excluded at 95.4% in {v['frac_M3V_excluded_95']*100:5.1f}%, at 99.73% in {v['frac_M3V_excluded_99p7']*100:5.1f}% "
          f"(median delta-chi2 at M3V {v['median_dchi_M3V']:.1f}); coverage of true position 95%: {v['coverage_true_in_95']:.3f}, 99.7%: {v['coverage_true_in_99p7']:.3f}; "
          f"by PA (99.7%): " + ' '.join(f'{k}:{x:.2f}' for k, x in fe.items()))
allt = np.concatenate(allt)
out['ring_threshold_quantiles'] = np.percentile(allt, [68.27, 95.45, 99.73]).tolist()
rs = np.array([1, 2, 3, 4, 6, 8.]); f95 = np.array([out['per_r'][str(r)]['frac_M3V_excluded_95'] for r in rs]); f99 = np.array([out['per_r'][str(r)]['frac_M3V_excluded_99p7'] for r in rs])
out['r50_95'] = float(np.interp(0.5, f95, rs)) if f95.max() >= 0.5 else None
out['r50_99p7'] = float(np.interp(0.5, f99, rs)) if f99.max() >= 0.5 else None
out['r90_99p7'] = float(np.interp(0.9, f99, rs)) if f99.max() >= 0.9 else None
print(f"{m}: offset at which the M3V is excluded in 50% of trials: {out['r50_95']}\" (95.4%), {out['r50_99p7']}\" (99.73%); 90% at 99.73%: {out['r90_99p7']}\"; "
      f"ring-calibrated thresholds {np.round(out['ring_threshold_quantiles'], 2).tolist()} vs M3V-calibrated {np.round(T, 2).tolist()}")
to_json(out, f'{BASE}/out/ring_{m}.json')
