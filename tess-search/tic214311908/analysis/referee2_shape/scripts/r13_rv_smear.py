"""r13: Gaia DR3 combined-CCF radial-velocity error under the white-dwarf scenario (semi-quantitative).
For grvs_mag > 12 Gaia DR3 derives the RV and its error from the sum of the per-transit CCFs (rv_method_used 2;
Katz+2023), the error following Zucker (2003): sigma^2 = -[N C''/C * C^2/(1 - C^2)]^-1 at the peak. A K1 orbit
sampled at the 21 RVS transits smears the combined CCF. For per-transit CCFs that are Gaussians of sigma w (12 and
16 km/s: RVS R = 11500 instrumental + M-dwarf lines, template-matched) with peak C0 (0.3 and 0.6), the ratio
sigma(K1)/sigma(0) is computed for random orbital phases (1000 draws, median and 5th percentile).
Output: out/rv_smear.json"""
import json, numpy as np
rng = np.random.default_rng(8)
v = np.linspace(-250, 250, 5001); dv = v[1] - v[0]
def zucker(C):
    i = np.argmax(C); i = min(max(i, 2), len(v) - 3)
    c = C[i]; c2 = (C[i+1] - 2*C[i] + C[i-1])/dv**2
    return np.sqrt(-1/(c2/c*c**2/(1 - c**2)))
out = {}
for w in (12., 16.):
    for C0 in (0.3, 0.6):
        base = zucker(C0*np.exp(-0.5*(v/w)**2))
        row = {}
        for M2, K in ((0.2, 34.7), (0.3, 46.7), (0.45, 61.3), (0.6, 73.0), (0.8, 85.9), (1.0, 96.7)):
            r = []
            for j in range(1000):
                ph = rng.uniform(0, 2*np.pi, 21); vi = K*np.sin(ph)
                C = C0*np.mean(np.exp(-0.5*((v[:, None] - vi[None, :])/w)**2), axis=1)
                r.append(zucker(C)/base)
            r = np.array(r)
            row[str(M2)] = dict(K1=K, ratio_median=float(np.median(r)), ratio_p05=float(np.percentile(r, 5)), pred_err_median_kms=float(2.217*np.median(r)), pred_err_p05_kms=float(2.217*np.percentile(r, 5)))
        out[f'w{int(w)}_C{C0}'] = row
        print(w, C0, {k: (round(x['ratio_p05'], 1), round(x['ratio_median'], 1)) for k, x in row.items()})
out['comparison'] = dict(target_rv_error=2.217, comp_median=2.47, comp_p90=4.24, comp_p99=7.75, note='comparison: 1642 M dwarfs of the same grvs_mag (../shape/out/wd_checks.json)')
json.dump(out, open('out/rv_smear.json', 'w'), indent=1)
# ---- likelihood ratio using the empirical RV-error distribution of comparison M dwarfs of the same grvs_mag
import pandas as pd
from scipy.stats import gaussian_kde
m = pd.read_csv('tic214311908/analysis/referee_star/data/mdw_cmd.csv')
sel = m[(m.grvs_mag > 12.3) & (m.grvs_mag < 12.8) & m.radial_velocity_error.notna()].radial_velocity_error.values
le = np.log(sel); kde = gaussian_kde(le)
obs = 2.217
out['empirical'] = dict(n=int(len(sel)), pct=[float(x) for x in np.percentile(sel, [0.5, 1, 5, 16, 50, 84, 99])], frac_below_obs=float(np.mean(sel < obs)))
out['m2lnL_wd_vs_planet'] = {}
for key in ('w16_C0.3', 'w12_C0.6'):
    row = {}
    for M2, x in out[key].items():
        # under the WD hypothesis log(err_obs) = log(err_intrinsic) + log(ratio), ratio drawn from the smearing simulation (approx. by its median, and by the 5th percentile = most favourable)
        res = {}
        for lab, rr in (('median', x['ratio_median']), ('p05', x['ratio_p05'])):
            Lw = float(kde(np.log(obs/rr))[0]); Lp = float(kde(np.log(obs))[0])
            res[lab] = dict(frac_intrinsic_below=float(np.mean(sel < obs/rr)), m2lnL=float(-2*np.log(max(Lw, 1e-300)/Lp)))
        row[M2] = res
    out['m2lnL_wd_vs_planet'][key] = row
    print(key, {k: (round(v['p05']['m2lnL'], 1), round(v['median']['m2lnL'], 1), v['p05']['frac_intrinsic_below']) for k, v in row.items()})
print(out['empirical'])
json.dump(out, open('out/rv_smear.json', 'w'), indent=1)
