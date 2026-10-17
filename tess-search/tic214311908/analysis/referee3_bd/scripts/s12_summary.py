"""s12: combine the referee-3 results into the final channel estimates.
(a) FPP / PTP / BD rows of s02 (sign-fixed r05, f_BD 1.3e-4, BD e laws) for the adopted and variant settings;
(b) Gaia RVS factor (s05, median and range of the mass-averaged LR over CCF width/peak/e/sampling, empirical BD masses)
    applied to the BD-on-M3V and BD-on-primary terms;
(c) shape-noise correction: ratio of my calibrated grazing-BD likelihood ratio on the cleaned curves (s11) to the
    exp(-dchi2/2) carried by the referee2_fpp r03 grids (r11: 0.05 REF, 0.11 INV) - applied as an alternative;
(d) my independent evidence-ratio estimate of BD-on-M3V (s06) with f_BD/f_N;
(e) P(e > 0.5) for the BD e laws.
Output: out/s12_summary.json"""
import sys; sys.path.insert(0, 'tic214311908/analysis/referee3_bd/scripts')
from rlib3 import *
import pandas as pd, os
rng = np.random.default_rng(1)
S2 = pd.read_csv(f'{R3}/out/s02_fpp_fixed.csv')
G5 = json.load(open(f'{R3}/out/s05_gaia_rv.json'))
rvf = G5['mass_averaged_LR_range']['empirical_Pl16']      # [min, median, max]
rvu = G5['mass_averaged_LR_range']['uniform']
S11 = json.load(open(f'{R3}/out/s11_calibration.json')) if os.path.exists(f'{R3}/out/s11_calibration.json') else {}
S6 = json.load(open(f'{R3}/out/s06_evidence.json'))['results'] if os.path.exists(f'{R3}/out/s06_evidence.json') else {}
out = {'rv_factor_empirical_masses_min_med_max': rvf, 'rv_factor_uniform_masses_min_med_max': rvu}
def row(tag):
    r = S2[S2.tag == tag]
    return r.iloc[0].to_dict() if len(r) else None
LRgrid = {'ref': 0.05, 'inv': 0.11}
tabs = []
for ds in ('ref_beta', 'inv_beta', 'ref_white', 'inv_white'):
    for pe in ('kipping', 'rayleigh19', 'vaneylen', 'circ'):
        for be in ('bdsp', 'bdlp'):
            r = row(f'{ds} planet-e {pe} BD-e {be} f_BD {1.3e-4:g}')
            if r is None: continue
            bdm = r['BD_on_M3V'] + r['BD_on_primary']
            d = dict(ds=ds, planet_e=pe, bd_e=be, FPP=r['FPP'], FPP_system=r['FPP_system'], P_planet_any_star=r['P_planet_any_star'], PTP=r['PTP'],
                     NTP=r['NTP_pairA'], BD_total_noRV=r['BD_total'], BD_total_RVmed=r['BD_total'] - bdm*(1 - rvf[1]),
                     BD_total_RVmin=r['BD_total'] - bdm*(1 - rvf[0]))
            k = ds.split('_')[0]
            if k in S11 and 'LR_calibrated_BD_over_planet' in S11[k] and ds.endswith('beta'):
                corr = S11[k]['LR_calibrated_BD_over_planet']/LRgrid[k]
                d['shape_corr_factor'] = corr
                d['BD_total_RVmed_shapecal'] = d['BD_total_RVmed']*min(corr, 1e3)
            tabs.append(d)
out['table'] = tabs
for d in tabs:
    print({k: (f'{v:.2e}' if isinstance(v, float) else v) for k, v in d.items()})
# e tails
E = {'bdsp': (0.426, 0.092, 2.354, 0.73), 'bdlp': (1.879, 0.92, 2.470, 1.25)}
out['P_e_gt_0p5'] = {}
for k, (a, sa, b, sb) in E.items():
    A = np.clip(rng.normal(a, sa, 400000), 0.05, None); Bb = np.clip(rng.normal(b, sb, 400000), 0.1, None)
    e = rng.beta(A, Bb); out['P_e_gt_0p5'][k] = float(np.mean(e > 0.5)); out.setdefault('P_e_gt_0p4', {})[k] = float(np.mean(e > 0.4))
out['P_e_gt_0p5']['kipping'] = float(np.mean(rng.beta(0.867, 3.03, 400000) > 0.5))
# independent evidence (s06): P(BD on M3V)/P(TP on M3V) = f_BD/f_N * ratio, f_N 0.003
ind = {}
for key, R in S6.items():
    for fbd in (1.3e-4, 1e-3):
        for be in ('bdsp', 'bdlp', 'kipping'):
            rr = R['ratio_BD_over_TP'][f'BD e {be} / TP e kipping']
            ind[f'{key}|f_BD {fbd:g}|BD e {be}'] = dict(ratio=rr, post_BD_on_M3V_approx=fbd/0.003*rr*0.96)
out['independent_BD_on_M3V'] = ind
for k, v in ind.items(): print(k, f"{v['ratio']:.2e}", f"{v['post_BD_on_M3V_approx']:.2e}")
jdump(dict(method=__doc__, **out), f'{R3}/out/s12_summary.json')
print(out['P_e_gt_0p5'], out['P_e_gt_0p4'])
