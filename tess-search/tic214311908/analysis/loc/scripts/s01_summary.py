"""Collect the key numbers of all steps into out/summary.json and print compact tables (for the README)."""
import sys, json, os
sys.path.insert(0, 'tic214311908/analysis/loc/scripts')
import numpy as np
from l0lib import BASE, to_json, SECTORS
J = lambda f: json.load(open(f'{BASE}/out/{f}')) if os.path.exists(f'{BASE}/out/{f}') else None
S = {}
w = J('wcs_cal.json'); mp = J('mira_pos.json')
S['wcs'] = {s: dict(offset_px=w[s]['offset_px'], err=w[s]['offset_err_px'], cal_median=w[s]['offset_px_calmedian']) for s in w}
S['mira_apparent'] = {s: dict(dEN=mp[s]['base']['dEN'], err=mp[s]['base']['err'], radial=mp[s]['base']['rad'], tangential=mp[s]['base']['tan'],
                              flux=mp[s]['base']['flux'], resid_rms_pct=100 * mp[s]['base']['resid_rms_frac'], predicted_colour_shift=mp[s]['predicted_colour_shift_rad'],
                              variant_range_E=mp[s]['variant_range_E'], variant_range_N=mp[s]['variant_range_N']) for s in mp}
for m in ('diff', 'amp'):
    P = J(f'post_{m}_prf_all.json')
    if P is None: continue
    S[f'{m}_dchi2_real'] = P['dchi2_real']; S[f'{m}_kappa'] = P['kappa']
    S[f'{m}_posterior_model'] = {k: v for k, v in P['posterior_real'].items() if not k.endswith('_top')}
    S[f'{m}_posterior_empirical'] = P['posterior_empirical_position']
    S[f'{m}_pvalues'] = P['pvalues_chi2_star_minus_M3V']; S[f'{m}_localisation'] = P['localisation']; S[f'{m}_loc_thr'] = P['loc_thresholds']
    S[f'{m}_uncatalogued'] = P['uncatalogued']; S[f'{m}_free_fit_per_sector'] = P['real_free_fit_offset_from_M3V_arcsec']
    S[f'{m}_confusion'] = P['injection_confusion_combined_k1']
    for f in (f'ring_{m}.json', f'mismatch_{m}.json'):
        x = J(f)
        if x: S[f.replace('.json', '')] = x
for f in ('joint_fits.json', 'variants.json', 'mira_emp_fit.json', 'shape_profile.json', 'colourmap.json', 'cubecheck.json'):
    x = J(f)
    if x: S[f.replace('.json', '')] = x
ph = J('physics.json')
if ph:
    S['physics_main'] = [r for r in ph['stars'] if r['name'] in ('M3V', 'pairA', 'pairB', 'Mira') or r['dist_M3V'] < 8 or r['sep_mira'] < 9]
    S['physics_counts'] = dict(n=len(ph['stars']), depth_ge_1=sum(r['depth_ref'] >= 1 for r in ph['stars']),
                               feasible=[r['name'] for r in ph['stars'] if r['dark_occultor_ok']])
to_json(S, f'{BASE}/out/summary.json')
# compact print
for m in ('diff', 'amp'):
    if f'{m}_dchi2_real' not in S: continue
    print(f'== {m.upper()}  chi2 - chi2(best catalogued) per sector')
    for k, v in S[f'{m}_dchi2_real'].items():
        print(f"  {k:>8s}: M3V {v['M3V']:6.1f} pairA {v['pairA']:6.1f} pairB {v['pairB']:6.1f} Mira {v['Mira']:6.1f} G17.5(3.2\") {v['5925373272827260544']:6.1f} G17.7(5.7\") {v['5925373272827262208']:6.1f}  best {v['best'][-10:]}")
    print('  posterior (kappa-calibrated, equal priors):')
    for k, v in S[f'{m}_posterior_model'].items():
        if 'kcal' in k: print(f"   {k:14s} " + ', '.join(f'{a}: {b:.3g}' for a, b in v.items()))
    print('  empirical position-likelihood posterior:')
    for k, v in S[f'{m}_posterior_empirical'].items():
        if not k.endswith('_top'): print(f"   {k:14s} " + ', '.join(f'{a}: {b:.3g}' for a, b in v.items()))
