"""s17: collects the scenario results into one table (out/scenario_table.csv, out/summary.json) and prints it as
markdown for the README. dchi2 values are relative to the free planet fit of the same light curve; 'total' adds the
non-light-curve terms each scenario must also satisfy (phase-0.5 secondary, stellar priors, the background star's
photometry and parallax relative to a single-star fit of the same photometry, s12b, and for the WD scenario the
Gaia astrometric likelihood -2 ln L, s16)."""
import json, glob, os, numpy as np, pandas as pd
SD = 'tic214311908/analysis/shape'
J = lambda f: json.load(open(f'{SD}/out/{f}'))
rows = []
plan = {ds: J(f'fit_{ds}_free.json')['chi2_best'] for ds in ['inv', 'ref']}
SEC = {'inv': (-0.16, 0.79), 'ref': (0.29, 0.77)}
bgs = J('bg_single.json')
for f in sorted(glob.glob(f'{SD}/out/bg_eb_*.json')):
    d = json.load(open(f)); ds = d['dataset']
    null = plan[ds] + (SEC[ds][0]/SEC[ds][1])**2 + (0 if d['lconly'] else bgs[d['source']]['total'])
    rows.append(dict(group='background EB', scenario=os.path.basename(f)[6:-5].replace(ds + '_', ''), dataset=ds,
                     dchi2_lc=d['dchi2_lc_vs_planet'], dchi2_total=d['total'] - null,
                     detail=f"A {d['A']['M']:.2f} Msun, {d['A']['R']:.3f} Rsun, {d['A']['T']:.0f} K; B {d['B']['M']:.2f} Msun, {d['B']['R']:.2f} Rsun; T14 {d['T14_h']:.2f} h; sec {d['sec_ppt']:.2f} ppt; phot chi2 {d['chi2_phot']:.1f}, plx {d['chi2_plx']:.1f}",
                     file=os.path.basename(f)))
for f in sorted(glob.glob(f'{SD}/out/compact_*.json')):
    d = json.load(open(f)); ds = d['dataset']
    rows.append(dict(group='compact star occulted (bound)', scenario=os.path.basename(f)[8:-5].replace(ds + '_', ''), dataset=ds, dchi2_lc=d['chi2_lc'] - plan[ds], dchi2_total=d['dchi2_vs_planet_free'],
                     detail=f"M_A {d['M_A']:.2f}, R_A {d['R_A']:.4f}, R_B {d['R_B']:.3f}, b {d['b_occ']:.2f}, T14 {d['T14_h']:.2f} h, ingress {d['ingress_min']:.1f} min, sec {d['secondary_ppt']:.2f} ppt, e {d['e']:.2f}",
                     file=os.path.basename(f)))
for ds in ['inv', 'ref']:
    fn = f'{SD}/out/bound_eb_{ds}.json'
    if os.path.exists(fn):
        d = json.load(open(fn))
        for lab, v in d['best_by_class'].items():
            if lab in ('P5.69', 'P11.38', 'occulter_BD', 'no_third_light'):
                rows.append(dict(group='bound hierarchical EB', scenario=f'best {lab}', dataset=ds, dchi2_lc=v['chi2_lc'] - plan[ds], dchi2_total=v['dchi2_vs_planet_free'],
                                 detail=f"m1 {v['m1']}, {v['typ2']} m2 {v['m2']}, k {v['k']:.2f}, m0 {v['m0']}, f_src {v['f_src']:.3f}, P {v['P_EB']:.2f}, b {v['b']:.2f}, T14 {v['T14_h']:.2f} h",
                                 file=os.path.basename(fn)))
D = pd.DataFrame(rows)
D.to_csv(f'{SD}/out/scenario_table.csv', index=False)
# WD mass profile combined with astrometry
wd = {}
fn = f'{SD}/out/wd_mass_profile_ref_cool.json'
if os.path.exists(fn):
    a = J('wd_checks.json')['astrometry']
    masses = np.array([0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0, 1.2])
    for sig in ('0.22', '0.35'):
        Ls = np.array([a[f'M{m}_sig{sig}']['L'] for m in masses])
        prof = []
        for r in J('wd_mass_profile_ref_cool.json')['rows']:
            L = float(np.exp(np.interp(r['M_WD'], masses, np.log(Ls))))
            prof.append(dict(M=r['M_WD'], dchi2_lc_sec=r['dchi2_vs_planet_free'], m2lnL_astrom=-2*np.log(L), total=r['dchi2_vs_planet_free'] - 2*np.log(L),
                             Teff_needed=r['Teff_WD_needed'], R=r['R_WD'], sec_ppt=r['secondary_ppt']))
        wd[f'sigma_CCD_{sig}'] = dict(profile=prof, min_total=min(p['total'] for p in prof), M_at_min=min(prof, key=lambda p: p['total'])['M'])
json.dump(dict(scenarios=rows, wd_with_astrometry=wd), open(f'{SD}/out/summary.json', 'w'), indent=1, default=float)
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 200)
print(D[['group', 'scenario', 'dataset', 'dchi2_lc', 'dchi2_total']].round(1).to_string())
for k, v in wd.items():
    print(k, 'min total', round(v['min_total'], 1), 'at M', v['M_at_min'], [(p['M'], round(p['total'], 1)) for p in v['profile']])
