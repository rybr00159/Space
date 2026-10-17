"""r15: collects the referee-2 key numbers into out/summary.json, incl. a 2-D (dWP, dVP) Gaussian likelihood ratio
planet/WD from the r07 shift injections (sample means and covariances of each truth)."""
import json, numpy as np, pandas as pd
RD = 'tic214311908/analysis/referee2_shape'
J = lambda f: json.load(open(f'{RD}/out/{f}'))
S = {}
for ds in ('inv', 'ref'):
    f = J(f'mcmc_{ds}_free_fake.json')['post']
    S[f'free_{ds}'] = {k: [round(f[k][3], 5), round(f[k][4] - f[k][3], 5), round(f[k][3] - f[k][2], 5)] for k in ('k', 'b', 'rho', 'T14_h', 'flat', 'ingress_min')}
    S[f'free_{ds}']['frac_grazing'] = J(f'mcmc_{ds}_free_fake.json')['frac_grazing']
    sc = J(f'scen_{ds}_fake.json')
    S[f'dchi2_{ds}'] = {k: round(sc[k]['dchi2_vs_planet_free'], 2) for k in ('planet_single', 'planet_twin', 'vshape', 'wd_cool', 'wd_hot', 'wd_twin')}
    try:
        R = pd.read_csv(f'{RD}/out/inject_{ds}_shift.csv'); real = J(f'inject_{ds}_shift.json')['real']
        x = np.array([real['dWP'], real['dVP']])
        L = {}
        for tr in ('planet', 'wd'):
            q = R[R.truth == tr][['dWP', 'dVP']].values; mu = q.mean(0); C = np.cov(q.T)
            d = x - mu; L[tr] = float(np.exp(-0.5*d @ np.linalg.solve(C, d))/np.sqrt(np.linalg.det(2*np.pi*C)))
        S[f'inject_{ds}'] = dict(real=real, LR_planet_over_wd_2d=L['planet']/L['wd'], LR_1d=J(f'inject_{ds}_shift.json').get('LR_planet_over_wd_dWP_gauss'),
                                 planet=J(f'inject_{ds}_shift.json')['planet'], wd=J(f'inject_{ds}_shift.json')['wd'])
    except FileNotFoundError:
        pass
    try: S[f'inject_vshape_{ds}'] = J(f'inject_{ds}_vshape.json')
    except FileNotFoundError: pass
    S[f'dilution_{ds}'] = [(r['f_src'], round(r['dT'], 2), round(r['transit']['dchi2'], 1), round(r.get('occultation', {}).get('dchi2', np.nan), 1)) for r in J(f'dilution_{ds}.json')['rows']]
S['derived'] = J('derived.json')
S['secondary'] = J('secondary.json')['phase05']
S['wd_combined_min'] = {ds: J('wd_combined.json')[ds]['min_own'] for ds in ('inv', 'ref')}
S['wd_combined_min_ref_r18_secondary'] = {ds: J('wd_combined.json')[ds]['min_ref_r18'] for ds in ('inv', 'ref')}
S['kprofile_tb'] = J('kprofile_tb.json')
S['per_sector'] = {ds: {k: v for k, v in J('per_sector.json')[ds].items() if k != 'duration_scale'} for ds in ('inv', 'ref')}
S['xp'] = J('xp_test.json')['consistent']
S['rv_smear'] = J('rv_smear.json')['m2lnL_wd_vs_planet']
S['ruwe_sigma_ccd'] = J('sigma_ccd_check.json')
json.dump(S, open(f'{RD}/out/summary.json', 'w'), indent=1, default=float)
print(json.dumps({k: v for k, v in S.items() if k.startswith(('inject', 'dchi2', 'wd_comb'))}, indent=1, default=float)[:6000])
