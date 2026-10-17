"""s11: calibrated shape likelihood ratio for the grazing BD on the M3V from the s04 fake-epoch injections.
Lambda_obs from s03 (cleaned curves, beta noise). Under BD truth: empirical exceedance P(Lambda >= Lambda_obs) and a
maximum-likelihood fit of Lambda ~ s * chi2_nu (s, nu free; zeros floored at 1e-4) for the tail and the density at
Lambda_obs. Under planet truth: Gaussian KDE density at Lambda_obs and the empirical quantile of Lambda_obs.
LR_cal = f_BD(Lambda_obs)/f_planet(Lambda_obs), compared with the formal exp(-Lambda_obs/2) and with the
best-grazing-BD dchi2 that the referee2_fpp r03 grids carry (r11: 5.8 REF beta, 4.3 INV beta).
Output: out/s11_calibration.json"""
import sys; sys.path.insert(0, 'tic214311908/analysis/referee3_bd/scripts')
from rlib3 import *
import pandas as pd, os
from scipy import stats
from scipy.optimize import minimize
F = json.load(open(f'{R3}/out/s03_fit_real.json'))['fits']
r11 = json.load(open(f'{B}/referee2_fpp/out/r11_bd_noise.json'))
out = {}
for ds in ('ref', 'inv'):
    Lobs = F[f'{ds}_beta']['Lambda_grazBD']
    o = dict(Lambda_obs=Lobs, formal_LR=float(np.exp(-Lobs/2)))
    fb, fp = f'{R3}/out/s04_inject_{ds}_grazBD.csv', f'{R3}/out/s04_inject_{ds}_planet.csv'
    if os.path.exists(fb):
        Lb = np.clip(pd.read_csv(fb).Lambda_G.values, 1e-4, None)
        o['BD_truth'] = dict(n=len(Lb), mean=float(Lb.mean()), p50=float(np.median(Lb)), p95=float(np.percentile(Lb, 95)), max=float(Lb.max()),
                             n_ge_obs=int(np.sum(Lb >= Lobs)), frac_ge_4=float(np.mean(Lb >= 4)))
        nll = lambda p: -np.sum(stats.chi2.logpdf(Lb/np.exp(p[0]), np.exp(p[1])) - p[0])
        r = minimize(nll, [0., 0.], method='Nelder-Mead')
        s, nu = float(np.exp(r.x[0])), float(np.exp(r.x[1]))
        o['BD_truth'].update(scaled_chi2_s=s, scaled_chi2_nu=nu, tail_P_ge_obs=float(stats.chi2.sf(Lobs/s, nu)),
                             density_at_obs=float(stats.chi2.pdf(Lobs/s, nu)/s),
                             chi2_1_tail=float(stats.chi2.sf(Lobs, 1)), chi2_1_density=float(stats.chi2.pdf(Lobs, 1)))
    if os.path.exists(fp):
        Lp = pd.read_csv(fp).Lambda_G.values
        kde = stats.gaussian_kde(Lp)
        o['planet_truth'] = dict(n=len(Lp), mean=float(Lp.mean()), sd=float(Lp.std()), p16=float(np.percentile(Lp, 16)), p50=float(np.median(Lp)),
                                 p84=float(np.percentile(Lp, 84)), quantile_of_obs=float(np.mean(Lp < Lobs)), density_at_obs=float(kde(Lobs)[0]),
                                 frac_lt_4=float(np.mean(Lp < 4)))
    if 'BD_truth' in o and 'planet_truth' in o:
        o['LR_calibrated_BD_over_planet'] = o['BD_truth']['density_at_obs']/o['planet_truth']['density_at_obs']
        o['LR_formal_over_calibrated'] = o['formal_LR']/o['LR_calibrated_BD_over_planet']
    out[ds] = o
    print(ds, json.dumps(o, default=float), flush=True)
out['r11_grid_dchi2'] = {k: v for k, v in r11.items() if 'dchi2' in k.lower() or 'best' in k.lower()} if isinstance(r11, dict) else None
jdump(dict(method=__doc__, **out), f'{R3}/out/s11_calibration.json')
