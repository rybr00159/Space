"""s03: own fits of the CLEANED REF and INV light curves (../clean) for: planet (D = 1, k 0.02-0.35, b free incl. grazing),
grazing BD on the M3V (D = 1, k >= 0.19 [bd's bound] and k >= 0.20 [R_BD >= 0.745 R_J for R* 0.383]), and a large body
on a diluted host (k 0.25-1.0, D 0.01-0.6). Lambda = chi2(scenario) - chi2(planet). Profile of the grazing-BD chi2 in
rho_c. Noise: flux_err x beta_transit (clean's fake-epoch calibration) and white (flux_err only).
Output: out/s03_fit_real.json"""
import sys; sys.path.insert(0, 'tic214311908/analysis/referee3_bd/scripts')
from lc3 import *
import pickle, time
t0 = time.time()
bank = TemplateBank(np.r_[np.linspace(0.09, 0.15, 7), [0.19, 0.20, 0.22, 0.25, 0.28, 0.32, 0.40, 0.50, 0.65]], nb=30)
print('bank', len(bank.items), round(time.time() - t0, 1), flush=True)
res = {}
for ds in ('ref', 'inv'):
    for nm in ('beta', 'white'):
        lc = load_ref(beta=(nm == 'beta')) if ds == 'ref' else load_inv(beta=(nm == 'beta'))
        R = {'N': len(lc.t), 'n_windows': int(lc.nw), 'chi2_null': float(np.sum(lc.w*lc.yp**2))}
        for sc in ('planet', 'grazBD', 'grazBD20', 'comp'):
            R[sc] = fit_scen(lc, bank, sc, nstart=4)
        for sc in ('grazBD', 'grazBD20', 'comp'):
            R[f'Lambda_{sc}'] = R[sc]['chi2'] - R['planet']['chi2']
        # profile in rho_c for the grazing BD (k >= 0.19)
        prof = []
        for rho in (0.8, 1.2, 1.5, 1.8, 2.2, 2.7, 3.5, 5.0):
            best = None
            for k0, bb0 in ((0.19, 0.83), (0.21, 0.82), (0.25, 0.80), (0.30, 0.78), (0.19, 0.80)):
                f = lambda x: chi2_par(lc, np.r_[x[0], x[1], x[2], np.log(rho), x[3], x[4]], 'grazBD')
                r = minimize(f, [R['grazBD']['x'][0], np.log(k0), bb0, U1, U2], method='Nelder-Mead', options=dict(maxfev=1500, xatol=1e-6, fatol=1e-4, adaptive=True))
                if best is None or r.fun < best.fun: best = r
            prof.append(dict(rho_c=rho, dchi2_vs_planet=float(best.fun - R['planet']['chi2']), k=float(np.exp(best.x[1])), b=float(best.x[2]*(1 + np.exp(best.x[1])))))
        R['profile_rho_grazBD'] = prof
        res[f'{ds}_{nm}'] = R
        print(ds, nm, 'N', R['N'], 'nw', R['n_windows'], {sc: (round(R[sc]['chi2'], 2), round(R[sc]['k'], 3), round(R[sc]['b'], 3), round(R[sc]['rho_c'], 2), round(R[sc]['D'], 3), round(R[sc]['T14_h'], 2), round(R[sc]['flat'], 2), round(R[sc]['dT0_min'], 2)) for sc in ('planet', 'grazBD', 'grazBD20', 'comp')},
              'Lambda G/G20/C', round(R['Lambda_grazBD'], 2), round(R['Lambda_grazBD20'], 2), round(R['Lambda_comp'], 2), flush=True)
        print('   profile', [(p['rho_c'], round(p['dchi2_vs_planet'], 1), round(p['k'], 3)) for p in prof], flush=True)
jdump(dict(method=__doc__, fits=res), f'{R3}/out/s03_fit_real.json')
print('done', round(time.time() - t0, 1))
