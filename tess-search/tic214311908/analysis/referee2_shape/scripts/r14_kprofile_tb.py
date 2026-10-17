"""r14: (a) chi2 cost of a large occulter on the M3V itself (brown dwarf / late-M grazing): planet model with k fixed
at 0.14, 0.17, 0.20, 0.25 (b free incl. grazing, LD prior, rho ~ N(9.41, 0.88)) vs the free planet fit (r03);
(b) TESS-band brightness-temperature limit of the occulter from the r04 secondary (3-sigma upper limit, equal
occulted areas, blackbodies with the TESS response; star 3350 K) - also with the analyst's 2.2 ppt.
Output: out/kprofile_tb.json"""
import sys, json, numpy as np
sys.path.insert(0, 'tic214311908/analysis/referee2_shape/scripts')
from rlib import *
from scipy.optimize import minimize
out = {'kprofile': {}}
for ds in ('inv', 'ref'):
    S = json.load(open(f'{RD}/out/scen_{ds}_fake.json')); ref = S['planet_free']['chi2']
    D = load(ds, 'fake'); lc = LC(D)
    row = {}
    for kf in (0.14, 0.17, 0.20, 0.25):
        def c2(x):
            T0, P, b, lr, u1, u2 = x
            if not (0 <= b < 1 + kf and u1 + u2 < 1 and u1 > 0 and u1 + 2*u2 > 0): return 1e9
            return lc.chi2(lc.transit(T0, P, kf, aR_from_rho(10**lr, P), b, u1, u2)) + ((u1 - 0.125)/0.1)**2 + ((u2 - 0.533)/0.1)**2 + ((10**lr - 9.41)/0.88)**2
        best = None
        for b0 in (0.5, 0.8, 0.95, 1.05, 1.1):
            if b0 >= 1 + kf: continue
            r = minimize(c2, [S['planet_free']['x'][0], S['planet_free']['x'][1], b0, np.log10(9.41), 0.125, 0.533], method='Nelder-Mead', options=dict(maxiter=6000, maxfev=6000, xatol=1e-8, fatol=1e-5))
            if best is None or r.fun < best.fun: best = r
        row[str(kf)] = dict(dchi2=float(best.fun - ref), b=float(best.x[2]))
    out['kprofile'][ds] = row; print(ds, row, flush=True)
# brightness temperature
filt = np.loadtxt(f'{B}/shape/data/filt_TESS.Red.dat'); lam = filt[:, 0]*1e-10; tr = filt[:, 1]
def S_(T): return np.trapezoid(tr*lam/(lam**5*(np.exp(1.4388e-2/(lam*T)) - 1)), lam)   # photon-weighted
sec = json.load(open(f'{RD}/out/secondary.json'))['phase05']['planet_white']
up = sec['depth_ppt'] + 3*max(sec['empirical_std_ppt'], sec['empirical_robust_ppt'])
Ts = np.linspace(1500, 3350, 1851); rat = np.array([S_(t) for t in Ts])/S_(3350)
for lab, u in (('own', up), ('analyst_2p2', 2.2)):
    out[f'Tb_{lab}'] = dict(upper_ppt=float(u), ratio=float(u/15.9), Tb_K=float(np.interp(u/15.9, rat, Ts)))
print({k: v for k, v in out.items() if k.startswith('Tb')})
json.dump(out, open(f'{RD}/out/kprofile_tb.json', 'w'), indent=1)
