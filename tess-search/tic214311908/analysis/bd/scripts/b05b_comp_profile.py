"""b05b: profile of chi2(companion-hosted body) - chi2(planet) versus the radius ratio k on the companion
(k = R_occulter / R_companion; D free, b free, rho_c free, LD prior), real REF and INV-clean data, beta noise.
Also stores the k = 0.40 best fit used as the 'compBD40' injection truth.
Outputs: out/b05b_comp_profile.json, out/b05_comp_k40.json"""
import sys; sys.path.insert(0, 'tic214311908/analysis/bd/scripts')
from blib import *
from lclib import *
import fitlib as F
from scipy.optimize import minimize
NP = {'ref': 2, 'inv': 3}
FJ = jload(f'{BD}/out/b05_fits.json')
prof = {}; k40 = {}
for ds in ('ref', 'inv'):
    bt = {int(s): v['beta'] for s, v in FJ['beta'][ds].items()}
    d = load_npz(f'{BD}/data/{ds}_real.npz', npoly=NP[ds], beta=bt)
    c0 = FJ['fits'][f'{ds}_beta']['planet']['chi2']; dT = FJ['fits'][f'{ds}_beta']['planet']['x'][0]
    rows = []
    for kk in (0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 1.0):
        F.LIM['fixk'] = dict(k=(kk*0.999, kk*1.001), D=(0.005, 1.0))
        f_ = lambda x: F.chi2(d, x, 'fixk')
        best = None
        for bb in (0.0, 0.3, 0.6):
            for r in (5., 10., 20.):
                s0 = [dT, np.log(kk), bb, np.log(r), U1, U2, np.log(min(0.9, 0.0157/kk**2*1.1))]
                r_ = minimize(f_, s0, method='Nelder-Mead', options=dict(maxiter=2500, maxfev=2500, xatol=1e-6, fatol=1e-4, adaptive=True))
                if best is None or r_.fun < best.fun: best = r_
        best = minimize(f_, best.x, method='Nelder-Mead', options=dict(maxiter=2500, maxfev=2500, xatol=1e-7, fatol=1e-5, adaptive=True))
        _, k_, bb_, rho_, u1, u2, D_ = F.unpack(best.x, 'fixk')
        rows.append(dict(k=kk, dchi2=float(best.fun - c0), D=float(D_), b=float(bb_*(1 + k_)), rho_c=float(rho_), x=[float(v) for v in best.x]))
        print(ds, rows[-1], flush=True)
        if abs(kk - 0.40) < 1e-6: k40[ds] = rows[-1]
    prof[ds] = rows
jdump(dict(method=__doc__, profile=prof), f'{BD}/out/b05b_comp_profile.json')
jdump(k40, f'{BD}/out/b05_comp_k40.json')
