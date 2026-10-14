"""Injection-recovery validation of the EB model: inject synthetic EBs (known parameters) into the REAL light curves at an
ephemeris shifted by +0.25 in phase (so the windows around the fake eclipses contain real noise, real spot modulation and
no real eclipse), then run the identical ML fit + L3 profile.  Tests code correctness and calibrates how well L3 / k / i
are identified with this data.  Usage: inj_recover.py <config>  (A|B|C)"""
import os
os.environ['NUMBA_NUM_THREADS'] = '1'; os.environ['OMP_NUM_THREADS'] = '1'; os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np, sys, json, time, pickle
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/scripts')
import fit3, fitlib as F, lclib as L, fit2
from ebmodel import eb_flux
from scipy.optimize import least_squares
from multiprocessing import get_context
CFG = sys.argv[1]
TSH = L.T00 + 0.25 * L.P0
base = dict(ecosw=0.0009, esinw=0.0, q1a=0.40, q2a=0.25, q1b=0.40, q2b=0.30, s95=1.0, s105=1.0, s106=1.0)
CONF = {'A': dict(rsum=0.1685, k=0.419, bn=0.267, J=0.502, L3=0.966, q1a=0.39, q2a=0.17, q1b=0.39, q2b=0.9, esinw=0.02),
        'B': dict(rsum=0.2277, k=0.129, bn=0.746, J=0.333, L3=0.5),
        'C': dict(rsum=0.2172, k=0.708, bn=0.78, J=0.484, L3=0.903)}
TRUE = dict(base); TRUE.update(CONF[CFG]); TRUE['T0'] = TSH; TRUE['P'] = L.P0; TRUE['cosi'] = TRUE['bn'] * TRUE['rsum']
path = L.HERE + f'/results/data_inj_{CFG}.pkl'
if not os.path.exists(path):
    lcs = L.default_set()
    for lc in lcs:
        u1a, u2a = fit2.q2u(TRUE['q1a'], TRUE['q2a']); u1b, u2b = fit2.q2u(TRUE['q1b'], TRUE['q2b'])
        m = eb_flux(lc['t'], TSH, L.P0, TRUE['rsum'], TRUE['k'], TRUE['cosi'], TRUE['J'], TRUE['L3'], TRUE['ecosw'], TRUE['esinw'],
                    u1a, u2a, u1b, u2b, lc['texp'], fit3.NSUB[lc['sector']], fit2.XG, fit2.WG)
        lc['f'] = lc['f'] * m
    dd = F.prepare(lcs, hw=0.20, P=L.P0, T0=TSH, verbose=False)
    pickle.dump(dd, open(path, 'wb'))
d = fit3.load_data(path)
SC = ['s95', 's105', 's106']
BASE = ['dT0', 'dP', 'rsum', 'k', 'bn', 'J', 'L3', 'ecosw', 'esinw', 'q1a', 'q2a', 'q1b', 'q2b'] + SC
LB = dict(dT0=-30, dP=-50, rsum=0.03, k=0.05, bn=0.0, J=0.01, L3=0.0, ecosw=-0.05, esinw=-0.3, q1a=1e-3, q2a=1e-3, q1b=1e-3, q2b=1e-3)
UB = dict(dT0=30, dP=50, rsum=0.7, k=2.0, bn=0.999, J=1.6, L3=0.998, ecosw=0.05, esinw=0.3, q1a=0.999, q2a=0.999, q1b=0.999, q2b=0.999)
XS = dict(dT0=0.1, dP=0.5, rsum=0.01, k=0.05, bn=0.02, J=0.02, L3=0.01, ecosw=1e-4, esinw=0.01, q1a=0.05, q2a=0.05, q1b=0.05, q2b=0.05)
for s in SC: LB[s] = 0.6; UB[s] = 1.4; XS[s] = 0.02
def to_p(names, x, fixed):
    p = dict(fixed); p.update(dict(zip(names, x)))
    p['T0'] = TSH + p['dT0'] / 1440.0; p['P'] = L.P0 + p['dP'] * 1e-6
    for s in fit3.SECS: p.setdefault(f's{s}', 1.0)
    p['cosi'] = p['bn'] * p['rsum']; return p
def resfun(x, names, fixed):
    p = to_p(names, x, fixed); return np.concatenate([fit3.resid(d, p), fit3.prior_resid(p, 'geo')])
def solve(start, fixed):
    names = [n for n in BASE if n not in fixed]
    lb = np.array([LB[n] for n in names]); ub = np.array([UB[n] for n in names])
    x0 = np.clip(np.array([start[n] for n in names]), lb + 1e-6, ub - 1e-6)
    r = least_squares(resfun, x0, args=(names, fixed), bounds=(lb, ub), x_scale=np.array([XS[n] for n in names]), method='trf',
                      diff_step=1e-4, max_nfev=400, ftol=1e-10, xtol=1e-10, gtol=1e-10)
    p = to_p(names, r.x, fixed); pr = fit3.prior_resid(p, 'geo')
    return dict(tot=float(2 * r.cost), chi2=fit3.chi2(d, p), **{k: float(v) for k, v in p.items()})
DEF = dict(dT0=0.0, dP=0.0, rsum=0.2, k=0.7, bn=0.75, J=0.45, L3=0.85, ecosw=0.0, esinw=0.0, q1a=0.4, q2a=0.3, q1b=0.4, q2b=0.3, s95=1.0, s105=1.0, s106=1.0)
GEOM = [(0.15, 0.6, 0.8), (0.2, 0.8, 0.85), (0.25, 1.0, 0.9), (0.12, 0.4, 0.75), (0.3, 0.7, 0.92), (0.2, 0.5, 0.6), (0.18, 0.9, 0.7)]
def task(args):
    L3, g = args
    st = dict(DEF); st.update(rsum=g[0], k=g[1], bn=g[2])
    if L3 is None: return dict(L3fixed=None, **solve(st, {}))
    st['L3'] = L3; return dict(L3fixed=L3, **solve(st, dict(L3=L3)))
if __name__ == '__main__':
    t0 = time.time()
    true_chi2 = fit3.chi2(d, dict(TRUE, T0=TSH))
    grid = [0.0, 0.3, 0.5, 0.75, 0.87, 0.93, 0.96, 0.975, 0.985, 0.993]
    jobs = [(None, g) for g in GEOM] + [(L3, g) for L3 in grid for g in GEOM[:5]]
    with get_context("spawn").Pool(4) as pool:
        res = pool.map(task, jobs, chunksize=1)
    best = {}
    for r in res:
        k = 'free' if r['L3fixed'] is None else r['L3fixed']
        if k not in best or r['tot'] < best[k]['tot']: best[k] = r
    print(f'config {CFG}: true params', {k: TRUE[k] for k in ('rsum', 'k', 'bn', 'J', 'L3')}, 'chi2 at truth %.2f' % true_chi2)
    for k in ['free'] + grid:
        r = best[k]
        print(f"  L3fix={k}: dtot={r['tot']-best['free']['tot']:.2f}  L3={r['L3']:.4f} rsum={r['rsum']:.4f} k={r['k']:.3f} bn={r['bn']:.3f} J={r['J']:.3f} esinw={r['esinw']:.4f}")
    json.dump(dict(true=TRUE, true_chi2=true_chi2, best={str(k): v for k, v in best.items()}, all=res), open(L.HERE + f'/results/inj_recover_{CFG}.json', 'w'), indent=1)
    print('time', time.time() - t0)
