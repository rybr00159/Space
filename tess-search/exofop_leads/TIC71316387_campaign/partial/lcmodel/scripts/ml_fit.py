"""Maximum-likelihood fits (scipy least_squares, trf) of the geometric EB model + profile likelihood in L3.
Usage: ml_fit.py <variant>   variant: A = SPOC sectors s=1, TC-AP sectors (95,105,106) free scales
                                     B = all sector scales free except S68 (=1)
                                     C = all scales fixed to 1"""
import os
os.environ['NUMBA_NUM_THREADS'] = '1'; os.environ['OMP_NUM_THREADS'] = '1'; os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np, sys, json, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/scripts')
import fit3, lclib as L
from scipy.optimize import least_squares
from multiprocessing import Pool
VAR = sys.argv[1] if len(sys.argv) > 1 else 'A'
d = fit3.load_data(L.HERE + '/results/data_default.pkl')
SC = {'A': ['s95', 's105', 's106'], 'B': ['s1', 's28', 's95', 's105', 's106'], 'C': []}[VAR]
BASE = ['dT0', 'dP', 'rsum', 'k', 'bn', 'J', 'L3', 'ecosw', 'esinw', 'q1a', 'q2a', 'q1b', 'q2b'] + SC
LB = dict(dT0=-30, dP=-50, rsum=0.03, k=0.1, bn=0.0, J=0.01, L3=0.0, ecosw=-0.05, esinw=-0.3, q1a=1e-3, q2a=1e-3, q1b=1e-3, q2b=1e-3)
UB = dict(dT0=30, dP=50, rsum=0.7, k=1.6, bn=0.999, J=1.6, L3=0.998, ecosw=0.05, esinw=0.3, q1a=0.999, q2a=0.999, q1b=0.999, q2b=0.999)
XS = dict(dT0=0.1, dP=0.5, rsum=0.01, k=0.05, bn=0.02, J=0.02, L3=0.01, ecosw=1e-4, esinw=0.01, q1a=0.05, q2a=0.05, q1b=0.05, q2b=0.05)
for s in SC: LB[s] = 0.6; UB[s] = 1.4; XS[s] = 0.02
def to_p(names, x, fixed):
    p = dict(fixed); p.update(dict(zip(names, x)))
    p['T0'] = L.T00 + p['dT0'] / 1440.0; p['P'] = L.P0 + p['dP'] * 1e-6
    for s in fit3.SECS: p.setdefault(f's{s}', 1.0)
    p['cosi'] = p['bn'] * p['rsum']
    return p
def resfun(x, names, fixed):
    p = to_p(names, x, fixed)
    return np.concatenate([fit3.resid(d, p), fit3.prior_resid(p, 'geo')])
def solve(start, fixed):
    names = [n for n in BASE if n not in fixed]
    x0 = np.array([start[n] for n in names])
    lb = np.array([LB[n] for n in names]); ub = np.array([UB[n] for n in names])
    x0 = np.clip(x0, lb + 1e-6, ub - 1e-6)
    r = least_squares(resfun, x0, args=(names, fixed), bounds=(lb, ub), x_scale=np.array([XS[n] for n in names]),
                      method='trf', diff_step=1e-4, max_nfev=400, ftol=1e-10, xtol=1e-10, gtol=1e-10)
    p = to_p(names, r.x, fixed)
    c2 = fit3.chi2(d, p); pr = fit3.prior_resid(p, 'geo')
    return dict(cost2=float(2 * r.cost), chi2=c2, prior2=float(pr @ pr), nfev=int(r.nfev), status=int(r.status), names=names,
                **{k: float(v) for k, v in p.items()})
DEF = dict(dT0=0.0, dP=0.0, rsum=0.2, k=0.7, bn=0.75, J=0.45, L3=0.85, ecosw=0.0, esinw=0.0, q1a=0.4, q2a=0.3, q1b=0.4, q2b=0.3,
           s1=1.0, s28=1.0, s95=1.0, s105=1.0, s106=1.0)
GEOM = [(0.15, 0.6, 0.8), (0.2, 0.8, 0.85), (0.25, 1.0, 0.9), (0.12, 0.4, 0.75), (0.3, 0.7, 0.92), (0.2, 0.5, 0.6), (0.18, 0.9, 0.7)]
def task(args):
    L3, geom = args
    st = dict(DEF); st.update(rsum=geom[0], k=geom[1], bn=geom[2])
    if L3 is None:
        st['L3'] = 0.85
        return dict(L3fixed=None, start=geom, **solve(st, {}))
    st['L3'] = L3
    return dict(L3fixed=L3, start=geom, **solve(st, dict(L3=L3)))
if __name__ == '__main__':
    t0 = time.time()
    grid = [0.0, 0.3, 0.5, 0.65, 0.75, 0.82, 0.87, 0.9, 0.92, 0.94, 0.95, 0.96, 0.97, 0.98, 0.985, 0.99, 0.995]
    jobs = [(None, g) for g in GEOM] + [(L3, g) for L3 in grid for g in GEOM]
    with Pool(4) as pool:
        res = []
        for i, r in enumerate(pool.imap_unordered(task, jobs)):
            res.append(r)
            print(f"[{time.time()-t0:6.0f}s] L3fix={r['L3fixed']} start={r['start']} chi2={r['chi2']:.2f} prior2={r['prior2']:.2f} L3={r['L3']:.4f} rsum={r['rsum']:.4f} k={r['k']:.3f} bn={r['bn']:.3f} J={r['J']:.3f} ecosw={r['ecosw']:.5f} esinw={r['esinw']:.4f} q={r['q1a']:.2f},{r['q2a']:.2f},{r['q1b']:.2f},{r['q2b']:.2f} s={[round(r.get('s%d'%s,1),3) for s in fit3.SECS]} nfev={r['nfev']}", flush=True)
    json.dump(res, open(L.HERE + f'/results/ml_profile_{VAR}.json', 'w'), indent=1)
    print('done', time.time() - t0)
