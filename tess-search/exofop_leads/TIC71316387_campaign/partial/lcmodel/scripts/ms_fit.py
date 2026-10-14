"""ML fits of the main-sequence-constrained EB model (Pecaut & Mamajek relations, distance independent):
free M1, q=M2/M1, fR (common radius inflation, prior N(1.04,0.05)), dlnJ (surface-brightness-ratio scatter, prior N(0,0.1)),
bn, L3, ecosw, esinw, LD (u1+u2 prior), T0, P, TC-AP sector scales.  Several starts on an (M1, q) grid; also a profile over
M1 (fixed) to show how the light curve constrains the primary mass.  Usage: ms_fit.py [sigJ]"""
import os
os.environ['NUMBA_NUM_THREADS'] = '1'; os.environ['OMP_NUM_THREADS'] = '1'
import numpy as np, sys, json, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/scripts')
import fit3, lclib as L, ms_table as MS, fit2
from scipy.optimize import least_squares
from multiprocessing import get_context
SIGJ = float(sys.argv[1]) if len(sys.argv) > 1 else 0.1
d = fit3.load_data(L.HERE + '/results/data_default.pkl')
SC = ['s95', 's105', 's106']
BASE = ['dT0', 'dP', 'M1', 'q', 'fR', 'dlnJ', 'bn', 'L3', 'ecosw', 'esinw', 'q1a', 'q2a', 'q1b', 'q2b'] + SC
LB = dict(dT0=-30, dP=-50, M1=0.1, q=0.1, fR=0.85, dlnJ=-0.5, bn=0.0, L3=0.0, ecosw=-0.05, esinw=-0.3, q1a=1e-3, q2a=1e-3, q1b=1e-3, q2b=1e-3)
UB = dict(dT0=30, dP=50, M1=1.4, q=1.0, fR=1.3, dlnJ=0.5, bn=0.999, L3=0.998, ecosw=0.05, esinw=0.3, q1a=0.999, q2a=0.999, q1b=0.999, q2b=0.999)
XS = dict(dT0=0.1, dP=0.5, M1=0.02, q=0.02, fR=0.01, dlnJ=0.02, bn=0.02, L3=0.002, ecosw=1e-4, esinw=0.01, q1a=0.05, q2a=0.05, q1b=0.05, q2b=0.05)
for s in SC: LB[s] = 0.6; UB[s] = 1.4; XS[s] = 0.02
def to_p(names, x, fixed):
    p = dict(fixed); p.update(dict(zip(names, x)))
    p['T0'] = L.T00 + p['dT0'] / 1440.0; p['P'] = L.P0 + p['dP'] * 1e-6
    p.update(fit3.ms_derived(p))
    for s in fit3.SECS: p.setdefault(f's{s}', 1.0)
    p['cosi'] = p['bn'] * p['rsum']; return p
def priors(p):
    pr = fit3.prior_resid(p, 'geo')
    return np.concatenate([pr, [p['dlnJ'] / SIGJ, (p['fR'] - 1.04) / 0.05]])
def resfun(x, names, fixed):
    p = to_p(names, x, fixed)
    if p['M1'] * p['q'] < 0.09: p['q'] = 0.09 / p['M1']
    return np.concatenate([fit3.resid(d, p), priors(p)])
def solve(start, fixed):
    names = [n for n in BASE if n not in fixed]
    lb = np.array([LB[n] for n in names]); ub = np.array([UB[n] for n in names])
    x0 = np.clip(np.array([start[n] for n in names]), lb + 1e-6, ub - 1e-6)
    r = least_squares(resfun, x0, args=(names, fixed), bounds=(lb, ub), x_scale=np.array([XS[n] for n in names]), method='trf',
                      diff_step=1e-4, max_nfev=500, ftol=1e-10, xtol=1e-10, gtol=1e-10)
    p = to_p(names, r.x, fixed); pr = priors(p)
    return dict(tot=float(2 * r.cost), chi2=fit3.chi2(d, p), prior2=float(pr @ pr), **{k: float(v) for k, v in p.items()})
DEF = dict(dT0=0.9, dP=1.0, M1=0.5, q=0.5, fR=1.04, dlnJ=0.0, bn=0.3, L3=0.97, ecosw=0.00088, esinw=0.02, q1a=0.39, q2a=0.17, q1b=0.39, q2b=0.9,
           s95=0.94, s105=0.93, s106=0.97)
def task(args):
    M1fix, M1, q, L3 = args
    st = dict(DEF); st.update(M1=M1, q=q, L3=L3)
    fixed = {} if M1fix is None else dict(M1=M1fix)
    if M1fix is not None: st['M1'] = M1fix
    return dict(M1fixed=M1fix, start=(M1, q, L3), **solve(st, fixed))
if __name__ == '__main__':
    t0 = time.time()
    jobs = [(None, M1, q, L3) for M1 in (0.3, 0.5, 0.7, 0.9) for q in (0.4, 0.7) for L3 in (0.95, 0.98)]
    jobs += [(M1f, M1f, q, 0.97) for M1f in (0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1) for q in (0.35, 0.6, 0.85)]
    with get_context('spawn').Pool(4) as pool:
        res = pool.map(task, jobs, chunksize=1)
    json.dump(res, open(L.HERE + f'/results/ms_fit_sigJ{SIGJ}.json', 'w'), indent=1)
    best = {}
    for r in res:
        k = 'free' if r['M1fixed'] is None else r['M1fixed']
        if k not in best or r['tot'] < best[k]['tot']: best[k] = r
    for k in ['free'] + [x for x in sorted(best, key=lambda s: -1 if s == 'free' else s) if x != 'free']:
        r = best[k]
        print(f"M1fix={k}: tot={r['tot']:.2f} chi2={r['chi2']:.2f} prior2={r['prior2']:.2f} M1={r['M1']:.3f} M2={r['M2']:.3f} fR={r['fR']:.3f} dlnJ={r['dlnJ']:.3f} L3={r['L3']:.4f} rsum={r['rsum']:.4f} k={r['k']:.3f} J={r['J']:.3f} bn={r['bn']:.3f} esinw={r['esinw']:.4f}", flush=True)
    print('time', time.time() - t0)
