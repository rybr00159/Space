"""Total-occultation (flat-bottomed secondary, k<1, small impact) branch versus grazing (V-shaped, k~0.6-1.6, bn>0.7) branch:
best chi2 of each branch for the real data and for the injected configurations A (flat truth), B (low-L3 small-k truth) and C (grazing truth).
Branch membership after convergence: 'flat' if star 2 is fully occulted at secondary (b1 + k <= 1 with b1 = cos i / r1), 'grazing' if b1 + k > 1 and b1 > 1 - k (partial), with k > 0.55."""
import os
os.environ['NUMBA_NUM_THREADS'] = '1'; os.environ['OMP_NUM_THREADS'] = '1'
import numpy as np, sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/scripts')
import fit3, lclib as L
from scipy.optimize import least_squares
from multiprocessing import get_context
SETS = {'real': (L.HERE + '/results/data_default.pkl', L.T00), 'injA': (L.HERE + '/results/data_inj_A.pkl', L.T00 + 0.25 * L.P0),
        'injB': (L.HERE + '/results/data_inj_B.pkl', L.T00 + 0.25 * L.P0), 'injC': (L.HERE + '/results/data_inj_C.pkl', L.T00 + 0.25 * L.P0)}
BASE = ['dT0', 'dP', 'rsum', 'k', 'bn', 'J', 'L3', 'ecosw', 'esinw', 'q1a', 'q2a', 'q1b', 'q2b', 's95', 's105', 's106']
LB = dict(dT0=-30, dP=-50, rsum=0.03, k=0.05, bn=0.0, J=0.01, L3=0.0, ecosw=-0.05, esinw=-0.3, q1a=1e-3, q2a=1e-3, q1b=1e-3, q2b=1e-3, s95=0.6, s105=0.6, s106=0.6)
UB = dict(dT0=30, dP=50, rsum=0.7, k=2.0, bn=0.999, J=1.6, L3=0.998, ecosw=0.05, esinw=0.3, q1a=0.999, q2a=0.999, q1b=0.999, q2b=0.999, s95=1.4, s105=1.4, s106=1.4)
XS = dict(dT0=0.1, dP=0.5, rsum=0.01, k=0.05, bn=0.02, J=0.02, L3=0.002, ecosw=1e-4, esinw=0.01, q1a=0.05, q2a=0.05, q1b=0.05, q2b=0.05, s95=0.02, s105=0.02, s106=0.02)
FLAT = [dict(rsum=0.168, k=0.42, bn=0.27, L3=0.966), dict(rsum=0.165, k=0.5, bn=0.05, L3=0.975), dict(rsum=0.19, k=0.25, bn=0.55, L3=0.90)]
GRAZ = [dict(rsum=0.217, k=0.71, bn=0.78, L3=0.903), dict(rsum=0.25, k=1.0, bn=0.85, L3=0.85), dict(rsum=0.3, k=0.8, bn=0.9, L3=0.7), dict(rsum=0.4, k=1.4, bn=0.92, L3=0.3)]
DEF = dict(dT0=0.0, dP=0.0, J=0.5, ecosw=0.0009, esinw=0.0, q1a=0.39, q2a=0.25, q1b=0.39, q2b=0.3, s95=1.0, s105=1.0, s106=1.0)
_c = {}
def task(args):
    sname, branch, i, esin0 = args
    path, T0 = SETS[sname]
    if sname not in _c: _c[sname] = fit3.load_data(path)
    d = _c[sname]
    fixed = dict(esinw=0.0) if esin0 else {}
    names = [n for n in BASE if n not in fixed]
    def to_p(x):
        p = dict(fixed); p.update(dict(zip(names, x))); p['T0'] = T0 + p['dT0'] / 1440; p['P'] = L.P0 + p['dP'] * 1e-6
        for s in fit3.SECS: p.setdefault(f's{s}', 1.0)
        p['cosi'] = p['bn'] * p['rsum']; return p
    def rf(x):
        p = to_p(x); return np.concatenate([fit3.resid(d, p), fit3.prior_resid(p, 'geo')])
    st = dict(DEF); st.update((FLAT if branch == 'flat' else GRAZ)[i])
    lb = np.array([LB[n] for n in names]); ub = np.array([UB[n] for n in names])
    # for the grazing branch keep k >= 0.55 and bn >= 0.6 so the optimiser stays on that branch
    if branch == 'graz':
        lb[names.index('k')] = 0.55; lb[names.index('bn')] = 0.6
    x0 = np.clip(np.array([st[n] for n in names]), lb + 1e-6, ub - 1e-6)
    r = least_squares(rf, x0, bounds=(lb, ub), x_scale=np.array([XS[n] for n in names]), method='trf', diff_step=1e-4, max_nfev=500, ftol=1e-10, xtol=1e-10, gtol=1e-10)
    p = to_p(r.x); r1 = p['rsum'] / (1 + p['k']); b1 = p['cosi'] / r1
    return dict(set=sname, branch=branch, esinw0=esin0, start=i, tot=float(2 * r.cost), b1=float(b1), total_occ=bool(b1 + p['k'] <= 1), **{k: float(v) for k, v in p.items()})
if __name__ == '__main__':
    jobs = [(s, 'flat', i, e) for s in SETS for i in range(len(FLAT)) for e in (False, True)] + [(s, 'graz', i, e) for s in SETS for i in range(len(GRAZ)) for e in (False, True)]
    with get_context('spawn').Pool(4) as pool:
        res = pool.map(task, jobs, chunksize=1)
    json.dump(res, open(L.HERE + '/results/branch_test.json', 'w'), indent=1)
    for s in SETS:
        for e in (False, True):
            f = min([r for r in res if r['set'] == s and r['branch'] == 'flat' and r['esinw0'] == e], key=lambda r: r['tot'])
            g = min([r for r in res if r['set'] == s and r['branch'] == 'graz' and r['esinw0'] == e], key=lambda r: r['tot'])
            print(f"{s:5s} esinw0={e!s:5s}: flat-branch tot={f['tot']:.2f} (L3={f['L3']:.3f} k={f['k']:.3f} b1={f['b1']:.2f} total_occ={f['total_occ']}) | grazing-branch tot={g['tot']:.2f} (L3={g['L3']:.3f} k={g['k']:.3f} bn={g['bn']:.2f} J={g['J']:.2f} esinw={g['esinw']:.3f})  => chi2(graz) - chi2(flat) = {g['tot']-f['tot']:+.1f}", flush=True)
