"""Robustness of the geometric ML solution to data products and modelling choices.  For each variant: least_squares from
5 starts (incl. a grazing start), report the best solution and the best grazing (k>0.6, bn>0.6) solution."""
import os
os.environ['NUMBA_NUM_THREADS'] = '1'; os.environ['OMP_NUM_THREADS'] = '1'
import numpy as np, sys, json, pickle, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/scripts')
import fit3, fitlib as F, lclib as L
from scipy.optimize import least_squares
from multiprocessing import get_context
def dataset(name):
    path = L.HERE + f'/results/data_{name}.pkl'
    if not os.path.exists(path):
        if name == 'tc': lcs = [L.load_tc(s) for s in (1, 28, 68, 95, 105, 106)]
        elif name == 'qlp': lcs = [L.load_spoc(1), L.load_spoc(28), L.load_spoc(68), L.load_qlp(95), L.load_qlp(105), L.load_tc(106)]
        elif name == 'spoc_only': lcs = [L.load_spoc(1), L.load_spoc(28), L.load_spoc(68)]
        elif name == 'recent_only': lcs = [L.load_tc(95), L.load_tc(105), L.load_tc(106)]
        d = F.prepare(lcs, hw=0.20, verbose=False); pickle.dump(d, open(path, 'wb'))
    return fit3.load_data(path)
VARIANTS = {
    'default_esinw0': ('default', ['s95', 's105', 's106'], dict(esinw=0.0)),
    'allscales': ('default', ['s1', 's28', 's95', 's105', 's106'], {}),
    'noscales': ('default', [], {}),
    'tc_all_s1': ('tc', [], {}),
    'qlp95_105': ('qlp', ['s95', 's105', 's106'], {}),
    'spoc_only': ('spoc_only', [], {}),
    'recent_only': ('recent_only', ['s105', 's106'], {}),
    'ld_fixed_claret_M2.5_M4.5': ('default', ['s95', 's105', 's106'], dict(q1a=(0.25 + 0.37) ** 2, q2a=0.25 / (2 * 0.62), q1b=(0.16 + 0.47) ** 2, q2b=0.16 / (2 * 0.63))),
}
BASE = ['dT0', 'dP', 'rsum', 'k', 'bn', 'J', 'L3', 'ecosw', 'esinw', 'q1a', 'q2a', 'q1b', 'q2b']
LB = dict(dT0=-30, dP=-50, rsum=0.03, k=0.05, bn=0.0, J=0.01, L3=0.0, ecosw=-0.05, esinw=-0.3, q1a=1e-3, q2a=1e-3, q1b=1e-3, q2b=1e-3)
UB = dict(dT0=30, dP=50, rsum=0.7, k=2.0, bn=0.999, J=1.6, L3=0.998, ecosw=0.05, esinw=0.3, q1a=0.999, q2a=0.999, q1b=0.999, q2b=0.999)
XS = dict(dT0=0.1, dP=0.5, rsum=0.01, k=0.05, bn=0.02, J=0.02, L3=0.002, ecosw=1e-4, esinw=0.01, q1a=0.05, q2a=0.05, q1b=0.05, q2b=0.05)
for s in fit3.SECS: LB[f's{s}'] = 0.6; UB[f's{s}'] = 1.4; XS[f's{s}'] = 0.02
STARTS = [dict(rsum=0.168, k=0.42, bn=0.27, L3=0.966), dict(rsum=0.165, k=0.5, bn=0.05, L3=0.975), dict(rsum=0.25, k=0.1, bn=0.8, L3=0.0),
          dict(rsum=0.22, k=0.71, bn=0.78, L3=0.90), dict(rsum=0.2, k=0.25, bn=0.55, L3=0.9)]
DEF = dict(dT0=0.9, dP=1.0, J=0.5, ecosw=0.00088, esinw=0.02, q1a=0.39, q2a=0.17, q1b=0.39, q2b=0.9, s1=1.0, s28=1.0, s95=0.94, s105=0.93, s106=0.97)
_cache = {}
def task(args):
    vname, si = args
    dname, sc, fixed = VARIANTS[vname]
    if dname not in _cache: _cache[dname] = dataset(dname)
    d = _cache[dname]
    names = [n for n in BASE + sc if n not in fixed]
    def to_p(x):
        p = dict(fixed); p.update(dict(zip(names, x)))
        p['T0'] = L.T00 + p['dT0'] / 1440; p['P'] = L.P0 + p['dP'] * 1e-6
        for s in fit3.SECS: p.setdefault(f's{s}', 1.0)
        p['cosi'] = p['bn'] * p['rsum']; return p
    def rf(x):
        p = to_p(x); return np.concatenate([fit3.resid(d, p), fit3.prior_resid(p, 'geo')])
    st = dict(DEF); st.update(STARTS[si])
    lb = np.array([LB[n] for n in names]); ub = np.array([UB[n] for n in names])
    x0 = np.clip(np.array([st[n] for n in names]), lb + 1e-6, ub - 1e-6)
    r = least_squares(rf, x0, bounds=(lb, ub), x_scale=np.array([XS[n] for n in names]), method='trf', diff_step=1e-4, max_nfev=500, ftol=1e-10, xtol=1e-10, gtol=1e-10)
    p = to_p(r.x)
    return dict(variant=vname, start=si, tot=float(2 * r.cost), **{k: float(v) for k, v in p.items()})
if __name__ == '__main__':
    for dn in ('tc', 'qlp', 'spoc_only', 'recent_only'): dataset(dn)
    jobs = [(v, i) for v in VARIANTS for i in range(len(STARTS))]
    with get_context('spawn').Pool(4) as pool:
        res = pool.map(task, jobs, chunksize=1)
    json.dump(res, open(L.HERE + '/results/variants.json', 'w'), indent=1)
    for v in VARIANTS:
        rr = [r for r in res if r['variant'] == v]
        b = min(rr, key=lambda r: r['tot'])
        gz = [r for r in rr if r['k'] > 0.6 and r['bn'] > 0.6]
        g = min(gz, key=lambda r: r['tot']) if gz else None
        lo = [r for r in rr if r['L3'] < 0.5]
        l0 = min(lo, key=lambda r: r['tot']) if lo else None
        print(f"{v:28s} best: L3={b['L3']:.4f} rsum={b['rsum']:.4f} k={b['k']:.3f} bn={b['bn']:.3f} J={b['J']:.3f} esinw={b['esinw']:.4f} ecosw={b['ecosw']:.5f} s={[round(b['s%d'%s],3) for s in fit3.SECS]}"
              + (f" | best grazing dchi2={g['tot']-b['tot']:.1f} (k={g['k']:.2f}, L3={g['L3']:.3f})" if g else ' | no grazing solution reached')
              + (f" | low-L3 dchi2={l0['tot']-b['tot']:.1f} (L3={l0['L3']:.2f}, k={l0['k']:.2f})" if l0 else ''), flush=True)
