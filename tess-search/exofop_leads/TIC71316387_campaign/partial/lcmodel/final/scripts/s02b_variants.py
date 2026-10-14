"""Robustness of the geometric ML solution to data products and modelling choices: for each variant, best total-branch and
best partial-branch solutions (multi-start) and the L3 = 0.90 / 0.95 profile points.
Variants: default; all sector scales free (S68 reference); own TESScut photometry for all sectors; QLP for S95/S105;
SPOC sectors only (S1/S28/S68); recent sectors only (S95/S105/S106); LD fixed at Claret (2017) PHOENIX values for
3500 K / 3100 K dwarfs; circular orbit (e=0); +-0.30 d windows."""
import os
os.environ['NUMBA_NUM_THREADS'] = '1'; os.environ['OMP_NUM_THREADS'] = '1'; os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np, sys, json, pickle, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/final/lib')
import fit3, fitlib as F, lclib as L, geofit as G
from multiprocessing import get_context
R = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/final/results'
def build(name):
    path = f'{R}/data_{name}.pkl'
    if not os.path.exists(path):
        if name == 'spoc_only': lcs = [L.load_spoc(s) for s in (1, 28, 68)]
        elif name == 'recent_only': lcs = [L.load_tc(s) for s in (95, 105, 106)]
        elif name == 'default_hw030': lcs = L.default_set()
        d = F.prepare(lcs, hw=0.30 if name == 'default_hw030' else 0.20, verbose=False); pickle.dump(d, open(path, 'wb'))
    return path
# Claret 2017 PHOENIX (LSM) TESS: 3500 K: u1=0.1604 u2=0.4325 ; 3100 K: u1=0.1606 u2=0.4974  -> Kipping q1=(u1+u2)^2, q2=u1/(2(u1+u2))
def kq(u1, u2): return (u1 + u2) ** 2, u1 / (2 * (u1 + u2))
q1a, q2a = kq(0.1604, 0.4325); q1b, q2b = kq(0.1606, 0.4974)
V = {'default': ('default', ('s95', 's105', 's106'), {}),
     'allscales': ('default', ('s1', 's28', 's95', 's105', 's106'), {}),
     'tc_all': ('tc', ('s1', 's28', 's95', 's105', 's106'), {}),
     'qlp95_105': ('qlp', ('s95', 's105', 's106'), {}),
     'spoc_only': ('spoc_only', (), {}),
     'recent_only': ('recent_only', ('s105', 's106'), {}),
     'ld_claret_fixed': ('default', ('s95', 's105', 's106'), dict(q1a=q1a, q2a=q2a, q1b=q1b, q2b=q2b)),
     'circular': ('default', ('s95', 's105', 's106'), dict(ecosw=0.0, esinw=0.0)),
     'windows_0.30d': ('default_hw030', ('s95', 's105', 's106'), {})}
_c = {}
def task(a):
    vn, branch, si, L3fix = a
    dn, sc, fixed = V[vn]
    if dn not in _c:
        _c.clear(); _c[dn] = fit3.load_data(build(dn) if dn in ('spoc_only', 'recent_only', 'default_hw030') else f'{R}/data_{dn}.pkl')
    fx = dict(fixed)
    if L3fix is not None: fx['L3'] = L3fix
    st = dict((G.STARTS_TOTAL if branch == 'total' else G.STARTS_PARTIAL)[si])
    if L3fix is not None: st['L3'] = L3fix
    o = G.solve(_c[dn], st, fx, branch=branch, scales=sc, ld_prior=not ('q1a' in fixed))
    o.update(variant=vn, L3fix=L3fix, start=si)
    return o
if __name__ == '__main__':
    t0 = time.time()
    for dn in ('spoc_only', 'recent_only', 'default_hw030'): build(dn)
    jobs = []
    for vn in V:
        jobs += [(vn, 'total', i, None) for i in (0, 1, 2)] + [(vn, 'partial', i, None) for i in (0, 1, 2, 3)]
        jobs += [(vn, 'total', i, v) for v in (0.90, 0.95) for i in (0, 2)]
    jobs.sort(key=lambda j: j[0])
    with get_context('spawn').Pool(int(os.environ.get('NPROC', '4'))) as pool:
        res = pool.map(task, jobs, chunksize=4)
    json.dump(res, open(f'{R}/variants.json', 'w'), indent=1)
    for vn in V:
        rr = [r for r in res if r['variant'] == vn]
        bt = min([r for r in rr if r['branch'] == 'total' and r['L3fix'] is None], key=lambda r: r['tot'])
        bp = min([r for r in rr if r['branch'] == 'partial' and r['L3fix'] is None], key=lambda r: r['tot'])
        p90 = min([r for r in rr if r['L3fix'] == 0.90], key=lambda r: r['tot']); p95 = min([r for r in rr if r['L3fix'] == 0.95], key=lambda r: r['tot'])
        print(f"{vn:16s} total: L3={bt['L3']:.4f} r1={bt['r1']:.4f} k={bt['k']:.3f} b1={bt['b1']:.2f} J={bt['J']:.3f} esinw={bt['esinw']:+.4f} | "
              f"partial-total dchi2={bp['tot']-bt['tot']:+.1f} (L3={bp['L3']:.2f},k={bp['k']:.2f}) | dchi2(L3=0.90)={p90['tot']-bt['tot']:+.1f} dchi2(L3=0.95)={p95['tot']-bt['tot']:+.1f}"
              + (f" | s={[round(bt.get('s%d'%s,1),3) for s in fit3.SECS]}"), flush=True)
    print('time %.0f' % (time.time() - t0))
