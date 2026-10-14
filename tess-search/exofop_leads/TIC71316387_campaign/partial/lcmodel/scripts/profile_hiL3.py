"""Continuation of the L3 profile on the small-companion branch for L3 >= 0.96 (starts with k=0.4-0.8, small impact)."""
import os
os.environ['NUMBA_NUM_THREADS'] = '1'; os.environ['OMP_NUM_THREADS'] = '1'
import sys, json, numpy as np
sys.argv = [sys.argv[0], 'A']
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/scripts')
import ml_fit as M
import lclib as L
from multiprocessing import get_context
GE = [(0.165, k, bn) for k in (0.42, 0.5, 0.6, 0.75) for bn in (0.05, 0.15, 0.25)]
def task(args):
    L3, g = args
    st = dict(M.DEF); st.update(rsum=g[0], k=g[1], bn=g[2], L3=L3)
    return dict(L3fixed=L3, start=g, **M.solve(st, dict(L3=L3)))
if __name__ == '__main__':
    grid = [0.963, 0.966, 0.969, 0.972, 0.975, 0.978, 0.981, 0.984, 0.987, 0.99]
    jobs = [(x, g) for x in grid for g in GE]
    with get_context('spawn').Pool(4) as pool:
        res = pool.map(task, jobs, chunksize=1)
    json.dump(res, open(L.HERE + '/results/ml_profile_A_hiL3.json', 'w'), indent=1)
    best = {}
    for r in res:
        tot = r['chi2'] + r['prior2']
        if r['L3fixed'] not in best or tot < best[r['L3fixed']][0]: best[r['L3fixed']] = (tot, r)
    for k in grid:
        t, r = best[k]
        print(f"L3={k:.3f} tot={t:.2f} rsum={r['rsum']:.4f} k={r['k']:.3f} bn={r['bn']:.3f} J={r['J']:.3f} esinw={r['esinw']:.4f}", flush=True)
