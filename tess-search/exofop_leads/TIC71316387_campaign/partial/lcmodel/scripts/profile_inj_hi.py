"""High-L3 continuation of the L3 profile for the injected configuration A (validates where the upper L3 wall falls)."""
import os
os.environ['NUMBA_NUM_THREADS'] = '1'; os.environ['OMP_NUM_THREADS'] = '1'
import sys, json, numpy as np
sys.argv = [sys.argv[0], 'A']
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/scripts')
import inj_recover as I
import lclib as L
from multiprocessing import get_context
GE = [(0.165, k, bn) for k in (0.42, 0.5, 0.6) for bn in (0.05, 0.2)]
def task(args):
    L3, g = args
    st = dict(I.DEF); st.update(rsum=g[0], k=g[1], bn=g[2], L3=L3)
    return dict(L3fixed=L3, **I.solve(st, dict(L3=L3)))
if __name__ == '__main__':
    grid = [0.96, 0.966, 0.972, 0.975, 0.978, 0.981, 0.984, 0.987]
    with get_context('spawn').Pool(4) as pool:
        res = pool.map(task, [(x, g) for x in grid for g in GE], chunksize=1)
    J = json.load(open(L.HERE + '/results/inj_recover_A.json'))
    J['all'] += res
    json.dump(J, open(L.HERE + '/results/inj_recover_A.json', 'w'), indent=1)
    free = min(r['tot'] for r in J['all'] if r['L3fixed'] is None)
    for x in grid:
        b = min([r for r in res if r['L3fixed'] == x], key=lambda r: r['tot'])
        print(f"L3={x:.3f} dtot={b['tot']-free:.2f} k={b['k']:.3f} bn={b['bn']:.3f} J={b['J']:.3f}", flush=True)
