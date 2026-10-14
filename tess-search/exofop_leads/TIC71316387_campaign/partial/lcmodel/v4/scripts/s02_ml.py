"""Maximum-likelihood geometric fits of the default data set:
 stage 1: multi-start fits restricted to the TOTAL/ANNULAR branch and to the PARTIAL (grazing) branch;
 stage 2: profile likelihoods in L3 (fixed grid) and in r1 = R1/a (fixed grid; r1 sets the primary's density) for each branch,
          started from the stage-1 best fit of that branch and from one alternative start.
Usage: s02_ml.py <dataset tag> [scale set: recent|all|none]"""
import os
os.environ['NUMBA_NUM_THREADS'] = '1'; os.environ['OMP_NUM_THREADS'] = '1'; os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np, sys, json, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4/lib')
import fit3, geofit as G, lclib as L
from multiprocessing import get_context
R = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4/results'
TAG = sys.argv[1] if len(sys.argv) > 1 else 'default'
SCN = sys.argv[2] if len(sys.argv) > 2 else 'recent'
SCALES = {'recent': ('s95', 's105', 's106'), 'all': ('s1', 's28', 's95', 's105', 's106'), 'none': ()}[SCN]
_d = {}
def data():
    if 'd' not in _d: _d['d'] = fit3.load_data(f'{R}/data_{TAG}.pkl')
    return _d['d']
def task(a):
    kind, branch, start, fixed, use_r1 = a
    o = G.solve(data(), start, fixed, branch=branch, scales=SCALES, use_r1=use_r1)
    o.update(kind=kind, fixed=fixed)
    return o
if __name__ == '__main__':
    t0 = time.time()
    jobs = [('branch', 'total', s, {}, False) for s in G.STARTS_TOTAL] + [('branch', 'partial', s, {}, False) for s in G.STARTS_PARTIAL]
    with get_context('spawn').Pool(int(os.environ.get('NPROC', '4'))) as pool:
        r1 = pool.map(task, jobs, chunksize=1)
        best = {b: min([r for r in r1 if r['branch'] == b], key=lambda r: r['tot']) for b in ('total', 'partial')}
        for b in ('total', 'partial'):
            o = best[b]
            print(f"stage1 {b:8s}: chi2 {o['tot']:.2f}  L3={o['L3']:.4f} rsum={o['rsum']:.4f} r1={o['r1']:.4f} k={o['k']:.3f} bn={o['bn']:.3f} b1={o['b1']:.3f} "
                  f"J={o['J']:.3f} ecosw={o['ecosw']:.5f} esinw={o['esinw']:.4f} s={[round(o['s%d'%s],3) for s in fit3.SECS]}", flush=True)
        keys = ['rsum', 'k', 'x', 'J', 'L3', 'ecosw', 'esinw', 'q1a', 'q2a', 'q1b', 'q2b', 'dT0', 'dP'] + list(SCALES)
        L3g = [0.0, 0.3, 0.5, 0.7, 0.8, 0.85, 0.9, 0.93, 0.95, 0.96, 0.965, 0.97, 0.975, 0.98, 0.985, 0.99]
        r1g = [0.08, 0.09, 0.10, 0.105, 0.11, 0.115, 0.12, 0.13, 0.14, 0.15, 0.165, 0.18, 0.20, 0.22, 0.25]
        alt = {'total': dict(G.STARTS_TOTAL[1]), 'partial': dict(G.STARTS_PARTIAL[1])}
        jobs = []
        for b in ('total', 'partial'):
            st0 = {k: best[b][k] for k in keys}
            for v in L3g:
                for st in (st0, alt[b]):
                    jobs.append(('L3prof', b, dict(st, L3=v), dict(L3=v), False))
            for v in r1g:
                for st in (st0, alt[b]):
                    s2 = dict(st); s2['r1'] = v
                    jobs.append(('r1prof', b, s2, dict(r1=v), True))
        r2 = pool.map(task, jobs, chunksize=1)
    allres = r1 + r2
    json.dump(allres, open(f'{R}/ml_{TAG}_{SCN}.json', 'w'), indent=1)
    gbest = min(best.values(), key=lambda r: r['tot'])['tot']
    for kind, key, grid in (('L3prof', 'L3', L3g), ('r1prof', 'r1', r1g)):
        print(f'--- profile in {key} (chi2 - global best {gbest:.2f})')
        for v in grid:
            line = f'{key}={v:6.3f}: '
            for b in ('total', 'partial'):
                rr = [r for r in r2 if r['kind'] == kind and r['branch'] == b and abs(r['fixed'][key] - v) < 1e-9]
                o = min(rr, key=lambda r: r['tot'])
                line += f"| {b:7s} d={o['tot']-gbest:7.2f} L3={o['L3']:.3f} r1={o['r1']:.4f} k={o['k']:.3f} b1={o['b1']:.2f} J={o['J']:.3f} esinw={o['esinw']:+.3f} "
            print(line, flush=True)
    print('time %.0f s' % (time.time() - t0))
