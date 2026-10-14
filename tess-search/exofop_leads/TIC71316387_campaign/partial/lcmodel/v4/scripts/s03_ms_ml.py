"""ML fits of the main-sequence-constrained model (lib/msfit.py; smooth Pecaut & Mamajek relations; distance independent):
multi-start free fits and a profile over the primary mass M1 (fixed grid), several (q, bn, L3) starts each.
Usage: s03_ms_ml.py <dataset tag> <SIGJ> <FR_MU> <FR_SD> <out tag>"""
import os
os.environ['NUMBA_NUM_THREADS'] = '1'; os.environ['OMP_NUM_THREADS'] = '1'; os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np, sys, json, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4/lib')
import fit3, msfit as M, lclib as L
from multiprocessing import get_context
R = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4/results'
TAG = sys.argv[1]; SIGJ = float(sys.argv[2]); FRMU = float(sys.argv[3]); FRSD = float(sys.argv[4]); OUT = sys.argv[5]
_d = {}
def task(a):
    M.CFG.update(SIGJ=SIGJ, FR_MU=FRMU, FR_SD=FRSD)
    if 'd' not in _d: _d['d'] = fit3.load_data(f'{R}/data_{TAG}.pkl')
    start, fixed = a
    o = M.solve(_d['d'], start, fixed)
    o['fixedM1'] = fixed.get('M1', None); o['start'] = start
    return o
if __name__ == '__main__':
    t0 = time.time()
    jobs = [(dict(M1=m, q=q, L3=l3, x=0.3), {}) for m in (0.3, 0.45, 0.6, 0.8, 1.0) for q in (0.3, 0.5, 0.8) for l3 in (0.9, 0.97)]
    grid = [0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1.0, 1.1]
    jobs += [(dict(M1=m, q=q, x=x, L3=(0.97 if q < 0.6 else 0.93)), dict(M1=m)) for m in grid for q in (0.35, 0.55, 0.8) for x in (0.2, 0.75)]
    with get_context('spawn').Pool(int(os.environ.get('NPROC', '4'))) as pool:
        res = pool.map(task, jobs, chunksize=1)
    json.dump(res, open(f'{R}/msml_{OUT}.json', 'w'), indent=1)
    free = [r for r in res if r['fixedM1'] is None]; b = min(free, key=lambda r: r['tot'])
    print(f"[{OUT}] SIGJ={SIGJ} fR~N({FRMU},{FRSD})")
    def show(r, lab):
        print(f"{lab}: tot={r['tot']:.2f} (d={r['tot']-b['tot']:+.2f}) chi2={r['chi2']:.2f} prior={r['prior2']:.2f} M1={r['M1']:.3f} M2={r['M2']:.3f} fR={r['fR']:.3f} "
              f"dlnJ={r['dlnJ']:+.3f} L3={r['L3']:.4f} r1={r['r1']:.4f} k={r['k']:.3f} J={r['J']:.3f} b1={r['b1']:.2f} total={r['total']} esinw={r['esinw']:+.4f}", flush=True)
    show(b, 'free best ')
    for m in grid:
        rr = [r for r in res if r['fixedM1'] == m]
        show(min(rr, key=lambda r: r['tot']), f'M1={m:4.2f}   ')
    print('time %.0f' % (time.time() - t0))
