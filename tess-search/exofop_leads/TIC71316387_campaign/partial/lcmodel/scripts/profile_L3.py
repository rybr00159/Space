"""Profile likelihood in L3 for the geometric model (TC-AP sector scales s95,s105,s106 free; SPOC sectors s=1).
For each L3 on a grid, optimise all other parameters from several starts (Nelder-Mead + Powell)."""
import os
os.environ['NUMBA_NUM_THREADS'] = '1'
import numpy as np, sys, pickle, json, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/scripts')
import fit2, lclib as L
from scipy.optimize import minimize
from multiprocessing import Pool
d = pickle.load(open(L.HERE + '/results/data_default.pkl', 'rb'))
P0 = dict(T0=L.T00, P=L.P0, rsum=0.2, k=0.7, cosi=0.12, J=0.45, L3=0.85, ecosw=0.0, esinw=0.0, q1a=0.4, q2a=0.3, q1b=0.4, q2b=0.3, s95=1.0, s105=1.0, s106=1.0)
FREE = [k for k in P0 if k != 'L3']
STARTS = [dict(rsum=0.15, k=0.6, cosi=0.11), dict(rsum=0.22, k=0.8, cosi=0.17), dict(rsum=0.30, k=1.0, cosi=0.25),
          dict(rsum=0.12, k=0.4, cosi=0.09), dict(rsum=0.4, k=0.7, cosi=0.36)]
def run(L3):
    best = None
    post = fit2.Post(d, FREE, dict(L3=L3))
    for s in STARTS:
        pp = dict(P0); pp.update(s)
        x0 = np.array([pp[n] for n in FREE])
        # pick cosi so the start is finite and depth roughly right: scan cosi
        cands = []
        for f in np.linspace(0.3, 0.98, 15):
            pp['cosi'] = f * pp['rsum']; x = np.array([pp[n] for n in FREE]); v = post(x)
            if np.isfinite(v): cands.append((v, x))
        if not cands: continue
        x0 = max(cands, key=lambda c: c[0])[1]
        nf = lambda x: -post(x) if np.isfinite(post(x)) else 1e30
        r = minimize(nf, x0, method='Nelder-Mead', options=dict(maxfev=12000, xatol=1e-9, fatol=1e-3, adaptive=True))
        r = minimize(nf, r.x, method='Powell', options=dict(maxfev=12000, xtol=1e-8, ftol=1e-9))
        r = minimize(nf, r.x, method='Nelder-Mead', options=dict(maxfev=8000, xatol=1e-9, fatol=1e-4, adaptive=True))
        if best is None or r.fun < best.fun: best = r
    p = post.full(best.x)
    c2 = fit2.chi2(d, p)
    out = dict(L3=L3, lnpost=float(-best.fun), chi2=float(c2), **{k: float(v) for k, v in p.items()})
    print(f"L3={L3:.3f} chi2={c2:.2f} lnpost={-best.fun:.2f} rsum={p['rsum']:.4f} k={p['k']:.3f} cosi={p['cosi']:.4f} b/rsum={p['cosi']/p['rsum']:.3f} J={p['J']:.3f} ecosw={p['ecosw']:.5f} esinw={p['esinw']:.4f} s95={p['s95']:.3f} s105={p['s105']:.3f} s106={p['s106']:.3f}", flush=True)
    return out
if __name__ == '__main__':
    grid = [0.0, 0.5, 0.7, 0.8, 0.85, 0.9, 0.92, 0.94, 0.96, 0.97, 0.98, 0.99]
    t0 = time.time()
    with Pool(4) as pool:
        res = pool.map(run, grid, chunksize=1)
    json.dump(res, open(L.HERE + '/results/profile_L3_geo.json', 'w'), indent=1)
    print('done', time.time() - t0)
