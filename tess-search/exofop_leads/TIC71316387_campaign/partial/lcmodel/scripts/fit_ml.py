"""Maximum-likelihood fits + profile likelihood in L3 (and timing of the likelihood)."""
import numpy as np, sys, json, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/scripts')
import lclib as L, fitlib as F
from scipy.optimize import minimize, differential_evolution
d = F.prepare(L.default_set())
np.save(L.HERE + '/results/prep_info.npy', d['info'], allow_pickle=True)
print('N points', d['t'].size, 'events', len(d['evtab']))
p0 = dict(T0=L.T00, P=L.P0, rsum=0.2, k=0.7, cosi=0.12, J=0.45, L3=0.85, ecosw=0.0, esinw=0.0, q1a=0.4, q2a=0.3, q1b=0.4, q2b=0.3,
          s95=1.0, s105=1.0, s106=1.0)
post = F.Posterior(d, fixed={})
x0 = np.array([p0[n] for n in post.free])
t0 = time.time(); v = post(x0); print('lnpost %.1f, %.4f s per call' % (v, time.time() - t0))
t0 = time.time(); v = post(x0); print('lnpost %.1f, %.4f s per call' % (v, time.time() - t0))

def opt(post, x0, n=3):
    best = None
    rng = np.random.default_rng(1)
    for k in range(n):
        xs = x0 if k == 0 else x0 * (1 + 0.02 * rng.normal(size=x0.size)) 
        r = minimize(lambda x: -post(x), xs, method='Nelder-Mead', options=dict(maxiter=20000, maxfev=20000, xatol=1e-8, fatol=1e-3, adaptive=True))
        r = minimize(lambda x: -post(x), r.x, method='Powell', options=dict(maxiter=20000, xtol=1e-7, ftol=1e-7))
        if best is None or r.fun < best.fun: best = r
    return best

# profile in L3: fix L3 at grid values, optimise everything else (several starts in k/cosi/rsum)
out = []
grid = [0.0, 0.3, 0.5, 0.7, 0.8, 0.85, 0.9, 0.93, 0.95, 0.97, 0.98, 0.99]
starts = [dict(rsum=0.15, k=0.5, cosi=0.13), dict(rsum=0.25, k=0.8, cosi=0.2), dict(rsum=0.2, k=1.0, cosi=0.15), dict(rsum=0.35, k=0.6, cosi=0.3)]
for L3 in grid:
    best = None
    for s in starts:
        pp = dict(p0); pp.update(s); pp['L3'] = L3
        # choose cosi so the eclipse is grazing-ish
        postL = F.Posterior(d, fixed=dict(L3=L3))
        x0 = np.array([pp[n] for n in postL.free])
        if not np.isfinite(postL(x0)):
            pp['cosi'] = 0.6 * pp['rsum']; x0 = np.array([pp[n] for n in postL.free])
        r = opt(postL, x0, n=1)
        if best is None or r.fun < best[0].fun: best = (r, postL)
    r, postL = best
    p = postL.full(r.x)
    chi2 = -2 * F.log_like(d, p)
    out.append(dict(L3=L3, lnpost=-r.fun, chi2=chi2, **{k: float(v) for k, v in p.items()}))
    print(f"L3={L3:.2f} chi2={chi2:.1f} lnpost={-r.fun:.1f} rsum={p['rsum']:.4f} k={p['k']:.3f} cosi={p['cosi']:.4f} b/rsum={p['cosi']/p['rsum']:.3f} J={p['J']:.3f} ecosw={p['ecosw']:.5f} esinw={p['esinw']:.4f} s95={p['s95']:.3f} s105={p['s105']:.3f} s106={p['s106']:.3f}", flush=True)
json.dump(out, open(L.HERE + '/results/profile_L3.json', 'w'), indent=1)
