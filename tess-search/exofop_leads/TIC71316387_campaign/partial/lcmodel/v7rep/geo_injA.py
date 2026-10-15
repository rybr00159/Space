"""Uniform-prior geometric emcee chain (identical set-up to v6 c03_mcmc.py: free = BASE + 5 depth scales, kmax 3,
DE/snooker moves, 56 walkers) on INJECTION A (truth L3 0.972, k 0.46, b0 0.10, rsum 0.162, J 0.52), started at the
injection-A profile minimum (L3 = 0.965) as the real-data chain was started at the real-data best fit.  Purpose:
does a true L3 = 0.972 give a geometric posterior pulled to L3 ~ 0.93 like the real data (prior-volume effect)?"""
import os
os.environ['OMP_NUM_THREADS'] = '1'
import numpy as np, sys, pickle, json, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6/lib')
import fit as Fi, emcee
from multiprocessing import Pool
V = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6'
R = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v7rep'
budget = float(sys.argv[1]); nproc = int(sys.argv[2])
skeys = ['SPOC1', 'SPOC28', 'TC95', 'TC105', 'TC106']
d = pickle.load(open(f'{V}/results/data_injA.pkl', 'rb'))
M = Fi.Model(d, scale_keys=skeys)
prof = json.load(open(f'{V}/results/profile_L3_injA.json'))['profile']
st = dict(min(prof, key=lambda r: r['chi2'])['p'])
for k in skeys: st.setdefault(f's_{k}', 1.0)
free = Fi.BASE + [f's_{k}' for k in skeys]
post = Fi.Posterior(M, st, free, kmax=3.0)
def lp(x): return post(x)
x0 = np.array([st[n] for n in free]); assert np.isfinite(lp(x0))
sc = dict(T0=2e-4, P=5e-7, rsum=0.005, k=0.02, b0=0.01, J=0.01, L3=0.005, ecosw=2e-4, esinw=0.005, usa=0.02, u1a=0.02, usb=0.02, u1b=0.02)
rng = np.random.default_rng(7); p0 = []
while len(p0) < 56:
    x = x0 + rng.normal(size=x0.size) * np.array([sc.get(n, 0.01) for n in free])
    if np.isfinite(lp(x)): p0.append(x)
fn = f'{R}/chain_geoinjA.h5'
if os.path.exists(fn): os.remove(fn)
be = emcee.backends.HDFBackend(fn); be.reset(56, len(free))
json.dump(dict(free=free, fixed=st, dataset='injA', skeys=skeys, kmax=3.0, chain=fn, truth=d['truth']), open(f'{R}/chain_geoinjA_meta.json', 'w'), indent=1)
t0 = time.time()
with Pool(nproc) as pool:
    sam = emcee.EnsembleSampler(56, len(free), lp, pool=pool, backend=be, moves=[(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)])
    for i, _ in enumerate(sam.sample(p0, iterations=100000)):
        if (i + 1) % 500 == 0:
            n = be.iteration; ch = be.get_chain(discard=n // 2, flat=True); j = free.index('L3')
            try: tau = be.get_autocorr_time(discard=n // 2, tol=0)
            except Exception: tau = np.full(len(free), np.nan)
            print(f'iter {n} t={time.time()-t0:.0f}s max tau={np.nanmax(tau):.0f} L3 med {np.median(ch[:, j]):.4f} [{np.percentile(ch[:, j], 2.5):.4f},{np.percentile(ch[:, j], 97.5):.4f}] '
                  f'frac(L3>=0.972) {np.mean(ch[:, j] >= 0.972):.4f}', flush=True)
        if time.time() - t0 > budget: break
print('stopped at iteration', be.iteration, flush=True)
