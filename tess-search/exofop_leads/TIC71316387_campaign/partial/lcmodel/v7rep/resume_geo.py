"""Resume the v6 geometric emcee chain (c03_mcmc.py setup: default data set, free = BASE + 5 depth scales,
kmax 3, uniform priors, DE/snooker moves) from its last state for a fixed wall-clock budget.  Works on a COPY
of v6/results/chain_geo.h5 (v7rep/chain_geoext.h5)."""
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
d = pickle.load(open(f'{V}/results/data_default.pkl', 'rb'))
M = Fi.Model(d, scale_keys=skeys)
st = json.load(open(f'{V}/results/best_default.json'))['p']
free = Fi.BASE + [f's_{k}' for k in skeys]
post = Fi.Posterior(M, st, free, kmax=3.0)
def lp(x): return post(x)
be = emcee.backends.HDFBackend(f'{R}/chain_geoext.h5')
print('resuming from iteration', be.iteration, flush=True)
json.dump(dict(free=free, fixed=st, dataset='default', skeys=skeys, kmax=3.0, note='v6 geo chain resumed by v7rep/resume_geo.py'),
          open(f'{R}/chain_geoext_meta.json', 'w'), indent=1)
t0 = time.time()
with Pool(nproc) as pool:
    sam = emcee.EnsembleSampler(be.shape[0], be.shape[1], lp, pool=pool, backend=be,
                                moves=[(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)])
    for i, _ in enumerate(sam.sample(be.get_last_sample(), iterations=100000)):
        if (i + 1) % 250 == 0:
            n = be.iteration
            ch = be.get_chain(discard=n // 2, flat=True); j = free.index('L3')
            try: tau = be.get_autocorr_time(discard=n // 2, tol=0)
            except Exception: tau = np.full(len(free), np.nan)
            print(f'iter {n} t={time.time()-t0:.0f}s max tau(2nd half)={np.nanmax(tau):.0f} L3 med {np.median(ch[:, j]):.4f} '
                  f'[{np.percentile(ch[:, j], 5):.4f},{np.percentile(ch[:, j], 95):.4f}] min {ch[:, j].min():.3f}', flush=True)
        if time.time() - t0 > budget: break
print('stopped at iteration', be.iteration, flush=True)
