"""emcee sampling of the geometric two-sphere model (uniform priors: L3 in [0,0.999), k in (0.08,kmax), b0 in [0,1),
J in (0,3), rsum in (0.03,0.8), |ecosw|<0.05, |esinw|<0.4; LD priors as in fit.py; free depth scales for the listed
product keys, uniform in (0.6,1.5)).
Usage: c03_mcmc.py <dataset> <start json> <tag> <nsteps> <nwalkers> <nproc> <kmax> <scale keys comma-separated|none>"""
import numpy as np, sys, pickle, json, time, os
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6/lib')
import fit as Fi
import emcee
from multiprocessing import Pool
V = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6'
dsname, startf, tag = sys.argv[1], sys.argv[2], sys.argv[3]
nsteps, nw, nproc = int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
kmax = float(sys.argv[7]); skeys = [] if sys.argv[8] == 'none' else sys.argv[8].split(',')
d = pickle.load(open(f'{V}/results/data_{dsname}.pkl', 'rb'))
M = Fi.Model(d, scale_keys=skeys)
st = json.load(open(startf))['p']
free = Fi.BASE + [f's_{k}' for k in skeys]
for k in skeys: st.setdefault(f's_{k}', 1.0)
post = Fi.Posterior(M, st, free, kmax=kmax)
def lp(x): return post(x)
x0 = np.array([st[n] for n in free])
assert np.isfinite(lp(x0)), 'start outside prior'
sc = dict(T0=2e-4, P=5e-7, rsum=0.005, k=0.02, b0=0.01, J=0.01, L3=0.005, ecosw=2e-4, esinw=0.005, usa=0.02, u1a=0.02, usb=0.02, u1b=0.02)
rng = np.random.default_rng(7)
p0 = []
while len(p0) < nw:
    x = x0 + rng.normal(size=x0.size) * np.array([sc.get(n, 0.01) for n in free])
    if np.isfinite(lp(x)): p0.append(x)
fn = f'{V}/results/chain_{tag}.h5'
if os.path.exists(fn): os.remove(fn)
backend = emcee.backends.HDFBackend(fn); backend.reset(nw, len(free))
t0 = time.time()
with Pool(nproc) as pool:
    sam = emcee.EnsembleSampler(nw, len(free), lp, pool=pool, backend=backend,
                                moves=[(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)])
    for i, _ in enumerate(sam.sample(p0, iterations=nsteps)):
        if (i + 1) % 500 == 0:
            try: tau = sam.get_autocorr_time(tol=0)
            except Exception: tau = np.full(len(free), np.nan)
            ch = sam.get_chain(discard=(i + 1) // 2, flat=True); lpb = sam.get_log_prob(discard=(i + 1) // 2, flat=True)
            j = free.index('L3')
            print(f'step {i+1} t={time.time()-t0:.0f}s acc={np.mean(sam.acceptance_fraction):.3f} max tau={np.nanmax(tau):.0f} '
                  f'L3 med {np.median(ch[:, j]):.3f} [{np.percentile(ch[:, j], 5):.3f},{np.percentile(ch[:, j], 95):.3f}] best lp {lpb.max():.2f}', flush=True)
json.dump(dict(free=free, fixed=st, dataset=dsname, skeys=skeys, kmax=kmax, nsteps=nsteps, nw=nw), open(f'{V}/results/chain_{tag}_meta.json', 'w'), indent=1)
print('done', time.time() - t0)
