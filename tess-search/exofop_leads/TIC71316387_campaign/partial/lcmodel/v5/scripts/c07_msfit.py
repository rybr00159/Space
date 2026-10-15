"""Main-sequence-constrained EB fit (lib/msfit.py): ML from several starts, then emcee.
mode 'free'  : L3 free (uniform) -> light-curve-only estimate of f_EB with MS-consistent stars; distance-free.
mode 'bound' : f_EB fixed by putting the MS pair at the G star's distance (dmag nuisance) -> bound hypothesis.
Usage: c07_msfit.py <dataset> <mode> <tag> <nsteps> <nwalkers> <nproc> [sig_r sig_j]"""
import os
os.environ['OMP_NUM_THREADS'] = '1'
import numpy as np, sys, pickle, json, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v5/lib')
import fit as Fi, opt, msfit, emcee
from multiprocessing import Pool
V = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v5'
dsname, mode, tag = sys.argv[1:4]; nsteps, nw, nproc = int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
sig_r = float(sys.argv[7]) if len(sys.argv) > 7 else msfit.SIG_R; sig_j = float(sys.argv[8]) if len(sys.argv) > 8 else msfit.SIG_J
d = pickle.load(open(f'{V}/results/data_{dsname}.pkl', 'rb'))
skeys = [k for k in d['info'] if k not in ('SPOC68',)] if dsname in ('default', 'qlp') else [k for k in d['info'] if not k.endswith('68')]
M = Fi.Model(d, scale_keys=skeys)
best = json.load(open(f'{V}/results/best_default.json'))['p'] if dsname == 'default' else d.get('truth', json.load(open(f'{V}/results/best_default.json'))['p'])
base = dict(T0=best['T0'], P=best['P'], ecosw=best['ecosw'], esinw=0.0, usa=0.64, u1a=0.3, usb=0.64, u1b=0.3, dR1=0.0, dR2=0.0, dJ=0.0)
for k in skeys: base[f's_{k}'] = best.get(f's_{k}', 1.0)
free = ['T0', 'P', 'M1', 'q', 'dR1', 'dR2', 'dJ', 'b0', 'L3' if mode == 'free' else 'dmag', 'ecosw', 'esinw', 'usa', 'u1a', 'usb', 'u1b'] + [f's_{k}' for k in skeys]
post = msfit.MSPosterior(M, base, free, mode=mode, sig_r=sig_r, sig_j=sig_j)
steps = dict(T0=2e-4, P=5e-7, M1=0.03, q=0.05, dR1=0.02, dR2=0.02, dJ=0.03, b0=0.03, L3=0.01, dmag=0.05, ecosw=3e-4, esinw=0.01, usa=0.03, u1a=0.05, usb=0.03, u1b=0.05)
starts = [dict(M1=0.6, q=0.65, b0=0.65, L3=0.95, dmag=0.0), dict(M1=0.45, q=0.75, b0=0.6, L3=0.97, dmag=0.0), dict(M1=0.75, q=0.55, b0=0.7, L3=0.93, dmag=0.0),
          dict(M1=0.3, q=0.8, b0=0.55, L3=0.985, dmag=0.0)]
mls = []
t0 = time.time()
mlfile = os.environ.get('MLSTART')
if mlfile:
    j = json.load(open(mlfile)); pj = dict(base); pj.update({k: v for k, v in j['p'].items() if k in free or k in base})
    if mode == 'bound' and 'dmag' not in j['p']: pj['dmag'] = 0.0
    starts = [{k: pj[k] for k in free}] + ([] if os.environ.get("MLONLY_SEED") else starts)
for s in starts:
    pp = dict(base); pp.update(s)
    x0 = [pp[n] for n in free]
    if not np.isfinite(post(x0)):
        # adjust b0 to make the start valid
        for b in np.linspace(0.3, 0.95, 14):
            pp['b0'] = b; x0 = [pp[n] for n in free]
            if np.isfinite(post(x0)): break
    if not np.isfinite(post(x0)): print('invalid start', s); continue
    x, fv = opt.optimise_lbfgs(post, x0, Fi.steps_for(free), Fi.bounds_for(free), n_cycles=2, nm_polish=1000)
    p = post.full(x); g = post.to_geom(p)
    mls.append((fv, p, g))
    print(f"start {s}: -lp {fv:.2f} chi2 {M.chi2(g):.2f} M1 {p['M1']:.3f} q {p['q']:.3f} b0 {p['b0']:.3f} L3 {g['L3']:.4f} rsum {g['rsum']:.4f} k {g['k']:.3f} J {g['J']:.3f} "
          f"dR {p['dR1']:+.3f} {p['dR2']:+.3f} dJ {p['dJ']:+.3f} esinw {p['esinw']:+.4f} t={time.time()-t0:.0f}s", flush=True)
mls.sort(key=lambda r: r[0])
fv, pbest, gbest = mls[0]
json.dump(dict(mode=mode, dataset=dsname, nlp=float(fv), chi2=float(M.chi2(gbest)), p={k: float(v) for k, v in pbest.items()}, geom={k: float(v) for k, v in gbest.items()},
               all=[dict(nlp=float(a), p={k: float(v) for k, v in b.items()}) for a, b, c in mls]), open(f'{V}/results/msml_{tag}.json', 'w'), indent=1)
if nsteps <= 0: sys.exit()
x0 = np.array([pbest[n] for n in free])
rng = np.random.default_rng(17)
p0 = []
while len(p0) < nw:
    x = x0 + rng.normal(size=x0.size) * np.array([steps.get(n, 0.01) * 0.2 for n in free])
    if np.isfinite(post(x)): p0.append(x)
fn = f'{V}/results/chain_{tag}.h5'
if os.path.exists(fn): os.remove(fn)
backend = emcee.backends.HDFBackend(fn); backend.reset(nw, len(free))
def lp(x): return post(x)
with Pool(nproc) as pool:
    sam = emcee.EnsembleSampler(nw, len(free), lp, pool=pool, backend=backend, moves=[(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)])
    for i, _ in enumerate(sam.sample(p0, iterations=nsteps)):
        if (i + 1) % 500 == 0:
            try: tau = sam.get_autocorr_time(tol=0)
            except Exception: tau = np.full(len(free), np.nan)
            ch = sam.get_chain(discard=(i + 1) // 2, flat=True)
            jm = free.index('M1'); jl = free.index('L3' if mode == 'free' else 'dmag')
            print(f'step {i+1} t={time.time()-t0:.0f}s acc={np.mean(sam.acceptance_fraction):.3f} max tau={np.nanmax(tau):.0f} M1 {np.median(ch[:, jm]):.3f} '
                  f'[{np.percentile(ch[:, jm], 5):.3f},{np.percentile(ch[:, jm], 95):.3f}] {free[jl]} {np.median(ch[:, jl]):.4f} [{np.percentile(ch[:, jl], 5):.4f},{np.percentile(ch[:, jl], 95):.4f}]', flush=True)
json.dump(dict(free=free, fixed=base, mode=mode, dataset=dsname, skeys=skeys, sig_r=sig_r, sig_j=sig_j), open(f'{V}/results/chain_{tag}_meta.json', 'w'), indent=1)
print('done')
