"""Profile likelihood (maximum posterior over all other parameters) in k, L3 or rsum, with the scaled L-BFGS-B
optimiser (lib/opt.py).  Seeds at every grid point: the previous grid point's solution (continuation, run in both
directions from the best-fit value) and a fixed list of good solutions from earlier fits (seed file).
Usage: c04_profile.py <dataset> <seed json list> <tag> <scale keys comma|none> <which: L3|k|rsum> [grid comma]"""
import os
os.environ['OMP_NUM_THREADS'] = '1'
import numpy as np, sys, pickle, json, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6/lib')
import fit as Fi, opt
V = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6'
name, seedf, tag, sk, which = sys.argv[1:6]
skeys = [] if sk == 'none' else sk.split(',')
d = pickle.load(open(f'{V}/results/data_{name}.pkl', 'rb'))
M = Fi.Model(d, scale_keys=skeys)
seeds0 = json.load(open(seedf))['seeds']
for s in seeds0:
    for k in skeys: s.setdefault(f's_{k}', 1.0)
allfree = Fi.BASE + [f's_{k}' for k in skeys]
grids = dict(L3=[0.0, 0.3, 0.5, 0.65, 0.75, 0.82, 0.87, 0.9, 0.92, 0.94, 0.955, 0.965, 0.975, 0.983, 0.99, 0.995],
             k=[0.10, 0.13, 0.16, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0],
             rsum=[0.12, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.22, 0.25, 0.28, 0.32])
grid = [float(x) for x in sys.argv[6].split(',')] if len(sys.argv) > 6 else grids[which]
free = [n for n in allfree if n != which]
lo, hi = Fi.BOUNDS[which]
def solve(g, seeds):
    best = None
    for s in seeds:
        s = dict(s); s[which] = min(max(g, lo), hi)
        post = Fi.Posterior(M, s, free)
        x0 = [s[n] for n in free]
        if not np.isfinite(post(x0)): continue
        try:
            x, fv = opt.optimise_lbfgs(post, x0, Fi.steps_for(free), Fi.bounds_for(free), n_cycles=2, nm_polish=800)
        except ValueError:
            continue
        if best is None or fv < best[0]: best = (fv, post.full(x))
    return best
bestseed = min(seeds0, key=lambda s: M.chi2(s))
order = np.argsort(np.abs(np.array(grid) - bestseed[which]))
i0 = int(order[0])
res = {}
t0 = time.time()
for direction in (range(i0, len(grid)), range(i0 - 1, -1, -1)):
    prev = None
    for i in direction:
        g = grid[i]
        extra = [prev] if prev is not None else ([res[grid[i0]]['p']] if grid[i0] in res else [])
        b = solve(g, list(seeds0) + extra)
        if b is None: print(f'{which}={g}: no valid seed', flush=True); continue
        fv, p = b; prev = p
        res[g] = dict(val=g, nlp=float(fv), chi2=float(M.chi2(p)), p={n: float(v) for n, v in p.items()})
        print(f"{which}={g}: nlp {fv:.2f} chi2 {res[g]['chi2']:.2f} L3 {p['L3']:.4f} rsum {p['rsum']:.4f} k {p['k']:.3f} b0 {p['b0']:.3f} J {p['J']:.3f} "
              f"esinw {p['esinw']:+.4f} u1a {p['u1a']:.2f} u1b {p['u1b']:.2f} t={time.time()-t0:.0f}s", flush=True)
        json.dump(dict(which=which, dataset=name, skeys=skeys, profile=[res[k] for k in sorted(res)]), open(f'{V}/results/profile_{which}_{tag}.json', 'w'), indent=1)
