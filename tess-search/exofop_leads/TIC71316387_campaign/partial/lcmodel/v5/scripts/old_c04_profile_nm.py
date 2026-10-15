"""Profile likelihood in third light: for each fixed L3 maximise the posterior over all other parameters, from
starts on both branches (k < 1: star 1 larger; k > 1: star 2 larger), each start seeded from the neighbouring
grid point's solution and from fresh guesses.  Also profiles in k and in rsum.
Usage: c04_profile.py <dataset> <start json> <tag> <scale keys comma|none> <which: L3|k|rsum>"""
import os
os.environ['OMP_NUM_THREADS'] = '1'
import numpy as np, sys, pickle, json, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v5/lib')
import fit as Fi, opt
V = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v5'
name, startf, tag, sk, which = sys.argv[1:6]
skeys = [] if sk == 'none' else sk.split(',')
d = pickle.load(open(f'{V}/results/data_{name}.pkl', 'rb'))
M = Fi.Model(d, scale_keys=skeys)
st = json.load(open(startf))
st = st['results'][0]['p'] if 'results' in st else st['p']
for k in skeys: st.setdefault(f's_{k}', 1.0)
allfree = Fi.BASE + [f's_{k}' for k in skeys]
steps = dict(T0=2e-4, P=5e-7, rsum=0.01, k=0.05, b0=0.03, J=0.03, L3=0.01, ecosw=0.0003, esinw=0.01, usa=0.03, u1a=0.05, usb=0.03, u1b=0.05)
grids = dict(L3=[0.0, 0.3, 0.5, 0.65, 0.75, 0.82, 0.87, 0.9, 0.92, 0.94, 0.955, 0.965, 0.975, 0.983, 0.99],
             k=[0.10, 0.13, 0.16, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0],
             rsum=[0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.25, 0.28, 0.32, 0.36, 0.40])
grid = grids[which]
free = [n for n in allfree if n != which]
fresh = [dict(k=0.5, b0=0.8, rsum=0.2, L3=0.9, J=0.47), dict(k=1.5, b0=0.7, rsum=0.2, L3=0.9, J=0.5)]
out = []
prev = {}
t0 = time.time()
for g in grid:
    cands = []
    seeds = ([dict(st)] + [dict(st, **f) for f in fresh]) if not prev else ([dict(st)] + list(prev.values()))
    for s in seeds:
        s = dict(s); s[which] = g
        if which == 'L3': pass
        # make the start valid
        post = Fi.Posterior(M, s, free)
        x0 = [s[n] for n in free]
        if not np.isfinite(post(x0)):
            continue
        try:
            x, fv = opt.optimise(post, x0, [steps.get(n, 0.02) for n in free], n_restart=2, maxiter=2500)
        except ValueError:
            continue
        p = post.full(x)
        cands.append((fv, p))
    cands.sort(key=lambda c: c[0])
    best_lo = [c for c in cands if c[1]['k'] < 1]; best_hi = [c for c in cands if c[1]['k'] >= 1]
    rec = dict(val=g)
    for lab, cc in (('klt1', best_lo), ('kge1', best_hi)):
        if cc:
            fv, p = cc[0]
            rec[lab] = dict(nlp=float(fv), chi2=float(M.chi2(p)), p={n: float(p[n]) for n in p})
            prev[lab] = p
    out.append(rec)
    msg = f'{which}={g}: ' + ' | '.join(f"{lab} nlp {rec[lab]['nlp']:.2f} L3 {rec[lab]['p']['L3']:.3f} rsum {rec[lab]['p']['rsum']:.4f} k {rec[lab]['p']['k']:.3f} b0 {rec[lab]['p']['b0']:.3f} J {rec[lab]['p']['J']:.3f}" for lab in ('klt1', 'kge1') if lab in rec)
    print(msg, f't={time.time()-t0:.0f}s', flush=True)
    json.dump(dict(which=which, dataset=name, skeys=skeys, profile=out), open(f'{V}/results/profile_{which}_{tag}.json', 'w'), indent=1)
