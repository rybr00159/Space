"""First maximum-posterior fits with free depth scales for every product except the reference.
Stage 1: LD, T0, P, ecosw fixed; stage 2: everything free.  Several starts on both k branches.
Usage: c02_first_fit.py <dataset> <reference key> <out tag>"""
import os
os.environ['OMP_NUM_THREADS'] = '1'
import numpy as np, sys, pickle, json, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v5/lib')
import fit as Fi, opt
V = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v5'
name = sys.argv[1]; ref = sys.argv[2]; tag = sys.argv[3]
d = pickle.load(open(f'{V}/results/data_{name}.pkl', 'rb'))
keys = [k for k in d['info'] if k != ref] if ref != 'none' else []
M = Fi.Model(d, scale_keys=keys)
p0 = dict(T0=3165.6697, P=1.3036712, rsum=0.20, k=0.6, b0=0.75, J=0.47, L3=0.90, ecosw=0.0, esinw=0.0, usa=0.64, u1a=0.30, usb=0.64, u1b=0.25)
for k in keys: p0[f's_{k}'] = 1.0
t = time.time(); [M.chi2(p0) for _ in range(10)]; print('chi2 call', (time.time() - t) / 10 * 1e3, 'ms; npts', d['t'].size, flush=True)
f1 = ['rsum', 'k', 'b0', 'J', 'L3', 'esinw'] + [f's_{k}' for k in keys]
f2 = Fi.BASE + [f's_{k}' for k in keys]
steps = dict(T0=3e-4, P=1e-6, rsum=0.02, k=0.08, b0=0.05, J=0.04, L3=0.02, ecosw=0.0005, esinw=0.02, usa=0.03, u1a=0.05, usb=0.03, u1b=0.05)
results = []
for L3s, ks, b0s, rs in ((0.9, 0.6, 0.75, 0.2), (0.6, 0.3, 0.9, 0.25), (0.97, 0.45, 0.3, 0.13), (0.95, 1.6, 0.6, 0.2)):
    pp = dict(p0); pp.update(L3=L3s, k=ks, b0=b0s, rsum=rs)
    post = Fi.Posterior(M, pp, f1)
    x, fval = opt.optimise(post, [pp[n] for n in f1], [steps.get(n, 0.02) for n in f1], n_restart=2, maxiter=3000)
    pp = post.full(x)
    post = Fi.Posterior(M, pp, f2)
    x, fval = opt.optimise(post, [pp[n] for n in f2], [steps.get(n, 0.02) for n in f2], n_restart=2, maxiter=5000)
    p = post.full(x)
    print(f'start L3={L3s} k={ks}: -logpost {fval:.2f} chi2 {M.chi2(p):.2f} ' + ' '.join(f'{n}={p[n]:.5g}' for n in f2), flush=True)
    results.append(dict(start=[L3s, ks, b0s, rs], nlp=float(fval), chi2=float(M.chi2(p)), p={n: float(p[n]) for n in p}))
results.sort(key=lambda r: r['nlp'])
json.dump(dict(dataset=name, ref=ref, npts=int(d['t'].size), results=results), open(f'{V}/results/first_fit_{tag}.json', 'w'), indent=1)
