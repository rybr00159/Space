"""Targeted check: L3 profile of the geometric model with e sin(omega) FIXED at 0 (physically expected for a 1.3-d
M-dwarf binary; e cos(omega) left free because it is set by the secondary-eclipse timing, which also contains the
~8 s Roemer delay).  Same likelihood, priors and optimiser as v6 c04_profile.py; seeds = the real-data free-e profile
solutions at the 3 nearest L3 grid values (esinw reset to 0) and the e0 variant.  Usage: prof_esin0.py <L3>"""
import os
os.environ['OMP_NUM_THREADS'] = '1'
import numpy as np, sys, pickle, json, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6/lib')
import fit as Fi, opt
V = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6'
R = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v7rep'
L3v = float(sys.argv[1])
skeys = ['SPOC1', 'SPOC28', 'TC95', 'TC105', 'TC106']
d = pickle.load(open(f'{V}/results/data_default.pkl', 'rb'))
M = Fi.Model(d, scale_keys=skeys)
prof = json.load(open(f'{V}/results/profile_L3_default.json'))['profile']
prof = sorted(prof, key=lambda r: abs(r['val'] - L3v))[:3]
seeds = [dict(r['p']) for r in prof] + [dict(json.load(open(f'{V}/results/variants_a.json'))['e0']['p'])]
free = [n for n in Fi.BASE + [f's_{k}' for k in skeys] if n not in ('L3', 'esinw')]
best = None; t0 = time.time()
for s in seeds:
    s = dict(s); s['L3'] = L3v; s['esinw'] = 0.0
    post = Fi.Posterior(M, s, free)
    x0 = [s[n] for n in free]
    if not np.isfinite(post(x0)): continue
    x, fv = opt.optimise_lbfgs(post, x0, Fi.steps_for(free), Fi.bounds_for(free), n_cycles=2, nm_polish=800)
    p = post.full(x); c = M.chi2(p)
    print(f'  seed -> nlp {fv:.2f} chi2 {c:.2f} k {p["k"]:.3f} b0 {p["b0"]:.3f} rsum {p["rsum"]:.4f} J {p["J"]:.3f} u1a {p["u1a"]:.2f} u1b {p["u1b"]:.2f} t={time.time()-t0:.0f}s', flush=True)
    if best is None or fv < best[0]: best = (fv, p, c)
fv, p, c = best
json.dump(dict(L3=L3v, nlp=float(fv), chi2=float(c), p={k: float(v) for k, v in p.items()}), open(f'{R}/prof_esin0_L3_{L3v}.json', 'w'), indent=1)
print(f'L3={L3v}: best nlp {fv:.2f} chi2 {c:.2f}', flush=True)
