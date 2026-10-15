"""Robustness variants of the geometric ML fit (scaled L-BFGS-B, seeds = best-fit families):
  e0        : circular orbit (ecosw = esinw = 0)
  esin0     : esinw = 0 only
  ldclaret  : LD fixed to Claret (2017) PHOENIX-COND TESS quadratic (LSM, logg 5), star 1 3550 K (0.164,0.426), star 2 3150 K (0.157,0.488)
  lin       : linear (not quadratic) per-event baselines
  noscale   : all depth scales fixed to 1
  tcall     : own TESScut photometry for all six sectors
  qlp       : QLP SAP for S95/S105 instead of own photometry
  early     : S1+S28 only ; late: S68+S95+S105+S106 only
Usage: c08_variants.py <variant> [variant ...]"""
import os
os.environ['OMP_NUM_THREADS'] = '1'
import numpy as np, sys, pickle, json, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6/lib')
import fit as Fi, opt
V = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6'
seeds = json.load(open(f'{V}/results/seeds_default.json'))['seeds']
out = {}
fn = os.environ.get('VFILE', f'{V}/results/variants.json')
if os.path.exists(fn): out = json.load(open(fn))
def subset(d, keys):
    m = np.isin(d['key'], keys)
    ev_keep = np.unique(d['ev'][m])
    dd = {k: v[m] for k, v in d.items() if isinstance(v, np.ndarray) and v.shape == d['t'].shape}
    # renumber events
    remap = {e: i for i, e in enumerate(ev_keep)}
    dd['ev'] = np.array([remap[e] for e in dd['ev']])
    dd['evtab'] = [d['evtab'][e] for e in ev_keep]
    dd['info'] = {k: v for k, v in d['info'].items() if k in keys}
    ev = dd['ev']; dd['starts'] = np.r_[0, np.where(np.diff(ev) != 0)[0] + 1]; dd['ends'] = np.r_[dd['starts'][1:], ev.size]
    return dd
for var in sys.argv[1:]:
    t0 = time.time()
    dsn = 'tcall' if var == 'tcall' else 'qlp' if var == 'qlp' else 'default'
    d = pickle.load(open(f'{V}/results/data_{dsn}.pkl', 'rb'))
    if var == 'early': d = subset(d, ['SPOC1', 'SPOC28'])
    if var == 'late': d = subset(d, ['SPOC68', 'TC95', 'TC105', 'TC106'])
    ref = {'tcall': 'TC68', 'early': 'SPOC28', 'late': 'SPOC68'}.get(var, 'SPOC68')
    skeys = [] if var == 'noscale' else [k for k in d['info'] if k != ref]
    M = Fi.Model(d, scale_keys=skeys, order=1 if var == 'lin' else 2)
    fixed_over = {}
    if var == 'e0': fixed_over = dict(ecosw=0.0, esinw=0.0)
    if var == 'esin0': fixed_over = dict(esinw=0.0)
    if var == 'ldclaret': fixed_over = dict(u1a=0.164, usa=0.590, u1b=0.157, usb=0.645)
    free = [n for n in Fi.BASE + [f's_{k}' for k in skeys] if n not in fixed_over]
    best = None
    for s in seeds:
        s = dict(s); s.update(fixed_over)
        for k in skeys: s.setdefault(f's_{k}', 1.0)
        if var in ('early',): s['s_SPOC1'] = 1.0
        post = Fi.Posterior(M, s, free)
        x0 = [s[n] for n in free]
        if not np.isfinite(post(x0)): continue
        x, fv = opt.optimise_lbfgs(post, x0, Fi.steps_for(free), Fi.bounds_for(free), n_cycles=2, nm_polish=800)
        if best is None or fv < best[0]: best = (fv, post.full(x))
    fv, p = best
    out[var] = dict(nlp=float(fv), chi2=float(M.chi2(p)), npts=int(d['t'].size), p={k: float(v) for k, v in p.items()})
    print(f"{var:9s} chi2 {out[var]['chi2']:.2f} (N={d['t'].size}) L3 {p['L3']:.4f} rsum {p['rsum']:.4f} k {p['k']:.3f} b0 {p['b0']:.3f} J {p['J']:.3f} ecosw {p['ecosw']:+.5f} esinw {p['esinw']:+.4f} "
          f"u1a {p['u1a']:.2f} u1b {p['u1b']:.2f} " + ' '.join(f"{k}={v:.3f}" for k, v in p.items() if k.startswith('s_')) + f" t={time.time()-t0:.0f}s", flush=True)
    json.dump(out, open(fn, 'w'), indent=1)
