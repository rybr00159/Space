"""Geometric EB fitting helpers (maximum likelihood with scipy least_squares).
Parameters: dT0 [min], dP [1e-6 d], rsum=(R1+R2)/a (or r1=R1/a when 'r1' is given/fixed), k=R2/R1, x (impact coordinate,
see BRANCH), J (disc-averaged TESS surface-brightness ratio star2/star1), L3, ecosw, esinw, LD (Kipping q1,q2 per star),
per-sector depth scales s<sec> (s=1 for sectors not listed).
Impact parameter bn = cos i / rsum (0..1).  The eclipse is TOTAL/ANNULAR when bn <= |1-k|/(1+k), PARTIAL otherwise.
  branch 'any'     : bn = x
  branch 'total'   : bn = x * |1-k|/(1+k)
  branch 'partial' : bn = |1-k|/(1+k) + x * (1 - |1-k|/(1+k))
Priors (as residuals): u1+u2 ~ N(0.63, 0.10) for each star (Claret 2017 TESS, 3000-6000 K dwarfs span 0.58-0.74)."""
import numpy as np, sys
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/final/lib')
import fit3, lclib as L
from scipy.optimize import least_squares
SECS = fit3.SECS
LB = dict(dT0=-30, dP=-50, rsum=0.03, r1=0.02, k=0.05, x=0.0, J=0.005, L3=0.0, ecosw=-0.05, esinw=-0.3, q1a=1e-3, q2a=1e-3, q1b=1e-3, q2b=1e-3)
UB = dict(dT0=30, dP=50, rsum=0.7, r1=0.6, k=2.5, x=0.9999, J=2.0, L3=0.998, ecosw=0.05, esinw=0.3, q1a=0.999, q2a=0.999, q1b=0.999, q2b=0.999)
XS = dict(dT0=0.1, dP=0.5, rsum=0.01, r1=0.005, k=0.05, x=0.02, J=0.02, L3=0.005, ecosw=1e-4, esinw=0.01, q1a=0.05, q2a=0.05, q1b=0.05, q2b=0.05)
for s in SECS:
    LB[f's{s}'] = 0.6; UB[f's{s}'] = 1.4; XS[f's{s}'] = 0.02
DEF = dict(dT0=0.0, dP=0.0, rsum=0.17, k=0.45, x=0.3, J=0.5, L3=0.97, ecosw=0.0009, esinw=0.0, q1a=0.39, q2a=0.25, q1b=0.39, q2b=0.3,
           s95=0.95, s105=0.94, s106=0.97, s1=1.0, s28=1.0, s68=1.0)

def bn_of(x, k, branch):
    c = abs(1 - k) / (1 + k)
    if branch == 'any': return x
    if branch == 'total': return x * c
    return c + x * (1 - c)

def to_p(names, xv, fixed, T0ref, branch):
    p = dict(fixed); p.update(dict(zip(names, xv)))
    p['T0'] = T0ref + p['dT0'] / 1440.0; p['P'] = L.P0 + p['dP'] * 1e-6
    if 'r1' in p: p['rsum'] = p['r1'] * (1 + p['k'])
    for s in SECS: p.setdefault(f's{s}', 1.0)
    p['bn'] = bn_of(p['x'], p['k'], branch)
    p['cosi'] = p['bn'] * p['rsum']
    return p

def derived(p):
    r1 = p['rsum'] / (1 + p['k']); r2 = p['rsum'] - r1
    b1 = p['cosi'] / r1
    return dict(r1=r1, r2=r2, b1=b1, inc=float(np.degrees(np.arccos(min(1, p['cosi'])))),
                total=bool(p['bn'] <= abs(1 - p['k']) / (1 + p['k'])), fEB=1 - p['L3'])

def solve(d, start, fixed=None, branch='any', T0ref=L.T00, scales=('s95', 's105', 's106'), use_r1=False, max_nfev=600, ld_prior=True):
    fixed = dict(fixed or {})
    base = ['dT0', 'dP', 'r1' if use_r1 else 'rsum', 'k', 'x', 'J', 'L3', 'ecosw', 'esinw', 'q1a', 'q2a', 'q1b', 'q2b'] + list(scales)
    names = [n for n in base if n not in fixed]
    st = dict(DEF); st.update(start)
    if use_r1 and 'r1' not in st: st['r1'] = st['rsum'] / (1 + st['k'])
    lb = np.array([LB[n] for n in names]); ub = np.array([UB[n] for n in names])
    x0 = np.clip(np.array([st[n] for n in names], float), lb + 1e-6, ub - 1e-6)
    def rf(xv):
        p = to_p(names, xv, fixed, T0ref, branch)
        pr = fit3.prior_resid(p, 'geo') if ld_prior else np.zeros(0)
        return np.concatenate([fit3.resid(d, p), pr])
    r = least_squares(rf, x0, bounds=(lb, ub), x_scale=np.array([XS[n] for n in names]), method='trf', diff_step=1e-4,
                      max_nfev=max_nfev, ftol=1e-10, xtol=1e-10, gtol=1e-10)
    p = to_p(names, r.x, fixed, T0ref, branch)
    out = dict(tot=float(2 * r.cost), chi2=fit3.chi2(d, p), branch=branch, nfev=int(r.nfev), **{k: float(v) for k, v in p.items()})
    out.update(derived(p))
    return out

STARTS_TOTAL = [dict(rsum=0.165, k=0.49, x=0.1, L3=0.975, J=0.5), dict(rsum=0.17, k=0.42, x=0.6, L3=0.965, J=0.5),
                dict(rsum=0.19, k=0.25, x=0.8, L3=0.90, J=0.45), dict(rsum=0.16, k=0.6, x=0.3, L3=0.985, J=0.5),
                dict(rsum=0.165, k=1.8, x=0.3, L3=0.975, J=0.25)]
STARTS_PARTIAL = [dict(rsum=0.25, k=0.1, x=0.5, L3=0.0, J=0.3), dict(rsum=0.22, k=0.71, x=0.4, L3=0.90, J=0.48),
                  dict(rsum=0.255, k=0.79, x=0.5, L3=0.91, J=0.42), dict(rsum=0.3, k=0.9, x=0.6, L3=0.7, J=0.5),
                  dict(rsum=0.21, k=0.3, x=0.3, L3=0.8, J=0.4), dict(rsum=0.4, k=1.0, x=0.8, L3=0.3, J=0.5)]
