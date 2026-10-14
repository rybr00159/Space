"""Main-sequence-constrained EB model (distance independent).
Free: M1, q = M2/M1, fR (common radius-inflation factor), dlnJ (log offset of the surface-brightness ratio), bn, L3,
ecosw, esinw, LD (q1,q2 per star), dT0, dP, sector depth scales.
Derived with the smooth Pecaut & Mamajek relations (ms_smooth): R_i = fR R_MS(M_i); a from Kepler's third law with M1+M2
and P; rsum = (R1+R2)/a; k = R2/R1; J = SB_T(M2)/SB_T(M1) exp(dlnJ).
Priors: fR ~ N(FR_MU, FR_SD) truncated to [0.9, 1.3]; dlnJ ~ N(0, SIGJ); LD sums as in the geometric fit; uniform in M1, M2."""
import numpy as np, sys
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4/lib')
import fit3, lclib as L, ms_smooth as MS
from scipy.optimize import least_squares
GM_SUN = 1.32712440018e20; RSUN = 6.957e8; DAY = 86400.0
def a_over_R(Mtot, P):
    return (GM_SUN * Mtot * (P * DAY) ** 2 / (4 * np.pi ** 2)) ** (1 / 3) / RSUN
CFG = dict(FR_MU=1.04, FR_SD=0.05, SIGJ=0.15)
def ms_derived(p):
    M1 = p['M1']; M2 = p['M1'] * p['q']
    R1 = MS.radius(M1) * p['fR']; R2 = MS.radius(M2) * p['fR']
    a = a_over_R(M1 + M2, p['P'])
    J = MS.SB_T(M2) / MS.SB_T(M1) * np.exp(p['dlnJ'])
    return dict(rsum=(R1 + R2) / a, k=R2 / R1, J=J, R1=R1, R2=R2, a=a, M2=M2)
LB = dict(dT0=-30, dP=-50, M1=0.09, q=0.05, fR=0.9, dlnJ=-0.8, x=0.0, L3=0.0, ecosw=-0.05, esinw=-0.3, q1a=1e-3, q2a=1e-3, q1b=1e-3, q2b=1e-3)
UB = dict(dT0=30, dP=50, M1=1.25, q=1.0, fR=1.3, dlnJ=0.8, x=0.9999, L3=0.998, ecosw=0.05, esinw=0.3, q1a=0.999, q2a=0.999, q1b=0.999, q2b=0.999)
XS = dict(dT0=0.1, dP=0.5, M1=0.02, q=0.02, fR=0.01, dlnJ=0.02, x=0.02, L3=0.003, ecosw=1e-4, esinw=0.01, q1a=0.05, q2a=0.05, q1b=0.05, q2b=0.05)
for s in fit3.SECS:
    LB[f's{s}'] = 0.6; UB[f's{s}'] = 1.4; XS[f's{s}'] = 0.02
DEF = dict(dT0=0.0, dP=0.0, M1=0.45, q=0.45, fR=1.04, dlnJ=0.0, x=0.3, L3=0.97, ecosw=0.0009, esinw=0.0, q1a=0.39, q2a=0.25, q1b=0.39, q2b=0.3,
           s95=0.95, s105=0.94, s106=0.97)
def to_p(names, xv, fixed, T0ref):
    p = dict(fixed); p.update(dict(zip(names, xv)))
    p['T0'] = T0ref + p['dT0'] / 1440.0; p['P'] = L.P0 + p['dP'] * 1e-6
    p.update(ms_derived(p))
    for s in fit3.SECS: p.setdefault(f's{s}', 1.0)
    p['bn'] = p['x']; p['cosi'] = p['bn'] * p['rsum']
    return p
def prior_resid(p):
    pr = list(fit3.prior_resid(p, 'geo'))
    pr += [(p['fR'] - CFG['FR_MU']) / CFG['FR_SD'], p['dlnJ'] / CFG['SIGJ']]
    return np.array(pr)
def solve(d, start, fixed=None, T0ref=L.T00, scales=('s95', 's105', 's106'), max_nfev=600):
    fixed = dict(fixed or {})
    base = ['dT0', 'dP', 'M1', 'q', 'fR', 'dlnJ', 'x', 'L3', 'ecosw', 'esinw', 'q1a', 'q2a', 'q1b', 'q2b'] + list(scales)
    names = [n for n in base if n not in fixed]
    st = dict(DEF); st.update(start)
    lb = np.array([LB[n] for n in names]); ub = np.array([UB[n] for n in names])
    if 'M1' in fixed:
        i = names.index('q'); lb[i] = max(lb[i], 0.09 / fixed['M1'])
    x0 = np.clip(np.array([st[n] for n in names], float), lb + 1e-6, ub - 1e-6)
    def rf(xv):
        p = to_p(names, xv, fixed, T0ref)
        pen = [min(0.0, p['M2'] - 0.08) / 0.002]   # soft floor at the hydrogen-burning limit
        return np.concatenate([fit3.resid(d, p), prior_resid(p), pen])
    r = least_squares(rf, x0, bounds=(lb, ub), x_scale=np.array([XS[n] for n in names]), method='trf', diff_step=1e-4,
                      max_nfev=max_nfev, ftol=1e-10, xtol=1e-10, gtol=1e-10)
    p = to_p(names, r.x, fixed, T0ref); pr = prior_resid(p)
    r1 = p['rsum'] / (1 + p['k'])
    return dict(tot=float(2 * r.cost), chi2=fit3.chi2(d, p), prior2=float(pr @ pr), r1=r1, b1=p['cosi'] / r1,
                total=bool(p['bn'] <= abs(1 - p['k']) / (1 + p['k'])), **{k: float(v) for k, v in p.items()})
