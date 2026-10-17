"""fitlib: constrained ML fits of the transit shape (bd/). Parameters x = [dT0 (d), ln k, beta_b (b = beta_b (1 + k)),
ln rho_c, u1, u2, ln D]. Scenario constraints:
  'planet'  : D = 1, k in [0.03, 0.30]                    (planet or any small body on the M3V)
  'grazBD'  : D = 1, k in [KG_LO, KG_HI] = [0.19, 0.33]   (brown dwarf on the M3V: R_BD 0.72-1.25 R_J, R* 0.383)
  'compBD'  : D in [0.008, 0.45], k in [0.25, 1.2]          (BD / giant planet on a bound 0.08-0.3 Msun companion)
rho_c in [0.3, 300] g/cm3 (eccentricity free); |dT0| < 0.01 d; Gaussian LD prior (U1, U2, SU) of the M3V."""
import numpy as np
from scipy.optimize import minimize
from blib import *
from lclib import cad_curves
KG_LO, KG_HI = 0.19, 0.33
LIM = dict(planet=dict(k=(0.03, 0.30), D=(1.0, 1.0)), grazBD=dict(k=(KG_LO, KG_HI), D=(1.0, 1.0)),
           compBD=dict(k=(0.25, 1.2), D=(0.008, 0.45)))

def unpack(x, scen):
    dT0, lk, bb, lr, u1, u2 = x[:6]
    lD = x[6] if len(x) > 6 else 0.
    return dT0, np.exp(lk), bb, np.exp(lr), u1, u2, np.exp(lD)

def model_s(d, x, scen):
    dT0, k, bb, rho, u1, u2, D = unpack(x, scen)
    b = bb*(1 + k)
    return D*d.template(cad_curves(k, b, float(aR_of_rho(rho)), u1, u2), dT0)

def chi2(d, x, scen):
    dT0, k, bb, rho, u1, u2, D = unpack(x, scen)
    L = LIM[scen]
    if not (L['k'][0] <= k <= L['k'][1]) or not (0 <= bb < 1) or not (0.3 <= rho <= 300) or abs(dT0) > 0.01: return 1e12
    if L['D'][0] < L['D'][1] and not (L['D'][0] <= D <= L['D'][1]): return 1e12
    if not (u1 + u2 < 1 and u1 > 0 and u1 + 2*u2 > 0): return 1e12
    s = model_s(d, x, scen)
    r = d.y - s
    return d.Q(r, r) + ((u1 - U1)/SU)**2 + ((u2 - U2)/SU)**2

def starts(scen, dT0=0.):
    if scen == 'planet':
        return [[dT0, np.log(0.12), 0.45, np.log(9.), U1, U2], [dT0, np.log(0.13), 0.75, np.log(4.), U1, U2],
                [dT0, np.log(0.115), 0.1, np.log(15.), U1, U2]]
    if scen == 'grazBD':
        return [[dT0, np.log(k), bb, np.log(r), U1, U2] for k, bb, r in ((0.21, 0.88, 2.), (0.26, 0.86, 1.5), (0.31, 0.85, 1.0), (0.23, 0.9, 4.))]
    if scen == 'compBD':
        return [[dT0, np.log(k), bb, np.log(r), U1, U2, np.log(0.0157/k**2*1.15)] for k, bb, r in ((0.3, 0.3, 8.), (0.4, 0.2, 10.), (0.55, 0.3, 15.), (0.8, 0.2, 20.), (0.35, 0.6, 5.))]

def fit(d, scen, dT0=0., extra_starts=()):
    f = lambda x: chi2(d, x, scen)
    best = None
    for s in list(starts(scen, dT0)) + list(extra_starts):
        r = minimize(f, np.array(s, float), method='Nelder-Mead', options=dict(maxiter=3000, maxfev=3000, xatol=1e-6, fatol=1e-4, adaptive=True))
        r = minimize(f, r.x, method='Nelder-Mead', options=dict(maxiter=2000, maxfev=2000, xatol=1e-7, fatol=1e-5, adaptive=True))
        if best is None or r.fun < best.fun: best = r
    dT0_, k, bb, rho, u1, u2, D = unpack(best.x, scen)
    b = bb*(1 + k); aR = float(aR_of_rho(rho))
    inc = np.arccos(min(b/aR, 1.0))
    x14 = np.sqrt(max((1 + k)**2 - b**2, 0))/(aR*np.sin(inc)); x23 = np.sqrt(max((1 - k)**2 - b**2, 0))/(aR*np.sin(inc))
    T14 = P/np.pi*np.arcsin(min(x14, 1)); T23 = P/np.pi*np.arcsin(min(x23, 1)) if b < 1 - k else 0.
    return dict(chi2=float(best.fun), x=[float(v) for v in best.x], dT0_min=dT0_*1440, k=k, b=b, rho_c=rho, u1=u1, u2=u2, D=D,
                T14_h=T14*24, flat=T23/T14 if T14 > 0 else 0., depth_ppt=float(D*k**2*1e3))
