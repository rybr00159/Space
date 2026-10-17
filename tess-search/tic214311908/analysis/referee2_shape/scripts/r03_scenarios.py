"""r03: maximum-likelihood comparison on the REAL light curves (both photometries) of
  planet_free   : batman planet, rho free, LD prior (as r02)
  planet_single : + rho ~ N(9.41, 0.88)
  vshape        : planet with b >= 1 - k (no flat part), rho free
  wd_cool / wd_hot : hot compact source (WD) totally occulted by the M3V (own overlap model, rwd.py), circular,
                  Kepler a from M_WD + M_M; R_M ~ N(0.381, 0.012), M_M ~ N(0.370, 0.011); R_WD = R_Nauenberg(M) x f,
                  f in [1, 1.3] (cool) or [1, 2.5] below 0.45 Msun (hot); light curve only (secondary added in r05)
  wd_profile    : wd_cool at fixed M_WD = 0.17 ... 1.2
  wd_twin       : WD occulted by one twin component (R 0.289 +- 0.009, M 0.264 +- 0.008), cool branch
All with per-window offsets profiled; the LD prior / stellar priors are included in chi2 as penalties so that
dchi2 values compare like with like (planet_free includes its LD-prior penalty).
usage: r03_scenarios.py DATASET [beta]  -> out/scen_{ds}_{beta}.json"""
import sys, json, numpy as np
sys.path.insert(0, 'tic214311908/analysis/referee2_shape/scripts')
from rlib import *
from rwd import *
from scipy.optimize import minimize
ds = sys.argv[1]; beta = sys.argv[2] if len(sys.argv) > 2 else 'fake'
D = load(ds, beta); lc = LC(D); sup = Super(lc)
U1, U2, SU = 0.125, 0.533, 0.10
def pen_ld(u1, u2):
    if not (u1 + u2 < 1 and u1 > 0 and u1 + 2*u2 > 0): return 1e9
    return ((u1 - U1)/SU)**2 + ((u2 - U2)/SU)**2
def c2_planet(x, rho_prior=None, vshape=False, dil=1.0):
    T0, P, k, bx, lr, u1, u2 = x
    if not (0.01 < k < 0.9 and -1 < lr < 2.5): return 1e9
    if vshape:
        b = (1 - k) + np.clip(bx, 0, 1)*2*k
    else:
        b = abs(bx)
        if b >= 1 + k: return 1e9
    m = lc.transit(T0, P, k, aR_from_rho(10**lr, P), b, u1, u2)
    m = 1 + (m - 1)/dil
    v = lc.chi2(m) + pen_ld(u1, u2)
    if rho_prior: v += ((10**lr - rho_prior[0])/rho_prior[1])**2
    return v
def multistart(f, starts, method='Nelder-Mead', opts=None):
    best = None
    for s in starts:
        r = minimize(f, s, method=method, options=opts or dict(maxiter=8000, maxfev=8000, xatol=1e-7, fatol=1e-5))
        r = minimize(f, r.x, method=method, options=opts or dict(maxiter=8000, maxfev=8000, xatol=1e-8, fatol=1e-6))
        if best is None or r.fun < best.fun: best = r
    return best
out = dict(dataset=ds, beta=beta, info=D['info'], n=lc.n)
st = [[T0REF, PREF, 0.12, b0, lr0, U1, U2] for b0 in (0.1, 0.45, 0.75) for lr0 in (0.6, 1.0, 1.2)]
r = multistart(c2_planet, st); out['planet_free'] = dict(chi2=r.fun, x=r.x.tolist())
P_ = r.x[1]; T0_ = r.x[0]
r = multistart(lambda x: c2_planet(x, (9.41, 0.88)), st); out['planet_single'] = dict(chi2=r.fun, x=r.x.tolist())
st2 = [[T0REF, PREF, 0.12*1.41, b0, np.log10(15.5), U1, U2] for b0 in (0.1, 0.4, 0.7)]
r = multistart(lambda x: c2_planet(x, (15.5, 1.45), dil=2.0), st2); out['planet_twin'] = dict(chi2=r.fun, x=r.x.tolist())
stv = [[T0REF, PREF, k0, bx0, lr0, U1, U2] for k0 in (0.14, 0.25, 0.4) for bx0 in (0.05, 0.3) for lr0 in (0.3, 0.8)]
r = multistart(lambda x: c2_planet(x, vshape=True), stv); out['vshape'] = dict(chi2=r.fun, x=r.x.tolist(), b=float((1 - r.x[2]) + np.clip(r.x[3], 0, 1)*2*r.x[2]))
# ---- WD occulted by an M dwarf
def wd_unpack(y, MA=None, branch='cool', host='single'):
    T0, P, lMA, xf, bocc, fA, RB, MB = y
    MA = 10**lMA if MA is None else MA
    fmax = 1.3 if (branch == 'cool' or MA >= 0.45) else 2.5
    RA = r_nauenberg(MA)*(1 + np.clip(xf, 0, 1)*(fmax - 1))
    return T0, P, MA, RA, abs(bocc), fA, RB, MB
HOST = dict(single=((0.381, 0.012), (0.370, 0.011)), twin=((0.289, 0.009), (0.264, 0.008)))
def c2_wd(y, MA=None, branch='cool', host='single'):
    T0, P, MA, RA, bocc, fA, RB, MB = wd_unpack(y, MA, branch, host)
    if not (0.15 <= MA <= 1.38 and 0 <= fA < 0.2 and 0 <= y[3] <= 1 and 0.05 < RB < 0.6 and 0.05 < MB < 0.6): return 1e9
    a = kepler_a_rsun(MA + MB, P)
    cosi = bocc*RB/a   # bocc = impact parameter in units of the occulter radius
    if cosi >= 1: return 1e9
    m = wd_occultation(sup, T0, P, a, cosi, RA, RB, fA)
    (rp, sr), (mp, sm) = HOST[host]
    return lc.chi2(m) + ((RB - rp)/sr)**2 + ((MB - mp)/sm)**2
def wd_fit(MA=None, branch='cool', host='single'):
    rp, mp = HOST[host][0][0], HOST[host][1][0]
    starts = []
    for lm in ((np.log10(0.25), np.log10(0.6), np.log10(0.9)) if MA is None else (0.,)):
        for bo in (0.05, 0.4, 0.7):
            for xf in (0.0, 0.5):
                starts.append([T0_, P_, lm, xf, bo, 0.0157, rp, mp])
    f = (lambda y: c2_wd(y, None, branch, host)) if MA is None else (lambda y: c2_wd(y, MA, branch, host))
    r = multistart(f, starts)
    T0, P, MA_, RA, bocc, fA, RB, MB = wd_unpack(r.x, MA, branch, host)
    a = kepler_a_rsun(MA_ + MB, P)
    sec, RE = secondary_depth(MA_, RA, MB, RB, P, bocc, fA, u=(0.125, 0.533))
    lcc = r.fun - ((RB - HOST[host][0][0])/HOST[host][0][1])**2 - ((MB - HOST[host][1][0])/HOST[host][1][1])**2
    # durations: occulter radius RB passes over the small source
    cosi = bocc*RB/a; si = np.sqrt(1 - cosi**2)
    T14 = P/np.pi*np.arcsin(min(1, np.sqrt(max(0, (RB + RA)**2 - (bocc*RB)**2))/(a*si)))
    T23 = P/np.pi*np.arcsin(min(1, np.sqrt(max(0, (RB - RA)**2 - (bocc*RB)**2))/(a*si)))
    return dict(chi2=float(r.fun), chi2_lc=float(lcc), M_WD=float(MA_), R_WD=float(RA), b_occ=float(bocc), fA=float(fA), R_occ=float(RB), M_occ=float(MB),
                a_rsun=float(a), T14_h=float(T14*24), ingress_min=float((T14 - T23)/2*1440), secondary_ppt=float(sec*1e3), R_E=float(RE), x=r.x.tolist())
out['wd_cool'] = wd_fit(None, 'cool'); print('wd_cool', out['wd_cool'], flush=True)
out['wd_hot'] = wd_fit(None, 'hot'); print('wd_hot', out['wd_hot'], flush=True)
out['wd_twin'] = wd_fit(None, 'cool', 'twin'); print('wd_twin', out['wd_twin'], flush=True)
out['wd_profile'] = []
for MA in [0.17, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2]:
    rr = wd_fit(MA, 'cool'); out['wd_profile'].append(rr); print('profile', MA, rr, flush=True)
ref = out['planet_free']['chi2']
for kk in ['planet_single', 'planet_twin', 'vshape', 'wd_cool', 'wd_hot', 'wd_twin']:
    out[kk]['dchi2_vs_planet_free'] = out[kk]['chi2'] - ref
    print(kk, round(out[kk]['dchi2_vs_planet_free'], 2))
json.dump(out, open(f'{RD}/out/scen_{ds}_{beta}.json', 'w'), indent=1, default=float)
