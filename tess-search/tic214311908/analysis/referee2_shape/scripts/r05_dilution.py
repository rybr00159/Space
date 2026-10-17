"""r05: independent check of the shape-only maximum-dilution profile (transit geometry, k < 1).
A source supplying a fraction f_src of the light is transited by an opaque body (k <= 1), b in [0, 1+k), a/R free,
LD of the source free (Kipping q1, q2 in [0, 1], no prior). For each f_src the minimum chi2 over many starts
(Nelder-Mead in transformed parameters) is compared with the planet_free ML chi2 of r03 (which includes its
LD-prior penalty, ~0.1-0.5). Also the same for a uniform-disk source occulted by a larger body (k in [1, 60]).
usage: r05_dilution.py DATASET -> out/dilution_{ds}.json"""
import sys, json, numpy as np
sys.path.insert(0, 'tic214311908/analysis/referee2_shape/scripts')
from rlib import *
from scipy.optimize import minimize
ds = sys.argv[1]
D = load(ds, 'fake'); lc = LC(D)
S = json.load(open(f'{RD}/out/scen_{ds}_fake.json'))
ref = S['planet_free']['chi2']; T0c, Pc = S['planet_free']['x'][:2]
def sig(x): return 1/(1 + np.exp(-x))
def unpack(y, kr):
    T0, P, zk, zb, laR, zq1, zq2 = y
    k = kr[0] + (kr[1] - kr[0])*sig(zk); b = sig(zb)*(1 + k); q1, q2 = sig(zq1), sig(zq2)
    u1 = 2*np.sqrt(q1)*q2; u2 = np.sqrt(q1)*(1 - 2*q2)
    return T0, P, k, b, 10**laR, u1, u2
def c2(y, fs, kr):
    T0, P, k, b, aR, u1, u2 = unpack(y, kr)
    if not (1.5 < aR < 500) or abs(T0 - T0c) > 0.01: return 1e9
    m = lc.transit(T0, P, k, aR, b, u1, u2)
    return lc.chi2(1 + fs*(m - 1))
def inv_sig(p): p = np.clip(p, 1e-6, 1 - 1e-6); return np.log(p/(1 - p))
def best(fs, kr):
    dneed = min(0.95, 0.016/fs)
    ks = [np.sqrt(dneed)*f for f in (0.9, 1.2, 1.7)] if kr[1] <= 1 else [1.5, 4., 15.]
    rbest = None
    for k0 in ks:
        k0 = float(np.clip(k0, kr[0]*1.001 + 1e-4, kr[1]*0.999))
        for bf in (0.1, 0.5, 0.8, 0.93):
            b0 = bf*(1 + k0); chord = np.sqrt(max((1 + k0)**2 - b0**2, 1e-4)); aR0 = np.clip(Pc/(np.pi*0.07)*chord, 2, 400)
            y0 = [T0c, Pc, inv_sig((k0 - kr[0])/(kr[1] - kr[0])), inv_sig(bf), np.log10(aR0), inv_sig(0.4), inv_sig(0.35)]
            r = minimize(c2, y0, args=(fs, kr), method='Nelder-Mead', options=dict(maxiter=5000, maxfev=5000, xatol=1e-7, fatol=1e-4))
            r = minimize(c2, r.x, args=(fs, kr), method='Nelder-Mead', options=dict(maxiter=4000, maxfev=4000, xatol=1e-8, fatol=1e-5))
            if rbest is None or r.fun < rbest.fun: rbest = r
    T0, P, k, b, aR, u1, u2 = unpack(rbest.x, kr)
    T14, T23 = t14(P, k, aR, b)
    return dict(chi2=float(rbest.fun), dchi2=float(rbest.fun - ref), k=float(k), b=float(b), aR=float(aR), u1=float(u1), u2=float(u2), T14_h=float(T14*24), flat=float(T23/T14) if T14 > 0 else 0.)
out = dict(dataset=ds, planet_free_chi2=ref, rows=[])
for fs in [1.0, 0.2, 0.09, 0.06, 0.045, 0.035, 0.031, 0.0225]:
    row = dict(f_src=fs, dT=float(-2.5*np.log10(fs)), transit=best(fs, (0.01, 1.0)))
    if fs <= 0.06: row['occultation'] = best(fs, (1.0, 60.0))
    out['rows'].append(row); print(json.dumps(row), flush=True)
json.dump(out, open(f'{RD}/out/dilution_{ds}.json', 'w'), indent=1)
