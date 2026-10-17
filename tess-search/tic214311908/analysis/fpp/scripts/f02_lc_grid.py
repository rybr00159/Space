"""f02: light-curve likelihood grid for 'planet-like' transits (dark body, radius ratio k) on a host that supplies a
fraction D of the (deblended) M3V flux, with circular-equivalent host density rho_circ and impact parameter b.
For each (D, rho_circ, b) the radius ratio k is optimised (bounded Brent in ln k) and the curvature d2chi2/dk2 is
stored for a Laplace approximation of the k-integral. Per-window constant offsets are profiled analytically.
Model: batman quadratic LD (0.125, 0.533; Claret 2017 TESS, 3350 K) with exposure integration (1800 s: 31, 600 s: 11,
200 s: 5 sub-samples), f = 1 + D (transit - 1); T0, P fixed (fpplib). The same grid serves every planet hypothesis:
TP (D = 1, rho = 9.4), TP-twin (D = 0.5, rho = 15.5), TP diluted by a companion (D = 1 - D_c), PTP on a companion
(D = D_c, rho = rho_c), and planets on background stars.
Data: INV = phot model B with the RRab leak removed (search/), REF = referee local-window photometry (shape/s01 sigma).
Usage: f02_lc_grid.py <inv|ref>   -> out/lc_grid_<ds>.npz"""
import numpy as np, sys, time, batman
from scipy.optimize import minimize_scalar
sys.path.insert(0, 'tic214311908/analysis/fpp/scripts')
from fpplib import *
ds = sys.argv[1]
d = inv_clean_windows(0.30) if ds == 'inv' else ref_windows()
t, f, sig, cad = d['t'], d['f'], d['sig'], d['cad']
win = d['sector']*100000 + d['epoch']
_, wi = np.unique(win, return_inverse=True); nw = wi.max() + 1
w = 1/sig**2; sw = np.bincount(wi, w, nw)
SS = {1800: 31, 600: 11, 200: 5}
groups = []
p0 = batman.TransitParams(); p0.t0 = T0; p0.per = P; p0.rp = 0.1; p0.a = 25.; p0.inc = 89.; p0.ecc = 0.; p0.w = 90.
p0.limb_dark = 'quadratic'; p0.u = list(LD_M3V)
for c in np.unique(cad):
    idx = np.where(cad == c)[0]
    groups.append((idx, batman.TransitModel(p0, t[idx], supersample_factor=SS[int(c)], exp_time=c/DAY)))
def aR_rho(rho): return (rho*G*(P*DAY)**2/(3*np.pi))**(1/3.)
def model(k, aR, b, D):
    cosi = b/aR
    if cosi >= 1: return np.ones(len(t))
    p0.rp = k; p0.a = aR; p0.inc = np.degrees(np.arccos(cosi))
    m = np.empty(len(t))
    for idx, tm in groups: m[idx] = tm.light_curve(p0)
    return 1 + D*(m - 1)
def chi2(mod):
    r = f - mod
    o = np.bincount(wi, r*w, nw)/sw
    r = r - o[wi]
    return float(np.sum(r*r*w))
chi2_null = chi2(np.ones(len(t)))
Ds = np.array([1.0, 0.97, 0.93, 0.86, 0.75, 0.62, 0.5, 0.39, 0.29, 0.2, 0.14, 0.09, 0.06, 0.04, 0.025])
rhos = np.logspace(np.log10(1.0), np.log10(300.), 26)
bs = np.round(np.arange(0, 1.901, 0.05), 3)
C = np.full((len(Ds), len(rhos), len(bs)), np.nan); K = np.full_like(C, np.nan); SK = np.full_like(C, np.nan)
t1 = time.time()
for i, D in enumerate(Ds):
    for j, rho in enumerate(rhos):
        aR = aR_rho(rho)
        for l, b in enumerate(bs):
            kmin = max(0.01, b - 0.999) if b > 1 else 0.01   # the body must reach the disc
            fun = lambda lk: chi2(model(np.exp(lk), aR, b, D))
            r = minimize_scalar(fun, bounds=(np.log(kmin + 1e-4), np.log(1.6)), method='bounded', options=dict(xatol=2e-3))
            k = float(np.exp(r.x)); c0 = float(r.fun)
            h = 0.01*k
            cp = chi2(model(k + h, aR, b, D)); cm = chi2(model(max(k - h, 1e-4), aR, b, D))
            d2 = (cp - 2*c0 + cm)/h**2
            C[i, j, l] = c0; K[i, j, l] = k; SK[i, j, l] = 1/np.sqrt(d2/2) if d2 > 0 else np.nan
    print(ds, 'D', D, 'done', round(time.time() - t1), 's; best chi2', np.nanmin(C[i]) - chi2_null, flush=True)
np.savez(f'{F}/out/lc_grid_{ds}.npz', D=Ds, rho=rhos, b=bs, chi2=C, k=K, sigk=SK, chi2_null=chi2_null, n=len(t),
         note='chi2 includes profiled per-window offsets; sigk = Laplace width of the k-likelihood (exp(-chi2/2))')
print('saved; global min dchi2 vs null', np.nanmin(C) - chi2_null)
