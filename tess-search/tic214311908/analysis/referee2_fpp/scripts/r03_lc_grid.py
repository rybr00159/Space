"""r03: light-curve likelihood grid for 'planet-like' transits (dark body, radius ratio k, quadratic LD) on a host that
supplies a fraction D of the deblended M3V-system flux, with circular-equivalent density rho_c and impact parameter b.
Own implementation (differences from fpp/f02):
  * data: REF = referee_phot/lc_M3V_referee_windows.csv (primary) and INV = phot/lc_M3V_all.csv (flux_norm rescaled by
    935/979, quality_ok, no flares); |dt| < 0.30 d from the phase-1 ephemeris (T0 3840.9141, P 5.6915087)
  * per-window QUADRATIC baseline (INV) / LINEAR baseline (REF) profiled analytically (the INV curve still contains the
    0.6485-d RRab leak, which a quadratic absorbs within +-0.3 d)
  * noise: per-sector multiplier = referee fake-window err_fake/err_white (REF, referee_phot/out/per_sector.csv) or
    scale x beta from 62 fake ephemerides (INV, referee2_shape/out/beta_inv.json); variant 'white' without beta
  * exposure integration by the midpoint rule (batman's supersample_factor uses end-point offsets -> 1.2x too long)
  * k integrated numerically on a ln k grid (13 points over +-4 local widths) instead of a Laplace approximation
Stores for each (D, rho_c, b): chi2_min, k_best, I = int exp(-(chi2 - chi2_ref)/2) dlnk.
Usage: r03_lc_grid.py <ref|inv> <beta|white>   (env R03_T0 overrides T0, R03_TAG suffixes the output name)"""
import numpy as np, pandas as pd, sys, time, batman, json
from scipy.optimize import minimize_scalar
sys.path.insert(0, 'tic214311908/analysis/referee2_fpp/scripts')
from rlib import *
ds, nm = sys.argv[1], sys.argv[2]
import os
T0 = float(os.environ.get('R03_T0', T0))   # optional ephemeris override (T0 sensitivity test)
TAG = os.environ.get('R03_TAG', '')
LD = (0.125, 0.533)
if ds == 'ref':
    d = pd.read_csv(f'{B}/referee_phot/lc_M3V_referee_windows.csv')
    d = d[np.isfinite(d.flux_norm) & np.isfinite(d.flux_err)]
    t, f, e, sec, cad = d.time_btjd_tdb.values, d.flux_norm.values, d.flux_err.values, d.sector.values, d.cadence_s.values
    ps = pd.read_csv(f'{B}/referee_phot/out/per_sector.csv').set_index('sector')
    mult = {s: (ps.loc[s, 'err_adopt']/ps.loc[s, 'err_white'] if nm == 'beta' else 1.0) for s in ps.index}
    npoly = 2
else:
    d = pd.read_csv(f'{B}/phot/lc_M3V_all.csv', usecols=['time_btjd_tdb', 'sector', 'cadence_s', 'quality_ok', 'flux_norm', 'flux_err', 'flare_flag'])
    d = d[(d.quality_ok == 1) & (d.flare_flag == 0) & np.isfinite(d.flux_norm) & np.isfinite(d.flux_err)]
    sc = 935./979.
    t, f, e, sec, cad = d.time_btjd_tdb.values, 1 + (d.flux_norm.values - 1)*sc, d.flux_err.values*sc, d.sector.values, d.cadence_s.values
    bj = jload('referee2_shape/out/beta_inv.json')
    mult = {int(s): (v['scale']*v['beta'] if nm == 'beta' else v['scale']) for s, v in bj.items() if isinstance(v, dict)}
    npoly = 3
ep = np.round((t - T0)/P).astype(int); dt = t - (T0 + ep*P)
m = np.abs(dt) < 0.30
t, f, e, sec, cad, ep, dt = t[m], f[m], e[m], sec[m], cad[m].astype(int), ep[m], dt[m]
sig = e*np.array([mult[int(s)] for s in sec])
w = 1/sig**2
win = sec*100000 + ep
_, wi = np.unique(win, return_inverse=True); nw = wi.max() + 1
# per-window polynomial baseline: normal matrices
X = np.vstack([dt**p for p in range(npoly)]).T
A = np.zeros((nw, npoly, npoly))
for a_ in range(npoly):
    for b_ in range(npoly):
        A[:, a_, b_] = np.bincount(wi, w*X[:, a_]*X[:, b_], nw)
Ainv = np.linalg.inv(A)
def chi2(mod):
    r = f - mod
    s = np.vstack([np.bincount(wi, w*r*X[:, p], nw) for p in range(npoly)]).T
    return float(np.sum(w*r*r) - np.einsum('wi,wij,wj->', s, Ainv, s))
# fine-dt model per cadence with midpoint exposure integration
TG = np.linspace(-0.31, 0.31, 1241)
pp = batman.TransitParams(); pp.t0 = 0.; pp.per = P; pp.rp = 0.1; pp.a = 25.; pp.inc = 89.; pp.ecc = 0.; pp.w = 90.
pp.limb_dark = 'quadratic'; pp.u = list(LD)
MODS = {}
for c in np.unique(cad):
    nsub = {1800: 25, 600: 9, 200: 5, 120: 3}.get(int(c), 5)
    off = (np.arange(nsub) + 0.5)/nsub - 0.5
    tt = (TG[:, None] + off[None, :]*c/DAY).ravel()
    MODS[int(c)] = (batman.TransitModel(pp, tt), nsub, np.where(cad == c)[0])
def model(k, aR, b, D):
    cosi = b/aR
    out = np.ones(len(t))
    if cosi >= 1: return out
    pp.rp = k; pp.a = aR; pp.inc = np.degrees(np.arccos(cosi))
    for c, (tm, nsub, idx) in MODS.items():
        lc = tm.light_curve(pp).reshape(len(TG), nsub).mean(axis=1)
        out[idx] = np.interp(dt[idx], TG, lc)
    return 1 + D*(out - 1)
def aR_of_rho(rho): return (rho*G*(P*DAY)**2/(3*np.pi))**(1/3.)
chi2_null = chi2(np.ones(len(t)))
Ds = np.array([1.0, 0.8, 0.65, 0.5, 0.38, 0.28, 0.2, 0.14, 0.1, 0.07, 0.05, 0.035])
rhos = np.logspace(0, np.log10(200.), 30)
bs = np.round(np.arange(0, 1.501, 0.05), 3)
CM = np.full((len(Ds), len(rhos), len(bs)), np.nan); KB = np.full_like(CM, np.nan); LI = np.full_like(CM, np.nan)
t1 = time.time()
for i, D in enumerate(Ds):
    for j, rho in enumerate(rhos):
        aR = aR_of_rho(rho)
        for l, b in enumerate(bs):
            lkmin = np.log(max(0.02, b - 0.98)) if b > 1 else np.log(0.02)
            fun = lambda lk: chi2(model(np.exp(lk), aR, b, D))
            r = minimize_scalar(fun, bounds=(lkmin, np.log(1.5)), method='bounded', options=dict(xatol=1e-3))
            lk0, c0 = float(r.x), float(r.fun)
            # local width from a finite difference, then numerical integral on 13 points over +-4 widths
            h = 0.02
            cp = fun(lk0 + h); cmn = fun(max(lk0 - h, lkmin))
            d2 = (cp - 2*c0 + cmn)/h**2
            wid = 1/np.sqrt(d2/2) if d2 > 0 else 0.3
            wid = min(max(wid, 0.005), 0.5)
            lg = lk0 + np.linspace(-4, 4, 13)*wid
            lg = lg[(lg >= lkmin) & (lg <= np.log(1.5))]
            cg = np.array([fun(x) for x in lg])
            CM[i, j, l] = min(c0, cg.min()); KB[i, j, l] = float(np.exp(lg[np.argmin(cg)] if cg.min() < c0 else lk0))
            LI[i, j, l] = float(np.trapezoid(np.exp(-0.5*(cg - c0)), lg)) * np.exp(-0.5*(c0 - CM[i, j, l]))
    print(ds, nm, 'D', D, round(time.time() - t1), 's; best dchi2 vs null', round(np.nanmin(CM[i]) - chi2_null, 2), flush=True)
cref = np.nanmin(CM)
# LI is relative to each cell's own minimum; convert to the global reference: I = LI * exp(-(CM - cref)/2)
I = LI*np.exp(-0.5*(CM - cref))
np.savez(f'{R2}/out/r03_grid_{ds}_{nm}{TAG}.npz', D=Ds, rho=rhos, b=bs, chi2=CM, k=KB, I=I, chi2_null=chi2_null, chi2_ref=cref, n=len(t))
print('saved; global best dchi2 vs null', cref - chi2_null, 'n', len(t))
