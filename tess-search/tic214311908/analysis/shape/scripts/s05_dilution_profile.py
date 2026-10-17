"""s05: maximum-dilution / blend test from the transit shape (generalised Vanderburg+19 / Seager & Mallen-Ornelas 03).
Model: a luminous source contributing a fraction f_src of the aperture flux is eclipsed by an opaque body of radius
ratio k (any k: k < 1 = transit/annular geometry; k > 1 = occultation of a smaller star by a larger body, which can be
total and flat-bottomed with 100% depth). Free: T0, k, b (0..1+k), a/R_src, Kipping q1,q2 (LD of the eclipsed source,
uniform, i.e. any LD). P fixed at the planet value. The occulter's own light only adds to the dilution, so this
geometric family contains every spherical-star blend at the primary eclipse (no ellipsoidal/reflection: those vary
on P/2 and are removed by the 0.75-d biweight anyway).
For a grid of f_src the minimum chi2 is found by multi-start bounded Powell (DE check at f_src = 1, 0.2, 0.031) (two regimes: k in [0.01,1], k in [1,60]).
Also: profiles of chi2 vs k for f_src = 1 (M3V host) and 0.5 (twin host), with and without the stellar-density
prior, which constrain grazing brown-dwarf / stellar companions of the M3V (s08 uses them).
usage: s05_dilution_profile.py DATASET   ->  out/dilution_profile_{DATASET}.json"""
import sys, json, numpy as np
sys.path.insert(0, 'tic214311908/analysis/shape/scripts')
from shapelib import *
from scipy.optimize import differential_evolution, minimize
ds = sys.argv[1]
lc = LC({'inv': f'{SD}/data/inv_win.npz', 'ref': f'{SD}/data/ref_win.npz'}[ds])
fit = json.load(open(f'{SD}/out/fit_{ds}_free.json'))
P = fit['post']['P']['med']; T0c = fit['post']['T0']['med']
def model(x, fsrc):
    T0, lk, bf, laR, q1, q2 = x
    k = 10**lk; b = bf*(1+k); aR = 10**laR
    u1 = 2*np.sqrt(q1)*q2; u2 = np.sqrt(q1)*(1-2*q2)
    m = lc.transit(T0, P, k, aR, b, u1, u2)
    return 1 + fsrc*(m-1)
def c2(x, fsrc, rho_prior=None):
    T0, lk, bf, laR, q1, q2 = x
    val = lc.chi2(model(x, fsrc))
    if rho_prior is not None:
        rho = rho_from_aR(10**laR, P)
        val += ((rho - rho_prior[0])/rho_prior[1])**2
    return val
def best(fsrc, kr, rho_prior=None, seed=1, laR_rng=(np.log10(3), np.log10(400)), use_de=False):
    """min chi2: multi-start bounded Powell (starts on a grid of b and k; k from the depth the source needs);
    optional differential evolution check (use_de)."""
    bounds = [(T0c-0.006, T0c+0.006), (np.log10(kr[0]), np.log10(kr[1])), (0, 1), laR_rng, (0, 1), (0, 1)]
    dneed = min(0.95, 0.0158/fsrc)
    starts = []
    if kr[1] <= 1.0001:
        kgrid = [np.clip(np.sqrt(dneed)*f, kr[0], kr[1]) for f in (0.9, 1.1, 1.5, 2.2)]
    else:
        kgrid = [np.clip(x, kr[0], kr[1]) for x in (1.5, 3, 8, 25)]
    for kk in sorted(set(np.round(kgrid, 5))):
        for bf in (0.05, 0.35, 0.6, 0.75, 0.85, 0.93, 0.98):
            b = bf*(1+kk); chord = np.sqrt(max((1+kk)**2 - b**2, 1e-4))
            aR = np.clip(P/(np.pi*0.069)*chord, 3.1, 390)
            starts.append([T0c, np.log10(kk), bf, np.log10(aR), 0.4, 0.3])
    bestr = None
    for x0 in starts:
        r = minimize(c2, x0, args=(fsrc, rho_prior), method='Powell', bounds=bounds, options=dict(xtol=1e-5, ftol=1e-6, maxfev=6000))
        if bestr is None or r.fun < bestr.fun: bestr = r
    if use_de:
        r = differential_evolution(c2, bounds, args=(fsrc, rho_prior), seed=seed, popsize=12, maxiter=200, tol=1e-7, polish=True, init='sobol')
        if r.fun < bestr.fun: bestr = r
    r = bestr
    x = r.x; k = 10**x[1]; b = x[2]*(1+k); aR = 10**x[3]
    T14, T23 = t14_t23(P, k, aR, b)
    m = model(x, fsrc); d_obs = 1 - m.min()
    und = 1 - lc.transit(x[0], P, k, aR, b, 2*np.sqrt(x[4])*x[5], np.sqrt(x[4])*(1-2*x[5])).min()
    return dict(chi2=float(r.fun), T0=float(x[0]), k=float(k), b=float(b), aR=float(aR), q1=float(x[4]), q2=float(x[5]),
                T14_h=float(T14*24), T23_h=float(T23*24), undiluted_depth_model=float(und), obs_depth_model=float(d_obs),
                rho_src_circ=float(rho_from_aR(aR, P)), rho_occulter_eff=float(rho_from_aR(aR/k, P)))
out = dict(dataset=ds, P=P, chi2_flat=float(lc.chi2(np.ones(lc.n))), planet_free_chi2=fit['chi2_best'])
fgrid = [1.0, 0.5, 0.3, 0.2, 0.13, 0.09, 0.06, 0.045, 0.035, 0.031, 0.0225, 0.018]
prof = []
for fs in fgrid:
    row = dict(f_src=fs, dmag=float(-2.5*np.log10(fs)))
    row['transit_k_lt_1'] = best(fs, (0.01, 1.0), use_de=(fs in (1.0, 0.2, 0.031)))
    row['occultation_k_gt_1'] = best(fs, (1.0, 60.0))
    prof.append(row)
    print(fs, round(row['transit_k_lt_1']['chi2'], 2), round(row['occultation_k_gt_1']['chi2'], 2), row['transit_k_lt_1']['k'], row['occultation_k_gt_1']['k'], flush=True)
out['profile'] = prof
# chi2 vs k for the M3V (f=1) and twin (f=0.5) hosts, rho free or with prior
kp = []
for kk in [0.10, 0.14, 0.17, 0.20, 0.23, 0.27, 0.32, 0.40, 0.55, 0.8]:
    row = dict(k=kk)
    for lab, fs, pri in [('M3V_free', 1.0, None), ('M3V_rho', 1.0, (9.41, 0.88)), ('twin_rho', 0.5, (15.5, 1.45))]:
        row[lab] = best(fs, (kk, kk*1.0001), pri)
    kp.append(row)
    print('k', kk, {l: round(row[l]['chi2'], 2) for l in ['M3V_free', 'M3V_rho', 'twin_rho']}, flush=True)
out['k_profile'] = kp
json.dump(out, open(f'{SD}/out/dilution_profile_{ds}.json', 'w'), indent=1)
