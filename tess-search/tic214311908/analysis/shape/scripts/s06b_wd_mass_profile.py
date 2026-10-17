"""s06b: chi2 profile of the 'WD totally occulted by the single M3V' scenario versus the WD mass (circular orbit,
WD radius tied to its mass as in s06, light-curve + phase-0.5 secondary + R/M priors of the M3V), combined with the
astrometric likelihood L(M2) of s16 (RUWE 1.076; sigma_CCD 0.22 and 0.35 mas).
Multi-start bounded Powell at each fixed mass. Two radius branches: 'hot' (R = R_Nauenberg x 1-2.5 below 0.45 Msun, x 1-1.3 above; the s06 default, appropriate for
young/hot He WDs) and 'cool' (x 1-1.3 at all masses, appropriate below ~8000 K). For each best fit the WD Teff
needed to supply f_A of the TESS flux is computed (Kurucz logg 5 surface fluxes, specs.py), so that the radius
branch can be checked for consistency.
usage: s06b_wd_mass_profile.py DATASET [hot|cool] -> out/wd_mass_profile_{DATASET}_{branch}.json"""
import sys, json, numpy as np
branch = sys.argv[2] if len(sys.argv) > 2 else 'hot'
sys.argv = [sys.argv[0], sys.argv[1], 'single']
ds = sys.argv[1]
src = open('tic214311908/analysis/shape/scripts/s06_compact.py').read()
src = src[:src.index('bounds = [')]
exec(src)
from scipy.optimize import minimize
from specs import S as SBF
if branch == 'cool':
    def unpack(x):
        T0, MA, lRA, bB, fA, RB, MB = x[:7]
        lRA = np.log10(r_nauenberg(MA)*(1 + lRA*0.3))
        return T0, MA, 10**lRA, bB, fA, RB, MB, 0., 90.
_Tg = np.exp(np.linspace(np.log(2500), np.log(50000), 500)); _sT = np.array([SBF(t, 'T', True) for t in _Tg])
def teff_needed(fA, RA, RB):
    r = fA/(1 - fA)*SBF(3350, 'T')*RB**2/RA**2
    return float(np.interp(r, _sT, _Tg)) if _sT.min() < r < _sT.max() else float('nan')
astro = json.load(open(f'{SD}/out/wd_checks.json'))['astrometry']
rows = []
for MA in [0.17, 0.2, 0.25, 0.3, 0.35, 0.45, 0.6, 0.8, 1.0, 1.2]:
    def f(y):
        x = [y[0], MA, y[1], y[2], y[3], y[4], y[5]]
        return c2(x)
    bnds = [(T0c-0.008, T0c+0.008), (0, 1), (0, 1.05), (0.005, 0.06), (0.37, 0.40), (0.34, 0.40)]
    best = None
    for xr in (0.0, 0.5, 1.0):
        for bb in (0.05, 0.4, 0.7, 0.9):
            r = minimize(f, [T0c, xr, bb, 0.0155, 0.381, 0.370], method='Powell', bounds=bnds, options=dict(xtol=1e-6, ftol=1e-7, maxfev=8000))
            if best is None or r.fun < best.fun: best = r
    x = [best.x[0], MA, *best.x[1:]]
    T0, MA_, RA, bB, fA, RB, MB, e, w = unpack(x)
    a, k, aRA = geom(MA_, RA, RB, MB)
    T14, T23 = t14_t23(P, k, aRA, bB*k)
    sec, _ = secondary(x)
    La = {s: min([v['L'] for kk, v in astro.items() if kk.startswith(f'M{m}_sig{s}')] or [np.nan]) for s in ('0.22', '0.35') for m in [min([0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0, 1.2], key=lambda q: abs(q - MA))]}
    row = dict(M_WD=MA, chi2=float(best.fun), chi2_lc=float(lc.chi2(model(x))), dchi2_vs_planet_free=float(best.fun - fit['chi2_best']),
               dchi2_lc_vs_planet_free=float(lc.chi2(model(x)) - fit['chi2_best']), R_WD=float(RA), b_occ=float(bB), f_A=float(fA), T14_h=float(T14*24),
               T23_h=float(T23*24), secondary_ppt=float(sec*1e3), Teff_WD_needed=teff_needed(fA, RA, RB), L_astrometry_nearest_mass=La)
    rows.append(row); print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
json.dump(dict(dataset=ds, branch=branch, planet_free_chi2=fit['chi2_best'], rows=rows), open(f'{SD}/out/wd_mass_profile_{ds}_{branch}.json', 'w'), indent=1)
