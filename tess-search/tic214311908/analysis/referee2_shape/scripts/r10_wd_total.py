"""r10: white dwarf occulted by the M3V - joint chi2 profile in M_WD of light curve + phase-0.5 secondary
(+ M-dwarf R, M priors), re-optimised with the secondary term included (WD radius free within the cool branch
R_Nauenberg x [1, 1.3]; the hot branch x [1, 2.5] below 0.45 Msun is also given), then combined with the r06 RUWE
likelihood (-2 ln L, density version, NSL scan pattern, edge-on orbit; sigma_eff = 0.18 (catalogue parallax error),
0.22 (mean of the five catalogue errors) and 0.25 mas (largest implied)).
Secondary measurement: r04 INV, WD-transit template: -0.20 +- 0.90 ppt (empirical); alternative REF r18 value
+0.29 +- 0.77 ppt (../referee_phot) for the REF light curve.
usage: r10_wd_total.py DATASET -> out/wd_total_{ds}.json"""
import sys, json, numpy as np
sys.path.insert(0, 'tic214311908/analysis/referee2_shape/scripts')
from rlib import *
from rwd import *
from scipy.optimize import minimize
ds = sys.argv[1]
D = load(ds, 'fake'); lc = LC(D); sup = Super(lc)
S = json.load(open(f'{RD}/out/scen_{ds}_fake.json'))
sec = json.load(open(f'{RD}/out/secondary.json'))['phase05']['wdtransit_white']
SEC = dict(own=(sec['depth_ppt']*1e-3, max(sec['empirical_std_ppt'], sec['empirical_robust_ppt'])*1e-3), ref_r18=(0.29e-3, 0.77e-3))
ref_chi2 = S['planet_free']['chi2']
def c2(y, MA, branch, secm):
    T0, P, xf, bocc, fA, RB, MB = y
    if not (0 <= xf <= 1 and 0 <= fA < 0.2 and 0.2 < RB < 0.6 and 0.2 < MB < 0.6 and 0 <= bocc < 1.2): return 1e9
    fmax = 1.3 if (branch == 'cool' or MA >= 0.45) else 2.5
    RA = r_nauenberg(MA)*(1 + xf*(fmax - 1))
    a = kepler_a_rsun(MA + MB, P); cosi = bocc*RB/a
    m = wd_occultation(sup, T0, P, a, cosi, RA, RB, fA)
    v = lc.chi2(m) + ((RB - 0.381)/0.012)**2 + ((MB - 0.370)/0.011)**2
    d2, _ = secondary_depth(MA, RA, MB, RB, P, bocc, fA)
    return v + ((d2 - SEC[secm][0])/SEC[secm][1])**2
rows = []
for MA in [0.17, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2]:
    pr = [r for r in S['wd_profile'] if abs(r['M_WD'] - MA) < 1e-6][0]
    T0, P = pr['x'][0], pr['x'][1]
    for branch in (('cool', 'hot') if MA < 0.45 else ('cool',)):
        for secm in ('own', 'ref_r18'):
            best = None
            for xf in (0.0, 0.5, 1.0):
                for bo in (0.0, 0.3, 0.6, 0.85):
                    y0 = [T0, P, xf, bo, pr['fA'], 0.381, 0.370]
                    r = minimize(c2, y0, args=(MA, branch, secm), method='Powell', bounds=[(T0 - 0.005, T0 + 0.005), (P - 2e-5, P + 2e-5), (0, 1), (0, 1.15), (0.005, 0.05), (0.33, 0.43), (0.33, 0.41)],
                                 options=dict(xtol=1e-6, ftol=1e-7, maxfev=6000))
                    if best is None or r.fun < best.fun: best = r
            T0b, Pb, xf, bocc, fA, RB, MB = best.x
            fmax = 1.3 if (branch == 'cool' or MA >= 0.45) else 2.5
            RA = r_nauenberg(MA)*(1 + xf*(fmax - 1))
            d2, RE = secondary_depth(MA, RA, MB, RB, Pb, bocc, fA)
            a = kepler_a_rsun(MA + MB, Pb); cosi = bocc*RB/a
            lcc = lc.chi2(wd_occultation(sup, T0b, Pb, a, cosi, RA, RB, fA))
            row = dict(M_WD=MA, branch=branch, sec_meas=secm, chi2=float(best.fun), dchi2_total_vs_planet=float(best.fun - ref_chi2), dchi2_lc_vs_planet=float(lcc - ref_chi2),
                       R_WD=float(RA), b_occ=float(bocc), fA=float(fA), secondary_ppt=float(d2*1e3), sec_term=float(((d2 - SEC[secm][0])/SEC[secm][1])**2))
            rows.append(row); print(json.dumps(row), flush=True)
out = dict(dataset=ds, planet_free_chi2=ref_chi2, secondary_used=SEC, rows=rows, ruwe={})
for sm in ('plx', '0.22', '0.25'):
    try:
        R = json.load(open(f'{RD}/out/ruwe_nsl_{sm}.json'))
    except FileNotFoundError:
        continue
    ms = sorted(float(k) for k in R['masses'])
    tab = []
    for row in rows:
        if row['branch'] != 'cool': continue
        mm = min(ms, key=lambda q: abs(q - row['M_WD']))
        ru = R['masses'][str(mm) if str(mm) in R['masses'] else repr(mm)]['nsl']['m2lnL_density']
        tab.append(dict(M_WD=row['M_WD'], sec_meas=row['sec_meas'], lc_sec=row['dchi2_total_vs_planet'], ruwe=ru, total=row['dchi2_total_vs_planet'] + ru))
    out['ruwe'][sm] = dict(table=tab, min_total_own=min(t['total'] for t in tab if t['sec_meas'] == 'own'), min_total_ref_r18=min(t['total'] for t in tab if t['sec_meas'] == 'ref_r18'))
    print(sm, out['ruwe'][sm]['min_total_own'], out['ruwe'][sm]['min_total_ref_r18'])
json.dump(out, open(f'{RD}/out/wd_total_{ds}.json', 'w'), indent=1)
