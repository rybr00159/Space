"""s08: bound (hierarchical) eclipsing-binary scenarios inside the unresolved M3V system at 64.3 pc.
The aperture light of the 'M3V' (M_T,tot = 13.09 - 4.04 = 9.05) is shared by an eclipsing pair (eclipsed star 1,
mass m1, and occulter 2, mass m2; main-sequence K/M dwarfs from the Mamajek sequence, or a brown dwarf with
L_T = 0 and R = 0.075-0.11 Rsun) and, if light is left over, a third star m0 (MS) that dilutes the eclipses.
Luminosity budget in Ks, where the M-dwarf relations are tight and metallicity-insensitive (the target is 0.75
mag over-luminous in the optical CMD, so optical MS colour relations cannot be used for it): component
MKs(M) and R(MKs) from Mann+19 / Mann+15 (they give the adopted 0.370 Msun / 0.381 Rsun for the single star at
MKs 6.44); combinations whose EB is brighter than the system in Ks by > 0.1 mag are rejected; the third light m0
takes the remaining Ks flux. TESS-band fluxes use T - Ks(M) from the Mamajek sequence, Teff from Mamajek. This covers: the M3V itself eclipsed by a small star/BD
(m0 = 0), a twin component with a companion (m0 = m1 ~ 0.26), and an M-dwarf EB diluted by a brighter primary
(including the open 'dT 3.5-4.5 companion' and '4-5 au companion' cases).
Orbits: P = 5.6915 d circular (secondary at phase 0.5 constrained by the measured phase-0.5 depth: INV
-0.16 +- 0.79 ppt (s10 empirical), REF +0.29 +- 0.77 ppt, added as a chi2 term) or P = 11.383 d (the two eclipses alternate in
the 5.69-d windows and are both fitted). Free: impact parameter b (grid + refinement), T0 fixed at the planet value.
usage: s08_bound_eb.py DATASET -> out/bound_eb_{DATASET}.json, out/bound_eb_{DATASET}.csv"""
import sys, json, numpy as np, pandas as pd
sys.path.insert(0, 'tic214311908/analysis/shape/scripts')
from shapelib import *
from eem import *
from scipy.optimize import minimize_scalar
ds = sys.argv[1]
lc = LC({'inv': f'{SD}/data/inv_win.npz', 'ref': f'{SD}/data/ref_win.npz'}[ds])
fit = json.load(open(f'{SD}/out/fit_{ds}_free.json'))
P = fit['post']['P']['med']; T0 = fit['post']['T0']['med']
SEC = {'inv': (-0.16e-3, 0.79e-3), 'ref': (0.29e-3, 0.77e-3)}[ds]  # s10 empirical (INV), referee r18 (REF)
MT_TOT = 13.09 - 5*np.log10(64.31/10)
BPRP_OBS = 2.756
def comp_bprp(parts):
    """parts: list of (M_G, Bp-Rp, G-Rp) -> composite Bp-Rp"""
    fbp = 0; frp = 0
    for mg, bprp, grp in parts:
        mrp = mg - grp; mbp = mrp + bprp
        fbp += 10**(-0.4*mbp); frp += 10**(-0.4*mrp)
    return -2.5*np.log10(fbp/frp)
LT_TOT = 10**(-0.4*MT_TOT)
def ld(teff):
    return (0.13, 0.53) if teff < 3900 else (0.30, 0.30)
m_ms = np.round(np.concatenate([np.arange(0.08, 0.20, 0.01), np.arange(0.20, 0.42, 0.02)]), 3)
bds = [('BD', 0.03, 0.105), ('BD', 0.05, 0.09), ('BD', 0.07, 0.08), ('BD', 0.07, 0.075)]
MKS_TOT = 6.4406; FK_TOT = 10**(-0.4*MKS_TOT)
def comp(m):
    mk = mann_MKs(m); R = mann_radius(mk); T = ms_Teff(m)
    fk = 10**(-0.4*mk); fT = 10**(-0.4*(mk + t_minus_ks(m)))
    return R, T, fk, fT
_, _, FK_ref, FT_ref = comp(0.370)
LT_TOT = FT_ref*FK_TOT/FK_ref  # TESS flux of the whole system scaled from the single-star solution
rows = []
for m1 in m_ms:
    R1, T1, K1f, L1 = comp(m1)
    cands = [('MS', m2, comp(m2)[0]) for m2 in m_ms if m2 <= m1] + bds
    for typ, m2, R2 in cands:
        if typ == 'MS':
            _, _, K2f, L2 = comp(m2)
        else:
            K2f, L2 = 0.0, 0.0
        K0 = FK_TOT - K1f - K2f
        if K0 < -0.0965*FK_TOT:
            continue
        m0 = np.nan; L0 = 0.0
        if K0 > 0.01*FK_TOT:
            m0 = float(mann_mass(-2.5*np.log10(K0)))
            L0 = comp(min(max(m0, 0.08), 0.69))[3]
        bprp = np.nan
        fsrc = min(1.0, (L1 + L2)/(L0 + L1 + L2))
        k = R2/R1
        T2 = comp(m2)[1] if typ == 'MS' else 1500.
        for Peb in (P, 2*P):
            a = (G*(m1+m2)*MSUN*(Peb*DAY)**2/(4*np.pi**2))**(1/3.)/RSUN
            aR = a/R1
            def c2(b):
                m = eb_flux(lc, T0, Peb, k, aR, b, ld(T1), ld(T2) if typ == 'MS' else (0.3, 0.3), L2/L1, fsrc)
                v = lc.chi2(m)
                return v
            bs = np.linspace(0, 1 + k, 45)
            cs = np.array([c2(b) for b in bs])
            i = int(np.argmin(cs))
            lo, hi = bs[max(0, i-1)], bs[min(len(bs)-1, i+1)]
            r = minimize_scalar(c2, bounds=(lo, hi), method='bounded', options=dict(xatol=1e-4))
            bbest = r.x if r.fun < cs[i] else bs[i]; cbest = min(r.fun, cs[i])
            # secondary eclipse depth (diluted) at phase 0.5 for P_EB = P: compute with batman at mid-eclipse
            sec = np.nan
            if abs(Peb - P) < 1e-6:
                if L2 > 0:
                    import batman
                    pp = batman.TransitParams(); pp.t0 = 0; pp.per = Peb; pp.rp = 1/k; pp.a = aR/k
                    pp.inc = np.degrees(np.arccos(min(1, bbest/aR))); pp.ecc = 0; pp.w = 90; pp.limb_dark = 'quadratic'; pp.u = list(ld(T2))
                    fB = batman.TransitModel(pp, np.array([0.0])).light_curve(pp)[0]
                    sec = fsrc*(L2/(L1+L2))*(1 - fB)
                else:
                    sec = 0.0
                cbest_tot = cbest + ((sec - SEC[0])/SEC[1])**2
            else:
                cbest_tot = cbest
            T14, T23 = t14_t23(Peb, k, aR, bbest)
            rows.append(dict(m1=m1, R1=R1, Teff1=T1, typ2=typ, m2=m2, R2=R2, k=k, m0=m0, f_src=fsrc, L2_L1=L2/L1, bprp=bprp, P_EB=Peb,
                             b=bbest, chi2_lc=cbest, sec_ppt=sec*1e3 if np.isfinite(sec) else np.nan, chi2=cbest_tot,
                             T14_h=T14*24, T23_h=T23*24, aR=aR))
    print(m1, len(rows), flush=True)
D = pd.DataFrame(rows)
D['dchi2_vs_planet_free'] = D.chi2 - fit['chi2_best']
D.to_csv(f'{SD}/out/bound_eb_{ds}.csv', index=False)
best = D.sort_values('chi2').head(15)
print(best.to_string())
out = dict(dataset=ds, planet_free_chi2=fit['chi2_best'], n_models=len(D), best=best.to_dict(orient='records'),
           best_by_class={})
for lab, m in [('P5.69', np.abs(D.P_EB - P) < 1e-6), ('P11.38', np.abs(D.P_EB - 2*P) < 1e-6),
               ('occulter_BD', D.typ2 == 'BD'), ('occulter_MS', D.typ2 == 'MS'),
               ('no_third_light', D.m0.isna()), ('with_third_light', D.m0.notna()),
               ('twin_like_m0_ge_0.2', D.m0 >= 0.2), ('dilute_fsrc_lt_0.1', D.f_src < 0.1)]:
    if m.sum():
        out['best_by_class'][lab] = D[m].sort_values('chi2').iloc[0].to_dict()
json.dump(out, open(f'{SD}/out/bound_eb_{ds}.json', 'w'), indent=1, default=float)
