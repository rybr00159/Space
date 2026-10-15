"""Predicted eclipse-depth ratios depth(band)/depth(TESS) for
 (a) a bound K/M pair at the G star's distance (hierarchical triple), pairs K5+K7 ... M2+M3, and
 (b) a background EB whose eclipsed star has G/F colour (and, generally, as a function of Teff).

Model
  * Spectra: BT-Settl CIFIST grid, warped (smooth multiplicative correction) so their synthetic B-V, V-Rc,
    V-Ic, V-Ks equal the Mamajek (2022.04) dwarf colours; validated against independent Gaia Bp-Rp, G-V
    and the Stassun TIC-8 T-G relation (differences <~0.04 mag, see sed3 check in the log).
  * Absolute scale: Mamajek M_V for each pair component; the G star = observed system (V = 11.857 APASS,
    d = 297.4 pc, A_V = 0.02) minus the pair, with a warped spectrum at Teff 5883 K.
  * Geometry: circular orbit, P = 1.3036712 d, a from Kepler with Mamajek masses, Mamajek radii,
    quadratic LD per band (Claret 2011/2017); inclination solved so that the TESS deep-eclipse
    template amplitude equals the observed one.  Light curves are computed in every band and fitted
    with the SAME TESS template as the data, so "depth" means the same thing for data and model.
"""
import sys, json
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/v3/scripts')
import lib3 as L
import sed3 as S

BANDS = ['TESS', 'WASP', 'WASP_flat', 'WASP_red', 'WASP_blue', 'V', 'B', 'g', 'R', 'I', 'Ks', 'W1', 'W2']
GSUN, RSUN, MSUN = 6.674e-8, 6.957e10, 1.989e33
V_SYS, DIST, AV = 11.857, 297.4, 0.02
T_G = 5883.0


def ld_for(teff, logg, b):
    return S.ld_quadratic(teff, logg, b)


def template_amp(tgrid, y, shp, which):
    s = L.template(tgrid, 0.0, shp[which]['T'], shp[which]['p'], shp[which]['q'], t0=0.0, per=L.P)
    return float(np.sum(y * s) / np.sum(s * s))


def lightcurves(R1, R2, a, inc, ld1, ld2, tgrid):
    """fraction of each star's own flux hidden during deep (2 in front of 1) and shallow (1 in front of 2)
    eclipses, vs time from conjunction"""
    ph = 2 * np.pi * tgrid / L.P
    d = a * np.sqrt(np.sin(ph) ** 2 + (np.cos(np.radians(inc)) * np.cos(ph)) ** 2)
    k1 = np.array([S.occulted_fraction(x, R1, R2, ld1, nr=1500) for x in d])
    k2 = np.array([S.occulted_fraction(x, R2, R1, ld2, nr=1500) for x in d])
    return k1, k2


def main():
    tab = S.mamajek()
    T = json.load(open(L.V3 + '/results/tess3.json'))
    shp = T['shapes']
    D_T, S_T = T['adopted_deep']['mean'], T['adopted_shallow']['mean']
    tgrid = np.linspace(-0.12, 0.12, 241)
    # ---- G star
    gm = S.star_model(teff=T_G, tab=tab)
    out = dict(TESS_deep=D_T, TESS_shallow=S_T, bands=BANDS, pairs={}, background={})
    Mv_sys = V_SYS - 5 * np.log10(DIST / 10) - AV
    FV_sys = S.ZP['V'] * 10 ** (-0.4 * Mv_sys)
    rel = lambda m, b: S.band_flux(m['spec'], b) / S.band_flux(m['spec'], 'V')
    cache = {}

    def comp(spt):
        if spt not in cache:
            cache[spt] = S.star_model(spt, tab=tab)
        return cache[spt]

    pairs = [('K5V', 'K7V'), ('K5V', 'M0V'), ('K7V', 'M0V'), ('K7V', 'M1V'), ('K7V', 'M2V'), ('K8V', 'M2V'),
             ('M0V', 'M1V'), ('M0V', 'M2V'), ('M0V', 'M3V'), ('M1V', 'M2V'), ('M1V', 'M3V'), ('M2V', 'M3V'),
             ('K5V', 'M1V'), ('K6V', 'M1V'), ('K5V', 'M2V')]
    print('%-8s %6s %6s | fT_1 fT_2 | inc   FWHM_h | sT_pred/obs | dKs  dI  dV | WASP_d WASP_s V_d V_s g_d Ks_d W1_d' % ('pair', 'Teff1', 'Teff2'))
    for s1, s2 in pairs:
        c1, c2 = comp(s1), comp(s2)
        r1, r2 = c1['row'], c2['row']
        F1V = S.ZP['V'] * 10 ** (-0.4 * r1['Mv']); F2V = S.ZP['V'] * 10 ** (-0.4 * r2['Mv'])
        FGV = FV_sys - F1V - F2V
        F = {}
        for b in BANDS:
            F[b] = (FGV * rel(gm, b), F1V * rel(c1, b), F2V * rel(c2, b))
        frac = {b: (F[b][1] / sum(F[b]), F[b][2] / sum(F[b])) for b in BANDS}
        Mtot = (r1['Msun'] + r2['Msun']) * MSUN
        a = (GSUN * Mtot * (L.P * 86400) ** 2 / (4 * np.pi ** 2)) ** (1 / 3) / RSUN
        R1, R2 = r1['R_Rsun'], r2['R_Rsun']
        lg1 = np.log10(GSUN * r1['Msun'] * MSUN / (R1 * RSUN) ** 2); lg2 = np.log10(GSUN * r2['Msun'] * MSUN / (R2 * RSUN) ** 2)
        ldT1, ldT2 = S.ld_quadratic(r1['Teff'], lg1, 'TESS'), S.ld_quadratic(r2['Teff'], lg2, 'TESS')

        def deepT(inc):
            k1, k2 = lightcurves(R1, R2, a, inc, ldT1, ldT2, tgrid)
            return template_amp(tgrid, k1 * frac['TESS'][0], shp, 'd')
        imin = np.degrees(np.arccos(min((R1 + R2) / a, 1.0)))
        incs = np.linspace(imin + 0.01, 89.99, 120)
        dv = np.array([deepT(i) for i in incs])
        row = dict(spt1=s1, spt2=s2, Teff1=r1['Teff'], Teff2=r2['Teff'], R1=R1, R2=R2, M1=r1['Msun'], M2=r2['Msun'], a_Rsun=a,
                   frac_TESS=list(frac['TESS']), frac_V=list(frac['V']), frac_WASP=list(frac['WASP']), frac_Ks=list(frac['Ks']))
        if dv.max() < D_T:
            row['feasible'] = False
            row['max_TESS_deep'] = float(dv.max())
            out['pairs'][s1 + '+' + s2] = row
            print('%-8s infeasible: max TESS deep %.0f ppm < %.0f' % (s1 + '+' + s2, dv.max() * 1e6, D_T * 1e6))
            continue
        # lowest inclination (most grazing) solution; partial eclipses are V-shaped as observed
        j = np.where(dv >= D_T)[0][0]
        inc = brentq(lambda i: deepT(i) - D_T, incs[max(j - 1, 0)], incs[j])
        row['feasible'] = True; row['inc'] = inc; row['b_impact'] = float(a * np.cos(np.radians(inc)) / (R1 + R2))
        amps = {}
        for b in BANDS:
            ld1 = ld_for(r1['Teff'], lg1, b); ld2 = ld_for(r2['Teff'], lg2, b)
            k1, k2 = lightcurves(R1, R2, a, inc, ld1, ld2, tgrid)
            amps[b] = (template_amp(tgrid, k1 * frac[b][0], shp, 'd'), template_amp(tgrid, k2 * frac[b][1], shp, 's'))
            if b == 'TESS':
                y = k1 * frac[b][0]
                hm = tgrid[y >= 0.5 * y.max()]
                row['FWHM_h_deep'] = float((hm.max() - hm.min()) * 24)
                row['undiluted_deep_TESS'] = float(k1.max())
        row['amps'] = {b: list(v) for b, v in amps.items()}
        row['ratio_deep'] = {b: amps[b][0] / amps['TESS'][0] for b in BANDS}
        row['ratio_shallow'] = {b: amps[b][1] / amps['TESS'][1] for b in BANDS}
        row['TESS_shallow_pred'] = amps['TESS'][1]
        # contrasts of the pair (combined) vs G star
        for b in ('Ks', 'I', 'V', 'TESS', 'RP', 'BP', 'G') if False else ('Ks', 'I', 'V', 'TESS'):
            row['d' + b] = float(-2.5 * np.log10((F[b][1] + F[b][2]) / F[b][0]))
        row['pair_frac_TESS'] = float(frac['TESS'][0] + frac['TESS'][1])
        out['pairs'][s1 + '+' + s2] = row
        print('%-8s %6.0f %6.0f | %.3f %.3f | %5.2f %.2f | %.0f/%.0f | %.2f %.2f %.2f | %.3f %.3f %.3f %.3f %.3f %.3f %.3f' % (
            s1 + '+' + s2, r1['Teff'], r2['Teff'], frac['TESS'][0], frac['TESS'][1], inc, row['FWHM_h_deep'],
            amps['TESS'][1] * 1e6, S_T * 1e6, row['dKs'], row['dI'], row['dV'],
            row['ratio_deep']['WASP'], row['ratio_shallow']['WASP'], row['ratio_deep']['V'], row['ratio_shallow']['V'],
            row['ratio_deep']['g'], row['ratio_deep']['Ks'], row['ratio_deep']['W1']), flush=True)

    # ---- (b) background / generic eclipsed star of temperature Teff1, small light fraction
    # ratio = (F1,b/F1,T)/(Ftot,b/Ftot,T) * kappa_b/kappa_T ; Ftot = system (G star dominated)
    print('\nGeneric eclipsed star (light fraction << 1), ratio depth(b)/depth(TESS) for grazing / mid / central geometry')
    geo = {'grazing': 0.95, 'mid': 0.6, 'deep': 0.25}     # d_min/(R1+R2), R2 = 0.8 R1
    teffs = [3400, 3600, 3850, 4100, 4440, 4830, 5270, 5600, 5900, 6200, 6600, 7000]
    for te in teffs:
        m = S.star_model(teff=te, tab=tab)
        lg = 4.5 if te > 4500 else 4.8
        row = {}
        for gname, x in geo.items():
            d = x * 1.8
            kap = {b: S.occulted_fraction(d, 1.0, 0.8, ld_for(te, lg, b)) for b in BANDS}
            row[gname] = {b: float(rel(m, b) / rel(gm, b) * rel(gm, 'TESS') / rel(m, 'TESS') * kap[b] / kap['TESS']) for b in BANDS}
        out['background'][str(te)] = row
        print('Teff %5d  WASP %.3f/%.3f/%.3f   V %.3f/%.3f/%.3f   g %.3f   Ks %.3f  W1 %.3f' % (
            te, row['grazing']['WASP'], row['mid']['WASP'], row['deep']['WASP'], row['grazing']['V'], row['mid']['V'], row['deep']['V'],
            row['mid']['g'], row['mid']['Ks'], row['mid']['W1']))
    # system colour check (bound K7+M1 case vs observed): V-Ks observed 11.857-10.341 = 1.516
    json.dump(out, open(L.V3 + '/results/predict3.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
