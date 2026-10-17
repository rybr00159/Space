"""r01b: diagnose why the (dMKs, r(G-Ks)) and (dMG, r(G-W1)) planes give different twin fractions among
over-luminous clean stars (r01). Uses the same local ridges as r01 (re-derived by importing r01's definitions).
Output: out/r01b_planes.json (printed summary)"""
import json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
H = 'tic214311908/analysis/host'
B = 'tic214311908/analysis/referee2_host'
R = {k: (v['slope'], v['zero_at_target_colour']) for k, v in json.load(open(f'{B}/out/r01_population.json'))['ridge_local'].items()}
C0 = 2.756149
F = pd.read_csv(f'{H}/data/field_mdwarfs.csv')
x = F.bp_rp; fx = 1.162004 + 0.011464 * x + 0.049255 * x ** 2 - 0.005879 * x ** 3
F['Cs'] = (F.phot_bp_rp_excess_factor - fx) / (0.0059898 + 8.817481e-12 * F.phot_g_mean_mag ** 7.618399)
dm = 5 * np.log10(F.parallax / 100)
F['MK'] = F.ks_m + dm; F['MG'] = F.phot_g_mean_mag + dm; F['GK'] = F.phot_g_mean_mag - F.ks_m; F['GW1'] = F.phot_g_mean_mag - F.w1mpro; F['KW1'] = F.ks_m - F.w1mpro
F['GJ'] = F.phot_g_mean_mag - F.j_m; F['GH'] = F.phot_g_mean_mag - F.h_m; F['JK'] = F.j_m - F.ks_m; F['HK'] = F.h_m - F.ks_m
for q in R: F['d' + q] = F[q] - (R[q][0] * (F.bp_rp - C0) + R[q][1])
clean = (F.ruwe < 1.1) & (F.ipd_frac_multi_peak == 0) & (F.ipd_gof_harmonic_amplitude < 0.03) & (F.parallax / F.parallax_error > 20) & \
        (F.ph_qual == 'AAA') & (F.ks_msigcom < 0.05) & (F.Cs.abs() < 3) & ((F.bp_rp - C0).abs() < 0.08)
wok = F.wise_cc.astype(str).str[:2].eq('00') & F.wise_ph_qual.astype(str).str[:2].eq('AA') & (F.wise_ext == 0)
D = F[clean].copy(); Dw = F[clean & wok].copy()
out = {}
for lab, d in (('K-plane window dMK -0.85..-0.65 (all clean)', D[D.dMK.between(-0.85, -0.65)]), ('K-plane window, WISE-ok subset', Dw[Dw.dMK.between(-0.85, -0.65)]),
               ('W-plane window dMG -0.785..-0.585 (WISE-ok)', Dw[Dw.dMG.between(-0.785, -0.585)]), ('dMG window, all clean', D[D.dMG.between(-0.785, -0.585)]),
               ('normal dMK -0.1..0.1', D[D.dMK.between(-0.1, 0.1)]), ('normal dMK, WISE-ok', Dw[Dw.dMK.between(-0.1, 0.1)])):
    s = {}
    for q in ('dGK', 'dGW1', 'dKW1', 'dGJ', 'dGH', 'dJK', 'dMK', 'dMG'):
        v = d[q].dropna().values
        s[q] = (round(float(np.median(v)), 3), round(float(1.4826 * np.median(np.abs(v - np.median(v)))), 3)) if len(v) else None
    s['n'] = len(d); s['frac_rGK_lt_0.03'] = round(float(np.mean(d.dGK < 0.03)), 3); s['frac_rGW1_lt_0.03'] = round(float(np.mean(d.dGW1.dropna() < 0.03)), 3)
    s['dMG_minus_dMK_median'] = round(float(np.median(d.dMG - d.dMK)), 3)
    out[lab] = s; print(lab, s)
# stars in the W-plane window with low r(G-W1): what are their r(G-Ks), dMK?
w = Dw[Dw.dMG.between(-0.785, -0.585)]
lo = w[w.dGW1 < 0.03]; hi = w[w.dGW1 >= 0.03]
print('W-window low-rGW1 stars:', len(lo), 'median dMK', np.median(lo.dMK), 'median rGK', np.median(lo.dGK), 'median rKW1', np.median(lo.dKW1))
print('W-window high-rGW1 stars:', len(hi), 'median dMK', np.median(hi.dMK), 'median rGK', np.median(hi.dGK), 'median rKW1', np.median(hi.dKW1))
out['W_window_low_rGW1'] = dict(n=len(lo), med_dMK=float(np.median(lo.dMK)), med_rGK=float(np.median(lo.dGK)), med_rKW1=float(np.median(lo.dKW1)), med_dMG=float(np.median(lo.dMG)))
out['W_window_high_rGW1'] = dict(n=len(hi), med_dMK=float(np.median(hi.dMK)), med_rGK=float(np.median(hi.dGK)), med_rKW1=float(np.median(hi.dKW1)), med_dMG=float(np.median(hi.dMG)))
# correlation of rGK and rGW1 among normal stars and over-luminous stars
for lab, d in (('normal', Dw[Dw.dMK.between(-0.2, 0.2)]), ('overlum', Dw[Dw.dMK.between(-0.95, -0.55)])):
    print(lab, 'corr(rGK, rGW1)=', np.corrcoef(d.dGK, d.dGW1)[0, 1], 'sd rKW1', 1.4826 * np.median(np.abs(d.dKW1 - np.median(d.dKW1))), 'median rKW1', np.median(d.dKW1))
    out[f'corr_{lab}'] = float(np.corrcoef(d.dGK, d.dGW1)[0, 1])
# 2D: dMG - dMK (= -r(G-Ks) + const) distribution: twins -> ~0, metal-rich -> positive
json.dump(out, open(f'{B}/out/r01b_planes.json', 'w'), indent=1)
