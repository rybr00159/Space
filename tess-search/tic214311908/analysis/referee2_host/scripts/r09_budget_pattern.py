"""r09: (a) companion light budget: how much of the 0.75-mag over-luminosity could an unresolved companion supply?
Empirical density ratio n(x_t + D, y_t) / n(x_t, y_t) in the (dMKs, r(G-Ks)) plane of clean field stars as a function
of the companion's contribution D (D = 0.753 -> twin; D -> 0 -> no companion), i.e. the relative likelihood that the
PRIMARY sits at dMKs = x_t + D with the observed colour residual (dy = 0 approximation; |dy| <= 0.015 for q 0.6-1).
Companion flux fraction in Ks = 1 - 10^(-0.4 D); TESS-band fraction approximated with the Mamajek MS flux ratios
(G/Ks average, as in h15).
(b) Wavelength pattern of the M3V's colour excess r(G-J):r(G-H):r(G-Ks):r(G-W1) compared with (i) the field
metallicity pattern (slopes of each residual against dMKs among clean stars, |dMK| < 1) and (ii) a contamination
pattern from the Mira's halo (excess flux proportional to Mira/M3V flux ratio x lambda; Mira magnitudes J 8.3, H 7.3,
Ks 6.68, W1 6.0 assumed - typical C/O Mira colours; only the shape matters): least-squares amplitude and chi2 for each.
Output: out/r09_budget_pattern.json"""
import json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
H = 'tic214311908/analysis/host'
B = 'tic214311908/analysis/referee2_host'
r01 = json.load(open(f'{B}/out/r01_population.json'))
R = {k: (v['slope'], v['zero_at_target_colour']) for k, v in r01['ridge_local'].items()}
C0 = 2.756149
F = pd.read_csv(f'{H}/data/field_mdwarfs.csv')
x = F.bp_rp; fx = 1.162004 + 0.011464 * x + 0.049255 * x ** 2 - 0.005879 * x ** 3
F['Cs'] = (F.phot_bp_rp_excess_factor - fx) / (0.0059898 + 8.817481e-12 * F.phot_g_mean_mag ** 7.618399)
dm = 5 * np.log10(F.parallax / 100)
F['MK'] = F.ks_m + dm; F['GK'] = F.phot_g_mean_mag - F.ks_m; F['GJ'] = F.phot_g_mean_mag - F.j_m; F['GH'] = F.phot_g_mean_mag - F.h_m; F['GW1'] = F.phot_g_mean_mag - F.w1mpro
for q in ('MK', 'GK', 'GJ', 'GH', 'GW1'): F['d' + q] = F[q] - (R[q][0] * (F.bp_rp - C0) + R[q][1])
clean = (F.ruwe < 1.1) & (F.ipd_frac_multi_peak == 0) & (F.ipd_gof_harmonic_amplitude < 0.03) & (F.parallax / F.parallax_error > 20) & \
        (F.ph_qual == 'AAA') & (F.ks_msigcom < 0.05) & (F.Cs.abs() < 3) & ((F.bp_rp - C0).abs() < 0.08)
wok = F.wise_cc.astype(str).str[:2].eq('00') & F.wise_ph_qual.astype(str).str[:2].eq('AA') & (F.wise_ext == 0)
D = F[clean & np.isfinite(F.dGK)]
T = r01['target']; tx = T['cmd']['raw_plx']['dMK']; ty = T['colour_residuals']['dGK']
X = D[['dMK', 'dGK']].values; bw = np.array([0.06, 0.012])
def dens(p): return np.mean(np.exp(-0.5 * np.sum(((X - p) / bw) ** 2, axis=1)))
d0 = dens(np.array([tx, ty])); bud = []
for Dd in np.arange(0.0, 0.76, 0.05):
    r = dens(np.array([tx + Dd, ty])) / d0
    bud.append(dict(D=float(Dd), ratio=float(r), companion_flux_frac_Ks=float(1 - 10 ** (-0.4 * Dd))))
    print(f'D {Dd:.2f}  companion Ks-flux fraction {1 - 10 ** (-0.4 * Dd):.2f}  density ratio {r:.3f}')
out = dict(budget=bud)
# (b) pattern
tv = np.array([T['colour_residuals']['dGJ'], T['colour_residuals']['dGH'], T['colour_residuals']['dGK'], T['colour_residuals']['dGW1']])
te = np.array([0.028, 0.031, 0.026, 0.03])
s = D.dMK.abs() < 1.0
slopes = []
for q in ('dGJ', 'dGH', 'dGK', 'dGW1'):
    k = s & np.isfinite(D[q]) & ((q != 'dGW1') | wok[D.index])
    slopes.append(np.polyfit(D.dMK[k], D[q][k], 1)[0])
pm = np.array(slopes)
lam = np.array([1.25, 1.65, 2.17, 3.4])
mira = np.array([8.3, 7.3, 6.68, 6.0]); m3v = np.array([11.363, 10.740, 10.482, 10.318])
pc = 10 ** (-0.4 * (mira - m3v)) * lam            # excess flux fraction pattern (arbitrary units)
pc_mag = pc / pc[2]                               # small contamination -> magnitude excess proportional to flux fraction
res = {}
for lab, p in (('metallicity (field slopes)', pm / pm[2]), ('Mira halo contamination', pc_mag), ('Mira, no lambda scaling', (10 ** (-0.4 * (mira - m3v))) / 10 ** (-0.4 * (mira[2] - m3v[2])))):
    A = np.sum(tv * p / te ** 2) / np.sum(p ** 2 / te ** 2); chi = float(np.sum(((tv - A * p) / te) ** 2))
    res[lab] = dict(pattern_rel_to_Ks=[float(v) for v in p], amplitude_at_Ks=float(A), chi2_4pts=chi)
    print(lab, res[lab])
out['pattern'] = dict(target_residuals_J_H_Ks_W1=tv.tolist(), errors=te.tolist(), field_slopes_vs_dMK=pm.tolist(), fits=res,
                      note='J/H/Ks share the 2MASS epoch; errors treated as independent; Mira magnitudes are assumed typical values (shape-only test)')
json.dump(out, open(f'{B}/out/r09_budget_pattern.json', 'w'), indent=1)
