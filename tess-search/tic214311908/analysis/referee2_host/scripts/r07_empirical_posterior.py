"""r07: P(twin) for the M3V without metallicity labels, without a fitted mixture and without a CMD scatter model.
P(twin | data) = f_T * n(x_t + dx, y_t + dy) / n(x_t, y_t)
  n(.)  = empirical KDE density of ALL clean field stars at the M3V colour (singles + every kind of binary), so the
          single-star density is not modelled at all;
  f_T   = literature prior fraction of systems with a q >= 0.8 companion that survives the M3V's Gaia data
          (RUWE, PM anomaly, RV error, IPD): 0.8%, 1.5%, 2.35% (= analyst's own h15 baseline, generous) and 3.5%
          (x1.5 triple enhancement);
  (dx, dy) = twin component position: dx = +0.753 (q=1; at fixed combined BP-RP), dy = 0 (q=1) or +0.015 (q=0.8).
Variables (r01 local ridges, raw parallaxes): x = dMKs, y = r(G-Ks) [2MASS plane]; x = dMG, y = r(G-W1) [WISE
plane, AllWISE-clean stars]; 3D (dMKs, r(G-Ks), r(G-W1)). XP-shape [Fe/H] (analyst's regression applied to field
stars, h12b_xp_pred_field.csv) is added as a conditional factor p(xp | r(G-Ks), twin) / p(xp | r(G-Ks), all at x_t)
with linear-Gaussian conditionals fitted on normal (|dMK|<0.15) and window (dMK -0.9..-0.6) stars.
Wide-pair factor (h14 likelihood at -0.11 vs +0.29): 0.36-0.45, applied optionally (biased high in [Fe/H] if A is
itself an unresolved binary: A has RUWE 13.7).
Photometric systematics: Ks 0.04 fainter; W1 0.03 fainter.
Bootstrap (200) over field stars for the KDE ratios.
Output: out/r07_empirical_posterior.json"""
import json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from scipy.stats import norm
H = 'tic214311908/analysis/host'
B = 'tic214311908/analysis/referee2_host'
r01 = json.load(open(f'{B}/out/r01_population.json'))
R = {k: (v['slope'], v['zero_at_target_colour']) for k, v in r01['ridge_local'].items()}
C0 = 2.756149
F = pd.read_csv(f'{H}/data/field_mdwarfs.csv')
x = F.bp_rp; fx = 1.162004 + 0.011464 * x + 0.049255 * x ** 2 - 0.005879 * x ** 3
F['Cs'] = (F.phot_bp_rp_excess_factor - fx) / (0.0059898 + 8.817481e-12 * F.phot_g_mean_mag ** 7.618399)
dm = 5 * np.log10(F.parallax / 100)
F['MK'] = F.ks_m + dm; F['MG'] = F.phot_g_mean_mag + dm; F['GK'] = F.phot_g_mean_mag - F.ks_m; F['GW1'] = F.phot_g_mean_mag - F.w1mpro
for q in ('MK', 'MG', 'GK', 'GW1'): F['d' + q] = F[q] - (R[q][0] * (F.bp_rp - C0) + R[q][1])
F['vtan'] = 4.74047 * np.hypot(F.pmra, F.pmdec) / F.parallax
clean = (F.ruwe < 1.1) & (F.ipd_frac_multi_peak == 0) & (F.ipd_gof_harmonic_amplitude < 0.03) & (F.parallax / F.parallax_error > 20) & \
        (F.ph_qual == 'AAA') & (F.ks_msigcom < 0.05) & (F.Cs.abs() < 3) & ((F.bp_rp - C0).abs() < 0.08)
wok = F.wise_cc.astype(str).str[:2].eq('00') & F.wise_ph_qual.astype(str).str[:2].eq('AA') & (F.wise_ext == 0)
DK = F[clean & np.isfinite(F.dGK)].reset_index(drop=True)
DW = F[clean & wok & np.isfinite(F.dGW1) & np.isfinite(F.dGK)].reset_index(drop=True)
T = r01['target']; tK = T['cmd']['raw_plx']['dMK']; tG = T['cmd']['raw_plx']['dMG']; yK = T['colour_residuals']['dGK']; yW = T['colour_residuals']['dGW1']
print(f'target dMK {tK:+.3f} dMG {tG:+.3f} r(G-Ks) {yK:+.3f} r(G-W1) {yW:+.3f}')
def kde_ratio(X, t, shift, bw):
    """ratio n(t+shift)/n(t) for a product-Gaussian KDE with bandwidths bw over the rows of X"""
    def dens(p): return np.mean(np.exp(-0.5 * np.sum(((X - p) / bw) ** 2, axis=1)))
    a = dens(np.asarray(t) + np.asarray(shift)); b = dens(np.asarray(t)); return a / b, a, b
def boot(X, t, shift, bw, n=200, seed=5):
    rng = np.random.default_rng(seed); r0 = kde_ratio(X, t, shift, bw)[0]; bs = []
    for _ in range(n): bs.append(kde_ratio(X[rng.integers(0, len(X), len(X))], t, shift, bw)[0])
    return r0, np.percentile(bs, [16, 50, 84])
out = dict(target=dict(dMK=tK, dMG=tG, rGK=yK, rGW1=yW))
FT = (0.008, 0.015, 0.0235, 0.035)
res = {}
def record(lab, X, t, shift, bw):
    r0, pc = boot(X, t, shift, bw)
    res[lab] = dict(ratio=float(r0), ratio_16_50_84=[float(v) for v in pc], n=int(len(X)), **{f'P_twin_fT{f}': float(f * r0 / (1 + 0 * r0)) for f in FT})
    print(f"{lab:62s} ratio {r0:8.4f} [{pc[0]:.4f},{pc[2]:.4f}]  P(twin) fT=0.8%: {0.008 * r0:.4f}  1.5%: {0.015 * r0:.4f}  2.35%: {0.0235 * r0:.4f}  3.5%: {0.035 * r0:.4f}", flush=True)
# 1D CMD only
record('CMD only dMKs', DK[['dMK']].values, [tK], [0.753], [0.06])
record('CMD only dMG', DW[['dMG']].values, [tG], [0.753], [0.06])
# 2D
for dy in (0.0, 0.015):
    record(f'dMKs + r(G-Ks), dy={dy}', DK[['dMK', 'dGK']].values, [tK, yK], [0.753, dy], [0.06, 0.012])
    record(f'dMG + r(G-W1), dy={dy}', DW[['dMG', 'dGW1']].values, [tG, yW], [0.753, dy], [0.06, 0.012])
    record(f'3D dMKs + r(G-Ks) + r(G-W1), dy={dy}', DW[['dMK', 'dGK', 'dGW1']].values, [tK, yK, yW], [0.753, dy, dy], [0.06, 0.012, 0.012])
# bandwidth sensitivity
record('dMKs + r(G-Ks), bw y 0.02', DK[['dMK', 'dGK']].values, [tK, yK], [0.753, 0.0], [0.08, 0.02])
record('3D, bw 0.08/0.02/0.02', DW[['dMK', 'dGK', 'dGW1']].values, [tK, yK, yW], [0.753, 0.0, 0.0], [0.08, 0.02, 0.02])
# photometric systematics
record('dMKs + r(G-Ks), Ks 0.04 fainter', DK[['dMK', 'dGK']].values, [tK + 0.04, yK - 0.04], [0.753, 0.0], [0.06, 0.012])
record('dMG + r(G-W1), W1 0.03 fainter', DW[['dMG', 'dGW1']].values, [tG, yW - 0.03], [0.753, 0.0], [0.06, 0.012])
record('3D, Ks 0.04 and W1 0.03 fainter', DW[['dMK', 'dGK', 'dGW1']].values, [tK + 0.04, yK - 0.04, yW - 0.03], [0.753, 0.0, 0.0], [0.06, 0.012, 0.012])
record('3D, Ks 0.03 fainter', DW[['dMK', 'dGK', 'dGW1']].values, [tK + 0.03, yK - 0.03, yW], [0.753, 0.0, 0.0], [0.06, 0.012, 0.012])
# kinematic subsample (v_tan 30-60, target 43.5)
kin = DW.vtan.between(30, 60)
record('3D, v_tan 30-60 subsample', DW[kin][['dMK', 'dGK', 'dGW1']].values, [tK, yK, yW], [0.753, 0.0, 0.0], [0.06, 0.012, 0.012])
kinK = DK.vtan.between(30, 60)
record('dMKs + r(G-Ks), v_tan 30-60', DK[kinK][['dMK', 'dGK']].values, [tK, yK], [0.753, 0.0], [0.06, 0.012])
out['kde_ratios'] = res
# ---------- XP conditional factor
P = pd.read_csv(f'{H}/out/h12b_xp_pred_field.csv').set_index('source_id')
Fi = F.set_index('source_id')
P = P.join(Fi[['dMK', 'dGK', 'ph_qual', 'ks_msigcom']], how='left', rsuffix='_r')
P = P[(P.ph_qual == 'AAA') & (P.ks_msigcom < 0.05) & np.isfinite(P.dGK) & np.isfinite(P.feh_xp)]
xp_t, exp_t = 0.12922249223373056, 0.0772
def condfit(d):
    A = np.vstack([np.ones(len(d)), d.dGK.values]).T; c, *_ = np.linalg.lstsq(A, d.feh_xp.values, rcond=None)
    r = d.feh_xp.values - A @ c; s = 1.4826 * np.median(np.abs(r - np.median(r))); return c, s
nor = P[P.dMK.abs() < 0.15]; win = P[P.dMK.between(-0.9, -0.6)]
xpres = {}
for dy in (0.0, 0.015):
    cn, sn = condfit(nor); cw, sw = condfit(win)
    pt = norm.pdf(xp_t, cn[0] + cn[1] * (yK + dy), sn); pa = norm.pdf(xp_t, cw[0] + cw[1] * yK, sw)
    rng = np.random.default_rng(9); bs = []
    for _ in range(300):
        nb = nor.iloc[rng.integers(0, len(nor), len(nor))]; wb = win.iloc[rng.integers(0, len(win), len(win))]
        a, sa = condfit(nb); b, sb = condfit(wb); bs.append(norm.pdf(xp_t, a[0] + a[1] * (yK + dy), sa) / norm.pdf(xp_t, b[0] + b[1] * yK, sb))
    xpres[f'dy{dy}'] = dict(n_normal=len(nor), n_window=len(win), normal_fit=[float(v) for v in cn], normal_sd=float(sn), window_fit=[float(v) for v in cw], window_sd=float(sw),
                            pred_twin=float(cn[0] + cn[1] * (yK + dy)), pred_all=float(cw[0] + cw[1] * yK), factor=float(pt / pa), factor_16_50_84=[float(v) for v in np.percentile(bs, [16, 50, 84])])
    print('XP conditional factor', dy, xpres[f'dy{dy}'])
# XP-only factor (not conditioned on colour) for comparison
xpres['unconditional'] = dict(factor=float(norm.pdf(xp_t, np.median(nor.feh_xp), 1.4826 * np.median(np.abs(nor.feh_xp - np.median(nor.feh_xp)))) /
                                           norm.pdf(xp_t, np.median(win.feh_xp), 1.4826 * np.median(np.abs(win.feh_xp - np.median(win.feh_xp))))))
print('XP unconditional factor', xpres['unconditional'])
out['xp_conditional'] = xpres
# ---------- combined
def comb(ratio, fT, xpf=1.0, wide=1.0):
    o = fT * ratio * xpf * wide; return o / (1 + o - fT * ratio * xpf * wide + fT * ratio * xpf * wide) if False else o
summ = {}
base3 = res['3D dMKs + r(G-Ks) + r(G-W1), dy=0.0']['ratio']; base2 = res['dMKs + r(G-Ks), dy=0.0']['ratio']
xpf = xpres['dy0.0']['factor']
for lab, rr in (('2D (dMKs, r(G-Ks))', base2), ('3D (dMKs, r(G-Ks), r(G-W1))', base3), ('3D, Ks+0.04 & W1+0.03', res['3D, Ks 0.04 and W1 0.03 fainter']['ratio']),
                ('3D, Ks+0.03', res['3D, Ks 0.03 fainter']['ratio']), ('3D, vtan 30-60', res['3D, v_tan 30-60 subsample']['ratio']), ('CMD only', res['CMD only dMKs']['ratio'])):
    for fT in FT:
        for xl, xf in (('noXP', 1.0), ('XP', xpf)):
            for wl, wf in (('nowide', 1.0), ('wide0.45', 0.45)):
                summ[f'{lab} | fT {fT} | {xl} | {wl}'] = float(fT * rr * xf * wf)
out['P_twin_combined'] = summ
for k, v in summ.items():
    if ('fT 0.015' in k or 'fT 0.0235' in k): print(f'{k:75s} P(twin) = {v:.4f}')
json.dump(out, open(f'{B}/out/r07_empirical_posterior.json', 'w'), indent=1)
