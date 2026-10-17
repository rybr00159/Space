"""f01: pixel-localisation likelihoods for the FPP calculation.
Inputs (read-only): referee_loc/out/inj_base.npz (1000 injections of the observed signal on each of 7 sources through
the referee's full image pipeline; real combined chi2 differences; real combined Delta-chi2 sky map on a 0.5" grid,
M-dwarf colour frame, co-moving with the M3V), referee_loc/out/summary.json (23 analysis variants),
loc/out/post_*.json (investigator's kappa-calibrated posteriors, DIFF and AMP).

1. Catalogued alternatives X (pair A, pair B, Mira, G17.5 star 3.1", G17.7 star 5.5", G17.0 star 12"): statistic
   s = chi2(X) - chi2(M3V). Likelihood ratio LR_X = p(s_obs | X true) / p(s_obs | M3V true), with each density
   estimated from the 1000 injections by (a) a Gaussian fit and (b) a Student-t fit (heavier tails; adopted, it is
   the conservative choice when s_obs lies beyond all injections of X). The empirical tail fraction P(s >= s_obs | X)
   is reported too (resolution 1e-3).
2. Uncatalogued source at sky offset x: L(x)/L(M3V) = exp(-Delta chi2_map(x) / (2 kappa)). kappa calibrates the
   heavy-tailed noise: the referee's injection thresholds give kappa = 4.59/2.30 = 2.0 (68%), 13.80/6.18 = 2.2 (95%) and
   56.6/11.83 = 4.8 (99.7%). Effective area A_eff = sum over the map of L(x)/L(M3V) * 0.25 arcsec^2, by annulus.
Output: out/loc_likelihoods.json"""
import numpy as np, json, sys
from scipy import stats
sys.path.insert(0, 'tic214311908/analysis/fpp/scripts')
from fpplib import *
z = np.load(f'{B}/referee_loc/out/inj_base.npz')
MAIN = ['M3V', 'pairA', 'pairB', 'Mira', 'G175', 'G177', 'G170']
real = z['real_comb_dchi2']
out = dict(method=__doc__.split('\n')[0], catalogued={}, uncatalogued={})
for j, X in enumerate(MAIN):
    if X == 'M3V': continue
    sM = z['M3V__comb_c'][:, j] - z['M3V__comb_c'][:, 0]
    sX = z[f'{X}__comb_c'][:, j] - z[f'{X}__comb_c'][:, 0]
    so = float(real[j])
    g_M = stats.norm(*stats.norm.fit(sM)); g_X = stats.norm(*stats.norm.fit(sX))
    tpX = stats.t.fit(sX); tpM = stats.t.fit(sM)
    # guard against nu -> inf / tiny nu
    t_M = stats.t(*tpM); t_X = stats.t(*tpX)
    lr_g = g_X.logpdf(so) - g_M.logpdf(so)
    lr_t = t_X.logpdf(so) - t_M.logpdf(so)
    out['catalogued'][X] = dict(s_obs=so, M3V_truth_quantiles=np.percentile(sM, [2.5, 50, 97.5]).tolist(),
                                X_truth_quantiles=np.percentile(sX, [2.5, 50, 97.5, 99.85]).tolist(),
                                frac_X_truth_ge_obs=float(np.mean(sX >= so)), frac_M3V_truth_le_obs=float(np.mean(sM <= so)),
                                log10LR_gauss=float(lr_g/np.log(10)), log10LR_t=float(lr_t/np.log(10)),
                                t_params_X=list(map(float, tpX)), t_params_M3V=list(map(float, tpM)))
    print(X, json.dumps(out['catalogued'][X])[:400])
# investigator / referee variants for the 3.1" star (G175)
summ = jload('referee_loc/out/summary.json')
out['catalogued_variants_note'] = ('referee_loc summary: Gaussian log10 L(M3V)/L(G17.5) = 3.9 at baseline, 1.4-4.2 over 23 variants '
                                   '(1.4 only if the star had BP-RP 2.76); loc AMP kappa-calibrated posterior P(G17.5) = 0.023 vs P(M3V) = 0.977 '
                                   '(log10 LR = 1.63); loc DIFF P(G17.5) = 8e-6')
out['G175_log10LR_range_for_sensitivity'] = [1.4, 1.63, 3.9, 4.2]
# ---- uncatalogued: effective area
GE, GN = z['GE'], z['GN']
mp = z['real_gsum'].reshape(len(GN), len(GE)).astype(float)   # chi2(x) - chi2(M3V hypothesis)
EE, NN = np.meshgrid(GE, GN); RR = np.hypot(EE, NN)
cell = 0.25
ann = [(0, 1), (1, 1.25), (1.25, 3), (3, 5), (5, 8), (8, 15), (15, 40)]
for kap in (1.0, 1.5, 2.0, 2.2, 3.0, 4.8):
    L = np.exp(-np.clip(mp, -50, None)/(2*kap))
    tot = float(L.sum()*cell)
    by = {f'{a}-{b}': float(L[(RR >= a) & (RR < b)].sum()*cell) for a, b in ann}
    out['uncatalogued'][f'kappa_{kap}'] = dict(A_eff_total_arcsec2=tot, A_eff_by_annulus=by)
    print('kappa', kap, 'A_eff', round(tot, 2), {k: round(v, 2) for k, v in by.items()})
out['map_value_at_M3V_position'] = float(mp[np.argmin(np.abs(GN)), np.argmin(np.abs(GE))])
out['map_min'] = float(mp.min())
out['kappa_baseline'] = 2.2
jdump(out, 'out/loc_likelihoods.json')
