"""r01: label-free, empirical test of 'metal-rich single' vs 'near-equal twin' for the M3V using the field population
(the analyst's raw Gaia-archive query host/data/field_mdwarfs.csv is re-used as raw data; all derived quantities are
recomputed here with an independent, LOCAL ridge definition, not the analyst's global 4th-order polynomials).

Ridge: stars within |BP-RP - 2.756| < 0.12 (clean: RUWE<1.1, ipd_frac_multi_peak==0, ipd_gof_harmonic_amplitude<0.03,
plx/err>20, 2MASS AAA, ks_msigcom<0.05, |C*|<3). For M_Ks and M_G: local linear slope vs BP-RP from an iteratively
clipped fit (|res|<0.25 core), zero point = KDE mode of colour-corrected values within |dBPRP|<0.04. Colours (G-J, G-H,
G-Ks, G-W1, Ks-W1, J-Ks, H-Ks): local linear fit + median zero point (WISE quantities only AllWISE cc '00..', ph 'AA..', ext 0).

Twin logic: an equal twin at fixed combined BP-RP has the colour residuals of its components; its components sit at
x_comp = x_obs + 0.753 (q=1). So p(y | twin, x_obs) is the empirical y-distribution of stars at x ~ x_obs + 0.753, and
the twin fraction among over-luminous stars can be estimated by a 2-component mixture of the y-values of stars at
x ~ x_obs: twin component = (empirical, fixed) normal-star y distribution; single component = free Gaussian.
Measurement-error convolution: intrinsic widths are deconvolved with the median field errors and re-convolved with the
target's errors when evaluating P(twin | target).
Binary control: over-luminous stars with RUWE 1.4-8 (astrometric binaries) and Gaia NSS stars: their y should lie on
the 'twin' (binary) line if the discriminant works.
Output: out/r01_population.json, figs/r01_population.png"""
import json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from scipy.stats import gaussian_kde, norm
from scipy.optimize import minimize
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
H = 'tic214311908/analysis/host'
B = 'tic214311908/analysis/referee2_host'
F = pd.read_csv(f'{H}/data/field_mdwarfs.csv')
x = F.bp_rp; fx = 1.162004 + 0.011464 * x + 0.049255 * x ** 2 - 0.005879 * x ** 3
F['Cs'] = (F.phot_bp_rp_excess_factor - fx) / (0.0059898 + 8.817481e-12 * F.phot_g_mean_mag ** 7.618399)
dm = 5 * np.log10(F.parallax / 100)
F['MK'] = F.ks_m + dm; F['MG'] = F.phot_g_mean_mag + dm
F['GK'] = F.phot_g_mean_mag - F.ks_m; F['GJ'] = F.phot_g_mean_mag - F.j_m; F['GH'] = F.phot_g_mean_mag - F.h_m
F['GW1'] = F.phot_g_mean_mag - F.w1mpro; F['KW1'] = F.ks_m - F.w1mpro; F['JK'] = F.j_m - F.ks_m; F['HK'] = F.h_m - F.ks_m
F['eGK'] = np.hypot(F.ks_msigcom, 0.003); F['eGW1'] = np.hypot(F.w1mpro_error, 0.003)
F['eMK'] = np.hypot(F.ks_msigcom, 5 / np.log(10) * F.parallax_error / F.parallax)
C0 = 2.756149
# target (Gaia DR3 raw values; both raw and ZP-corrected parallax)
T = dict(bp_rp=C0, G=14.281043, J=11.363, H=10.740, K=10.482, W1=10.332 - 0.014, plx_raw=15.4973, plx_zp=15.549)
eT = dict(K=0.026, J=0.028, H=0.031, W1=0.03)   # W1: unWISE error + AllWISE-scale offset error (analyst: +-0.013) -> 0.03 adopted
clean = (F.ruwe < 1.1) & (F.ipd_frac_multi_peak == 0) & (F.ipd_gof_harmonic_amplitude < 0.03) & (F.parallax / F.parallax_error > 20) & \
        (F.ph_qual == 'AAA') & (F.ks_msigcom < 0.05) & (F.Cs.abs() < 3)
wok = F.wise_cc.astype(str).str[:2].eq('00') & F.wise_ph_qual.astype(str).str[:2].eq('AA') & (F.wise_ext == 0)
win = (F.bp_rp - C0).abs() < 0.12
out = dict()
# ---------------- local ridges
def mode(v):
    v = v[np.isfinite(v)]; kde = gaussian_kde(v, bw_method=0.2); g = np.linspace(np.percentile(v, 2), np.percentile(v, 98), 2000); return g[np.argmax(kde(g))]
R = {}
for q in ('MK', 'MG'):
    s = clean & win; xx = F.bp_rp[s].values - C0; yy = F[q][s].values
    k = np.ones(len(xx), bool); c = np.polyfit(xx, yy, 1)
    for _ in range(10):
        zp = mode(yy - c[0] * xx); r = yy - c[0] * xx - zp; k = np.abs(r) < 0.25; c = np.polyfit(xx[k], yy[k], 1)
    zp = mode((yy - c[0] * xx)[np.abs(xx) < 0.04]); R[q] = (c[0], zp)
for q in ('GK', 'GJ', 'GH', 'GW1', 'KW1', 'JK', 'HK'):
    s = clean & win & (wok if 'W1' in q else True) & np.isfinite(F[q]); xx = F.bp_rp[s].values - C0; yy = F[q][s].values
    k = np.ones(len(xx), bool)
    for _ in range(10):
        c = np.polyfit(xx[k], yy[k], 1); r = yy - np.polyval(c, xx); sd = 1.4826 * np.median(np.abs(r[k] - np.median(r[k]))); k = np.abs(r) < 4 * sd
    zp = np.median((yy - c[0] * xx)[np.abs(xx) < 0.04]); R[q] = (c[0], zp)
out['ridge_local'] = {k: dict(slope=float(v[0]), zero_at_target_colour=float(v[1])) for k, v in R.items()}
def resid(df, q, col=None):
    return df[q] - (R[q][0] * (df.bp_rp - C0) + R[q][1])
for q in R: F['d' + q] = resid(F, q)
Tt = {}
for lab, plx in (('raw_plx', T['plx_raw']), ('zp_plx', T['plx_zp'])):
    d = 5 * np.log10(plx / 100)
    Tt[lab] = dict(dMK=T['K'] + d - R['MK'][1], dMG=T['G'] + d - R['MG'][1])
Tc = dict(dGK=T['G'] - T['K'] - R['GK'][1], dGJ=T['G'] - T['J'] - R['GJ'][1], dGH=T['G'] - T['H'] - R['GH'][1], dGW1=T['G'] - T['W1'] - R['GW1'][1],
          dKW1=T['K'] - T['W1'] - R['KW1'][1], dJK=T['J'] - T['K'] - R['JK'][1], dHK=T['H'] - T['K'] - R['HK'][1])
out['target'] = dict(cmd=Tt, colour_residuals=Tc)
print('target', Tt, {k: round(v, 3) for k, v in Tc.items()})
tx = Tt['raw_plx']['dMK']; txg = Tt['raw_plx']['dMG']
# fraction of clean stars at the target colour at least as over-luminous
s = clean & ((F.bp_rp - C0).abs() < 0.08)
out['frac_as_bright'] = dict(dMK=float(np.mean(F.dMK[s] <= tx)), dMG=float(np.mean(F.dMG[s] <= txg)), n=int(s.sum()))
print('fraction as bright', out['frac_as_bright'])
# ---------------- empirical y distributions vs x (clean sample)
s = clean & ((F.bp_rp - C0).abs() < 0.08)
D = F[s].copy()
bins = np.arange(-1.2, 0.61, 0.1)
prof = []
for lo in bins[:-1]:
    k = D.dMK.between(lo, lo + 0.1)
    if k.sum() >= 8:
        prof.append(dict(x=float(lo + 0.05), n=int(k.sum()), med_rGK=float(np.median(D.dGK[k])), sd_rGK=float(1.4826 * np.median(np.abs(D.dGK[k] - np.median(D.dGK[k]))))))
out['profile_rGK_vs_dMK_clean'] = prof
# ---------------- conditional mixture
def cond_mixture(D, xcol, ycol, ecol, tx, ty, tey, xshift=0.753, hw=0.1, nboot=300, label='', seed=1):
    """twin template: y of stars with x in [tx+xshift-hw, tx+xshift+hw]; target window: x in [tx-hw, tx+hw]"""
    Dn = D[(D[xcol] - (tx + xshift)).abs() < hw]; Do = D[(D[xcol] - tx).abs() < hw]
    yn = Dn[ycol].values; en = np.median(Dn[ecol]); yo = Do[ycol].values; eo = Do[ecol].values
    # twin template: robust Gaussian (median, MAD), intrinsic width deconvolved
    def tmpl(y):
        mu = np.median(y); sd = 1.4826 * np.median(np.abs(y - mu)); return mu, sd
    def fit(yo, eo, muT, sdT_obs):
        # twin component for field stars: N(muT, sdT_obs) ; single: N(muS, sqrt(sS^2 + eo^2)) ; weights fT
        def nll(p):
            fT = 1 / (1 + np.exp(-p[0])); muS = p[1]; sS = np.exp(p[2])
            L = fT * norm.pdf(yo, muT, sdT_obs) + (1 - fT) * norm.pdf(yo, muS, np.sqrt(sS ** 2 + eo ** 2))
            return -np.sum(np.log(L + 1e-300))
        best = None
        for f0 in (-3, -1, 0.5):
            r = minimize(nll, [f0, np.median(yo), np.log(0.03)], method='Nelder-Mead', options=dict(maxiter=4000, xatol=1e-6, fatol=1e-8))
            if best is None or r.fun < best.fun: best = r
        p = best.x; return 1 / (1 + np.exp(-p[0])), p[1], np.exp(p[2]), best.fun
    muT, sdTo = tmpl(yn); sdTi = np.sqrt(max(sdTo ** 2 - en ** 2, 1e-6))
    fT, muS, sS, nll0 = fit(yo, eo, muT, sdTo)
    def post(fT, muT, sdTi, muS, sS):
        a = fT * norm.pdf(ty, muT, np.hypot(sdTi, tey)); b = (1 - fT) * norm.pdf(ty, muS, np.hypot(sS, tey)); return a / (a + b), a / max(b, 1e-300) * (1 - fT) / max(fT, 1e-12)
    P, LR = post(fT, muT, sdTi, muS, sS)
    # likelihood-ratio profile for fT=0 (is a twin component required?)
    rng = np.random.default_rng(seed); bs = []
    for _ in range(nboot):
        yn_b = rng.choice(yn, len(yn)); ii = rng.integers(0, len(yo), len(yo))
        m_, s_ = tmpl(yn_b); si = np.sqrt(max(s_ ** 2 - en ** 2, 1e-6))
        f_, mS_, sS_, _ = fit(yo[ii], eo[ii], m_, s_)
        bs.append((f_, post(f_, m_, si, mS_, sS_)[0], post(f_, m_, si, mS_, sS_)[1]))
    bs = np.array(bs)
    res = dict(label=label, n_template=int(len(yn)), n_window=int(len(yo)), twin_template_mu=float(muT), twin_template_sd_obs=float(sdTo), twin_template_sd_intr=float(sdTi),
               median_err_template=float(en), median_err_window=float(np.median(eo)), window_median_y=float(np.median(yo)),
               frac_window_y_ge_target=float(np.mean(yo >= ty)), frac_template_y_ge_target=float(np.mean(yn >= ty)),
               fT_twin_fraction_in_window=float(fT), fT_16_50_84=[float(v) for v in np.percentile(bs[:, 0], [16, 50, 84])],
               single_mu=float(muS), single_sd_intr=float(sS), target_y=float(ty), target_err=float(tey),
               LR_twin_over_single_y_only=float(LR), LR_16_84=[float(v) for v in np.percentile(bs[:, 2], [16, 84])],
               P_twin=float(P), P_twin_16_50_84=[float(v) for v in np.percentile(bs[:, 1], [16, 50, 84])],
               P_twin_95pct_upper=float(np.percentile(bs[:, 1], 95)))
    print(f"[{label}] n_tmpl={len(yn)} n_win={len(yo)} tmpl mu={muT:+.3f} sd_obs={sdTo:.3f} sd_int={sdTi:.3f} | window median y={np.median(yo):+.3f} "
          f"frac>=target {res['frac_window_y_ge_target']:.2f} (template {res['frac_template_y_ge_target']:.3f}) | fT={fT:.3f} {res['fT_16_50_84']} single mu={muS:+.3f} sd={sS:.3f} | "
          f"LR={LR:.3f} P(twin)={P:.4f} {res['P_twin_16_50_84']} 95%<{res['P_twin_95pct_upper']:.3f}", flush=True)
    return res
runs = {}
eTK = np.hypot(eT['K'], 0.003); eTW = np.hypot(eT['W1'], 0.003)
ty = Tc['dGK']; tyw = Tc['dGW1']
Dk = D[np.isfinite(D.dGK)]; Dw = D[wok[s].values & np.isfinite(D.dGW1)]
runs['K_clean'] = cond_mixture(Dk, 'dMK', 'dGK', 'eGK', tx, ty, eTK, label='dMK-r(G-Ks) clean')
runs['K_clean_hw0.15'] = cond_mixture(Dk, 'dMK', 'dGK', 'eGK', tx, ty, eTK, hw=0.15, label='dMK-r(G-Ks) clean, window +-0.15')
runs['K_clean_Ks_fainter_0.04'] = cond_mixture(Dk, 'dMK', 'dGK', 'eGK', tx + 0.04, ty - 0.04, eTK, label='Ks 0.04 fainter')
runs['W_clean'] = cond_mixture(Dw, 'dMG', 'dGW1', 'eGW1', txg, tyw, eTW, label='dMG-r(G-W1) clean')
runs['W_clean_hw0.15'] = cond_mixture(Dw, 'dMG', 'dGW1', 'eGW1', txg, tyw, eTW, hw=0.15, label='dMG-r(G-W1) clean, window +-0.15')
runs['W_clean_W1_fainter_0.03'] = cond_mixture(Dw, 'dMG', 'dGW1', 'eGW1', txg, tyw - 0.03, eTW, label='W1 0.03 fainter')
# distance-matched subsample (40-100 pc) and kinematic subsample
vt = 4.74047 * np.hypot(D.pmra, D.pmdec) / D.parallax
runs['K_plx10_25'] = cond_mixture(Dk[Dk.parallax.between(10, 25)], 'dMK', 'dGK', 'eGK', tx, ty, eTK, label='dMK-r(G-Ks) plx 10-25')
runs['K_vtan_gt40'] = cond_mixture(Dk[(4.74047 * np.hypot(Dk.pmra, Dk.pmdec) / Dk.parallax) > 40], 'dMK', 'dGK', 'eGK', tx, ty, eTK, label='dMK-r(G-Ks) vtan>40')
runs['W_vtan_gt40'] = cond_mixture(Dw[(4.74047 * np.hypot(Dw.pmra, Dw.pmdec) / Dw.parallax) > 40], 'dMG', 'dGW1', 'eGW1', txg, tyw, eTW, label='dMG-r(G-W1) vtan>40')
out['conditional_mixture'] = runs
# ---------------- binary control: astrometric binaries (RUWE 1.4-8, IPD-clean) and NSS, over-luminous
sb = (F.ruwe.between(1.4, 8)) & (F.ipd_frac_multi_peak <= 2) & (F.parallax / F.parallax_error > 20) & (F.ph_qual == 'AAA') & (F.ks_msigcom < 0.05) & \
     (F.Cs.abs() < 3) & ((F.bp_rp - C0).abs() < 0.08)
ctrl = {}
for lab, sel in (('RUWE1.4-8', sb), ('NSS>0', (F.non_single_star > 0) & (F.parallax / F.parallax_error > 20) & (F.ph_qual == 'AAA') & (F.Cs.abs() < 3) & ((F.bp_rp - C0).abs() < 0.08))):
    for lo, hi in ((-0.95, -0.55), (-0.55, -0.25), (-0.15, 0.15)):
        k = sel & F.dMK.between(lo, hi) & np.isfinite(F.dGK)
        kc = clean & ((F.bp_rp - C0).abs() < 0.08) & F.dMK.between(lo, hi)
        ctrl[f'{lab} dMK {lo}..{hi}'] = dict(n=int(k.sum()), median_rGK=float(np.median(F.dGK[k])) if k.sum() else None,
                                              sd=float(1.4826 * np.median(np.abs(F.dGK[k] - np.median(F.dGK[k])))) if k.sum() else None,
                                              clean_same_dMK_median_rGK=float(np.median(F.dGK[kc])), clean_n=int(kc.sum()),
                                              median_rGW1=float(np.median(F.dGW1[k & wok])) if (k & wok).sum() else None,
                                              clean_same_dMG_median_rGW1=float(np.median(F.dGW1[clean & ((F.bp_rp - C0).abs() < 0.08) & F.dMG.between(lo, hi) & wok])))
        print('control', lab, lo, hi, ctrl[f'{lab} dMK {lo}..{hi}'])
out['binary_control'] = ctrl
json.dump(out, open(f'{B}/out/r01_population.json', 'w'), indent=1)
# ---------------- figure
fig, ax = plt.subplots(1, 3, figsize=(16, 5))
ax[0].scatter(D.dMK, D.dGK, s=2, c='0.6', label='clean field, |dBP-RP|<0.08')
k = sb & np.isfinite(F.dGK); ax[0].scatter(F.dMK[k], F.dGK[k], s=4, c='tab:orange', label='RUWE 1.4-8')
ax[0].plot(tx, ty, 'r*', ms=16, label='M3V'); ax[0].set_xlim(0.8, -1.4); ax[0].set_ylim(-0.2, 0.3); ax[0].set_xlabel('dMKs'); ax[0].set_ylabel('r(G-Ks)'); ax[0].legend(fontsize=7)
r = runs['K_clean']; yy = np.linspace(-0.15, 0.25, 300)
Do = Dk[(Dk.dMK - tx).abs() < 0.1]; Dn = Dk[(Dk.dMK - (tx + 0.753)).abs() < 0.1]
ax[1].hist(Dn.dGK, bins=40, range=(-0.15, 0.25), density=True, alpha=0.5, label=f'stars at dMKs={tx + 0.753:+.2f} (twin template)')
ax[1].hist(Do.dGK, bins=20, range=(-0.15, 0.25), density=True, alpha=0.5, label=f'stars at dMKs={tx:+.2f} (n={len(Do)})')
ax[1].axvline(ty, c='r'); ax[1].set_xlabel('r(G-Ks)'); ax[1].legend(fontsize=7)
Do = Dw[(Dw.dMG - txg).abs() < 0.1]; Dn = Dw[(Dw.dMG - (txg + 0.753)).abs() < 0.1]
ax[2].hist(Dn.dGW1, bins=40, range=(-0.15, 0.25), density=True, alpha=0.5, label='twin template')
ax[2].hist(Do.dGW1, bins=20, range=(-0.15, 0.25), density=True, alpha=0.5, label=f'stars at dMG={txg:+.2f} (n={len(Do)})')
ax[2].axvline(tyw, c='r'); ax[2].set_xlabel('r(G-W1)'); ax[2].legend(fontsize=7)
plt.tight_layout(); plt.savefig(f'{B}/figs/r01_population.png', dpi=110)
