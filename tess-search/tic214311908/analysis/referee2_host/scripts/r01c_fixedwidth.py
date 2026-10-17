"""r01c: identifiability check of the conditional (window) mixture of r01. The free width of the 'metal-rich single'
component trades off against the twin fraction fT. Here the single component's INTRINSIC width is fixed to
k x the twin template's intrinsic width (k = 1, 1.25, 1.5, 2) - the same device the analyst used in h12c - and fT,
the target's y-only likelihood ratio and P(twin | x, y) are re-derived in both planes, with bootstrap.
Twin template = stars at x_obs + shift (q=1: 0.753; also 'binary-like' q~0.8-1 template with shift 0.65-0.75).
Output: out/r01c_fixedwidth.json"""
import json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from scipy.stats import norm
from scipy.optimize import minimize_scalar
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
F['eGK'] = np.hypot(F.ks_msigcom, 0.003); F['eGW1'] = np.hypot(F.w1mpro_error, 0.003)
clean = (F.ruwe < 1.1) & (F.ipd_frac_multi_peak == 0) & (F.ipd_gof_harmonic_amplitude < 0.03) & (F.parallax / F.parallax_error > 20) & \
        (F.ph_qual == 'AAA') & (F.ks_msigcom < 0.05) & (F.Cs.abs() < 3) & ((F.bp_rp - C0).abs() < 0.08)
wok = F.wise_cc.astype(str).str[:2].eq('00') & F.wise_ph_qual.astype(str).str[:2].eq('AA') & (F.wise_ext == 0)
Dk = F[clean & np.isfinite(F.dGK)]; Dw = F[clean & wok & np.isfinite(F.dGW1)]
T = r01['target']; tx = T['cmd']['raw_plx']['dMK']; txg = T['cmd']['raw_plx']['dMG']; ty = T['colour_residuals']['dGK']; tyw = T['colour_residuals']['dGW1']
eTK, eTW = np.hypot(0.026, 0.003), np.hypot(0.03, 0.003)
def run(D, xc, yc, ec, tx, ty, tey, kw, shifts=(0.753,), hw=0.1, nboot=400, seed=2):
    rng = np.random.default_rng(seed)
    yn = np.concatenate([D[(D[xc] - (tx + s)).abs() < hw][yc].values for s in shifts]); en = np.median(D[(D[xc] - (tx + shifts[0])).abs() < hw][ec])
    Do = D[(D[xc] - tx).abs() < hw]; yo = Do[yc].values; eo = Do[ec].values
    def one(yn, yo, eo):
        muT = np.median(yn); sT = 1.4826 * np.median(np.abs(yn - muT)); sTi = np.sqrt(max(sT ** 2 - en ** 2, 1e-6)); sSi = kw * sTi
        def nll_f(fT, muS):
            L = fT * norm.pdf(yo, muT, sT) + (1 - fT) * norm.pdf(yo, muS, np.sqrt(sSi ** 2 + eo ** 2)); return -np.sum(np.log(L + 1e-300))
        best = (np.inf, None, None)
        for fT in np.linspace(0, 0.9, 91):
            r = minimize_scalar(lambda m: nll_f(fT, m), bounds=(-0.05, 0.2), method='bounded')
            if r.fun < best[0]: best = (r.fun, fT, r.x)
        _, fT, muS = best
        a = norm.pdf(ty, muT, np.hypot(sTi, tey)); b = norm.pdf(ty, muS, np.hypot(sSi, tey))
        LR = a / b; P = fT * a / (fT * a + (1 - fT) * b)
        return fT, LR, P, muS
    fT, LR, P, muS = one(yn, yo, eo)
    bs = np.array([one(rng.choice(yn, len(yn)), *(lambda ii: (yo[ii], eo[ii]))(rng.integers(0, len(yo), len(yo))))[:3] for _ in range(nboot)])
    return dict(k_width=kw, n_window=int(len(yo)), fT=float(fT), fT_16_84=[float(np.percentile(bs[:, 0], 16)), float(np.percentile(bs[:, 0], 84))], single_mu=float(muS),
                LR_y=float(LR), P_twin=float(P), P_twin_16_50_84=[float(v) for v in np.percentile(bs[:, 2], [16, 50, 84])], P_twin_95=float(np.percentile(bs[:, 2], 95)))
out = {}
for lab, D, xc, yc, ec, x0, y0, e0 in (('K', Dk, 'dMK', 'dGK', 'eGK', tx, ty, eTK), ('W', Dw, 'dMG', 'dGW1', 'eGW1', txg, tyw, eTW),
                                        ('K_Ks+0.04', Dk, 'dMK', 'dGK', 'eGK', tx + 0.04, ty - 0.04, eTK), ('W_W1+0.03', Dw, 'dMG', 'dGW1', 'eGW1', txg, tyw - 0.03, eTW)):
    for kw in (1.0, 1.25, 1.5, 2.0):
        for sh_lab, sh in (('q1', (0.753,)), ('q0.8-1', (0.65, 0.70, 0.753))):
            r = run(D, xc, yc, ec, x0, y0, e0, kw, sh, nboot=150)
            out[f'{lab}_k{kw}_{sh_lab}'] = r
            print(f'{lab:10s} k={kw:4.2f} {sh_lab:7s} n={r["n_window"]} fT={r["fT"]:.2f} [{r["fT_16_84"][0]:.2f},{r["fT_16_84"][1]:.2f}] muS={r["single_mu"]:+.3f} LR_y={r["LR_y"]:.3f} '
                  f'P(twin)={r["P_twin"]:.4f} [{r["P_twin_16_50_84"][0]:.4f},{r["P_twin_16_50_84"][2]:.4f}] 95%<{r["P_twin_95"]:.4f}', flush=True)
json.dump(out, open(f'{B}/out/r01c_fixedwidth.json', 'w'), indent=1)
