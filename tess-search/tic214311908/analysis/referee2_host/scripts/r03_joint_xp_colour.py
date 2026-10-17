"""r03: joint (r(G-Ks), XP-[Fe/H]) test, so that the two metallicity indicators are not multiplied as if independent.
Inputs: analyst's XP predictions for field stars (host/out/h12b_xp_pred_field.csv: same PCA+ridge model applied to
normal and over-luminous field M dwarfs, RUWE<1.2, G<15) joined with r(G-Ks), r(G-W1), dMKs from r01's local ridges.
(1) correlation of r(G-Ks) and XP-[Fe/H] among normal and over-luminous stars (are they independent given H?);
(2) 2-component 2D mixture among over-luminous stars (-1.0 < dMK < -0.5): twin component = 2D Gaussian of normal
    stars (|dMK| < 0.15, shifted to the twin's component position), single = 2D Gaussian with free mean and covariance
    fixed to k x the normal-star covariance (k = 1, 1.5); target's errors added (r: 0.026; XP: 0.077 MC noise);
(3) P(twin | target) jointly vs the product of the 1D likelihood ratios.
Output: out/r03_joint.json"""
import json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from scipy.stats import multivariate_normal as mvn
from scipy.optimize import minimize
H = 'tic214311908/analysis/host'
B = 'tic214311908/analysis/referee2_host'
r01 = json.load(open(f'{B}/out/r01_population.json'))
R = {k: (v['slope'], v['zero_at_target_colour']) for k, v in r01['ridge_local'].items()}
C0 = 2.756149
F = pd.read_csv(f'{H}/data/field_mdwarfs.csv').set_index('source_id')
P = pd.read_csv(f'{H}/out/h12b_xp_pred_field.csv').set_index('source_id')
P = P.join(F[['ks_m', 'phot_g_mean_mag', 'parallax', 'w1mpro', 'ph_qual', 'ks_msigcom', 'wise_cc', 'wise_ph_qual', 'wise_ext', 'ipd_frac_multi_peak']], how='left')
P['dMK'] = P.ks_m + 5 * np.log10(P.parallax / 100) - (R['MK'][0] * (P.bp_rp - C0) + R['MK'][1])
P['rGK'] = P.phot_g_mean_mag - P.ks_m - (R['GK'][0] * (P.bp_rp - C0) + R['GK'][1])
P = P[(P.ph_qual == 'AAA') & (P.ks_msigcom < 0.05) & np.isfinite(P.rGK) & np.isfinite(P.feh_xp)]
out = {}
tr, tx = r01['target']['colour_residuals']['dGK'], 0.12922249223373056
er, ex = 0.026, 0.0772
nor = P[P.dMK.abs() < 0.15]; ovl = P[P.dMK.between(-1.0, -0.5)]; win = P[P.dMK.between(-0.9, -0.6)]
for lab, d in (('normal |dMK|<0.15', nor), ('overlum -1.0..-0.5', ovl), ('window -0.9..-0.6', win)):
    out[lab] = dict(n=len(d), corr=float(np.corrcoef(d.rGK, d.feh_xp)[0, 1]), med_rGK=float(np.median(d.rGK)), med_xp=float(np.median(d.feh_xp)),
                    sd_rGK=float(1.4826 * np.median(np.abs(d.rGK - np.median(d.rGK)))), sd_xp=float(1.4826 * np.median(np.abs(d.feh_xp - np.median(d.feh_xp)))))
    print(lab, out[lab])
def robust_cov(X):
    med = np.median(X, axis=0); sd = 1.4826 * np.median(np.abs(X - med), axis=0)
    k = (np.abs((X - med) / sd) < 3).all(axis=1); C = np.cov(X[k].T); return X[k].mean(axis=0), C
muT, CT = robust_cov(nor[['rGK', 'feh_xp']].values)
Y = win[['rGK', 'feh_xp']].values
res = {}
for kk in (1.0, 1.5, 2.0):
    CS = CT * kk
    def nll(p):
        f = 1 / (1 + np.exp(-p[0])); m = p[1:3]
        L = f * mvn.pdf(Y, muT, CT) + (1 - f) * mvn.pdf(Y, m, CS); return -np.sum(np.log(L + 1e-300))
    best = min((minimize(nll, [f0, 0.07, 0.12], method='Nelder-Mead') for f0 in (-3, -1, 0.5)), key=lambda r: r.fun)
    f = 1 / (1 + np.exp(-best.x[0])); mS = best.x[1:3]
    E = np.diag([er ** 2, ex ** 2])
    # field-star measurement noise is already in CT (normal stars have similar G<15 S/N) -> add the target's extra error only approx: use CT + E/2
    a = mvn.pdf([tr, tx], muT, CT + E / 2); b = mvn.pdf([tr, tx], mS, CS + E / 2)
    a_r = mvn.pdf(tr, muT[0], CT[0, 0] + er ** 2 / 2); b_r = mvn.pdf(tr, mS[0], CS[0, 0] + er ** 2 / 2)
    a_x = mvn.pdf(tx, muT[1], CT[1, 1] + ex ** 2 / 2); b_x = mvn.pdf(tx, mS[1], CS[1, 1] + ex ** 2 / 2)
    LRj = a / b; LRp = (a_r / b_r) * (a_x / b_x)
    # bootstrap
    rng = np.random.default_rng(1); bsP = []
    for _ in range(200):
        Yb = Y[rng.integers(0, len(Y), len(Y))]; Nb = nor[['rGK', 'feh_xp']].values[rng.integers(0, len(nor), len(nor))]
        mT_, CT_ = robust_cov(Nb)
        def nl(p):
            ff = 1 / (1 + np.exp(-p[0])); return -np.sum(np.log(ff * mvn.pdf(Yb, mT_, CT_) + (1 - ff) * mvn.pdf(Yb, p[1:3], CT_ * kk) + 1e-300))
        rb = minimize(nl, best.x, method='Nelder-Mead'); ff = 1 / (1 + np.exp(-rb.x[0]))
        aa = mvn.pdf([tr, tx], mT_, CT_ + E / 2); bb = mvn.pdf([tr, tx], rb.x[1:3], CT_ * kk + E / 2); bsP.append(ff * aa / (ff * aa + (1 - ff) * bb))
    res[f'k{kk}'] = dict(n_window=len(Y), fT=float(f), single_mean=mS.tolist(), twin_mean=muT.tolist(), twin_cov=CT.tolist(), LR_joint=float(LRj), LR_product_of_marginals=float(LRp),
                         LR_r_only=float(a_r / b_r), LR_xp_only=float(a_x / b_x), P_twin_joint=float(f * LRj / (f * LRj + 1 - f)),
                         P_twin_16_50_84=[float(v) for v in np.percentile(bsP, [16, 50, 84])], P_twin_95=float(np.percentile(bsP, 95)))
    print(kk, {k: (np.round(v, 4) if isinstance(v, float) else v) for k, v in res[f'k{kk}'].items() if k != 'twin_cov'})
out['mixture'] = res
json.dump(out, open(f'{B}/out/r03_joint.json', 'w'), indent=1)
