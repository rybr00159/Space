"""h13: luminosity-independent (colour-only) metallicity of the M3V, and the [Fe/H]-dependence of the CMD offset.
Calibrators: data/calibrators_gaia.csv (h04); main scale = Mann+13-type ([Fe/H] from Mann+15, Mann+13 wide binaries,
Newton+14 (+cpm), Rojas-Ayala+12, Terrien+15, Gaidos+14, Birky+20 APOGEE-Cannon); LAMOST Ding+22 is kept separate.
Duplicates: inverse-variance mean per Gaia source. Cuts: -1.3 < dMKs < 1.3 (removes giants, pre-MS stars and
cross-match failures that contaminate Birky+20),  2.40 < BP-RP < 3.15, parallax_over_error > 20, 2MASS AAA,
ks_msigcom < 0.05, |C*| < 3 sigma (Riello+21), ipd_frac_multi_peak <= 2, RUWE < 1.2 (and a RUWE < 1.1 subset).
Ridges at fixed BP-RP from the field sample (h03; RUWE < 1.2, AAA): KDE mode for M_Ks and M_G (the single-star
sequence) and running medians for the colours G-Ks, G-J, G-H, J-H, H-Ks, J-Ks, G-RP, BP-G, Ks-W1, W1-W2 (AllWISE).
For each quantity: residual = a + b [Fe/H] fitted to calibrators by iterative 3-sigma clipping (MAD scatter; the
clipped fraction removes unresolved binaries/bad photometry); slope, scatter (includes [Fe/H] label noise).
Colour-only [Fe/H] estimator: ridge regression of [Fe/H] on the colour residuals (no absolute magnitude), 10-fold CV
-> calibrated scatter; applied to the M3V (2MASS JHKs; W1/W2 from unWISE put on the AllWISE scale, h06b).
Twin test: an equal-mass twin has the colours of a single star of the same type, so the colour-only [Fe/H] measures
the composition; compare with the [Fe/H] needed to explain dMKs = -0.75 if single.
Output: out/h13_colour_feh.json, figs/h13_colour_feh.png, data/calib_sel.csv"""
import json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from scipy.stats import gaussian_kde
from scipy.optimize import minimize
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
B = 'tic214311908/analysis/host'
def Cstar(df):
    x = df.bp_rp; f = 1.162004 + 0.011464 * x + 0.049255 * x ** 2 - 0.005879 * x ** 3
    return (df.phot_bp_rp_excess_factor - f) / (0.0059898 + 8.817481e-12 * df.phot_g_mean_mag ** 7.618399)
def add_cols(d):
    d['MK'] = d.ks_m + 5 * np.log10(d.parallax / 100); d['MG'] = d.phot_g_mean_mag + 5 * np.log10(d.parallax / 100)
    d['GK'] = d.phot_g_mean_mag - d.ks_m; d['GJ'] = d.phot_g_mean_mag - d.j_m; d['GH'] = d.phot_g_mean_mag - d.h_m
    d['JH'] = d.j_m - d.h_m; d['HK'] = d.h_m - d.ks_m; d['JK'] = d.j_m - d.ks_m; d['GRP'] = d.phot_g_mean_mag - d.phot_rp_mean_mag
    d['BPG'] = d.phot_bp_mean_mag - d.phot_g_mean_mag; d['KW1'] = d.ks_m - d.w1mpro; d['W1W2'] = d.w1mpro - d.w2mpro; d['GW1'] = d.phot_g_mean_mag - d.w1mpro
    d['Cs'] = Cstar(d); return d
F = add_cols(pd.read_csv(f'{B}/data/field_mdwarfs.csv'))
fq = (F.ruwe < 1.2) & (F.ph_qual == 'AAA') & (F.ks_msigcom < 0.05) & (F.Cs.abs() < 3) & (F.ipd_frac_multi_peak <= 2)
Fs = F[fq & F.bp_rp.between(2.3, 3.25)]
wq = Fs.wise_ph_qual.astype(str).str[:2].eq('AA') & Fs.wise_cc.astype(str).str[:2].eq('00') & (Fs.wise_ext == 0)
x0s = np.arange(2.36, 3.20, 0.02)
def mode_ridge(x, y):
    pts = []
    for x0 in x0s:
        k = np.abs(x - x0) < 0.025
        if k.sum() < 60: continue
        yy = y[k]; kde = gaussian_kde(yy, bw_method=0.25); gr = np.linspace(np.percentile(yy, 1), np.percentile(yy, 99), 500)
        pts.append((x0, gr[np.argmax(kde(gr))]))
    pts = np.array(pts); return np.polyfit(pts[:, 0], pts[:, 1], 4)
def med_ridge(x, y):
    pts = []
    for x0 in x0s:
        k = np.abs(x - x0) < 0.025
        if k.sum() < 40: continue
        pts.append((x0, np.median(y[k])))
    pts = np.array(pts); return np.polyfit(pts[:, 0], pts[:, 1], 4)
ridge = {'MK': mode_ridge(Fs.bp_rp.values, Fs.MK.values), 'MG': mode_ridge(Fs.bp_rp.values, Fs.MG.values)}
COLS = ['GK', 'GJ', 'GH', 'JH', 'HK', 'JK', 'GRP', 'BPG', 'KW1', 'W1W2', 'GW1']
for c in COLS:
    sel = wq if c in ('KW1', 'W1W2', 'GW1') else np.ones(len(Fs), bool)
    ridge[c] = med_ridge(Fs.bp_rp.values[sel], Fs[c].values[sel])
# ---- calibrators
C = pd.read_csv(f'{B}/data/calibrators_gaia.csv')
C = C[np.isfinite(C.feh) & np.isfinite(C.bp_rp)]
C['w'] = 1 / np.clip(C.efeh.fillna(0.12), 0.05, 0.5) ** 2
C['is_lamost'] = C.scale.eq('LAMOST')
def agg(g):
    w = g.w.values; return pd.Series(dict(feh=np.sum(w * g.feh) / np.sum(w), efeh=1 / np.sqrt(np.sum(w)), src='+'.join(sorted(set(g.src))), n=len(g)))
M13 = C[~C.is_lamost].groupby('source_id').apply(agg).reset_index()
LAM = C[C.is_lamost].groupby('source_id').apply(agg).reset_index()
gcols = [c for c in C.columns if c not in ('src', 'name', 'ra', 'dec', 'feh', 'efeh', 'scale', 'row', 'd00', 'w', 'is_lamost')]
gaia = C.drop_duplicates('source_id')[gcols]
M13 = add_cols(M13.merge(gaia, on='source_id')); LAM = add_cols(LAM.merge(gaia, on='source_id'))
def qual(d, ruwe=1.2):
    return d.bp_rp.between(2.40, 3.15) & (d.parallax / d.parallax_error > 20) & (d.ph_qual == 'AAA') & (d.ks_msigcom < 0.05) & \
        (d.Cs.abs() < 3) & (d.ipd_frac_multi_peak <= 2) & (d.ruwe < ruwe)
for d in (M13, LAM):
    d['dMK'] = d.MK - np.polyval(ridge['MK'], d.bp_rp); d['dMG'] = d.MG - np.polyval(ridge['MG'], d.bp_rp)
    for c in COLS:
        d['r' + c] = d[c] - np.polyval(ridge[c], d.bp_rp)
    d['wiseok'] = d.wise_ph_qual.astype(str).str[:2].eq('AA') & d.wise_cc.astype(str).str[:2].eq('00') & (d.wise_ext == 0)
cal = M13[qual(M13) & M13.dMK.between(-1.3, 1.3)].copy(); lam = LAM[qual(LAM) & LAM.dMK.between(-1.3, 1.3)].copy()   # CMD sanity cut: removes giants/pre-MS/mismatches (|dMKs| > 1.3)
print('M13-scale calibrators', len(cal), cal.src.str.split('+').explode().value_counts().to_dict(), '; LAMOST', len(lam))
# LAMOST vs M13 offset on common stars
com = M13[['source_id', 'feh']].merge(LAM[['source_id', 'feh']], on='source_id', suffixes=('_m13', '_lam'))
print('common M13/LAMOST:', len(com), 'median diff (LAMOST - M13)', float(np.median(com.feh_lam - com.feh_m13)) if len(com) else None,
      'corr', float(np.corrcoef(com.feh_lam, com.feh_m13)[0, 1]) if len(com) > 3 else None)
cal.to_csv(f'{B}/data/calib_sel.csv', index=False)
# ---- target
plx = 15.549
T = dict(bp_rp=2.756149, phot_g_mean_mag=14.281043, phot_bp_mean_mag=15.830741, phot_rp_mean_mag=13.074592, j_m=11.363, h_m=10.740, ks_m=10.482,
         w1mpro=10.332 - 0.014, w2mpro=10.229 - 0.030, parallax=plx, phot_bp_rp_excess_factor=1.452786)
Td = add_cols(pd.DataFrame([T])); Td['Cs'] = 0.71
Td['dMK'] = Td.MK - np.polyval(ridge['MK'], Td.bp_rp); Td['dMG'] = Td.MG - np.polyval(ridge['MG'], Td.bp_rp)
for c in COLS:
    Td['r' + c] = Td[c] - np.polyval(ridge[c], Td.bp_rp)
terr = dict(rGW1=0.03, rGK=0.027, rGJ=0.029, rGH=0.032, rJH=0.042, rHK=0.040, rJK=0.038, rGRP=0.004, rBPG=0.006, rKW1=0.035, rW1W2=0.035, dMK=0.027, dMG=0.006)
out = dict(target={k: float(Td[k].iloc[0]) for k in ['dMK', 'dMG'] + ['r' + c for c in COLS]})
print('target residuals', {k: round(v, 3) for k, v in out['target'].items()})
# ---- mixture fit residual = a + b feh
def mixfit(x, y, ey):
    """robust linear fit residual = a + b [Fe/H]: iterative 3-sigma clipping with MAD scatter (replaces the earlier
    2-component mixture, which degenerated once the CMD sanity cut removed the giant/mismatch outliers); bootstrap errors."""
    def clipfit(xx, yy):
        k = np.ones(len(xx), bool)
        for _ in range(15):
            c = np.polyfit(xx[k], yy[k], 1); r = yy - np.polyval(c, xx); s_ = 1.4826 * np.median(np.abs(r[k] - np.median(r[k])))
            k_new = np.abs(r) < 3 * s_
            if (k_new == k).all(): break
            k = k_new
        return c, s_, k
    c, s_, k = clipfit(x, y)
    rng = np.random.default_rng(1); bs = []
    for _ in range(200):
        ii = rng.integers(0, len(x), len(x)); cb, _, _ = clipfit(x[ii], y[ii]); bs.append(cb)
    bs = np.array(bs)
    return dict(a=float(c[1]), b=float(c[0]), b_err=float(np.std(bs[:, 0])), a_err=float(np.std(bs[:, 1])), scatter=float(s_),
                f_out=float(1 - k.mean()), n=int(len(x)), method='3-sigma clipped linear fit, MAD scatter')
fits = {}
for q in ['dMK', 'dMG'] + ['r' + c for c in COLS]:
    for lab, sel in (('ruwe1.2', np.ones(len(cal), bool)), ('ruwe1.1', (cal.ruwe < 1.1).values)):
        k = sel & np.isfinite(cal[q].values)
        if q in ('rKW1', 'rW1W2', 'rGW1'): k &= cal.wiseok.values
        x = cal.feh.values[k]; y = cal[q].values[k]; ey = np.full(k.sum(), 0.03)
        fits[f'{q}_{lab}'] = mixfit(x, y, ey)
    f = fits[f'{q}_ruwe1.1']
    print(f'{q:6s}: slope {f["b"]:+.3f} +- {f["b_err"]:.3f} mag/dex, zero {f["a"]:+.3f}, scatter {f["scatter"]:.3f}, f_out {f["f_out"]:.2f} (n={f["n"]}); '
          f'target {out["target"][q]:+.3f} -> [Fe/H] {((out["target"][q] - f["a"]) / f["b"]):+.2f} (if slope != 0)')
out['fits'] = fits
# ---- colour-only regression ([Fe/H] on colour residuals), CV
feats = ['rGK', 'rGJ', 'rGH', 'rJH', 'rHK', 'rGRP', 'rBPG', 'rKW1', 'rW1W2', 'bp_rp']
k = (cal.ruwe < 1.2) & cal.wiseok & np.isfinite(cal[feats]).all(axis=1)
X = cal.loc[k, feats].values; y = cal.feh.values[k]
def ridge_fit(X, y, lam):
    mu = X.mean(0); sd = X.std(0); Z = (X - mu) / sd
    A = Z.T @ Z + lam * np.eye(Z.shape[1]); w = np.linalg.solve(A, Z.T @ (y - y.mean()))
    return lambda Xn: y.mean() + ((Xn - mu) / sd) @ w, w
rng = np.random.default_rng(5)
best = None
for lam in (0.1, 1, 3, 10, 30, 100):
    fold = rng.integers(0, 10, len(y)); pred = np.zeros(len(y))
    for f_ in range(10):
        tr = fold != f_; mdl, _ = ridge_fit(X[tr], y[tr], lam); pred[~tr] = mdl(X[~tr])
    res = y - pred; sd = 1.4826 * np.median(np.abs(res - np.median(res)))
    print(f'colour-only ridge lam={lam}: CV robust scatter {sd:.3f} dex, rms {np.std(res):.3f}, corr {np.corrcoef(pred, y)[0,1]:.2f}')
    if best is None or sd < best[1]: best = (lam, sd, np.std(res), pred.copy())
lam, sdcv, rmscv, predcv = best
mdl, wts = ridge_fit(X, y, lam)
xt = np.array([[out['target'][f] if f != 'bp_rp' else 2.756149 for f in feats]])
fe_col = float(mdl(xt)[0])
# measurement error of target features -> MC
mc = []
for _ in range(2000):
    xx = xt.copy()
    for j, f in enumerate(feats):
        if f in terr: xx[0, j] += rng.normal(0, terr[f])
    mc.append(mdl(xx)[0])
fe_col_err = float(np.hypot(np.std(mc), sdcv))
# bias check: CV prediction vs truth in the metal-rich tail; and vs dMK at fixed feh (binarity leakage)
hi = y > 0.25
out['colour_only'] = dict(features=feats, lam=lam, n=int(len(y)), cv_robust_scatter=sdcv, cv_rms=rmscv, weights=dict(zip(feats, wts.tolist())),
                          target_feh=fe_col, target_feh_err=fe_col_err, target_meas_err=float(np.std(mc)),
                          cv_bias_metalrich=float(np.median(predcv[hi] - y[hi])) if hi.sum() else None,
                          cv_pred_of_feh_gt_0p3=dict(n=int((y > 0.3).sum()), median_pred=float(np.median(predcv[y > 0.3])) if (y > 0.3).sum() else None))
print('colour-only [Fe/H] of the M3V:', fe_col, '+-', fe_col_err, ' (CV scatter', sdcv, ')', out['colour_only']['cv_pred_of_feh_gt_0p3'])
# dependence of colour-only prediction on dMK at fixed feh (do overluminous calibrators get biased predictions?)
dmk = cal.dMK.values[k]
rr = predcv - y
for lo, hi_ in ((-0.3, 0.3), (-1.0, -0.4)):
    s = (dmk > lo) & (dmk < hi_)
    print(f'   CV residual for calibrators with dMK in ({lo},{hi_}): median {np.median(rr[s]):+.3f} (n={s.sum()})')
# ---- the [Fe/H] that would explain dMK if single, from the dMK-[Fe/H] fit (core component)
f = fits['dMK_ruwe1.1']
out['feh_needed_single'] = dict(value=float((out['target']['dMK'] - f['a']) / f['b']),
                                err=float(abs(f['scatter'] / f['b'])), note='core dMK = a + b [Fe/H], scatter = intrinsic single-star scatter')
fg = fits['dMG_ruwe1.1']
out['feh_needed_single_MG'] = dict(value=float((out['target']['dMG'] - fg['a']) / fg['b']), err=float(abs(fg['scatter'] / fg['b'])))
print('feh needed (single) from dMK:', out['feh_needed_single'], ' from dMG:', out['feh_needed_single_MG'])
# the over-luminosity ratio dMG/dMK: twin -> ~1 (colour unchanged); metallicity -> b_MG/b_MK
out['dMG_over_dMK'] = dict(target=float(out['target']['dMG'] / out['target']['dMK']), metallicity_prediction=float(fg['b'] / f['b']), twin_prediction=1.0)
print('dMG/dMK target', out['dMG_over_dMK'])
out['ridges'] = {k: v.tolist() for k, v in ridge.items()}
# forward calibration of the colour-only regression (CV prediction = A0 * true + A1) for use as a likelihood
A = np.polyfit(y, predcv, 1); out['colour_only']['cv_slope'] = float(A[0]); out['colour_only']['cv_intercept'] = float(A[1])
out['colour_only']['cv_scatter_about_forward'] = float(np.std(predcv - np.polyval(A, y)))
# residual correlation of dMK and dMG among calibrators (core: |res| < 2.5 sigma)
k2 = (cal.ruwe < 1.1).values
rk = cal.dMK.values[k2] - (fits['dMK_ruwe1.1']['a'] + fits['dMK_ruwe1.1']['b'] * cal.feh.values[k2])
rg = cal.dMG.values[k2] - (fits['dMG_ruwe1.1']['a'] + fits['dMG_ruwe1.1']['b'] * cal.feh.values[k2])
cc_ = (np.abs(rk) < 2.5 * fits['dMK_ruwe1.1']['scatter']) & (np.abs(rg) < 2.5 * fits['dMG_ruwe1.1']['scatter'])
out['core_resid_corr_KG'] = float(np.corrcoef(rk[cc_], rg[cc_])[0, 1])
# high-precision (non-Birky) calibrators only: dMK slope and scatter
hp = (cal.ruwe < 1.1) & ~cal.src.eq('Birky20')
if hp.sum() > 30:
    out['fits_highprec'] = {}
    for q in ('dMK', 'dMG', 'rGK', 'rGW1'):
        kk = hp.values & np.isfinite(cal[q].values) & np.isfinite(cal.feh.values)
        if q == 'rGW1': kk &= cal.wiseok.values
        if kk.sum() > 15: out['fits_highprec'][q] = mixfit(cal.feh.values[kk], cal[q].values[kk], np.full(kk.sum(), 0.03))
    out['fits_highprec']['n'] = int(hp.sum())
    print('high-precision calibrators (non-Birky):', {q: (round(v['b'], 3), round(v['a'], 3), round(v['scatter'], 3), v['n']) for q, v in out['fits_highprec'].items() if q != 'n'}, hp.sum())
print('colour-only forward calibration', A, 'core resid corr K-G', out['core_resid_corr_KG'])
json.dump(out, open(f'{B}/out/h13_colour_feh.json', 'w'), indent=1)
# ---- figure
fig, ax = plt.subplots(2, 3, figsize=(15, 8)); ax = ax.ravel()
for i, q in enumerate(['dMK', 'dMG', 'rGK', 'rJK', 'rKW1']):
    kk = cal.ruwe < 1.1
    if q in ('rKW1',): kk &= cal.wiseok
    ax[i].scatter(cal.feh[kk], cal[q][kk], s=5, alpha=0.4)
    f = fits[f'{q}_ruwe1.1']; xx = np.linspace(-0.8, 0.6, 10); ax[i].plot(xx, f['a'] + f['b'] * xx, 'k-')
    ax[i].axhline(out['target'][q], c='r', ls='--', label='M3V'); ax[i].set_xlabel('[Fe/H] (spectroscopic)'); ax[i].set_ylabel(q); ax[i].legend()
    if q in ('dMK', 'dMG'): ax[i].invert_yaxis()
ax[5].scatter(y, predcv, s=5, alpha=0.4); ax[5].plot([-1, 0.6], [-1, 0.6], 'k-'); ax[5].axhline(fe_col, c='r', ls='--')
ax[5].set_xlabel('[Fe/H] spec'); ax[5].set_ylabel('colour-only [Fe/H] (CV)')
plt.tight_layout(); plt.savefig(f'{B}/figs/h13_colour_feh.png', dpi=110)
