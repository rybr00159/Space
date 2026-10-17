"""h12b: what the Gaia XP spectrum of the M3V says about contamination, [Fe/H] and a twin.
Inputs: data/xp/xps_chunk_*.npz (h11b; official sampled XP of calibrators, field M dwarfs, the M3V, A(+B), the Mira),
out/xp_target_sampled.csv (h01, consistency check), data/calib_sel.csv-like info from calibrators_gaia.csv and field_mdwarfs.csv.
Shape vector: log10 of the sampled flux on 420-1000 nm (2 nm) normalised by its mean over 620-1000 nm.
(1) Contamination: residual of the M3V's shape from the median of colour-matched (|dBP-RP| < 0.03) normal field stars,
    in units of the star-to-star scatter; compared with the distribution of the same statistic for the field stars
    themselves; fitted fraction of the Mira's spectrum (scaled by the Mira/M3V flux ratio) in the M3V's spectrum.
(2) [Fe/H] from XP shape (no parallax, no absolute magnitude): PCA (20-40 comps) of the shape vectors of Mann-scale
    calibrators with 2.45 < BP-RP < 3.10 and |dMKs| < 1.3 (main-sequence sanity cut, not used as a feature), then ridge regression on PCs + BP-RP; 10-fold CV picks n_comp and penalty;
    CV scatter, bias in the metal-rich tail; applied to the M3V with Monte Carlo over its XP flux errors.
    Same model applied to the LAMOST (Ding+22) stars as an external check, and to the field 'normal' vs
    'over-luminous' samples (the empirical twin null: over-luminous stars that are twins should look solar).
(3) Composite fits: M3V shape vs single templates and vs two-star composites (q = 0.6-1.0) built from the colour-binned
    median spectra of normal field stars, each weighted by its absolute flux (M_G from the field CMD ridge).
Output: out/h12b_xp.json, figs/h12b_xp_*.png, out/h12b_xp_pred_field.csv"""
import glob, json, numpy as np, pandas as pd, warnings, sys
warnings.filterwarnings('ignore')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
B = 'tic214311908/analysis/host'
sys.path.insert(0, f'{B}/scripts')
rng = np.random.default_rng(8)
wl = np.arange(336, 1021, 2.0)
win = (wl >= 420) & (wl <= 1000); wn = (wl >= 620) & (wl <= 1000)
Z = [np.load(f) for f in sorted(glob.glob(f'{B}/data/xp/xps_chunk_*.npz'))]
sid = np.concatenate([z['source_id'] for z in Z]); FL = np.vstack([z['flux'] for z in Z if len(z['source_id'])]).astype(float); FE = np.vstack([z['flux_error'] for z in Z if len(z['source_id'])]).astype(float)
_, ui = np.unique(sid, return_index=True); ui = np.sort(ui); sid, FL, FE = sid[ui], FL[ui], FE[ui]
ids = pd.read_csv(f'{B}/data/xp/xp_ids.csv')
TID, AID, MID = 5925373272850784768, 5925373272850782208, 5925373268528684288
def shape(f):
    f = np.clip(f, 1e-22, None)
    return np.log10(f[..., win] / f[..., wn].mean(axis=-1, keepdims=True))
S = shape(FL)
# target
it = np.where(sid == TID)[0][0]; tf, te = FL[it].copy(), FE[it].copy()
# consistency with the GaiaXPy calibration of the continuous spectrum (h01)
tg = pd.read_csv(f'{B}/out/xp_target_sampled.csv'); tg = tg[tg.source_id == TID].sort_values('wl_nm')
print('VizieR sampled vs GaiaXPy-calibrated M3V spectrum: median ratio', float(np.median(tf / tg.flux.values)), 'rms', float(np.std(tf / tg.flux.values)))
im = np.where(sid == MID)[0]; mf_viz = FL[im[0]] if len(im) else None
keepm = ~np.isin(sid, [TID, AID, MID]); sid, FL, FE = sid[keepm], FL[keepm], FE[keepm]; S = shape(FL)
St = shape(tf)
# metadata
C = pd.read_csv(f'{B}/data/calibrators_gaia.csv')
C = C[np.isfinite(C.feh)]
C['w'] = 1 / np.clip(C.efeh.fillna(0.12), 0.05, 0.5) ** 2
lam = C[C.scale == 'LAMOST'].groupby('source_id').feh.mean()
m13 = C[C.scale != 'LAMOST'].groupby('source_id').apply(lambda d: np.sum(d.w * d.feh) / np.sum(d.w))
srcs = C[C.scale != 'LAMOST'].groupby('source_id').src.apply(lambda s: '+'.join(sorted(set(s))))
gC = C.drop_duplicates('source_id').set_index('source_id')
F = pd.read_csv(f'{B}/data/field_mdwarfs.csv').set_index('source_id')
meta = pd.DataFrame(index=sid)
meta['set'] = ids.set_index('source_id').set.reindex(sid).values
meta['bp_rp'] = gC.bp_rp.reindex(sid).values
meta.loc[meta.bp_rp.isna(), 'bp_rp'] = F.bp_rp.reindex(sid[meta.bp_rp.isna().values]).values
meta['G'] = gC.phot_g_mean_mag.reindex(sid).values
meta.loc[meta.G.isna(), 'G'] = F.phot_g_mean_mag.reindex(sid[meta.G.isna().values]).values
meta['ruwe'] = gC.ruwe.reindex(sid).values
meta.loc[meta.ruwe.isna(), 'ruwe'] = F.ruwe.reindex(sid[meta.ruwe.isna().values]).values
meta['feh'] = m13.reindex(sid).values; meta['feh_lam'] = lam.reindex(sid).values; meta['src'] = srcs.reindex(sid).values
plx = gC.parallax.reindex(sid).values.copy(); plx[np.isnan(plx)] = F.parallax.reindex(sid[np.isnan(plx)]).values
ks = gC.ks_m.reindex(sid).values.copy(); ks[np.isnan(ks)] = F.ks_m.reindex(sid[np.isnan(ks)]).values
meta['MK'] = ks + 5 * np.log10(plx / 100); meta['MG'] = meta.G + 5 * np.log10(plx / 100)
ok = np.isfinite(S).all(1) & np.isfinite(meta.bp_rp.values)
print('XP spectra loaded', len(sid), 'usable', ok.sum(), meta.set.value_counts().to_dict())
out = {}
# ---------------- (1) contamination / peculiarity
nf = ok & (meta.set.values == 'field_normal')
cm = nf & (np.abs(meta.bp_rp.values - 2.756) < 0.03)
med = np.median(S[cm], axis=0); sd = 1.4826 * np.median(np.abs(S[cm] - med), axis=0)
def chi(s): return np.sqrt(np.mean(((s - med) / sd) ** 2))
chis = np.array([chi(S[i]) for i in np.where(cm)[0]])
out['contamination'] = dict(n_colour_matched_field=int(cm.sum()), target_rms_z=float(chi(St)), field_rms_z_median=float(np.median(chis)),
                            field_rms_z_p90=float(np.percentile(chis, 90)), target_percentile=float(np.mean(chis <= chi(St))))
print('peculiarity (rms z vs colour-matched normal field stars):', out['contamination'])
# Mira fraction fit
mf = mf_viz
# model: target = (1-e) * template(med) + e * Mira  (in flux, normalised at the same window)
tmpl = 10 ** med
tnorm = tf[win] / tf[wn].mean(); mnorm = mf[win] / tf[wn].mean()
best = []
for e in np.linspace(-0.05, 0.10, 301):
    mod = tmpl * (1 - 0) + e * mnorm
    mod = mod / (mod[wn[win]].mean() if False else 1.0)
    r = (np.log10(np.clip(tnorm, 1e-22, None)) - np.log10(np.clip(mod, 1e-22, None))) / sd
    best.append((e, np.sum(r ** 2)))
best = np.array(best); eb = best[np.argmin(best[:, 1]), 0]
dchi = best[:, 1] - best[:, 1].min()
# uncertainty: rescale chi2 so that min chi2/dof = 1
dof = win.sum() - 1; scale = best[:, 1].min() / dof
okk = dchi / scale < 1
out['mira_fraction'] = dict(best_fraction_of_mira_flux=float(eb), range_1sigma=[float(best[okk, 0].min()), float(best[okk, 0].max())],
                            note='fraction of the Mira XP flux (absolute) added to a colour-matched template, chi2 rescaled; target-relative contribution at 800 nm = e * Mira/M3V flux ratio')
ratio800 = mf[np.argmin(np.abs(wl - 800))] / tf[np.argmin(np.abs(wl - 800))]
out['mira_fraction']['mira_over_m3v_flux_800nm'] = float(ratio800)
print('Mira contamination fit:', out['mira_fraction'])
# ---------------- (2) XP [Fe/H] regression
dmk_all = meta.MK.values - np.polyval(np.polyfit(meta.bp_rp.values[nf], meta.MK.values[nf], 2), meta.bp_rp.values)
cal = ok & np.isfinite(meta.feh.values) & meta.set.isin(['cal_m13', 'cal_birky']).values & (meta.bp_rp.values > 2.45) & (meta.bp_rp.values < 3.10) & (meta.ruwe.values < 1.4) & (np.abs(dmk_all) < 1.3)
Xc = S[cal]; yc = meta.feh.values[cal]; cc = meta.bp_rp.values[cal]
mu = Xc.mean(0); U, sv, Vt = np.linalg.svd(Xc - mu, full_matrices=False)
def feats(Sx, bprp, nc):
    P = (Sx - mu) @ Vt[:nc].T
    return np.column_stack([P, bprp - 2.75, (bprp - 2.75) ** 2])
def fit(Xf, y, lam_):
    m_ = Xf.mean(0); s_ = Xf.std(0); Zf = (Xf - m_) / s_
    w = np.linalg.solve(Zf.T @ Zf + lam_ * np.eye(Zf.shape[1]), Zf.T @ (y - y.mean()))
    return lambda Xn: y.mean() + ((Xn - m_) / s_) @ w
fold = rng.integers(0, 10, len(yc)); bestcv = None
for nc in (10, 20, 30, 40):
    for lam_ in (1, 10, 30, 100, 300):
        pred = np.zeros(len(yc))
        for f_ in range(10):
            tr = fold != f_; mdl = fit(feats(Xc[tr], cc[tr], nc), yc[tr], lam_); pred[~tr] = mdl(feats(Xc[~tr], cc[~tr], nc))
        r = pred - yc; rs = 1.4826 * np.median(np.abs(r - np.median(r)))
        if bestcv is None or rs < bestcv[2]: bestcv = (nc, lam_, rs, np.std(r), pred.copy())
nc, lam_, rs, rms, predcv = bestcv
print(f'XP [Fe/H]: n_cal={len(yc)}, best n_comp={nc}, lambda={lam_}: CV robust scatter {rs:.3f}, rms {rms:.3f}, corr {np.corrcoef(predcv, yc)[0,1]:.2f}')
mdl = fit(feats(Xc, cc, nc), yc, lam_)
fe_t = float(mdl(feats(St[None, :], np.array([2.756149]), nc))[0])
mc = []
for _ in range(500):
    ff = tf + rng.normal(0, 1, len(tf)) * te
    mc.append(mdl(feats(shape(ff)[None, :], np.array([2.756149]), nc))[0])
# calibration of the regression (CV pred vs truth slope) -> de-bias
A = np.polyfit(yc, predcv, 1)       # pred = A0*true + A1
hi = yc > 0.25
# colour-local check: calibrators within 0.08 of the target colour
loc = np.abs(cc - 2.756) < 0.08
out['xp_feh'] = dict(n_cal=int(len(yc)), n_comp=nc, ridge_lambda=lam_, cv_robust_scatter=float(rs), cv_rms=float(rms),
                     cv_robust_scatter_local_colour=float(1.4826 * np.median(np.abs((predcv - yc)[loc] - np.median((predcv - yc)[loc])))),
                     n_local=int(loc.sum()), target_pred=fe_t, target_pred_mc_sd=float(np.std(mc)),
                     cv_pred_vs_true_slope=float(A[0]), cv_pred_vs_true_intercept=float(A[1]),
                     target_debiased=float((fe_t - A[1]) / A[0]),
                     cv_median_pred_true_gt_0p25=float(np.median(predcv[hi])) if hi.sum() else None, n_true_gt_0p25=int(hi.sum()),
                     cv_median_pred_true_m0p1_0p1=float(np.median(predcv[np.abs(yc) < 0.1])))
print('XP [Fe/H] of the M3V:', out['xp_feh'])
# probability calibration: among calibrators near the target's colour with CV-pred within +-0.05 of the target's prediction, true [Fe/H] distribution
near = loc & (np.abs(predcv - fe_t) < 0.07)
out['xp_feh']['true_feh_of_cal_with_similar_pred'] = dict(n=int(near.sum()), median=float(np.median(yc[near])) if near.sum() else None,
                                                          p16=float(np.percentile(yc[near], 16)) if near.sum() > 3 else None,
                                                          p84=float(np.percentile(yc[near], 84)) if near.sum() > 3 else None,
                                                          frac_gt_0p3=float(np.mean(yc[near] > 0.3)) if near.sum() else None)
print('  calibrators with similar XP prediction:', out['xp_feh']['true_feh_of_cal_with_similar_pred'])
# external check on LAMOST
lm = ok & np.isfinite(meta.feh_lam.values) & (meta.set.values == 'cal_lamost') & (meta.bp_rp.values > 2.45) & (meta.bp_rp.values < 3.10)
pl = mdl(feats(S[lm], meta.bp_rp.values[lm], nc)); yl = meta.feh_lam.values[lm]
out['xp_feh']['lamost_check'] = dict(n=int(lm.sum()), corr=float(np.corrcoef(pl, yl)[0, 1]) if lm.sum() > 3 else None,
                                     median_diff_pred_minus_lamost=float(np.median(pl - yl)) if lm.sum() else None,
                                     robust_sd=float(1.4826 * np.median(np.abs(pl - yl - np.median(pl - yl)))) if lm.sum() else None)
print('  LAMOST check', out['xp_feh']['lamost_check'])
# leakage of over-luminosity into the XP prediction: CV residual vs dMK among calibrators (binaries)
dmk_c = meta.MK.values[cal] - np.polyval(np.polyfit(meta.bp_rp.values[nf], meta.MK.values[nf], 2), cc)
for lo, hi_ in ((-0.3, 0.3), (-1.0, -0.45)):
    s_ = (dmk_c > lo) & (dmk_c < hi_)
    out['xp_feh'][f'cv_resid_dMK_{lo}_{hi_}'] = dict(n=int(s_.sum()), median=float(np.median((predcv - yc)[s_])) if s_.sum() else None)
print('  CV residual vs dMK:', {k: v for k, v in out['xp_feh'].items() if k.startswith('cv_resid_dMK')})
# ---------------- empirical twin null: field normal vs over-luminous at the target colour
ridgeK = np.polyfit(meta.bp_rp.values[nf], meta.MK.values[nf], 2)
fld = ok & meta.set.isin(['field_normal', 'field_overlum']).values
pf = mdl(feats(S[fld], meta.bp_rp.values[fld], nc))
dmk_f = meta.MK.values[fld] - np.polyval(ridgeK, meta.bp_rp.values[fld])
pd.DataFrame(dict(source_id=sid[fld], set=meta.set.values[fld], bp_rp=meta.bp_rp.values[fld], ruwe=meta.ruwe.values[fld], dMK_local=dmk_f, feh_xp=pf)).to_csv(f'{B}/out/h12b_xp_pred_field.csv', index=False)
for st in ('field_normal', 'field_overlum'):
    k = meta.set.values[fld] == st
    out[f'xp_feh_{st}'] = dict(n=int(k.sum()), median=float(np.median(pf[k])), p16=float(np.percentile(pf[k], 16)), p84=float(np.percentile(pf[k], 84)),
                               frac_gt_target=float(np.mean(pf[k] >= fe_t)))
    print(st, out[f'xp_feh_{st}'])
dmk_t = 6.4406 - np.polyval(ridgeK, 2.756149)
out['target_dMK_local'] = float(dmk_t)
# ---------------- (3) composite fits
nfm = nf & (meta.ruwe.values < 1.2)
bins = np.arange(2.30, 3.40, 0.04)
tm, tMG = [], []
F2 = pd.read_csv(f'{B}/data/field_mdwarfs.csv')
F2 = F2[(F2.ruwe < 1.2) & (F2.parallax / F2.parallax_error > 20)]
F2['MG'] = F2.phot_g_mean_mag + 5 * np.log10(F2.parallax / 100)
rMG = np.polyfit(F2.bp_rp, F2.MG, 3)
# templates: colour-binned medians of all available XP spectra with normal CMD position (field normal + calibrators with |dMK|<0.25)
allS_ok = ok & (np.abs(meta.MK.values - np.polyval(np.polyfit(meta.bp_rp.values[nf], meta.MK.values[nf], 2), meta.bp_rp.values)) < 0.25)
absF = FL / FL[:, wn].mean(1, keepdims=True)
tc = []
for b0 in bins:
    k = allS_ok & (np.abs(meta.bp_rp.values - b0) < 0.03)
    if k.sum() >= 8:
        tc.append((b0, np.median(np.log10(np.clip(absF[k][:, win], 1e-30, None)), axis=0), k.sum()))
tcol = np.array([t[0] for t in tc]); tsp = np.array([t[1] for t in tc])
def template(col):
    return np.array([np.interp(col, tcol, tsp[:, j]) for j in range(tsp.shape[1])])
# relative absolute flux in the normalisation window for a star of colour col: from M_G and the G-normalisation
def comp_shape(c1, c2, MG1, MG2):
    s1 = 10 ** template(c1); s2 = 10 ** template(c2)
    w1 = 10 ** (-0.4 * MG1); w2 = 10 ** (-0.4 * MG2)
    # template normalised to mean over 620-1000 nm; convert to G-flux scaling using the mean ratio of G-band to window flux ~ same for both: approximate
    f = w1 * s1 + w2 * s2
    return np.log10(f / f[wn[win]].mean())
cols = np.linspace(2.35, 3.35, 201)
res_q = []
for c1 in np.linspace(2.45, 2.80, 36):
    for c2 in np.linspace(c1, 3.35, 40):
        MG1 = np.polyval(rMG, c1); MG2 = np.polyval(rMG, c2)
        sh = comp_shape(c1, c2, MG1, MG2)
        r = (St - sh) / sd
        res_q.append((c1, c2, MG1, MG2, np.sum(r ** 2)))
res_q = np.array(res_q)
single = np.array([(c, np.sum(((St - template(c)) / sd) ** 2)) for c in cols])
bs = single[np.argmin(single[:, 1])]; bc = res_q[np.argmin(res_q[:, 4])]
out['composite'] = dict(best_single_colour=float(bs[0]), chi2_single=float(bs[1]), best_comp_c1=float(bc[0]), best_comp_c2=float(bc[1]),
                        dMG_comp=float(bc[3] - bc[2]), chi2_comp=float(bc[4]), n_points=int(win.sum()),
                        note='chi2 uses the star-to-star scatter of colour-matched normal field stars as noise; templates are solar-neighbourhood medians')
print('composite vs single:', out['composite'])
json.dump(out, open(f'{B}/out/h12b_xp.json', 'w'), indent=1)
# ---------------- figures
fig, ax = plt.subplots(3, 1, figsize=(10, 11))
ax[0].plot(wl[win], med, 'k-', label='median normal field (|dBP-RP|<0.03)'); ax[0].fill_between(wl[win], med - sd, med + sd, color='gray', alpha=0.3)
ax[0].plot(wl[win], St, 'r-', label='M3V'); ax[0].set_ylabel('log flux (norm)'); ax[0].legend()
ax[1].plot(wl[win], (St - med) / sd, 'r-', label='M3V'); ax[1].axhline(0, c='k')
mr = cal & False
for lab_, sel_, colr in (('cal [Fe/H]>+0.25', (np.abs(cc - 2.756) < 0.05) & (yc > 0.25), 'b'), ('cal [Fe/H]<-0.25', (np.abs(cc - 2.756) < 0.05) & (yc < -0.25), 'g')):
    if sel_.sum() > 2:
        ax[1].plot(wl[win], (np.median(Xc[sel_], 0) - med) / sd, colr + '-', label=f'{lab_} (n={sel_.sum()})')
ax[1].set_ylabel('residual / scatter'); ax[1].legend(fontsize=8)
ax[2].scatter(yc, predcv, s=4, alpha=0.4, label='calibrators (CV)'); ax[2].plot([-1, 0.6], [-1, 0.6], 'k-'); ax[2].axhline(fe_t, c='r', ls='--', label='M3V XP prediction')
ax[2].set_xlabel('[Fe/H] spectroscopic'); ax[2].set_ylabel('[Fe/H] from XP shape'); ax[2].legend()
plt.tight_layout(); plt.savefig(f'{B}/figs/h12b_xp_main.png', dpi=110)
fig, ax = plt.subplots(figsize=(7, 4.5))
for st, cl in (('field_normal', 'k'), ('field_overlum', 'm')):
    k = meta.set.values[fld] == st
    ax.hist(pf[k], bins=np.arange(-0.8, 0.8, 0.05), histtype='step', color=cl, density=True, label=f'{st} (n={k.sum()})')
ax.axvline(fe_t, c='r', ls='--', label='M3V'); ax.set_xlabel('[Fe/H] from XP shape'); ax.legend()
plt.tight_layout(); plt.savefig(f'{B}/figs/h12b_xp_twin_null.png', dpi=110)
