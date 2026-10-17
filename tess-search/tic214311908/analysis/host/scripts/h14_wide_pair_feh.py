"""h14: the wide companions as a metallicity probe. TIC 214311910 A (G 14.654) and B (G 15.741, 0.83" from A) share
one blended Gaia BP/RP measurement (the A entry; excess factor consistent with the A+B G sum), so the blend has
BP-RP = 2.878 at M_G(A+B) = 10.231, while A and B have separate G magnitudes (M_G 10.613, 11.700) at the M3V's
parallax (bound system, El-Badry+21 R = 2e-5). If A and B are single main-sequence stars of the system's [Fe/H],
their individual BP-RP follow from M_G and [Fe/H]; the predicted blend colour is compared with the observed one.
Empirical single-star model from Mann-scale calibrators (h13 calib set, RUWE < 1.1, 2.2 < BP-RP < 3.5, plx/err > 20):
  M_G = poly3(BP-RP) + (c0 + c1 (BP-RP - 2.75)) [Fe/H], intrinsic scatter fitted;
  BP-G and G-RP as functions of BP-RP from the field sample (|C*| < 3; cubic fit to 0.02-mag binned medians).
For each [Fe/H] on a grid: Monte Carlo over the intrinsic scatter (each star's M_G offset at fixed colour drawn from
the fitted scatter), G-magnitude errors of A and B (0.02 / 0.05 mag, blending), parallax; likelihood of the observed
blend BP-RP (error 0.01 mag + 0.02 mag for the window losing part of B's light). Also the same model applied to the
M3V: [Fe/H] needed if single, and twin case. Caveats: A or B may themselves be unresolved binaries (A's RUWE 13.7 is
attributed to B); B's 2-parameter G may be biased.
Output: out/h14_wide_pair_feh.json, figs/h14_wide_pair_feh.png"""
import json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from scipy.optimize import minimize
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
B_ = 'tic214311908/analysis/host'
C = pd.read_csv(f'{B_}/data/calibrators_gaia.csv')
def Cs(d):
    x = d.bp_rp; f = 1.162004 + 0.011464 * x + 0.049255 * x ** 2 - 0.005879 * x ** 3
    return (d.phot_bp_rp_excess_factor - f) / (0.0059898 + 8.817481e-12 * d.phot_g_mean_mag ** 7.618399)
C = C[(C.scale != 'LAMOST') & np.isfinite(C.feh)]
C['w'] = 1 / np.clip(C.efeh.fillna(0.12), 0.05, 0.5) ** 2
g = C.groupby('source_id').apply(lambda d: pd.Series(dict(feh=np.sum(d.w * d.feh) / np.sum(d.w)))).reset_index()
C = g.merge(C.drop_duplicates('source_id').drop(columns=['feh']), on='source_id')
C['Cs'] = Cs(C); C['MG'] = C.phot_g_mean_mag + 5 * np.log10(C.parallax / 100)
C['MKx'] = C.ks_m + 5 * np.log10(C.parallax / 100)
_m = C.bp_rp.between(2.2, 3.5) & np.isfinite(C.MKx) & np.isfinite(C.bp_rp)
_rr = np.polyfit(C.bp_rp[_m], C.MKx[_m], 2)
k = C.bp_rp.between(2.2, 3.5) & (C.parallax / C.parallax_error > 20) & (C.ruwe < 1.1) & (C.Cs.abs() < 3) & (C.ipd_frac_multi_peak <= 2) & (np.abs(C.MKx - np.polyval(_rr, C.bp_rp)) < 1.5)
c = C[k]
x = c.bp_rp.values; y = c.MG.values; fe = c.feh.values
def model(p, x, fe):
    return np.polyval(p[:4], x - 2.75) + (p[4] + p[5] * (x - 2.75)) * fe
def nll(p):
    s = np.exp(p[6]); so = np.exp(p[7]); f = 1 / (1 + np.exp(-p[8]))
    r = y - model(p, x, fe)
    core = np.exp(-0.5 * r ** 2 / s ** 2) / s; outl = np.exp(-0.5 * (r + 0.6) ** 2 / so ** 2) / so   # outliers: over-luminous binaries
    return -np.sum(np.log((1 - f) * core + f * outl + 1e-300))
p0 = np.r_[np.polyfit(x - 2.75, y, 3), -1.0, 0.0, np.log(0.2), np.log(0.4), -2]
r = minimize(nll, p0, method='Nelder-Mead', options=dict(maxiter=40000, maxfev=40000, xatol=1e-7, fatol=1e-9))
r = minimize(nll, r.x, method='Nelder-Mead', options=dict(maxiter=40000, maxfev=40000))
p = r.x; sig = float(np.exp(p[6]))
print(f'calibrators n={len(c)}; dMG/d[Fe/H] at BP-RP 2.75: {p[4]:+.3f} (slope with colour {p[5]:+.3f}); intrinsic scatter {sig:.3f} mag; outlier frac {1/(1+np.exp(-p[8])):.2f}')
F = pd.read_csv(f'{B_}/data/field_mdwarfs.csv')
F = F[(F.ruwe < 1.2) & (F.parallax / F.parallax_error > 20)]
F['Cs'] = Cs(F); F = F[F.Cs.abs() < 3]
F['BPG'] = F.phot_bp_mean_mag - F.phot_g_mean_mag; F['GRP'] = F.phot_g_mean_mag - F.phot_rp_mean_mag
bc = np.arange(2.32, 3.24, 0.02)
mb = np.array([np.median(F.BPG[np.abs(F.bp_rp - b0) < 0.01]) for b0 in bc]); mr_ = np.array([np.median(F.GRP[np.abs(F.bp_rp - b0) < 0.01]) for b0 in bc])
cb = np.polyfit(bc, mb, 3); cr = np.polyfit(bc, mr_, 3)
print('BP-G, G-RP poly residual sd', np.std(F.BPG - np.polyval(cb, F.bp_rp)), np.std(F.GRP - np.polyval(cr, F.bp_rp)))
xs = np.linspace(2.0, 3.8, 3601)
def colour_for(MG, feh, dM=0.0):
    # solve model(x, feh) = MG - dM (dM = star's intrinsic offset at fixed colour: brighter star -> negative)
    v = model(p, xs, feh) + dM
    return xs[np.argmin(np.abs(v - MG))]
plx = 15.549; DM = 5 * np.log10(100 / plx)
GA, GB = 14.654037, 15.740896
obs_bprp = 15.977035 - 13.099223
rng = np.random.default_rng(4)
fgrid = np.round(np.arange(-0.6, 0.61, 0.05), 3)
res = []
for f_ in fgrid:
    pred = []
    for _ in range(600):
        MGA = GA + rng.normal(0, 0.02) - DM - rng.normal(0, 0.036 * 5 / np.log(10) / plx)
        MGB = GB + rng.normal(0, 0.05) - DM
        xa = colour_for(MGA, f_, rng.normal(0, sig)); xb = colour_for(MGB, f_, rng.normal(0, sig))
        GAa, GBb = MGA + DM, MGB + DM
        BPa = GAa + np.polyval(cb, xa); RPa = GAa - np.polyval(cr, xa); BPb = GBb + np.polyval(cb, xb); RPb = GBb - np.polyval(cr, xb)
        BP = -2.5 * np.log10(10 ** (-0.4 * BPa) + 10 ** (-0.4 * BPb)); RP = -2.5 * np.log10(10 ** (-0.4 * RPa) + 10 ** (-0.4 * RPb))
        pred.append((BP - RP, xa, xb))
    pred = np.array(pred)
    e_obs = np.hypot(0.01, 0.02)
    L = np.mean(np.exp(-0.5 * (pred[:, 0] - obs_bprp) ** 2 / e_obs ** 2))
    res.append(dict(feh=float(f_), pred_bprp_med=float(np.median(pred[:, 0])), pred_bprp_sd=float(np.std(pred[:, 0])), xa=float(np.median(pred[:, 1])),
                    xb=float(np.median(pred[:, 2])), L=float(L)))
    print(f'[Fe/H] {f_:+.2f}: blend BP-RP {np.median(pred[:,0]):.3f} +- {np.std(pred[:,0]):.3f} (A {np.median(pred[:,1]):.3f}, B {np.median(pred[:,2]):.3f}) obs {obs_bprp:.3f} L {L:.3f}')
Ls = np.array([r_['L'] for r_ in res]); post = Ls / Ls.sum()
mean = float(np.sum(post * fgrid)); sd = float(np.sqrt(np.sum(post * (fgrid - mean) ** 2)))
# M3V in the same model
MG3 = 14.281043 - DM
fe_single = float(fgrid[np.argmin([abs(model(p, np.array([2.756149]), f_)[0] - MG3) for f_ in fgrid])])
fs = np.linspace(-1, 1, 2001); vv = np.array([model(p, np.array([2.756149]), f_)[0] for f_ in fs])
fe_single = float(fs[np.argmin(np.abs(vv - MG3))]); fe_twin = float(fs[np.argmin(np.abs(vv - (MG3 + 0.753)))])
out = dict(model_params=p.tolist(), intrinsic_scatter=sig, n_cal=int(len(c)), dMG_dfeh_at_2p75=float(p[4]), obs_blend_bprp=obs_bprp,
           grid=res, posterior_mean_feh=mean, posterior_sd_feh=sd, m3v_feh_if_single=fe_single, m3v_feh_if_twin=fe_twin,
           m3v_feh_single_err=float(sig / abs(p[4])))
print(f'Wide pair A+B posterior (flat prior): [Fe/H] = {mean:+.2f} +- {sd:.2f}; M3V needs {fe_single:+.2f} (single) or {fe_twin:+.2f} (twin) in the same model (+-{sig/abs(p[4]):.2f})')
json.dump(out, open(f'{B_}/out/h14_wide_pair_feh.json', 'w'), indent=1)
fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
ax[0].plot(fgrid, [r_['pred_bprp_med'] for r_ in res], 'k-'); ax[0].fill_between(fgrid, [r_['pred_bprp_med'] - r_['pred_bprp_sd'] for r_ in res], [r_['pred_bprp_med'] + r_['pred_bprp_sd'] for r_ in res], alpha=0.3)
ax[0].axhline(obs_bprp, c='r'); ax[0].set_xlabel('[Fe/H] of the system'); ax[0].set_ylabel('predicted A+B blend BP-RP')
ax[1].plot(fgrid, post, 'k-'); ax[1].axvline(fe_single, c='b', ls='--', label='M3V single needs'); ax[1].axvline(fe_twin, c='g', ls='--', label='M3V twin needs')
ax[1].set_xlabel('[Fe/H]'); ax[1].set_ylabel('posterior (flat prior)'); ax[1].legend()
plt.tight_layout(); plt.savefig(f'{B_}/figs/h14_wide_pair_feh.png', dpi=110)
