"""Step 4: combine sectors, compare with the injection calibrations, likelihood ratios, figures.

usage: python s04_combine.py <real_tag> <inj_tag e.g. prf_prf> [label]
Outputs: combine_<real_tag>__<inj_tag>.json, fig_injcal_<...>.png
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.stats import norm, gaussian_kde, chi2 as chi2dist
from scipy.optimize import least_squares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pipeline import *

real_tag = sys.argv[1] if len(sys.argv) > 1 else 'prf_main'
inj_tag = sys.argv[2] if len(sys.argv) > 2 else 'prf_prf'
label = sys.argv[3] if len(sys.argv) > 3 else f'{real_tag}__{inj_tag}'
rng = np.random.default_rng(11)

real = json.load(open(f'{HERE}/diff_{real_tag}.json'))
inj = np.load(f'{HERE}/inj_{inj_tag}.npz')
keys = list(inj['keys'])
ki = {k: i for i, k in enumerate(keys)}


def gauss_lnL(x, v):
    v = v[np.isfinite(v)]
    return norm.logpdf(x, np.mean(v), np.std(v, ddof=1))


def gauss2_lnL(x, V):
    V = V[np.all(np.isfinite(V), axis=1)]
    mu = V.mean(0); C = np.cov(V.T)
    d = x - mu
    return float(-0.5 * d @ np.linalg.solve(C, d) - 0.5 * np.log(np.linalg.det(2 * np.pi * C))), float(np.sqrt(d @ np.linalg.solve(C, d)))


out = dict(real_tag=real_tag, inj_tag=inj_tag, sectors={})
tot = dict(lnLR_spar=0.0, lnLR_EN=0.0, lnLR_dchi2=0.0)
fig, axes = plt.subplots(2, 6, figsize=(21, 7))
for j, sec in enumerate(SECTORS):
    r = real[str(sec)]['stacked']
    T = inj[f's{sec}_resampled_T']; N = inj[f's{sec}_resampled_N']
    o = {}
    for stat, key in [('s_par', 's_par'), ('dchi2', 'dchi2_TminusN'), ('f_N', 'f_N')]:
        x = r[key]
        vT = T[:, ki[key]]; vN = N[:, ki[key]]
        vT = vT[np.isfinite(vT)]; vN = vN[np.isfinite(vN)]
        lT, lN = gauss_lnL(x, vT), gauss_lnL(x, vN)
        zT = (x - np.mean(vT)) / np.std(vT, ddof=1); zN = (x - np.mean(vN)) / np.std(vN, ddof=1)
        # empirical tail fractions (two-sided about the median)
        pT = float(np.mean(np.abs(vT - np.median(vT)) >= abs(x - np.median(vT))))
        pN = float(np.mean(np.abs(vN - np.median(vN)) >= abs(x - np.median(vN))))
        # KDE
        try:
            kT = gaussian_kde(vT)(x)[0]; kN = gaussian_kde(vN)(x)[0]
            lkde = float(np.log(max(kT, 1e-300)) - np.log(max(kN, 1e-300)))
        except Exception:
            lkde = np.nan
        o[stat] = dict(real=x, T_mean=float(np.mean(vT)), T_std=float(np.std(vT, ddof=1)), N_mean=float(np.mean(vN)),
                       N_std=float(np.std(vN, ddof=1)), z_T=float(zT), z_N=float(zN), p_emp_T=pT, p_emp_N=pN,
                       lnLR_TN_gauss=float(lT - lN), lnLR_TN_kde=lkde, nsim=int(len(vT)))
    xEN = np.array([r['offT_E'], r['offT_N']])
    lT2, mT = gauss2_lnL(xEN, T[:, [ki['offT_E'], ki['offT_N']]])
    lN2, mN = gauss2_lnL(xEN, N[:, [ki['offT_E'], ki['offT_N']]])
    o['EN2d'] = dict(real=xEN.tolist(), lnLR_TN=lT2 - lN2, mahal_T=mT, mahal_N=mN)
    tot['lnLR_spar'] += o['s_par']['lnLR_TN_gauss']
    tot['lnLR_EN'] += lT2 - lN2
    tot['lnLR_dchi2'] += o['dchi2']['lnLR_TN_gauss']
    out['sectors'][str(sec)] = o
    # figure: s_par and dchi2 histograms
    for row, (key, lab) in enumerate([('s_par', 's_par (0 = target, 1 = neighbour)'), ('dchi2_TminusN', r'$\Delta\chi^2$ = $\chi^2_T-\chi^2_N$')]):
        ax = axes[row, j]
        vT = T[:, ki[key]]; vN = N[:, ki[key]]
        lo, hi = np.nanpercentile(np.r_[vT, vN], [0.5, 99.5])
        lo = min(lo, r[key]); hi = max(hi, r[key])
        bins = np.linspace(lo - 0.05 * (hi - lo), hi + 0.05 * (hi - lo), 50)
        ax.hist(vT, bins, color='tab:blue', alpha=0.55, label='injected on target')
        ax.hist(vN, bins, color='tab:orange', alpha=0.55, label='injected on neighbour')
        ax.axvline(r[key], color='k', lw=2, label='real')
        ax.set_title(f'S{sec}', fontsize=10); ax.set_xlabel(lab, fontsize=8)
        if j == 0:
            ax.legend(fontsize=7)

# combined statistic: inverse-variance weighted mean of s_par (weights from the H_N sims) and sum of dchi2
w = np.array([1 / out['sectors'][str(s)]['s_par']['N_std'] ** 2 for s in SECTORS])
w /= w.sum()
xr = sum(w[k] * real[str(s)]['stacked']['s_par'] for k, s in enumerate(SECTORS))
dr = sum(real[str(s)]['stacked']['dchi2_TminusN'] for s in SECTORS)
NS = 20000
combT = np.zeros(NS); combN = np.zeros(NS); dT = np.zeros(NS); dN = np.zeros(NS)
for k, s in enumerate(SECTORS):
    T = inj[f's{s}_resampled_T']; N = inj[f's{s}_resampled_N']
    iT = rng.integers(0, len(T), NS); iN = rng.integers(0, len(N), NS)
    combT += w[k] * T[iT, ki['s_par']]; combN += w[k] * N[iN, ki['s_par']]
    dT += T[iT, ki['dchi2_TminusN']]; dN += N[iN, ki['dchi2_TminusN']]
comb = dict(weights=w.tolist(), s_par_real=float(xr),
            s_par_T_mean=float(np.nanmean(combT)), s_par_T_std=float(np.nanstd(combT)),
            s_par_N_mean=float(np.nanmean(combN)), s_par_N_std=float(np.nanstd(combN)),
            z_T=float((xr - np.nanmean(combT)) / np.nanstd(combT)), z_N=float((xr - np.nanmean(combN)) / np.nanstd(combN)),
            frac_T_sims_ge_real=float(np.mean(combT >= xr)), frac_N_sims_le_real=float(np.mean(combN <= xr)),
            dchi2_sum_real=float(dr), dchi2_T_mean=float(dT.mean()), dchi2_T_std=float(dT.std()),
            dchi2_N_mean=float(dN.mean()), dchi2_N_std=float(dN.std()),
            z_dchi2_T=float((dr - dT.mean()) / dT.std()), z_dchi2_N=float((dr - dN.mean()) / dN.std()),
            frac_T_sims_dchi2_ge_real=float(np.mean(dT >= dr)))
comb['lnLR_TN_combined_spar_gauss'] = float(norm.logpdf(xr, comb['s_par_T_mean'], comb['s_par_T_std']) - norm.logpdf(xr, comb['s_par_N_mean'], comb['s_par_N_std']))
comb['lnLR_TN_combined_dchi2_gauss'] = float(norm.logpdf(dr, comb['dchi2_T_mean'], comb['dchi2_T_std']) - norm.logpdf(dr, comb['dchi2_N_mean'], comb['dchi2_N_std']))
comb.update({f'sum_{k}': v for k, v in tot.items()})
out['combined'] = comb
ax = axes[0, 5]
bins = np.linspace(min(combT.min(), xr) - 0.1, max(combN.max(), xr) + 0.1, 60)
ax.hist(combT, bins, color='tab:blue', alpha=0.55); ax.hist(combN, bins, color='tab:orange', alpha=0.55); ax.axvline(xr, color='k', lw=2)
ax.set_title('all sectors: weighted s_par', fontsize=10); ax.set_xlabel('combined s_par', fontsize=8)
ax = axes[1, 5]
bins = np.linspace(min(dT.min(), dr) - 10, max(dN.max(), dr) + 10, 60)
ax.hist(dT, bins, color='tab:blue', alpha=0.55); ax.hist(dN, bins, color='tab:orange', alpha=0.55); ax.axvline(dr, color='k', lw=2)
ax.set_title(r'all sectors: $\Sigma\Delta\chi^2$', fontsize=10); ax.set_xlabel(r'$\Sigma(\chi^2_T-\chi^2_N)$', fontsize=8)
fig.suptitle(f'Injection-recovery calibration ({inj_tag}: inject-PSF_analysis-PSF) vs real ({real_tag})', fontsize=11)
plt.tight_layout(); plt.savefig(f'{HERE}/fig_injcal_{label}.png', dpi=105); plt.close()
to_json(out, f'{HERE}/combine_{label}.json')
print(json.dumps(comb, indent=1))
for s in SECTORS:
    o = out['sectors'][str(s)]
    print(f"S{s}: s_par real {o['s_par']['real']:+.3f} | H_T {o['s_par']['T_mean']:+.3f}+-{o['s_par']['T_std']:.3f} (z={o['s_par']['z_T']:+.1f}, p_emp={o['s_par']['p_emp_T']:.4f}) | "
          f"H_N {o['s_par']['N_mean']:+.3f}+-{o['s_par']['N_std']:.3f} (z={o['s_par']['z_N']:+.1f}) | lnLR_TN gauss {o['s_par']['lnLR_TN_gauss']:+.1f} kde {o['s_par']['lnLR_TN_kde']:+.1f}; "
          f"EN2d lnLR {o['EN2d']['lnLR_TN']:+.1f} (mahal T {o['EN2d']['mahal_T']:.1f}, N {o['EN2d']['mahal_N']:.1f}); dchi2 lnLR {o['dchi2']['lnLR_TN_gauss']:+.1f}")
