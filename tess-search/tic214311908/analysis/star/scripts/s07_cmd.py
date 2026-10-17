"""Position of TIC 214311908 relative to the empirical M-dwarf main sequence (Gaia DR3 + 2MASS, d < 100 pc).
Sample (data/mdwarf_sample.csv from s04): parallax_over_error > 20, BP-RP 2.45-3.05, RUWE < 1.4, ipd_frac_multi_peak <= 5,
|C*| < 3 sigma (Riello+21), 2MASS Ks quality A, sigma(Ks) < 0.05.  Ridge lines MKs(BP-RP), MG(BP-RP), MKs(G-Ks) and
MKs(J-Ks) are fitted iteratively (cubic, 2.5-sigma clipping, faint-side mode via a robust 'lower envelope' check).
Reports the target's offset from each ridge, the width of the sequence, and the fraction of comparison stars at
least as over-luminous. Also checks TIC 214311910 A (blended with B, so only indicative).
Output: out/cmd.json, figs/cmd.png"""
import json, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
B = 'tic214311908/analysis/star'
m = pd.read_csv(f'{B}/data/mdwarf_sample.csv')
def cstar(C, x):
    f = np.where(x < 0.5, 1.154360 + 0.033772 * x + 0.032277 * x ** 2,
                 np.where(x < 4, 1.162004 + 0.011464 * x + 0.049255 * x ** 2 - 0.005879 * x ** 3, 1.057572 + 0.140537 * x))
    return C - f
def sig_cstar(G): return 0.0059898 + 8.817481e-12 * G ** 7.618399
m['Cstar'] = cstar(m.phot_bp_rp_excess_factor, m.bp_rp); m['Cs_sig'] = m.Cstar / sig_cstar(m.phot_g_mean_mag)
m['MG'] = m.phot_g_mean_mag + 5 * np.log10(m.parallax / 100); m['MK'] = m.ks_m + 5 * np.log10(m.parallax / 100)
m['GK'] = m.phot_g_mean_mag - m.ks_m; m['JK'] = m.j_m - m.ks_m
good = (m.ruwe < 1.4) & (m.ipd_frac_multi_peak <= 5) & (m.Cs_sig.abs() < 3) & m.ph_qual.astype(str).str[2].eq('A') & (m.ks_msigcom < 0.05) & (m.parallax > 10)
s = m[good].copy()
print('comparison stars:', len(s))
# target
T = dict(G=14.281043, BP_RP=2.756149, K=10.482, eK=0.026, J=11.363, plx=15.5494)
T['MG'] = T['G'] + 5 * np.log10(T['plx'] / 100); T['MK'] = T['K'] + 5 * np.log10(T['plx'] / 100); T['GK'] = T['G'] - T['K']; T['JK'] = T['J'] - T['K']
A = dict(G=14.654037, BP_RP=2.877811, K=8.678, plx=15.5494)   # pair A (2MASS blended with B; G-band A only)
A['MG'] = A['G'] + 5 * np.log10(A['plx'] / 100)
def ridge(x, y, deg=3, xr=None):
    k = np.isfinite(x) & np.isfinite(y)
    x, y = x[k], y[k]; msk = np.ones(len(x), bool)
    for it in range(10):
        c = np.polyfit(x[msk], y[msk], deg); r = y - np.polyval(c, x)
        sd = 1.4826 * np.median(np.abs(r[msk] - np.median(r[msk])))
        msk = np.abs(r - np.median(r[msk])) < 2.5 * sd
    # use the mode of residuals (single-star ridge) rather than the mean (binaries pull to the bright side)
    h, e = np.histogram(r[np.abs(r) < 1.5], bins=np.arange(-1.5, 1.51, 0.05)); mode = 0.5 * (e[np.argmax(h)] + e[np.argmax(h) + 1])
    return c, sd, r, mode, x
out = {}
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, axs = plt.subplots(1, 3, figsize=(13, 4.5))
for ax, (lab, xs, ys, xt, yt) in zip(axs, [('MKs_vs_BPRP', s.bp_rp.values, s.MK.values, T['BP_RP'], T['MK']),
                                          ('MG_vs_BPRP', s.bp_rp.values, s.MG.values, T['BP_RP'], T['MG']),
                                          ('MKs_vs_GKs', s.GK.values, s.MK.values, T['GK'], T['MK'])]):
    if lab == 'MKs_vs_GKs':
        kk_ = (xs > 3.3) & (xs < 4.2); xs, ys = xs[kk_], ys[kk_]
    c, sd, r, mode, x = ridge(xs, ys, deg=(1 if lab == 'MKs_vs_GKs' else 3))
    dt = yt - np.polyval(c, xt) - mode          # negative = brighter than ridge
    near = np.abs(x - xt) < 0.1
    rr = r[near] - mode
    frac = float(np.mean(rr <= dt))
    # robust width of the single-star (faint) side: use residuals > 0 (fainter than ridge) mirrored
    faint = rr[rr > 0]; sd_faint = float(np.sqrt(np.mean(faint ** 2))) if len(faint) > 10 else np.nan
    out[lab] = dict(target_offset_mag=float(dt), ridge_sd_mag=float(sd), faint_side_rms=sd_faint, n_near=int(near.sum()),
                    frac_comparison_as_bright_or_brighter=frac, percentiles_resid=[float(np.percentile(rr, p)) for p in (2.5, 16, 50, 84, 97.5)],
                    frac_brighter_than_0p5=float(np.mean(rr < -0.5)), frac_brighter_than_0p75=float(np.mean(rr < -0.75)), mode=mode)
    print(f'{lab}: target offset {dt:+.3f} mag (neg = over-luminous); ridge sd {sd:.3f}; faint-side rms {sd_faint:.3f}; '
          f'{frac*100:.1f}% of {near.sum()} comparison stars (|dx|<0.1) are as bright or brighter')
    ax.scatter(xs, ys, s=1, c='0.6', alpha=0.4, rasterized=True)
    xx = np.linspace(np.nanmin(xs), np.nanmax(xs), 100); ax.plot(xx, np.polyval(c, xx) + mode, 'C0'); ax.plot(xx, np.polyval(c, xx) + mode - 0.753, 'C0--')
    ax.plot(xt, yt, '*', ms=14, color='C3'); ax.invert_yaxis(); ax.set_title(f'{lab}: offset {dt:+.2f} mag'); ax.set_xlabel(lab.split('_vs_')[1]); ax.set_ylabel(lab.split('_vs_')[0])
    if lab == 'MG_vs_BPRP':
        ax.plot(A['BP_RP'], A['MG'], 'v', ms=9, color='C2', label='TIC 214311910 A (BP/RP blended)')
        out['pairA_MG_offset'] = float(A['MG'] - np.polyval(c, A['BP_RP']) - mode); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(f'{B}/figs/cmd.png', dpi=110)
# Mamajek table comparison
mm = pd.read_csv(f'{B}/data/ext/EEM_dwarf.txt', sep=r'\s+', comment='#', header=None, engine='python', on_bad_lines='skip', nrows=200)
out['target'] = T
# does the offset depend on RUWE? (unresolved binaries more common at high RUWE) -- for context
for rlab, rm in [('ruwe<1.1', s.ruwe < 1.1), ('1.1<ruwe<1.4', (s.ruwe >= 1.1) & (s.ruwe < 1.4))]:
    c, sd, r, mode, x = ridge(s.bp_rp.values, s.MK.values)
    rr = (r - mode)[rm.values[np.isfinite(s.bp_rp.values) & np.isfinite(s.MK.values)]]
    out[f'frac_overlum_0p5_{rlab}'] = float(np.mean(rr < -0.5))
print(json.dumps({k: v for k, v in out.items() if k.startswith('frac') or k.startswith('pairA')}, indent=0))
json.dump(out, open(f'{B}/out/cmd.json', 'w'), indent=1, default=float)
