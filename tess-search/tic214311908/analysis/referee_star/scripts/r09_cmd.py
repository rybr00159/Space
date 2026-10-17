"""Referee: is TIC 214311908 over-luminous for its colour, and are its colours internally consistent?
Sample: data/mdw_cmd.csv (own Gaia archive query: plx > 12 mas, plx/err > 20, BP-RP 2.55-2.95, 2MASS best neighbour).
Cuts: RUWE < 1.4, 2MASS ph_qual AAA, ks_msigcom < 0.05, |C*| < 3 sigma (Riello+21).
Ridge = KDE mode of the absolute magnitude in sliding colour bins (width 0.04), smoothed with a quadratic fit;
target offset relative to the ridge; fraction of comparison stars within +-0.03 in colour that are at least as bright
(relative to the ridge).  Also colour-colour residuals (G-Ks, G-RP, J-Ks vs BP-RP) to test whether BP-RP or Ks
of the target is anomalous (e.g. contamination by the Mira 17.8" away).
Output: out/r09_cmd.json, figs/r09_cmd.png"""
import json, numpy as np, pandas as pd
from scipy.stats import gaussian_kde
R = 'tic214311908/analysis/referee_star'
m = pd.read_csv(f'{R}/data/mdw_cmd.csv')
x = m.bp_rp
f = 1.162004 + 0.011464 * x + 0.049255 * x ** 2 - 0.005879 * x ** 3
m['Cs'] = (m.phot_bp_rp_excess_factor - f) / (0.0059898 + 8.817481e-12 * m.phot_g_mean_mag ** 7.618399)
m['MG'] = m.phot_g_mean_mag + 5 * np.log10(m.parallax / 100); m['MK'] = m.ks_m + 5 * np.log10(m.parallax / 100)
m['GK'] = m.phot_g_mean_mag - m.ks_m; m['GRP'] = m.phot_g_mean_mag - m.phot_rp_mean_mag; m['JK'] = m.j_m - m.ks_m
good = (m.ruwe < 1.4) & (m.ph_qual == 'AAA') & (m.ks_msigcom < 0.05) & (m.Cs.abs() < 3)
s = m[good].copy(); print('comparison stars', len(s), 'of', len(m))
plx = 15.4973 + 0.0521
T = dict(bp_rp=2.756149, G=14.281043, RP=13.074592, K=10.482, J=11.363)
T['MG'] = T['G'] + 5 * np.log10(plx / 100); T['MK'] = T['K'] + 5 * np.log10(plx / 100); T['GK'] = T['G'] - T['K']; T['GRP'] = T['G'] - T['RP']; T['JK'] = T['J'] - T['K']
def mode_ridge(xc, y, x0s, w=0.04):
    pts = []
    for x0 in x0s:
        k = np.abs(xc - x0) < w / 2
        if k.sum() < 50: continue
        yy = y[k]; kde = gaussian_kde(yy, bw_method=0.25); grid = np.linspace(np.percentile(yy, 1), np.percentile(yy, 99), 400)
        pts.append((x0, grid[np.argmax(kde(grid))]))
    pts = np.array(pts); c = np.polyfit(pts[:, 0], pts[:, 1], 2); return c, pts
out = {}
x0s = np.arange(2.60, 2.91, 0.01)
for lab, ycol in (('MKs', 'MK'), ('MG', 'MG')):
    c, pts = mode_ridge(s.bp_rp.values, s[ycol].values, x0s)
    res = s[ycol].values - np.polyval(c, s.bp_rp.values)
    dt = T[ycol] - np.polyval(c, T['bp_rp'])
    near = np.abs(s.bp_rp.values - T['bp_rp']) < 0.03
    frac = float(np.mean(res[near] <= dt))
    # robust width of the faint (single) side
    rr = res[near]; faint = rr[rr > 0]; wid = float(np.sqrt(np.mean(faint ** 2)))
    out[f'{lab}_vs_BPRP'] = dict(offset=float(dt), frac_as_bright=frac, n_near=int(near.sum()), faint_side_rms=wid,
                                 frac_brighter_than_0p5=float(np.mean(rr < -0.5)), frac_in_twin_band=float(np.mean((rr < -0.5) & (rr > -1.0))))
    print(f'{lab} vs BP-RP: target offset {dt:+.3f} mag; {frac*100:.1f}% of {near.sum()} stars (|dBP-RP|<0.03) as bright; faint-side rms {wid:.3f}; '
          f'frac < -0.5 mag: {np.mean(rr < -0.5):.3f}')
# colour-colour consistency at the target's BP-RP
for lab in ('GK', 'GRP', 'JK'):
    near = np.abs(s.bp_rp - T['bp_rp']) < 0.03
    cc = np.polyfit(s.bp_rp[np.abs(s.bp_rp - T['bp_rp']) < 0.15], s[lab][np.abs(s.bp_rp - T['bp_rp']) < 0.15], 1)
    r_ = s[lab] - np.polyval(cc, s.bp_rp); sd = 1.4826 * np.median(np.abs(r_[near] - np.median(r_[near])))
    dt = T[lab] - np.polyval(cc, T['bp_rp']) - np.median(r_[near])
    out[f'{lab}_resid'] = dict(target_resid=float(dt), robust_sd=float(sd), z=float(dt / sd), pct=float(np.mean(r_[near] - np.median(r_[near]) <= dt)))
    print(f'{lab} at fixed BP-RP: target residual {dt:+.3f} mag (robust sd {sd:.3f}; z {dt/sd:+.2f}; percentile {out[lab+"_resid"]["pct"]*100:.1f})')
# Do over-luminous stars (offset < -0.5) show a G-Ks residual like the target's? (binaries vs metal-rich)
c, _ = mode_ridge(s.bp_rp.values, s.MK.values, x0s); s['dMK'] = s.MK - np.polyval(c, s.bp_rp)
cc = np.polyfit(s.bp_rp, s.GK, 2); s['rGK'] = s.GK - np.polyval(cc, s.bp_rp)
for lo, hi in ((-0.2, 0.2), (-0.5, -0.2), (-0.9, -0.55), (-1.2, -0.9)):
    k = (s.dMK > lo) & (s.dMK < hi) & (np.abs(s.bp_rp - 2.756) < 0.1)
    print(f'  dMKs in ({lo},{hi}): n={k.sum()}, median G-Ks residual {np.median(s.rGK[k]):+.3f}, RUWE median {np.median(s.ruwe[k]):.3f}')
    out[f'GKres_dMK_{lo}_{hi}'] = dict(n=int(k.sum()), med_rGK=float(np.median(s.rGK[k])), med_ruwe=float(np.median(s.ruwe[k])))
out['target_rGK'] = float(T['GK'] - np.polyval(cc, T['bp_rp']))
print('target G-Ks residual (quadratic over full sample)', round(out['target_rGK'], 3))
# fraction of 'as over-luminous' stars with low RUWE like the target (RUWE<1.1)
k = (np.abs(s.bp_rp - 2.756) < 0.05)
for rl, rh in ((0, 1.1), (1.1, 1.4)):
    kk = k & (s.ruwe >= rl) & (s.ruwe < rh)
    print(f'  RUWE {rl}-{rh}: n={kk.sum()}, frac dMKs <= {out["MKs_vs_BPRP"]["offset"]:.2f}: {np.mean(s.dMK[kk] <= out["MKs_vs_BPRP"]["offset"]):.3f}')
    out[f'frac_overlum_ruwe_{rl}_{rh}'] = float(np.mean(s.dMK[kk] <= out['MKs_vs_BPRP']['offset']))
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(10, 4.5))
ax[0].scatter(s.bp_rp, s.MK, s=1, c='0.6', alpha=0.4, rasterized=True); xx = np.linspace(2.55, 2.95, 50)
ax[0].plot(xx, np.polyval(c, xx), 'C0'); ax[0].plot(xx, np.polyval(c, xx) - 0.753, 'C0--'); ax[0].plot(T['bp_rp'], T['MK'], '*', ms=14, c='C3')
ax[0].invert_yaxis(); ax[0].set_xlabel('BP-RP'); ax[0].set_ylabel('M_Ks')
ax[1].scatter(s.bp_rp, s.GK, s=1, c='0.6', alpha=0.4, rasterized=True); ax[1].plot(T['bp_rp'], T['GK'], '*', ms=14, c='C3'); ax[1].set_xlabel('BP-RP'); ax[1].set_ylabel('G-Ks')
fig.tight_layout(); import os; os.makedirs(f'{R}/figs', exist_ok=True); fig.savefig(f'{R}/figs/r09_cmd.png', dpi=100)
json.dump(out, open(f'{R}/out/r09_cmd.json', 'w'), indent=1)
