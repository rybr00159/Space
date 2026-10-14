"""Summary figures: (1) TESS localisation with the lateral-colour model, (2) Gaia IPD speckle calibration."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from loclib import R, SECTORS

BLUE, ORANGE, INK, MUTED, GRID = '#2a78d6', '#eb6834', '#2b2b2a', '#6f6e69', '#e4e3df'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED, 'ytick.color': MUTED,
                     'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.facecolor': '#fcfcfb', 'axes.facecolor': '#fcfcfb', 'legend.frameon': False})
comb = json.load(open(f'{R}/results/s05_combine.json'))
bs = json.load(open(f'{R}/results/x02_boresight.json'))
real = json.load(open(f'{R}/results/real_fast.json'))
chi = np.load(f'{R}/results/s06_chi2map.npz')
x12 = json.load(open(f'{R}/results/x12_red_end.json'))
x13 = json.load(open(f'{R}/results/x13_dwarf_giant.json'))
fig, ax = plt.subplots(2, 2, figsize=(11, 10))
# (a) raw sky offsets
a = ax[0, 0]
for t, col, mk in [('deep', BLUE, 'o'), ('shallow', ORANGE, 's')]:
    for i, s in enumerate(SECTORS):
        lab = f'S{s}_{t}'
        k = comb['labels'].index(lab)
        o = np.array(comb['obs'][k]); C = np.array(comb['cov'][k])
        a.errorbar(o[0], o[1], xerr=np.sqrt(C[0, 0]), yerr=np.sqrt(C[1, 1]), fmt=mk, color=col, ms=6, lw=1, capsize=0,
                   label=f'{t} eclipses' if i == 0 else None)
        a.annotate(f'S{s}', (o[0], o[1]), xytext=(4, 3), textcoords='offset points', fontsize=7, color=MUTED)
for s in SECTORS:
    pa = np.radians(bs[str(s)]['pa_to_boresight'])
    a.annotate('', xy=(-6.5 * np.sin(pa) * 0 + 1.2 * np.sin(pa), 1.2 * np.cos(pa)), xytext=(0, 0),
               arrowprops=dict(arrowstyle='->', color=MUTED, lw=0.8))
a.plot(0, 0, '*', color=INK, ms=12, label='G star (Gaia DR3)')
a.set_xlim(7, -7); a.set_ylim(-7, 7); a.set_aspect('equal')
a.set_xlabel('offset East (arcsec)'); a.set_ylabel('offset North (arcsec)')
a.set_title('(a) Measured difference-image offsets (bias-corrected)\ngrey arrows: direction to each sector\'s camera boresight', fontsize=9, color=INK)
a.legend(loc='lower left', fontsize=8)
# (b) boresight frame
a = ax[0, 1]
for t, col, mk, dx in [('deep', BLUE, 'o', -0.12), ('shallow', ORANGE, 's', 0.12)]:
    for i, s in enumerate(SECTORS):
        k = comb['labels'].index(f'S{s}_{t}')
        o = np.array(comb['obs'][k]); C = np.array(comb['cov'][k])
        pa = np.radians(bs[str(s)]['pa_to_boresight'])
        u = np.array([np.sin(pa), np.cos(pa)]); tn = np.array([u[1], -u[0]])
        rad = o @ u; tan = o @ tn
        a.errorbar(i + dx, rad, yerr=np.sqrt(u @ C @ u), fmt=mk, color=col, ms=6, lw=1, label=f'{t}: radial' if i == 0 else None)
        a.errorbar(i + dx, tan, yerr=np.sqrt(tn @ C @ tn), fmt=mk, mfc='none', color=col, ms=6, lw=1, label=f'{t}: tangential' if i == 0 else None)
a.axhline(0, color=MUTED, lw=0.8)
a.set_xticks(range(len(SECTORS))); a.set_xticklabels([f'S{s}\ncam{bs[str(s)]["cam"]}' for s in SECTORS])
a.set_ylabel('offset component (arcsec)\n(radial: + = towards boresight)')
a.set_title('(b) Same offsets in the field frame: all radially outward,\ntangential components ~0', fontsize=9, color=INK)
a.legend(fontsize=7, ncol=2, loc='lower center')
# (c) delta contours
a = ax[1, 0]
g = chi['g']; M = chi['M']
cs = a.contour(g, g, M, levels=[2.30, 6.18, 11.83, 28.74], colors=[BLUE, BLUE, BLUE, MUTED], linewidths=[1.6, 1.2, 0.9, 0.8],
               linestyles=['-', '--', ':', '-'])
a.clabel(cs, fmt={2.30: '1σ', 6.18: '2σ', 11.83: '3σ', 28.74: '5σ'}, fontsize=7)
d = comb['M2_delta']
a.plot(d[0], d[1], 'o', color=BLUE, ms=6, label=f'best fit ({d[0]:+.2f}, {d[1]:+.2f})"')
a.plot(0, 0, '*', color=INK, ms=12, label='G star')
for r in [1, 2, 3]:
    th = np.linspace(0, 2 * np.pi, 200); a.plot(r * np.cos(th), r * np.sin(th), color=GRID, lw=0.8)
a.set_xlim(3.2, -3.2); a.set_ylim(-3.2, 3.2); a.set_aspect('equal')
a.set_xlabel('source offset East (arcsec)'); a.set_ylabel('source offset North (arcsec)')
a.set_title('(c) Source position, 6 sectors jointly, colour levers profiled\n(contours: profile Δχ², injection-calibrated errors)', fontsize=9, color=INK)
a.legend(fontsize=8, loc='lower left')
# (d) lateral colour curve
a = ax[1, 1]
sc = 20.4
kref12 = x12['kref']
curve = x12['curve']
C = [c['C_med'] for c in curve]; y = [c['rad_mpx'] / 1e3 * sc for c in curve]; e = [c['err_mpx'] / 1e3 * sc for c in curve]
a.errorbar(C, y, yerr=e, fmt='o', color=INK, ms=4, lw=1, label='field stars (all), median per bin')
dw = [q for q in x13 if q['cls'] == 'dwarf' and q['n_unique'] >= 1 and q['lo'] >= 2.0]
if dw:
    a.errorbar([q['C_med'] for q in dw], [q['rad_mpx'] / 1e3 * sc for q in dw], yerr=[q['err_mpx'] / 1e3 * sc for q in dw], fmt='D',
               mfc='none', color=MUTED, ms=6, label='M dwarfs (2 stars)')
xx = np.linspace(0.3, 4.2, 10); a.plot(xx, kref12 * (xx - 0.786) * sc, color=MUTED, lw=0.8, ls='--', label='linear fit, BP-RP<1.7')
kref11 = comb['kref_mpx_per_mag']
for t, col, cval, csd in [('deep', BLUE, comb['models']['M2']['p'][2], comb['M2_c_emp_sd'][0]),
                          ('shallow', ORANGE, comb['models']['M2']['p'][3], comb['M2_c_emp_sd'][1])]:
    v = cval * kref12 / kref11
    a.axhspan(v - csd, v + csd, color=col, alpha=0.18, lw=0)
    a.axhline(v, color=col, lw=1.4, label=f'{t}-eclipse light (fit): {v:.1f}"')
a.axvline(0.786, color=MUTED, lw=0.6); a.text(0.82, 0.5, 'G star', fontsize=7, color=MUTED)
a.set_xlabel('Gaia BP-RP'); a.set_ylabel('radial PRF-fit displacement (arcsec)\nat the target field angle (+ = towards boresight)')
a.set_title('(d) TESS lateral colour measured with field stars (12 FFIs)\nvs the radial shift fitted to the eclipse light', fontsize=9, color=INK)
a.legend(fontsize=7, loc='lower left')
plt.tight_layout()
plt.savefig(f'{R}/plots/fig1_tess_localisation.png', dpi=130)
plt.close()

# ---------------- Figure 2: Gaia calibration
S = pd.read_csv(f'{R}/speckle/calib_sample_v06.csv')
B = S[(S.kind == 'binary') & (S.Gmag > 9.5) & (S.Gmag < 13.5) & ~S.gaia_resolved.astype(bool)]
Cn = S[(S.kind == 'single') & (S.Gmag > 9.5) & (S.Gmag < 13.5) & (S.dm1 >= 3.5)]
v07 = json.load(open(f'{R}/results/v07_gaia_likelihood.json'))
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
a = ax[0]
ramp = ['#9ec5f4', '#5598e7', '#2a78d6', '#1c5cab', '#0d366b']
bins = [(0, 1.5), (1.5, 2.5), (2.5, 3.5), (3.5, 4.5), (4.5, 8)]
for (lo, hi), col in zip(bins, ramp):
    s = B[(B.dm >= lo) & (B.dm < hi)]
    a.scatter(s.sep, s.IPDfmp, s=9, color=col, label=f'Δm_I {lo}-{hi}', lw=0)
m = B[(B.IPDfmp >= 28) & (B.IPDfmp <= 56) & (B.RUWE < 1.4) & (B.bfrac >= 0.8)]
a.scatter(m.sep, m.IPDfmp, s=60, facecolors='none', edgecolors=ORANGE, lw=1.4, label='matches target signature')
a.axhline(41, color=INK, lw=0.8, ls='--'); a.text(2.2, 43, 'target: 41', fontsize=8, color=INK)
a.set_xscale('log'); a.set_xlim(0.03, 3); a.set_xlabel('speckle separation (arcsec)'); a.set_ylabel('Gaia ipd_frac_multi_peak (%)')
a.set_title(f'(a) SOAR speckle binaries, G 9.5-13.5, Gaia-unresolved (n={len(B)})', fontsize=9, color=INK)
a.legend(fontsize=7, loc='upper left')
a = ax[1]
bb = np.arange(0, 102, 5)
a.hist(Cn.IPDfmp, bins=bb, color=MUTED, alpha=0.9, label=f'speckle-unresolved stars (n={len(Cn)})', log=True)
s = B[(B.sep > 0.5) & (B.sep < 1.3)]
a.hist(s.IPDfmp, bins=bb, histtype='step', color=BLUE, lw=1.6, label=f'binaries 0.5-1.3" (n={len(s)})', log=True)
a.axvline(41, color=INK, lw=0.8, ls='--')
a.set_xlabel('ipd_frac_multi_peak (%)'); a.set_ylabel('number of stars')
a.set_title('(b) Single stars almost never reach 41', fontsize=9, color=INK)
a.legend(fontsize=7)
a = ax[2]
grid = v07['fmp+ruwe+blend']['grid']
seps = sorted({tuple(gq['sep']) for gq in grid}); dms = sorted({tuple(gq['dm']) for gq in grid})
Z = np.full((len(dms), len(seps)), np.nan)
for gq in grid:
    if gq['n']:
        Z[dms.index(tuple(gq['dm'])), seps.index(tuple(gq['sep']))] = gq['frac']
from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list('b', ['#f3f7fd', '#9ec5f4', '#2a78d6', '#0d366b'])
im = a.imshow(Z, origin='lower', aspect='auto', cmap=cmap, vmin=0, vmax=0.4)
for gq in grid:
    if gq['n']:
        a.text(seps.index(tuple(gq['sep'])), dms.index(tuple(gq['dm'])), f"{gq['k']}/{gq['n']}", ha='center', va='center', fontsize=7,
               color='white' if (gq['frac'] or 0) > 0.2 else INK)
a.set_xticks(range(len(seps))); a.set_xticklabels([f'{s0:.2f}-{s1:.2f}' for s0, s1 in seps], rotation=45, fontsize=7)
a.set_yticks(range(len(dms))); a.set_yticklabels([f'{d0}-{d1}' for d0, d1 in dms])
a.set_xlabel('speckle separation (arcsec)'); a.set_ylabel('Δm_I (mag)')
sr = v07['fmp+ruwe+blend']
a.set_title(f'(c) Fraction matching IPDfmp 28-56 & RUWE<1.4 & BP-blend>=0.8\nsingles: {sr["singles_k"]}/{sr["singles_n"]}', fontsize=9, color=INK)
plt.colorbar(im, ax=a, fraction=0.04, label='fraction matching')
plt.tight_layout()
plt.savefig(f'{R}/plots/fig2_gaia_calibration.png', dpi=130)
print('ok')
