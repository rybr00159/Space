"""v4 summary figure: (A) per-sector bias-corrected offsets with injection 1-sigma ellipses and the away-from-boresight
directions, (B) calibrated confidence regions of delta (chromatic model M2) with the VHS feature direction and the Gaia
IPD-favoured annulus, (C) field-star lateral-colour curve with the fitted colour levers of the two eclipse types."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from v4lib import V4
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Wedge

LAB = os.environ.get('F10_LABEL', 'baseline')
C = json.load(open(f'{V4}/results/s05_{LAB}.json'))
Mz = np.load(f'{V4}/results/s05_{LAB}_map.npz')
F3 = json.load(open(f'{V4}/results/f03_theta_red.json'))
INK, INK2, GRID, SURF = '#0b0b0b', '#52514e', '#d9d8d4', '#fcfcfb'
COL = {'deep': '#2a78d6', 'shallow': '#eb6834'}
plt.rcParams.update({'font.size': 10, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
                     'figure.facecolor': SURF, 'axes.facecolor': SURF})
fig, axs = plt.subplots(1, 3, figsize=(17, 5.8))
# (A)
ax = axs[0]
obs = np.array(C['obs_bias_corrected']); cov = np.array(C['cov']); U = np.array(C['u_to_boresight'])
for o, c, l, u in zip(obs, cov, C['labels'], U):
    t = l.split('_')[1]
    w, v = np.linalg.eigh(c)
    ang = np.degrees(np.arctan2(v[1, 1], v[0, 1]))
    ax.add_patch(Ellipse((o[0], o[1]), 2 * np.sqrt(w[1]), 2 * np.sqrt(w[0]), angle=ang, fc='none', ec=COL[t], lw=1.2, alpha=0.8))
    ax.plot(o[0], o[1], 'o', color=COL[t], ms=7, mec=SURF, mew=1.5)
    if t == 'deep':
        ax.annotate(l.split('_')[0], (o[0], o[1]), xytext=(6, 4), textcoords='offset points', color=INK2, fontsize=8)
        ax.plot([0, -6.5 * u[0]], [0, -6.5 * u[1]], color=GRID, lw=1, zorder=0)
ax.plot(0, 0, '*', color=INK, ms=14)
ax.plot([], [], 'o', color=COL['deep'], label='deep eclipse'); ax.plot([], [], 'o', color=COL['shallow'], label='shallow eclipse')
ax.plot([], [], '-', color=GRID, label='away from boresight')
ax.set_xlim(7.5, -7.5); ax.set_ylim(-7.5, 7.5); ax.set_aspect('equal'); ax.grid(color=GRID, lw=0.5)
ax.set_xlabel('dE (arcsec)  [east to the left]'); ax.set_ylabel('dN (arcsec)')
ax.set_title('A. Per-sector offsets of the eclipse PSF from the G star', fontsize=10, color=INK, loc='left')
ax.legend(loc='lower left', fontsize=8, frameon=False)
# (B)
ax = axs[1]
g = Mz['g']; M = Mz['M']
thr = [C['thresholds'][k] for k in ('68.27', '95.45', '99.73')]
ax.contourf(g, g, M, levels=[-1, thr[0], thr[1], thr[2]], colors=['#2a78d6', '#7fb0ea', '#d3e4f8'], alpha=0.9)
ax.contour(g, g, M, levels=thr, colors=[INK2], linewidths=0.6)
ax.add_patch(Wedge((0, 0), 1.2, 0, 360, width=0.7, fc='none', ec='#1baf7a', lw=1.5, ls='--'))
pa = np.radians(345)
ax.plot([0, 1.1 * np.sin(pa)], [0, 1.1 * np.cos(pa)], color='#1baf7a', lw=2.5)
ax.plot(*C['delta'], 'x', color=INK, ms=9, mew=2)
ax.plot(0, 0, '*', color=INK, ms=14)
ax.text(1.1 * np.sin(pa) - 0.05, 1.1 * np.cos(pa) + 0.12, 'VHS residual\nPA ~345, <=1.1"', color=INK2, fontsize=8, ha='left')
ax.text(-1.25, -1.35, 'Gaia IPD-favoured\n0.5-1.2" (dashed)', color=INK2, fontsize=8, ha='right')
ax.set_xlim(2.0, -2.0); ax.set_ylim(-2.0, 2.0); ax.set_aspect('equal'); ax.grid(color=GRID, lw=0.5)
ax.set_xlabel('dE (arcsec)'); ax.set_ylabel('dN (arcsec)')
ax.set_title('B. Eclipsing source position (68/95/99.7% calibrated)', fontsize=10, color=INK, loc='left')
# (C)
ax = axs[2]
cv = F3['curve']
Cc = np.array([c['C'] for c in cv]); R = -np.array([c['rad'] for c in cv]); E = np.array([c['err'] for c in cv])
ax.errorbar(Cc, R, yerr=E, fmt='s', color=INK2, ms=5, capsize=0, label='field stars, binned')
rd = F3['reds']
ax.errorbar([r['C'] for r in rd], [-r['rad'] for r in rd], yerr=[r['err'] for r in rd], fmt='o', color=GRID, mec=INK2, ms=5, capsize=0, label='individual stars, BP-RP>2')
xx = np.linspace(0.3, 4.2, 10); ax.plot(xx, -F3['k_mean'] * (xx - 0.786), color=INK2, lw=1, ls=':', label='linear fit (C<1.75)')
names = C['types']
for i, t in enumerate(names):
    c, s = C['c'][i], C['c_sd_emp'][i]
    ax.axhspan(c - s, c + s, color=COL[t], alpha=0.25, lw=0)
    ax.axhline(c, color=COL[t], lw=2, label=f'{t}-eclipse lever {c:.1f} +- {s:.1f}"')
ax.axvline(0.786, color=INK, lw=0.8); ax.text(0.81, 0.3, 'G star', color=INK2, fontsize=8)
ax.set_xlim(0.3, 4.2); ax.set_ylim(-1, 7.5); ax.grid(color=GRID, lw=0.5)
ax.set_xlabel('Gaia BP-RP'); ax.set_ylabel('PRF-fit shift away from boresight rel. to G star (arcsec)')
ax.set_title('C. TESS lateral colour at the target (12 FFIs)', fontsize=10, color=INK, loc='left')
ax.legend(fontsize=7.5, frameon=False, loc='upper left')
plt.tight_layout()
plt.savefig(f'{V4}/plots/fig_localisation_summary_v4.png', dpi=110)
print('saved')
