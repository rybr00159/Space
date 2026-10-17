"""f08: figures. Palette: reference categorical slots (blue = planet on the M3V system, orange = false positive,
aqua = resolved-neighbour false positive), text in neutral ink, light surface.
figs/scenario_posteriors.png, figs/fpp_distribution.png, figs/variants.png, figs/ptp_map.png, figs/tri_runs.png"""
import numpy as np, json, sys, pandas as pd, os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, 'tic214311908/analysis/fpp/scripts')
from fpplib import *
BLUE, ORANGE, AQUA = '#2a78d6', '#eb6834', '#1baf7a'
SURF, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
plt.rcParams.update({'figure.facecolor': SURF, 'axes.facecolor': SURF, 'axes.edgecolor': INK2, 'axes.labelcolor': INK,
                     'xtick.color': INK2, 'ytick.color': INK2, 'text.color': INK, 'font.size': 9, 'axes.spines.top': False,
                     'axes.spines.right': False, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6})
base = jload(f'{F}/out/fpp_baseline.json')['result']
# 1. scenario posteriors
post = base['post']
names = sorted(post, key=lambda k: -post[k])
vals = np.array([max(post[k], 1e-20) for k in names])
cols = [BLUE if k.startswith('TP') else (AQUA if k.startswith('N') else ORANGE) for k in names]
fig, ax = plt.subplots(figsize=(7.5, 5.2))
ax.barh(range(len(names)), vals, color=cols, height=0.62, left=1e-20)
ax.set_xscale('log'); ax.set_xlim(1e-17, 2)
ax.set_yticks(range(len(names))); ax.set_yticklabels([n.replace('_', ' ') for n in names]); ax.invert_yaxis()
for i, v in enumerate(vals):
    if v > 1e-17: ax.text(v*1.3, i, f'{v:.1e}', va='center', fontsize=7.5, color=INK2)
ax.set_xlabel('posterior probability (baseline priors, INV light curve)')
ax.set_title(f"Scenario posteriors - FPP {base['FPP']:.4f}, NFPP {base['NFPP']:.1e}  (Mira: physically excluded, 0)", fontsize=9.5)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=BLUE, label='planet on the M3V or its twin'), Patch(color=ORANGE, label='false positive (bound / background)'),
                   Patch(color=AQUA, label='false positive on a resolved neighbour')], loc='lower right', frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(f'{F}/figs/scenario_posteriors.png', dpi=150); plt.close(fig)
# 2. FPP / NFPP distributions over the random prior draws
R = pd.read_csv(f'{F}/out/fpp_random.csv')
tri = pd.read_csv(f'{F}/out/tri_runs.csv') if os.path.exists(f'{F}/out/tri_runs.csv') else None
fig, axs = plt.subplots(1, 2, figsize=(9, 3.4))
for ax, col, thr, lab in ((axs[0], 'FPP', 0.015, 'FPP'), (axs[1], 'NFPP', 1e-3, 'NFPP')):
    x = np.log10(np.clip(R[col].values, 1e-9, None))
    ax.hist(x, bins=50, color=BLUE, edgecolor=SURF, linewidth=0.6)
    ax.axvline(np.log10(thr), color=INK, lw=1.2, ls='--'); ax.text(np.log10(thr), ax.get_ylim()[1]*0.93, f' threshold {thr:g}', fontsize=8)
    ax.axvline(np.log10(base[col]), color=ORANGE, lw=2); ax.text(np.log10(base[col]), ax.get_ylim()[1]*0.80, ' baseline', color=INK2, fontsize=8)
    ax.set_xlabel(f'log10 {lab}'); ax.set_ylabel('prior draws')
    f_ok = np.mean(R[col] < thr)
    ax.set_title(f'{lab} over {len(R)} joint prior draws: {100*f_ok:.0f}% below {thr:g}', fontsize=9)
    if tri is not None and col == 'FPP':
        a = tri[tri.group == 'A'][col].values
        ax.plot(np.log10(a), np.full(len(a), ax.get_ylim()[1]*0.05), '|', color=INK, ms=12, mew=1.5)
        ax.text(-4.4, ax.get_ylim()[1]*0.05, 'TRICERATOPS main runs ->', fontsize=7.5, color=INK2, va='center')
fig.tight_layout(); fig.savefig(f'{F}/figs/fpp_distribution.png', dpi=150); plt.close(fig)
# 3. named variants
V = pd.read_csv(f'{F}/out/fpp_variants.csv')
fig, ax = plt.subplots(figsize=(8, 15))
y = np.arange(len(V))
ax.plot(V.FPP, y, 'o', color=BLUE, ms=4.5, label='FPP (task definition)')
ax.plot(V.FPP_system, y, 's', color=AQUA, ms=3.5, label='FPP_system (any planet in the bound system)')
ax.plot(V.NFPP, y, 'd', color=ORANGE, ms=3.5, label='NFPP')
ax.axvline(0.015, color=INK, ls='--', lw=1); ax.axvline(0.001, color=INK2, ls=':', lw=1)
ax.set_xscale('log'); ax.set_xlim(1e-8, 1)
ax.set_yticks(y); ax.set_yticklabels(V.variant, fontsize=6.5); ax.invert_yaxis()
ax.set_xlabel('probability'); ax.legend(loc='lower left', frameon=False, fontsize=8)
ax.set_title('Named sensitivity variants (dashed: FPP 0.015; dotted: NFPP 0.001)', fontsize=9.5)
fig.tight_layout(); fig.savefig(f'{F}/figs/variants.png', dpi=150); plt.close(fig)
# 4. PTP posterior map (where the residual planet-on-companion probability sits) with the constraints
cz = np.load(f'{F}/out/companions_mc.npz')
# recompute baseline PTP weights cheaply from saved terms is complex; use f04's stored arrays if present
if os.path.exists(f'{F}/out/ptp_weights.npz'):
    w = np.load(f'{F}/out/ptp_weights.npz')
    fig, ax = plt.subplots(figsize=(6.5, 4.4))
    sep, dI, wp = cz['sep'], cz['dI'], w['ptp']
    m = (sep > 0.005) & (sep < 5)
    H, xe, ye = np.histogram2d(np.log10(sep[m]), dI[m], bins=[np.linspace(-2.3, 0.7, 46), np.linspace(0, 5, 41)], weights=wp[m])
    H = H/H.sum()
    im = ax.pcolormesh(xe, ye, H.T, cmap='Blues', shading='flat')
    plt.colorbar(im, ax=ax, label='fraction of the PTP posterior')
    ipd = jload('star/out/ipd_limits.json')['cells']
    for c in ipd:
        if c['n'] >= 5 and c['LR'] == c['LR'] and c['LR'] < 0.2:
            ax.add_patch(plt.Rectangle((np.log10(c['sep'][0]), c['dm'][0]), np.log10(c['sep'][1]) - np.log10(c['sep'][0]), c['dm'][1] - c['dm'][0],
                                       fill=False, ec=INK2, lw=0.5, hatch='//', alpha=0.5))
    ds_ = np.array([0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0]); dv_ = np.array([2.0, 4.24, 5.07, 5.80, 5.58, 6.17, 6.89])
    ax.plot(np.log10(ds_), dv_, color=ORANGE, lw=2, label='DECaPS2 5-sigma (referee-corrected)')
    hs = np.array([0.03, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0]); hv = np.array([2.0, 3.0, 4.0, 4.5, 5.0, 5.5, 6.0])
    ax.plot(np.log10(hs), hv, color=INK, lw=1.5, ls='--', label='hypothetical 8-m speckle 5-sigma')
    ax.axvline(np.log10(0.08), color=INK2, lw=0.8, ls=':'); ax.text(np.log10(0.08), 0.25, ' IPD calibration starts (0.08")', fontsize=7, color=INK2)
    ax.set_xlabel('log10 projected separation (arcsec)'); ax.set_ylabel('companion contrast Delta I (mag)'); ax.set_ylim(5, 0)
    ax.legend(loc='lower left', frameon=False, fontsize=7.5)
    ax.set_title('Where the residual "planet on a bound companion" (PTP) probability lies\nhatched: Gaia IPD cells with LR < 0.2', fontsize=9)
    fig.tight_layout(); fig.savefig(f'{F}/figs/ptp_map.png', dpi=150); plt.close(fig)
# 5. TRICERATOPS runs
if tri is not None:
    order = [g for g in ['A', 'F', 'C', 'G', 'D', 'E1', 'E2', 'E3', 'E4', 'E5'] if g in set(tri.group)]
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for ax, col, thr in ((axs[0], 'FPP', 0.015), (axs[1], 'NFPP', 1e-3)):
        for i, g in enumerate(order):
            v = tri[tri.group == g][col].values
            ax.plot(np.clip(v, 1e-12, None), np.full(len(v), i) + np.linspace(-0.15, 0.15, len(v)), 'o', color=BLUE, ms=4, alpha=0.85)
            if col == 'FPP':
                v2 = tri[tri.group == g]['FPP_no_x2P'].values
                ax.plot(v2, np.full(len(v2), i) + 0.3, '|', color=ORANGE, ms=6)
        ax.axvline(thr, color=INK, ls='--', lw=1); ax.set_xscale('log'); ax.set_xlabel(col)
        ax.set_yticks(range(len(order))); ax.set_yticklabels(order); ax.invert_yaxis()
    axs[0].set_title('TRICERATOPS FPP per run (dots); orange ticks: x2P scenarios removed', fontsize=9)
    axs[1].set_title('TRICERATOPS NFPP per run (D, E5: all TIC neighbours kept)', fontsize=9)
    fig.tight_layout(); fig.savefig(f'{F}/figs/tri_runs.png', dpi=150); plt.close(fig)
print('figures written')
