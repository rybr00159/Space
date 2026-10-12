"""Figures for the referee audit (static PNG; palette: target-host = blue #2a78d6, neighbour-host = orange #eb6834,
observed = ink #0b0b0b; sequential maps in one hue)."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

OUT = '/home/user/work/tic87/referee_audit'
BLUE, ORANGE, INK, MUTED, GRID = '#2a78d6', '#eb6834', '#0b0b0b', '#52514e', '#d9d8d4'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.spines.top': False, 'axes.spines.right': False})
SECT = ['13', '39', '66', '93', '104']

# ---------------- Figure 1: pixel test
r = json.load(open(f'{OUT}/r3_pixel_ratio_none_plane.json'))
fig, axs = plt.subplots(2, 5, figsize=(15, 7.0), gridspec_kw=dict(height_ratios=[1.15, 1]))
fig.subplots_adjust(right=0.9, hspace=0.35)
for k, s in enumerate(SECT):
    rs = r[s]; ax = axs[0, k]
    delta = np.array(rs['maps']['delta'], float) * 1e3
    xt, yt = rs['target_px']; xn, yn = rs['neighbour_px']
    ix, iy = int(round(xt)), int(round(yt))
    sub = delta[iy - 3:iy + 4, ix - 3:ix + 4]
    im = ax.imshow(sub, origin='lower', cmap='Blues', vmin=0, vmax=6, extent=(ix - 3.5, ix + 3.5, iy - 3.5, iy + 3.5))
    for (gx, gy), col, lab in ((rs['groups']['a'][0], ORANGE, 'a'), (rs['groups']['b'][0], BLUE, 'b'),
                               (rs['groups']['t'][0], MUTED, 't')):
        ax.add_patch(Rectangle((gx - 0.5, gy - 0.5), 1, 1, fill=False, ec=col, lw=2))
        ax.text(gx - 0.42, gy + 0.22, lab, color=col, fontsize=9, fontweight='bold')
    ax.plot(xt, yt, 'x', color=BLUE, ms=9, mew=2.2)
    ax.plot(xn, yn, '+', color=ORANGE, ms=11, mew=2.2)
    for yy_ in range(iy - 3, iy + 4):
        for xx_ in range(ix - 3, ix + 4):
            v = delta[yy_, xx_]
            if np.isfinite(v) and abs(xx_ - ix) <= 2 and abs(yy_ - iy) <= 2:
                ax.text(xx_ + 0.42, yy_ - 0.42, f'{v:.1f}', fontsize=6, ha='right', va='bottom',
                        color='white' if v > 3.5 else INK)
    ax.set_title(f'S{s}: relative depth D/O (ppt)', fontsize=9, color=INK)
    ax.set_xlabel('column (px)'); ax.set_ylabel('row (px)') if k == 0 else None
cax = fig.add_axes([0.915, 0.56, 0.008, 0.3]); cb = fig.colorbar(im, cax=cax); cb.set_label('ppt')
pairs = ['a/b', 'a/t', 'Ga/Gb', 'n/t']
for k, s in enumerate(SECT):
    ax = axs[1, k]; rs = r[s]
    for j, p in enumerate(pairs):
        if p not in rs['pairs']:
            continue
        pr = rs['pairs'][p]
        qt = pr['model_context']['scene']['Q_T_range']; qn = pr['model_context']['scene']['Q_N_range']
        ax.plot([j - 0.18] * 2, [qt[0], qt[2]], color=BLUE, lw=5, solid_capstyle='round', alpha=0.9)
        ax.plot([j + 0.18] * 2, [qn[0], qn[2]], color=ORANGE, lw=5, solid_capstyle='round', alpha=0.9)
        lo, mid, hi = pr['Q_null_noise_16_50_84']
        q = pr['Q']
        lo_ = max(lo, 0.05); hi_ = min(hi, 60)
        ax.plot([j, j], [lo_, hi_], color=INK, lw=1.2)
        ax.plot(j, min(max(q, 0.06), 60), 'D', color=INK, ms=6)
    ax.set_yscale('log'); ax.set_ylim(0.05, 70)
    ax.set_xticks(range(len(pairs))); ax.set_xticklabels(pairs)
    ax.axhline(1, color=GRID, lw=1, zorder=0)
    ax.set_title(f'S{s}: Q = delta(x)/delta(y)', fontsize=9, color=INK)
    if k == 0:
        ax.set_ylabel('depth ratio Q (log)')
axs[1, 0].plot([], [], color=BLUE, lw=5, label='predicted if target hosts (range)')
axs[1, 0].plot([], [], color=ORANGE, lw=5, label='predicted if neighbour hosts (range)')
axs[1, 0].plot([], [], 'D-', color=INK, ms=5, lw=1.2, label='observed, 16-84 % noise interval')
h_, l_ = axs[1, 0].get_legend_handles_labels()
fig.legend(h_, l_, loc='lower center', ncol=3, fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.03))
fig.suptitle('Model-free pixel test: relative transit depth on the neighbour side (a) vs far side (b) / target pixel (t). '
             'x = target, + = neighbour (Gaia, epoch-propagated)', fontsize=10, color=INK)
fig.savefig(f'{OUT}/fig_pixel_test.png', dpi=130, bbox_inches='tight'); plt.close(fig)

# ---------------- Figure 2: centroid.py injections
b = json.load(open(f'{OUT}/b1_centroid_audit.json'))
fig, axs = plt.subplots(2, 5, figsize=(15, 5.6), sharey='row')
for k, s in enumerate(SECT):
    bs = b[s]
    for row, key in ((0, 's_par'), (1, 'offset_sigma')):
        ax = axs[row, k]
        for j, (nm, col) in enumerate((('prf_target', BLUE), ('gauss0.6_target', BLUE),
                                       ('prf_neighbour', ORANGE), ('gauss0.6_neighbour', ORANGE))):
            vals = [o[key] for o in bs['injections'][nm] if o]
            xj = j + np.linspace(-0.15, 0.15, len(vals))
            ax.plot(xj, vals, 'o', color=col, ms=4.5, alpha=0.85, mec='white', mew=0.6)
        if bs['real'] and np.isfinite(bs['real'][key] if bs['real'][key] is not None else np.nan):
            ax.axhline(bs['real'][key], color=INK, lw=1.4, ls='--')
            ax.text(3.45, bs['real'][key], 'real', fontsize=7, color=INK, va='bottom', ha='right')
        elif row == 0:
            ax.text(1.5, 1.3, 'real: NaN', fontsize=8, color=INK, ha='center')
        ax.set_xticks(range(4)); ax.set_xticklabels(['T prf', 'T gau', 'N prf', 'N gau'], fontsize=8)
        if row == 0:
            ax.axhline(0, color=GRID, lw=1, zorder=0); ax.axhline(1, color=GRID, lw=1, zorder=0)
            ax.set_title(f'S{s}', color=INK)
        else:
            ax.axhline(3, color=GRID, lw=1, zorder=0)
axs[0, 0].set_ylabel('centroid.py s_par\n(0 = target, 1 = neighbour)')
axs[1, 0].set_ylabel("centroid.py 'offset_sigma'")
axs[0, 0].plot([], [], 'o', color=BLUE, label='injected on target'); axs[0, 0].plot([], [], 'o', color=ORANGE, label='injected on neighbour')
h_, l_ = axs[0, 0].get_legend_handles_labels()
fig.legend(h_, l_, loc='lower center', ncol=2, fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.04))
fig.suptitle('src/centroid.py core logic on transits injected into the real cutouts at fake epochs (11 phases per sector)',
             fontsize=10, color=INK)
fig.savefig(f'{OUT}/fig_centroidpy_injections.png', dpi=130, bbox_inches='tight'); plt.close(fig)
print('ok')
