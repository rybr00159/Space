"""Figures: per sector, the real DIFF image (S/N), and the residual S/N after fitting a*PSF+b at the M3V, the pair
and the Mira (measured apparent position), with the chi2 of each; stars marked at their apparent positions.
Output: figs/diff_images.png"""
import sys
sys.path.insert(0, 'tic214311908/analysis/loc/scripts')
from l2fit import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, axs = plt.subplots(len(SECTORS), 4, figsize=(12, 3.1 * len(SECTORS)))
for r, sec in enumerate(SECTORS):
    M = SectorModel(sec, 'DIFF', grid_half=5.0, grid_step=5.0)
    S = M.S
    full = np.full(S.ny * S.nx, np.nan); full[M.mask] = M.d_real / M.sig
    ext = None
    mx, my = M.mira_xy
    def show(ax, v, title):
        img = np.full(S.ny * S.nx, np.nan); img[M.mask] = v; img = img.reshape(S.ny, S.nx)
        h = ax.imshow(img, origin='lower', cmap='RdBu_r', vmin=-5, vmax=5)
        ax.set_xlim(mx - 5.6, mx + 5.6); ax.set_ylim(my - 5.6, my + 5.6); ax.set_title(title, fontsize=8)
        for nm, c in (('M3V', 'c'), ('pairA', 'm'), ('Mira', 'k'), ('5925373272827260544', 'y')):
            i = int(np.where(M.hyp.name == nm)[0][0]); ax.plot(M.hyp.x[i], M.hyp.y[i], '+', color=c, ms=9, mew=1.8)
        return h
    show(axs[r, 0], M.d_real / M.sig, f'S{sec} DIFF S/N (n={len(M.rv)} transits)')
    c2, a, b, sa = M.fit(M.d_real, 'H')
    for k, nm in enumerate(('M3V', 'pairA', 'Mira')):
        i = int(np.where(M.hyp.name == nm)[0][0])
        res = (M.d_real - a[i] * M.H[i] - b[i]) / M.sig
        show(axs[r, k + 1], res, f'resid after {nm}: chi2 {c2[i]:.0f}/{M.npix - 2}, a={a[i]:.1f}')
for a_ in axs.ravel(): a_.set_xticks([]); a_.set_yticks([])
fig.text(0.5, 0.003, 'cyan + M3V, magenta + pair, black + Mira (measured), yellow + G=17.5 star 3.2" from M3V; colour scale -5..+5 sigma (red = dims)', ha='center', fontsize=9)
fig.tight_layout(rect=(0, 0.01, 1, 1)); fig.savefig(f'{BASE}/figs/diff_images.png', dpi=75)
