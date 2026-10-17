"""Figures from the simulations: (top) combined delta-chi2 maps (DIFF, AMP) with the injection-calibrated
68/95/99.7% contours and the catalogued stars at the mean epoch of the data (co-moving frame of the M3V; positions
are true sky positions for a source of the M3V's colour); (bottom) per-sector free-fit positions of the real signal
(stars) over the clouds of free-fit positions recovered from injections on the M3V, the pair, the Mira and the
G = 17.5 star 3.2" from the M3V.
Output: figs/localisation.png"""
import sys, json
sys.path.insert(0, 'tic214311908/analysis/loc/scripts')
import numpy as np
from l0lib import BASE, IDS, stars
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
tags = [t for t in ('diff_prf_all', 'amp_prf_all')]
g = stars()
r3 = g[g.source_id == IDS['M3V']].iloc[0]
fig, axs = plt.subplots(2, 6, figsize=(22, 8.5), gridspec_kw=dict(height_ratios=[1.35, 1]))
for c, tag in enumerate(tags):
    try:
        Z = np.load(f'{BASE}/out/sim_{tag}.npz', allow_pickle=True); P = json.load(open(f'{BASE}/out/post_{tag}.json'))
    except FileNotFoundError:
        continue
    ax = plt.subplot2grid((2, 6), (0, 3 * c), colspan=3, fig=fig)
    axv = Z['grid_axis']; n = len(axv); gs = Z['real__gsum']; M = (gs - gs.min()).reshape(n, n)
    T = P['loc_thresholds_used']
    h = ax.imshow(np.sqrt(np.clip(M, 0, 400)), origin='lower', extent=[axv[0] - .5, axv[-1] + .5, axv[0] - .5, axv[-1] + .5], cmap='magma')
    plt.colorbar(h, ax=ax, label='sqrt(delta chi2)')
    ax.contour(axv, axv, M, levels=T, colors=['cyan', 'lime', 'w'], linewidths=1.3)
    ep = 2023.6 - 2016.0
    for _, s in g[(np.hypot(g.dE - r3.dE, g.dN - r3.dN) < 45) & (g['T'] < 17.6)].iterrows():
        xE = (s.dE - r3.dE) + (np.nan_to_num(s.pmra) - r3.pmra) * ep / 1000; xN = (s.dN - r3.dN) + (np.nan_to_num(s.pmdec) - r3.pmdec) * ep / 1000
        ax.plot(xE, xN, 'o', mfc='none', mec='w', ms=max(4, 16 - 1.2 * (s['T'] - 11)), mew=0.8)
    for nm, col, lab in (('M3V', 'cyan', 'M3V'), ('pairA', 'm', 'pair A/B'), ('Mira', 'r', 'Mira (Gaia)'), ):
        s = g[g.source_id == IDS[nm]].iloc[0]
        xE = (s.dE - r3.dE) + (s.pmra - r3.pmra) * ep / 1000; xN = (s.dN - r3.dN) + (s.pmdec - r3.pmdec) * ep / 1000
        ax.plot(xE, xN, '+', color=col, ms=14, mew=2, label=lab)
    for sid, lab in ((5925373272827260544, 'G17.5, 3.2"'), (5925373272827262208, 'G17.7, 5.7"')):
        s = g[g.source_id == sid].iloc[0]
        xE = (s.dE - r3.dE) + (s.pmra - r3.pmra) * ep / 1000; xN = (s.dN - r3.dN) + (s.pmdec - r3.pmdec) * ep / 1000
        ax.plot(xE, xN, 'x', color='y', ms=9, mew=1.5); ax.annotate(lab, (xE, xN), color='y', fontsize=7, xytext=(3, 3), textcoords='offset points')
    ax.set_xlim(30, -30); ax.set_ylim(-22, 22); ax.set_xlabel('dE from M3V (arcsec)'); ax.set_ylabel('dN (arcsec)')
    ax.legend(loc='lower left', fontsize=7)
    ax.set_title(f"{tag.split('_')[0].upper()}: combined 6-sector delta-chi2 map; contours 68/95/99.7% (injection-calibrated: "
                 f"{T[0]:.1f}/{T[1]:.1f}/{T[2]:.1f})", fontsize=8)
    # per-sector free fits
    if c == 0:
        names = list(Z['hyp_name']); secs = list(Z['sectors'])
        rf = P['real_free_fit_offset_from_M3V_arcsec']
        for k, s in enumerate(secs):
            a = axs[1, k]
            for nm, col in (('M3V', 'c'), ('pairA', 'm'), ('Mira', 'r'), ('5925373272827260544', 'y')):
                F = Z[f'inj_{nm}__free'][:, k]
                a.plot(F[:300, 0], F[:300, 1], '.', color=col, ms=1.5, alpha=0.5)
            a.plot(*rf[str(s)], '*', color='k', ms=14, mec='w')
            a.set_xlim(40, -40); a.set_ylim(-30, 30); a.set_title(f'S{s} DIFF free fit (black star) vs injections', fontsize=8)
            a.set_aspect('equal')
        axs[1, 0].text(38, -28, 'cyan M3V, magenta pair, red Mira, yellow G17.5 3.2" star', fontsize=7)
for a in axs[0]: a.set_visible(False)
fig.tight_layout(); fig.savefig(f'{BASE}/figs/localisation.png', dpi=75)
