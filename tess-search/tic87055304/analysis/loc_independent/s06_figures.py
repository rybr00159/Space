"""Figures: per-sector direct/difference images with fitted positions and Gaia stars, and
per-transit offsets in the sky plane.

usage: python s06_figures.py <real_tag> [prf|gauss] [reg_variant]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from pipeline import *

tag = sys.argv[1] if len(sys.argv) > 1 else 'prf_main'
psfkind = sys.argv[2] if len(sys.argv) > 2 else 'prf'
regv = sys.argv[3] if len(sys.argv) > 3 else 'local_r4'
arr = np.load(f'{HERE}/diff_{tag}.npz')
real = json.load(open(f'{HERE}/diff_{tag}.json'))
joint = json.load(open(f'{HERE}/joint_{tag}.json')) if os.path.exists(f'{HERE}/joint_{tag}.json') else None


def ellipse(ax, xy, cov, **kw):
    w, v = np.linalg.eigh(cov)
    ang = np.degrees(np.arctan2(v[1, 1], v[0, 1]))
    ax.add_patch(Ellipse(xy, 2 * np.sqrt(w[1]), 2 * np.sqrt(w[0]), angle=ang, fill=False, **kw))


fig, axes = plt.subplots(5, 5, figsize=(19, 19))
for r_, sec in enumerate(SECTORS):
    ctx = SectorContext(sec, psfkind=psfkind, reg_variant=regv, rfit=3.0)
    S = ctx.S; g = S.gaia
    D = arr[f's{sec}_D']; sg = arr[f's{sec}_sig']
    st = real[str(sec)]['stacked']
    mimg = S.flux[S.good].mean(0)
    # residuals after single-source fits
    m = ctx.mask
    def resid_single(i):
        c, cov, chi2, full = fit_fixed_positions(D, sg, ctx.psf, [ctx.x[i]], [ctx.y[i]], m)
        return D - full, chi2
    rT, cT = resid_single(ctx.it)
    rN, cN = resid_single(ctx.inb)
    cx, cy = ctx.pT
    lim = (cx - 4.5, cx + 4.5, cy - 4.5, cy + 4.5)
    ext = (-0.5, S.nx - 0.5, -0.5, S.ny - 0.5)
    vmax = np.nanmax(np.abs(D[m]))
    panels = [(np.log10(np.clip(mimg, 10, None)), 'viridis', None, f'S{sec} mean image (log e/s)'),
              (D, 'RdBu_r', (-vmax, vmax), 'difference image (OOT - in), e/s'),
              (D / sg, 'RdBu_r', (-8, 8), 'difference S/N'),
              (rT, 'RdBu_r', (-vmax, vmax), f'resid. after TARGET PSF  chi2={cT:.0f}'),
              (rN, 'RdBu_r', (-vmax, vmax), f'resid. after NEIGHBOUR PSF  chi2={cN:.0f}')]
    for c_, (im, cm, vr, title) in enumerate(panels):
        ax = axes[r_, c_]
        kw = dict(origin='lower', extent=ext, cmap=cm)
        if vr is not None:
            kw.update(vmin=vr[0], vmax=vr[1])
        ax.imshow(im, **kw)
        sel = (g.sep_t < 75) & (g.tmag < 17.5)
        ax.scatter(ctx.x[sel], ctx.y[sel], s=np.clip(5 * (18 - g.tmag[sel]) ** 1.5, 3, 250), facecolor='none', edgecolor='k' if c_ else 'w', lw=0.6)
        ax.plot(*ctx.pT, '+', color='lime' if c_ == 0 else 'g', ms=14, mew=2)
        ax.plot(*ctx.pN, 'x', color='m', ms=11, mew=2)
        ax.plot(ctx.x[ctx.ic3], ctx.y[ctx.ic3], 'D', color='c', ms=5, mfc='none')
        if c_ in (1, 2):
            ax.plot(st['free_x'], st['free_y'], 'o', color='orange', mec='k', ms=7)
            bc = np.array(real[str(sec)]['boot_cov_EN'])
            pc = ctx.jinv @ bc @ ctx.jinv.T
            ellipse(ax, (st['free_x'], st['free_y']), pc, color='orange', lw=1.5)
            ax.add_patch(plt.Circle(tuple(ctx.pT), ctx.rfit, fill=False, ls=':', color='k'))
        ax.set_xlim(lim[0], lim[1]); ax.set_ylim(lim[2], lim[3])
        ax.set_title(title, fontsize=9)
        if c_ == 0:
            # compass: 20" north and east arrows
            p0 = np.array([cx - 3.5, cy - 3.5])
            for vec, lab in [((0, 20.0), 'N'), ((20.0, 0), 'E')]:
                dp = ctx.jinv @ np.array(vec)
                ax.annotate('', xy=p0 + dp, xytext=p0, arrowprops=dict(arrowstyle='->', color='w', lw=1.5))
                ax.text(*(p0 + 1.25 * dp), lab, color='w', fontsize=9, ha='center', va='center')
            ax.set_ylabel(f'S{sec}  cam{S.camera} ccd{S.ccd}\nrow (px)')
fig.suptitle(f'Difference images ({psfkind}). + target, x neighbour TIC 87055315, diamond Gaia ...830336; '
             'orange: free single-PSF position with bootstrap 1-sigma ellipse; dotted: fit region', fontsize=12)
plt.tight_layout(rect=(0, 0, 1, 0.98)); plt.savefig(f'{HERE}/fig_diffimg_{tag}.png', dpi=90); plt.close()

# ---------------------------------------------------------------- per-transit offsets
fig, ax = plt.subplots(1, 2, figsize=(15, 7))
cols = dict(zip(SECTORS, ['tab:blue', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown']))
g = SectorContext(66, psfkind=psfkind, reg_variant=regv).S.gaia
for sec in SECTORS:
    v = real[str(sec)]
    pts = np.array([[p['offT_E'], p['offT_N']] for p in v['per_transit']])
    errs = [np.sqrt(np.diag(np.array(p['free_cov_sky']))) for p in v['per_transit']]
    for p_, e_ in zip(pts, errs):
        ax[0].errorbar(p_[0], p_[1], xerr=e_[0], yerr=e_[1], fmt='.', color=cols[sec], alpha=0.5, lw=0.8)
    ax[0].plot([], [], '.', color=cols[sec], label=f'S{sec} per transit')
    st = v['stacked']
    ax[0].plot(st['offT_E'], st['offT_N'], 's', color=cols[sec], mec='k', ms=9, label=f'S{sec} stacked')
    ellipse(ax[0], (st['offT_E'], st['offT_N']), np.array(v['boot_cov_EN']), color=cols[sec], lw=1.5)
    # per-transit s_par
    sp = [p['s_par'] for p in v['per_transit']]
    ax[1].plot(np.full(len(sp), SECTORS.index(sec)) + np.linspace(-0.15, 0.15, len(sp)), sp, 'o', color=cols[sec], alpha=0.6)
    b = v['bootstrap']['s_par']
    ax[1].errorbar(SECTORS.index(sec) + 0.3, st['s_par'], yerr=[[st['s_par'] - b['p16']], [b['p84'] - st['s_par']]], fmt='s', color=cols[sec], mec='k', ms=9, capsize=4)
sel = (g.sep_t < 40) & (g.phot_g_mean_mag < 20.5)
ax[0].scatter(g.de[sel], g.dn[sel], s=np.clip(4 * (21 - g.phot_g_mean_mag[sel]) ** 1.6, 4, 400), facecolor='none', edgecolor='gray', lw=0.8)
ax[0].plot(0, 0, 'r+', ms=18, mew=2.5, label='target')
nE, nN = g.de[g.source_id == NEIGH_ID].values[0], g.dn[g.source_id == NEIGH_ID].values[0]
ax[0].plot(nE, nN, 'mx', ms=14, mew=2.5, label='neighbour TIC 87055315')
c3 = g[g.source_id == C3_ID]
ax[0].plot(c3.de.values[0], c3.dn.values[0], 'cD', ms=8, mfc='none', mew=1.5, label='Gaia ...830336')
if joint is not None:
    ax[0].plot(*joint['best_EN'], '*', color='gold', mec='k', ms=20, label='joint 5-sector fit')
    ellipse(ax[0], joint['best_EN'], np.array(joint['boot_cov']), color='k', lw=2)
ax[0].set_xlim(25, -15); ax[0].set_ylim(-12, 20); ax[0].set_aspect('equal')
ax[0].set_xlabel('offset East of target (arcsec)'); ax[0].set_ylabel('offset North (arcsec)')
ax[0].legend(fontsize=7, ncol=2, loc='lower left'); ax[0].set_title('Free single-PSF positions: per transit (formal errors) and stacked (bootstrap 1-sigma)')
ax[1].axhline(0, color='r', label='target'); ax[1].axhline(1, color='m', label='neighbour')
ax[1].axhline(c3.de.values[0] * nE / (nE**2 + nN**2) + c3.dn.values[0] * nN / (nE**2 + nN**2), color='c', ls='--', label='...830336 (projected)')
ax[1].set_xticks(range(5)); ax[1].set_xticklabels([f'S{s}' for s in SECTORS])
ax[1].set_ylabel('s_par = projection on target->neighbour axis'); ax[1].legend(fontsize=8)
ax[1].set_title('per-transit (circles) and stacked (squares, bootstrap 68%)')
plt.tight_layout(); plt.savefig(f'{HERE}/fig_per_transit_offsets_{tag}.png', dpi=110); plt.close()
print('done')
