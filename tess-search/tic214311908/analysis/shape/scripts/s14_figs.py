"""s14: figures for the README (matplotlib PNG; palette slots 1-3 of the dataviz reference palette, data in gray).
figs/fold_models.png      : INV and REF folded, 10-min bins, with the planet ('single'), white-dwarf occultation,
                            best background-EB (bg3) and V-shape models; residuals below.
figs/rho_posteriors.png   : rho_circ posteriors (free fits) vs the single-star and twin densities; b vs rho.
figs/dilution_profile.png : min chi2 vs source flux fraction (transit and occultation geometries).
figs/per_sector.png       : per-sector depth and duration scale.
figs/secondary_scan.png   : occultation depth vs phase (INV).
"""
import sys, json, os, numpy as np, pandas as pd
sys.path.insert(0, 'tic214311908/analysis/shape/scripts')
from shapelib import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
C1, C2, C3, C4, GR = '#2a78d6', '#eb6834', '#1baf7a', '#6250d6', '#8a8984'
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True, 'grid.color': '#e6e5e0', 'grid.linewidth': 0.6})
os.makedirs(f'{SD}/figs', exist_ok=True)
def binned(x, y, w, edges):
    keep = (x >= edges[0]) & (x < edges[-1]); x, y, w = x[keep], y[keep], w[keep]
    idx = np.digitize(x, edges) - 1; n = len(edges) - 1
    sw = np.bincount(idx, w, n+1)[:n]; sy = np.bincount(idx, w*y, n+1)[:n]
    ok = sw > 0
    return 0.5*(edges[1:]+edges[:-1])[ok], (sy[ok]/sw[ok]), 1/np.sqrt(sw[ok])
# ------------- fold with models
fig, axes = plt.subplots(2, 2, figsize=(10, 6.2), gridspec_kw=dict(height_ratios=[3, 1.3]), sharex=True)
for j, ds in enumerate(['inv', 'ref']):
    lc = LC(f'{SD}/data/{ds}_win.npz')
    models = {}
    fsj = json.load(open(f'{SD}/out/fit_{ds}_single.json'))['best']
    models['planet (single host, rho prior)'] = (lc.transit(fsj['T0'], fsj['P'], fsj['k'], aR_from_rho(10**fsj['log10rho'], fsj['P']), fsj['b'], fsj['u1'], fsj['u2']), C1, '-')
    T0, P = fsj['T0'], fsj['P']
    fn = f'{SD}/out/compact_{ds}_single.json'
    if os.path.exists(fn):
        w = json.load(open(fn))
        models[f'WD occulted by the M3V (dchi2 {w["dchi2_vs_planet_free"]:+.1f})'] = (1 + w['f_A']*(lc.transit(w['T0'], P, w['k'], w['aR_A'], w['b_occ']*w['k'], 0.1, 0.2) - 1), C2, '-')
    fn = f'{SD}/out/dilution_profile_{ds}.json'
    if os.path.exists(fn):
        g = json.load(open(fn))
        row = [r for r in g['profile'] if abs(r['f_src'] - 0.031) < 1e-6]
        if row:
            r = min([row[0]['transit_k_lt_1'], row[0]['occultation_k_gt_1']], key=lambda q: q['chi2'])
            u1 = 2*np.sqrt(r['q1'])*r['q2']; u2 = np.sqrt(r['q1'])*(1-2*r['q2'])
            models['best geometric blend, f_src 0.031 (dchi2 %+.1f)' % (r['chi2'] - g['planet_free_chi2'])] = (1 + 0.031*(lc.transit(r['T0'], P, r['k'], r['aR'], r['b'], u1, u2) - 1), C3, '--')
    _, off = lc.chi2(models['planet (single host, rho prior)'][0], return_offsets=True)
    fcor = lc.f - off[lc.wi]
    dt = (lc.t - T0 - np.round((lc.t - T0)/P)*P)*24
    edges = np.arange(-6, 6.01, 1/6.)
    x, y, e = binned(dt, fcor, lc.w, edges)
    ax = axes[0, j]; axr = axes[1, j]
    ax.errorbar(x, (y-1)*1e3, e*1e3, fmt='o', ms=3, color=GR, ecolor='#c9c8c2', elinewidth=0.8, label='data, 10-min bins')
    for lab, (m, col, ls) in models.items():
        xm, ym, _ = binned(dt, m, lc.w, edges)
        ax.plot(xm, (ym-1)*1e3, color=col, lw=1.8, ls=ls, label=lab)
        if lab.startswith('planet'): base = (xm, ym)
    for lab, (m, col, ls) in models.items():
        xm, ym, _ = binned(dt, m, lc.w, edges)
        axr.plot(xm, (ym - base[1])*1e3, color=col, lw=1.6, ls=ls)
    axr.errorbar(x, (y - base[1])*1e3, e*1e3, fmt='o', ms=2.5, color=GR, ecolor='#c9c8c2', elinewidth=0.8)
    ax.set_title({'inv': 'Investigator photometry (model B), 979 e-/s scale', 'ref': 'Referee local-window photometry'}[ds], fontsize=9.5)
    ax.set_xlim(-3, 3); ax.set_ylabel('flux - 1 (ppt)'); axr.set_ylabel('minus planet (ppt)'); axr.set_xlabel('hours from mid-transit')
    ax.legend(fontsize=7, loc='lower left', frameon=False)
fig.tight_layout(); fig.savefig(f'{SD}/figs/fold_models.png', dpi=150); plt.close(fig)
# ------------- rho posteriors
fig, axs = plt.subplots(1, 2, figsize=(10, 3.6))
for ds, col in [('inv', C1), ('ref', C2)]:
    z = np.load(f'{SD}/out/fit_{ds}_free.npz')
    axs[0].hist(z['rho_circ'], bins=np.linspace(0, 30, 61), histtype='step', color=col, lw=1.6, density=True, label=f'{ds.upper()} free fit')
    axs[1].plot(z['rho_circ'][::4], z['b'][::4], '.', ms=1.5, color=col, alpha=0.4, label=f'{ds.upper()}')
for v, s, lab, ls in [(9.41, 0.88, 'single star 9.4', '-'), (15.5, 1.45, 'twin 15.5', '--')]:
    axs[0].axvspan(v-s, v+s, color='#d9d8d2', alpha=0.6); axs[0].axvline(v, color='#52514e', ls=ls, lw=1.2, label=lab)
axs[0].set_xlabel('rho_circ (g/cm3)'); axs[0].set_ylabel('posterior density'); axs[0].legend(fontsize=7.5, frameon=False)
axs[1].set_xlabel('rho_circ (g/cm3)'); axs[1].set_ylabel('impact parameter b'); axs[1].set_xlim(0, 30); axs[1].legend(fontsize=7.5, frameon=False, markerscale=6)
fig.tight_layout(); fig.savefig(f'{SD}/figs/rho_posteriors.png', dpi=150); plt.close(fig)
# ------------- dilution profile
if all(os.path.exists(f'{SD}/out/dilution_profile_{d}.json') for d in ['inv', 'ref']):
    fig, ax = plt.subplots(figsize=(7, 3.8))
    for ds, col in [('inv', C1), ('ref', C2)]:
        g = json.load(open(f'{SD}/out/dilution_profile_{ds}.json'))
        f = np.array([r['f_src'] for r in g['profile']]); dm = -2.5*np.log10(f)
        ct = np.array([r['transit_k_lt_1']['chi2'] for r in g['profile']]) - g['planet_free_chi2']
        co = np.array([r['occultation_k_gt_1']['chi2'] for r in g['profile']]) - g['planet_free_chi2']
        ax.plot(dm, ct, '-o', color=col, ms=3.5, lw=1.6, label=f'{ds.upper()}: occulter smaller than source (k<1)')
        ax.plot(dm, co, '--s', color=col, ms=3.5, lw=1.4, label=f'{ds.upper()}: occulter larger (k>1)')
    for x, lab in [(3.73, 'bg3 3.1"'), (4.08, 'bg5 5.5"')]:
        ax.axvline(x, color='#52514e', lw=0.9, ls=':'); ax.text(x+0.03, ax.get_ylim()[1]*0.9 if ax.get_ylim()[1] > 0 else 5, lab, fontsize=7.5, rotation=90, va='top')
    ax.axvspan(3.5, 4.5, color='#e6e5e0', alpha=0.6, lw=0)
    ax.axhline(9, color='#52514e', lw=0.8); ax.text(0.05, 9.5, 'dchi2 = 9', fontsize=7.5)
    ax.set_xlabel('source magnitude difference dT = -2.5 log10(f_src)'); ax.set_ylabel('min chi2 - planet chi2'); ax.set_yscale('symlog', linthresh=10)
    ax.legend(fontsize=7, frameon=False, loc='center left', bbox_to_anchor=(0.0, 0.55)); fig.tight_layout(); fig.savefig(f'{SD}/figs/dilution_profile.png', dpi=150); plt.close(fig)
# ------------- per sector
if os.path.exists(f'{SD}/out/per_sector.csv'):
    D = pd.read_csv(f'{SD}/out/per_sector.csv'); D = D[D.subset.str.isnumeric()]
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.4))
    for i, (ds, col) in enumerate([('inv', C1), ('invA', C3), ('ref', C2)]):
        q = D[D.dataset == ds]; x = np.arange(len(q)) + (i-1)*0.15
        axs[0].errorbar(x, q.depth_ppt, q.depth_err_ppt, fmt='o', color=col, ms=4, label=ds.upper())
        axs[1].errorbar(x, q.dur_scale, q.dur_scale_err, fmt='o', color=col, ms=4, label=ds.upper())
    for ax in axs: ax.set_xticks(range(len(q))); ax.set_xticklabels([f'S{s}' for s in q.subset]); ax.legend(fontsize=7.5, frameon=False)
    axs[0].set_ylabel('depth (ppt, 979 e-/s scale)'); axs[1].set_ylabel('duration scale (1 = global)')
    fig.tight_layout(); fig.savefig(f'{SD}/figs/per_sector.png', dpi=150); plt.close(fig)
# ------------- secondary scan
if os.path.exists(f'{SD}/out/secondary_scan.csv'):
    S_ = pd.read_csv(f'{SD}/out/secondary_scan.csv'); J = json.load(open(f'{SD}/out/secondary_phase.json'))
    fig, ax = plt.subplots(figsize=(7, 3))
    ax.plot(S_.phase, S_.depth*1e3, color=C1, lw=1.4, label='INV occultation depth')
    e = J['secondary']['inv_err_empirical_ppt']
    ax.axhspan(-3*e, 3*e, color='#e6e5e0', alpha=0.7, lw=0, label='+-3 sigma (empirical)')
    ax.axvline(0.5, color='#52514e', lw=0.8, ls=':')
    ax.set_xlabel('phase'); ax.set_ylabel('depth (ppt)'); ax.legend(fontsize=7.5, frameon=False)
    fig.tight_layout(); fig.savefig(f'{SD}/figs/secondary_scan.png', dpi=150); plt.close(fig)
print('figures written')
