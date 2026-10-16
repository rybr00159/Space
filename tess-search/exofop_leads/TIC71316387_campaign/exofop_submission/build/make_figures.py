"""Figures for the TIC 71316387 supporting report (CTOI 71316387.01).
Run with the work venv (numba model code):  /home/user/work/venv/bin/python make_figures.py
Inputs are the campaign's own result files; nothing is refitted here except per-sector depth scaling for plotting."""
import sys, json, pickle, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

SCR = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad'
V6 = f'{SCR}/tic713/lcmodel/v6'
L3D = f'{SCR}/tic713/localise/v3'
OUT = f'{SCR}/obsnote/report/figs'
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, f'{V6}/lib')
import fit as Fi  # noqa: E402

for f in ('LiberationSans-Regular', 'LiberationSans-Bold', 'LiberationSans-Italic'):
    font_manager.fontManager.addfont(f'/usr/share/fonts/truetype/liberation/{f}.ttf')

BLUE, ORANGE, AQUA = '#2a78d6', '#eb6834', '#1baf7a'
INK, INK2, GRID, MUTE = '#0b0b0b', '#52514e', '#e4e3df', '#9a9893'
plt.rcParams.update({
    'font.family': ['Liberation Sans', 'DejaVu Sans'], 'font.size': 8, 'axes.titlesize': 8.5,
    'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 6.8,
    'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white',
    'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
    'text.color': INK, 'axes.linewidth': 0.7, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
    'axes.grid': False, 'legend.frameon': False, 'lines.linewidth': 1.2, 'mathtext.default': 'regular'})


def tag(ax, s):
    ax.text(0.015, 0.975, s, transform=ax.transAxes, ha='left', va='top', fontsize=8.5, fontweight='bold', color=INK)


def binned(x, y, w):
    edges = np.arange(x.min(), x.max() + w, w)
    i = np.digitize(x, edges)
    xb, yb, eb = [], [], []
    for k in np.unique(i):
        m = i == k
        if m.sum() < 6:
            continue
        xb.append(x[m].mean()); yb.append(y[m].mean()); eb.append(y[m].std(ddof=1) / np.sqrt(m.sum()))
    return np.array(xb), np.array(yb), np.array(eb)


# ------------------------------------------------------------------ light curve data and adopted model
d = pickle.load(open(f'{V6}/results/data_default.pkl', 'rb'))
sk = ['SPOC1', 'SPOC28', 'TC95', 'TC105', 'TC106']
M = Fi.Model(d, scale_keys=sk)
ms = json.load(open(f'{V6}/results/msml_msfree_ml.json'))
p = {k: ms['geom'][k] for k in Fi.BASE}
for k in sk:
    p[f's_{k}'] = ms['geom'][f's_{k}']
chi, coefs, _ = M.chi2(p, return_coefs=True)
base = np.empty(d['t'].size)
for e, (a, b) in enumerate(zip(d['starts'], d['ends'])):
    dt = d['dt'][a:b]
    base[a:b] = coefs[e, 0] + coefs[e, 1] * dt + coefs[e, 2] * dt ** 2
y = d['y'] / base
scale = np.ones(d['t'].size)
for k in sk:
    scale[d['key'] == k] = p[f's_{k}']
y = 1 - (1 - y) / scale          # all products on the SPOC S68 depth scale
P, T0 = p['P'], p['T0']
fast = d['texp'] < 0.003         # 200-s sectors: S68, S95, S105, S106


def model_curve(ph, hw=0.17, texp=200 / 86400):
    tt = T0 + ph * P + np.linspace(-hw, hw, 3000)
    f = M.eb(p, t=tt, texp=np.full(tt.size, texp))
    return (tt - T0 - ph * P) * 24, 1 - (1 - p['L3']) * (1 - f)


xm0, m0 = model_curve(0.0)
xm5, m5 = model_curve(0.5)
dep_p = (1 - m0.min()) * 1e6
dep_s = (1 - m5.min()) * 1e6

# ------------------------------------------------------------------ FIGURE 1: light curve
fig = plt.figure(figsize=(6.9, 5.0))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.1], hspace=0.42, wspace=0.26, left=0.085, right=0.985, top=0.965, bottom=0.085)

# (a) fold at the CTOI period, odd vs even events
ax = fig.add_subplot(gs[0, 0])
for ph, c, lab in ((0.0, BLUE, 'odd events'), (0.5, ORANGE, 'even events')):
    s = (d['ph'] == ph) & fast
    x = (((d['t'][s] - T0) / P - ph + 0.5) % 1 - 0.5) * P * 24
    xb, yb, eb = binned(x, (y[s] - 1) * 1e3, 0.1)
    ax.errorbar(xb, yb, eb, fmt='o', ms=2.4, color=c, ecolor=c, elinewidth=0.6, alpha=0.95, label=lab, zorder=3)
ax.axhline(0, color=GRID, lw=0.7, zorder=1)
ax.set_xlim(-2.4, 2.4); ax.set_ylim(-6.9, 1.2)
ax.set_xlabel('Hours from mid-eclipse (folded at the CTOI period, 0.65184 d)')
ax.set_ylabel('Relative flux − 1 (ppt)')
ax.legend(loc='lower left', handletextpad=0.3, borderaxespad=0.3)
tag(ax, 'a')

# (b) per-sector depths
ax = fig.add_subplot(gs[0, 1])
dep = json.load(open(f'{V6}/results/depths_by_product.json'))
keys = ['SPOC1', 'SPOC28', 'SPOC68', 'TC95', 'TC105', 'TC106']
labs = ['S1\n2018', 'S28\n2020', 'S68\n2023', 'S95\n2025', 'S105\n2026', 'S106\n2026']
tp = dep['p']
tmpl = dict(tp)
mt0 = M.eb(tmpl, t=np.array([T0]), texp=np.array([200 / 86400]))
mt5 = M.eb(tmpl, t=np.array([T0 + 0.5 * P]), texp=np.array([200 / 86400]))
td_p = (1 - tmpl['L3']) * (1 - mt0[0]) * 1e6
td_s = (1 - tmpl['L3']) * (1 - mt5[0]) * 1e6
xi = np.arange(len(keys))
for kk, (typ, c, td, off) in enumerate((('prim', BLUE, td_p, -0.09), ('sec', ORANGE, td_s, 0.09))):
    v = np.array([dep['depths'][k][typ]['s'] for k in keys]) * td
    e = np.array([dep['depths'][k][typ]['err_adopted'] for k in keys]) * td
    ax.errorbar(xi + off, v, e, fmt='o', ms=3.4, color=c, ecolor=c, elinewidth=0.9, capsize=0,
                label=('deeper (primary) eclipse' if typ == 'prim' else 'shallower (secondary) eclipse'))
ax.axhspan(5770 - 100, 5770 + 100, color=BLUE, alpha=0.12, lw=0)
ax.axhspan(2600 - 55, 2600 + 55, color=ORANGE, alpha=0.14, lw=0)
ax.axhline(5770, color=BLUE, lw=0.7, ls='--'); ax.axhline(2600, color=ORANGE, lw=0.7, ls='--')
ax.set_xticks(xi); ax.set_xticklabels(labs, fontsize=6.6)
ax.set_xlim(-0.5, len(keys) - 0.5); ax.set_ylim(0, 7600)
ax.set_ylabel('Eclipse depth (ppm, diluted)')
ax.legend(loc='center right', handletextpad=0.3, bbox_to_anchor=(1.0, 0.63))
tag(ax, 'b')

# (c),(d) eclipses at P_EB with the adopted model, residuals below
for col, (ph, c, xm, mm, ttl) in enumerate(((0.0, BLUE, xm0, m0, 'Primary eclipse'), (0.5, ORANGE, xm5, m5, 'Secondary eclipse'))):
    sub = gs[1, col].subgridspec(2, 1, height_ratios=[3, 1], hspace=0.0)
    ax = fig.add_subplot(sub[0]); axr = fig.add_subplot(sub[1], sharex=ax)
    s = (d['ph'] == ph) & fast
    x = (((d['t'][s] - T0) / P - ph + 0.5) % 1 - 0.5) * P * 24
    xb, yb, eb = binned(x, (y[s] - 1) * 1e3, 0.1)
    ax.errorbar(xb, yb, eb, fmt='o', ms=2.4, color=INK2, ecolor=MUTE, elinewidth=0.6, zorder=2, label='TESS (6-min bins)')
    ax.plot(xm, (mm - 1) * 1e3, color=c, lw=1.4, zorder=3, label='adopted model')
    ax.axhline(0, color=GRID, lw=0.7, zorder=1)
    ax.set_xlim(-2.4, 2.4)
    ax.set_ylim((-6.9, 1.3) if ph == 0 else (-3.3, 0.95))
    ax.set_ylabel('Relative flux − 1 (ppt)')
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.text(0.98, 0.955, ttl + (' (phase 0)' if ph == 0 else ' (phase 0.5)'), transform=ax.transAxes, ha='right', va='top', fontsize=7.2, color=INK2)
    if col == 1:
        ax.annotate('model: flat bottom', xy=(0.12, -2.62), xytext=(1.05, -2.45), fontsize=6.6, color=INK2, va='center',
                    arrowprops=dict(arrowstyle='-', color=MUTE, lw=0.6, shrinkA=1, shrinkB=1))
    ax.legend(loc='lower left', handletextpad=0.3, borderaxespad=0.2)
    tag(ax, 'c' if col == 0 else 'd')
    r = yb - (np.interp(xb, xm, mm) - 1) * 1e3
    axr.errorbar(xb, r * 1e3, eb * 1e3, fmt='o', ms=2.0, color=INK2, ecolor=MUTE, elinewidth=0.5)
    axr.axhline(0, color=c, lw=1.0)
    axr.set_ylim(-650, 650); axr.set_yticks([-400, 0, 400])
    axr.set_ylabel('O−C\n(ppm)', fontsize=7)
    axr.set_xlabel('Hours from conjunction ($P_{EB}$ = 1.30367 d)')
fig.savefig(f'{OUT}/fig1_lightcurve.png', dpi=300)
fig.savefig(f'{OUT}/fig1_lightcurve.pdf')
plt.close(fig)

# ------------------------------------------------------------------ FIGURE 2: localisation
r = json.load(open(f'{L3D}/results/s05_baseline.json'))
z = np.load(f'{L3D}/results/s05_baseline_map.npz')
fig, axs = plt.subplots(1, 2, figsize=(6.9, 3.35), gridspec_kw=dict(width_ratios=[1, 1], wspace=0.28, left=0.085, right=0.985, top=0.97, bottom=0.13))
ax = axs[0]
for i, lab in enumerate(r['labels']):
    o = np.array(r['obs'][i]); u = np.array(r['u_to_boresight'][i]); c = np.array(r['cov'][i])
    deep = 'deep' in lab
    if deep:
        ax.annotate('', xy=(-1.55 * u[0], -1.55 * u[1]), xytext=(0, 0),
                    arrowprops=dict(arrowstyle='-|>', color=MUTE, lw=0.8, mutation_scale=7, shrinkA=0, shrinkB=0), zorder=1)
        ax.errorbar(o[0], o[1], xerr=np.sqrt(c[0, 0]), yerr=np.sqrt(c[1, 1]), fmt='o', ms=4.2, color=BLUE, mec='white', mew=0.7,
                    elinewidth=0.8, capsize=0, zorder=3)
        sec = lab.split('_')[0]
        dx, dy = {'S1': (0.55, 0.55), 'S28': (0.0, -0.78), 'S68': (-0.78, 0.0), 'S95': (0.05, -0.78), 'S105': (0.0, 0.8), 'S106': (0.0, -0.85)}.get(sec, (0.3, 0.3))
        ax.text(o[0] + dx, o[1] + dy, sec, fontsize=6.6, color=INK, ha='center', va='center')
ax.plot(0, 0, marker='*', ms=10, color=INK, mec='white', mew=0.5, zorder=4, ls='none')
ax.plot([], [], 'o', ms=4, color=BLUE, label='deeper eclipse, per sector (raw)')
ax.plot([], [], '-', color=MUTE, lw=0.8, label='direction away from camera boresight')
ax.plot([], [], '*', ms=8, color=INK, ls='none', label='G star (Gaia DR3)')
ax.set_xlim(5.0, -5.0); ax.set_ylim(-5.0, 5.0); ax.set_aspect('equal')
ax.axhline(0, color=GRID, lw=0.5, zorder=0); ax.axvline(0, color=GRID, lw=0.5, zorder=0)
ax.set_xlabel('ΔE (arcsec; east to the left)'); ax.set_ylabel('ΔN (arcsec)')
ax.legend(loc='lower left', handletextpad=0.3, borderaxespad=0.2, fontsize=6.4)
tag(ax, 'a')

ax = axs[1]
g = z['g']; Mp = z['M']
lv = [r['thresholds'][k] for k in ('68.27', '95.45', '99.73')]
ax.contourf(g, g, Mp, levels=[-1e-9] + lv, colors=['#2a78d6', '#86b4ea', '#d3e3f7'], zorder=1)
ax.contour(g, g, Mp, levels=lv, colors=[BLUE], linewidths=0.6, zorder=2)
t = np.linspace(0, 2 * np.pi, 400)
ax.plot(0.6 * np.sin(t), 0.6 * np.cos(t), '--', color=ORANGE, lw=1.0, zorder=3)
ax.plot(1.0 * np.sin(t), 1.0 * np.cos(t), ':', color=ORANGE, lw=1.2, zorder=3)
pa = np.radians(348)
ax.plot(1.13 * np.sin(pa), 1.13 * np.cos(pa), 'D', ms=4.2, color=ORANGE, mec='white', mew=0.6, zorder=4)
dl = r['M2_delta']
ax.plot(dl[0], dl[1], 'x', ms=6, mew=1.4, color=INK, zorder=5)
ax.plot(0, 0, marker='*', ms=10, color=INK, mec='white', mew=0.5, zorder=5, ls='none')
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
h = [Patch(color='#2a78d6', label='68%'), Patch(color='#86b4ea', label='95%'), Patch(color='#d3e3f7', label='99.7%'),
     Line2D([], [], marker='x', ls='none', color=INK, mew=1.4, ms=5, label='best fit (colour-corrected)'),
     Line2D([], [], ls='--', color=ORANGE, lw=1.0, label='VHS: Δ$K_s$ ≤ 2.5 excluded beyond 0.6 arcsec'),
     Line2D([], [], ls=':', color=ORANGE, lw=1.2, label='VHS: Δ$K_s$ ≤ 3.0 excluded beyond 1.0 arcsec'),
     Line2D([], [], marker='D', ls='none', color=ORANGE, mec='white', ms=4, label='unexplained VHS residual (Δ$K_s$ ≈ 4.0–4.4)')]
ax.legend(handles=h, loc='lower left', fontsize=6.0, handletextpad=0.35, borderaxespad=0.2, ncol=1, labelspacing=0.25)
ax.set_xlim(1.9, -1.9); ax.set_ylim(-1.9, 1.9); ax.set_aspect('equal')
ax.set_xlabel('ΔE (arcsec; east to the left)'); ax.set_ylabel('ΔN (arcsec)')
tag(ax, 'b')
fig.savefig(f'{OUT}/fig2_localisation.png', dpi=300)
fig.savefig(f'{OUT}/fig2_localisation.pdf')
plt.close(fig)

# ------------------------------------------------------------------ FIGURE 3: colour of the eclipsed light
cal = json.load(open(f'{L3D}/results/f02_colourw_theta.json'))   # wide calibration sample (the BP-RP 2.31 read-off)
pts = np.load(f'{L3D}/results/f02_pointsw.npz')
fig, axs = plt.subplots(1, 2, figsize=(6.9, 3.0), gridspec_kw=dict(wspace=0.27, left=0.085, right=0.985, top=0.965, bottom=0.15))
ax = axs[0]
ok = pts['ok']
ax.plot(pts['C'][ok], -pts['rad'][ok], '.', ms=1.6, color=MUTE, alpha=0.45, zorder=1, rasterized=True)
cv = cal['curve']
ax.errorbar([q['C'] for q in cv], [-q['rad'] for q in cv], [q['err'] for q in cv], fmt='s', ms=3.2, color=INK, ecolor=INK,
            elinewidth=0.8, zorder=3, label='field stars, binned (5,560 fits, 6 sectors)')
ax.axhspan(2.94 - 0.28, 2.94 + 0.28, color=BLUE, alpha=0.15, lw=0, zorder=0)
ax.axhline(2.94, color=BLUE, lw=1.0, zorder=2, label='eclipsing light: shift 2.94 ± 0.28 arcsec')
ax.axvspan(2.22, 2.40, color=BLUE, alpha=0.15, lw=0, zorder=0)
ax.axvline(0.786, color=INK2, lw=0.7, ls=':', zorder=2)
ax.text(0.83, 4.15, 'G star (BP−RP 0.79)', fontsize=6.4, color=INK2, va='center')
ax.text(2.45, 4.3, 'BP−RP ≈ 2.3', fontsize=6.8, color=BLUE, ha='left', va='center')
ax.set_xlim(0.3, 3.2); ax.set_ylim(-1.2, 4.6)
ax.set_xlabel('Gaia BP−RP colour (mag)')
ax.set_ylabel('Centroid shift away from boresight (arcsec)')
ax.legend(loc='lower right', handletextpad=0.3, borderaxespad=0.3, fontsize=6.4)
tag(ax, 'a')

ax = axs[1]
lam = {'g′': 0.475, 'V': 0.551, 'r′': 0.622, 'G': 0.64, 'i′': 0.763, 'TESS': 0.80, 'z′': 0.905, 'J': 1.25, 'H': 1.65, 'Ks': 2.15, 'W1': 3.35}
pred = {'g′': 1150, 'V': 1900, 'r′': 3075, 'G': 3600, 'i′': 6150, 'TESS': 5770, 'z′': 7200, 'J': 12100, 'H': 15800, 'Ks': 18000, 'W1': 18900}
ordk = ['g′', 'V', 'r′', 'TESS', 'z′', 'J', 'H', 'Ks', 'W1']
pred['V'] = 0.305 * 5770   # V/TESS 0.28-0.33 (model-dependent); band spans the full range
xl = np.array([lam[k] for k in ordk]); yl = np.array([pred[k] / 5770 for k in ordk])
fb = np.where(xl > 1.0, 0.20, 0.15)
lo_, hi_ = yl * (1 - fb), yl * (1 + fb)
lo_[ordk.index('V')] = 0.28 * 0.85; hi_[ordk.index('V')] = 0.33 * 1.15
ax.fill_between(xl, lo_, hi_, color=BLUE, alpha=0.13, lw=0)
ax.plot(xl, yl, '-', color=BLUE, lw=1.2, label='predicted for the adopted M2V + M4–4.5V pair')
for xg, lo, hi in ((0.551, 0.82, 1.31), (3.35, 0.70, 1.35)):
    ax.add_patch(plt.Rectangle((xg / 1.09, lo), xg * 1.09 - xg / 1.09, hi - lo, color=MUTE, alpha=0.35, lw=0))
ax.add_patch(plt.Rectangle((0, 0), 0, 0, color=MUTE, alpha=0.35, lw=0, label='expected if a G/F star were eclipsed'))
ax.errorbar([0.551 / 1.03], [0.418], [0.145], fmt='o', ms=4, color=ORANGE, ecolor=ORANGE, elinewidth=1.0, capsize=0, zorder=4,
            label='measured (SuperWASP, NEOWISE W1)')
ax.errorbar([0.551 * 1.03], [0.36], [0.11], fmt='o', ms=4, mfc='white', color=ORANGE, ecolor=ORANGE, elinewidth=1.0, capsize=0, zorder=4)
ax.errorbar([3.35], [3.67], [1.14], fmt='o', ms=4, color=ORANGE, ecolor=ORANGE, elinewidth=1.0, capsize=0, zorder=4)
ax.plot([0.80], [1.0], 'D', ms=4, color=INK, zorder=5, label='TESS (reference)')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(0.42, 4.2); ax.set_ylim(0.13, 7)
ax.set_xticks([0.5, 0.7, 1, 2, 3]); ax.set_xticklabels(['0.5', '0.7', '1', '2', '3'])
ax.set_yticks([0.2, 0.5, 1, 2, 5]); ax.set_yticklabels(['0.2', '0.5', '1', '2', '5'])
ax.minorticks_off()
ax.set_xlabel('Wavelength (µm)')
ax.set_ylabel('Deeper-eclipse depth / TESS depth')
ax.legend(loc='lower right', handletextpad=0.35, borderaxespad=0.3, fontsize=6.2, labelspacing=0.3)
ax.text(0.615, 0.40, 'SuperWASP\n(2 reductions)', fontsize=6.2, color=ORANGE, ha='left', va='top')
ax.text(3.1, 2.25, 'W1', fontsize=6.4, color=ORANGE, ha='right', va='center')
ax.text(0.551, 1.42, 'G/F', fontsize=6.2, color=INK2, ha='center', va='bottom')
ax.text(3.35, 1.42, 'G/F', fontsize=6.2, color=INK2, ha='center', va='bottom')
tag(ax, 'b')
fig.savefig(f'{OUT}/fig3_colour.png', dpi=300)
fig.savefig(f'{OUT}/fig3_colour.pdf')
plt.close(fig)

summary = dict(model='msml_msfree_ml.json geom', model_L3=p['L3'], model_depth_primary_ppm=dep_p, model_depth_secondary_ppm=dep_s,
               template_depth_primary_ppm=td_p, template_depth_secondary_ppm=td_s, chi2=float(chi))
json.dump(summary, open(f'{OUT}/figure_numbers.json', 'w'), indent=1)
print(json.dumps(summary, indent=1))
