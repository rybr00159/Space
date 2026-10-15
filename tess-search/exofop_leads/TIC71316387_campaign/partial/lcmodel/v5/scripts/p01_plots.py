"""Diagnostic figures (matplotlib PNG): folded eclipses with alternative models, profile likelihoods, O-C diagram,
out-of-eclipse harmonic fits, companion exclusion."""
import numpy as np, sys, json, pickle
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v5/lib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import fit as Fi, data as D
V = D.V
C1, C2, C3, C4, C5 = '#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#6250d6'
INK, INK2, GRID, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
plt.rcParams.update({'figure.facecolor': SURF, 'axes.facecolor': SURF, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2,
                     'ytick.color': INK2, 'text.color': INK, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'font.size': 9,
                     'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False})
d = pickle.load(open(f'{V}/results/data_default.pkl', 'rb'))
sk = ['SPOC1', 'SPOC28', 'TC95', 'TC105', 'TC106']
M = Fi.Model(d, scale_keys=sk)
pk = {r['val']: r for r in json.load(open(f'{V}/results/profile_k_default.json'))['profile']}
pl = {r['val']: r for r in json.load(open(f'{V}/results/profile_L3_default.json'))['profile']}
best = min(list(pk.values()) + list(pl.values()), key=lambda r: r['nlp'])
models = [('best geometric (k=%.2f, L3=%.3f)' % (best['p']['k'], best['p']['L3']), best['p'], C1)]
if 0.45 in pk: models.append(('k=0.45, L3=%.3f (MS-like)' % pk[0.45]['p']['L3'], pk[0.45]['p'], C3))
low = min([r for r in pl.values() if r['p']['L3'] < 0.7], key=lambda r: r['p']['L3'], default=None)
if low: models.append(('L3=%.2f, k=%.2f' % (low['p']['L3'], low['p']['k']), low['p'], C5))
if 0.6 in pk: models.append(('k=0.60 (dchi2=%.0f)' % (pk[0.6]['chi2'] - best['chi2']), pk[0.6]['p'], C2))
# detrend data with the best model's per-event baselines
chi, coefs, m = M.chi2(best['p'], return_coefs=True)
base = np.empty(d['t'].size)
for e, (a, b) in enumerate(zip(d['starts'], d['ends'])):
    dt = d['dt'][a:b]; base[a:b] = coefs[e, 0] + coefs[e, 1] * dt + coefs[e, 2] * dt ** 2
y = d['y'] / base
# undo product scale so all products are on the SPOC68 depth scale
scale = np.ones(d['t'].size)
for k in sk: scale[d['key'] == k] = best['p'][f's_{k}']
y = 1 - (1 - y) / scale
P = best['p']['P']; T0 = best['p']['T0']
fig, axs = plt.subplots(2, 2, figsize=(10, 6.4), gridspec_kw=dict(height_ratios=[3, 1.2]), sharex='col')
for col, (ph, lab) in enumerate(((0.0, 'primary'), (0.5, 'secondary'))):
    sel = (d['ph'] == ph) & (d['texp'] < 0.01)   # 200-s and 600-s data only for the binned plot
    x = (((d['t'][sel] - T0) / P - ph + 0.5) % 1 - 0.5) * P * 24
    bins = np.arange(-4.8, 4.81, 0.1); idx = np.digitize(x, bins)
    xb = np.array([x[idx == i].mean() for i in range(1, bins.size) if (idx == i).sum() > 5])
    yb = np.array([y[sel][idx == i].mean() for i in range(1, bins.size) if (idx == i).sum() > 5])
    eb = np.array([y[sel][idx == i].std() / np.sqrt((idx == i).sum()) for i in range(1, bins.size) if (idx == i).sum() > 5])
    ax = axs[0, col]; ax.errorbar(xb, (yb - 1) * 1e3, eb * 1e3, fmt='o', ms=2.5, color=INK2, ecolor=GRID, lw=0.8, label='S28+S68+S95+S105+S106, 6-min bins')
    tt = T0 + ph * P + np.linspace(-0.2, 0.2, 2000)
    for name, p, c in models:
        f = M.eb(p, t=tt, texp=np.full(tt.size, 200 / 86400)); mm = 1 - (1 - p['L3']) * (1 - f)
        ax.plot((tt - T0 - ph * P) * 24, (mm - 1) * 1e3, color=c, lw=1.6, label=name)
        if name.startswith('best'): mref = mm
    ax.set_title(f'{lab} eclipse', loc='left', fontsize=10); ax.set_ylabel('flux - 1 (ppt)')
    ax.set_xlim(-3.0, 3.0)
    if col == 0: ax.legend(fontsize=7, loc='lower left')
    ax = axs[1, col]
    for name, p, c in models:
        f = M.eb(p, t=tt, texp=np.full(tt.size, 200 / 86400)); mm = 1 - (1 - p['L3']) * (1 - f)
        ax.plot((tt - T0 - ph * P) * 24, (mm - mref) * 1e6, color=c, lw=1.4)
    ax.errorbar(xb, (yb - np.interp(xb, (tt - T0 - ph * P) * 24, mref)) * 1e6, eb * 1e6, fmt='o', ms=2.5, color=INK2, ecolor=GRID, lw=0.8)
    ax.set_ylabel('minus best (ppm)'); ax.set_xlabel('hours from conjunction'); ax.set_ylim(-900, 900)
fig.tight_layout(); fig.savefig(f'{V}/plots/fold_models.png', dpi=150); plt.close(fig)
# profiles
fig, axs = plt.subplots(1, 2, figsize=(9, 3.4))
for ax, prof, nm in ((axs[0], pl, 'L3 (third-light fraction)'), (axs[1], pk, 'k = R2/R1')):
    xs = sorted(prof); ys = [prof[x]['chi2'] - best['chi2'] for x in xs]
    ax.plot(xs, ys, 'o-', color=C1, ms=4, lw=1.6); ax.axhline(9, color=C2, lw=1, ls='--'); ax.axhline(4, color=C4, lw=1, ls=':')
    ax.set_ylim(-1, 60); ax.set_xlabel(nm); ax.set_ylabel('profile delta chi2')
axs[0].text(0.02, 9.6, 'delta chi2 = 9', color=INK2, fontsize=7); fig.tight_layout(); fig.savefig(f'{V}/plots/profiles.png', dpi=150); plt.close(fig)
# O-C
F = json.load(open(f'{V}/results/etvfit_default.json'))
ev = F['events_used']
fig, ax = plt.subplots(figsize=(9, 3.6))
for ph, c, lab in ((0.0, C1, 'primary'), (0.5, C2, 'secondary')):
    e = [x for x in ev if x['ph'] == ph]
    tt = np.array([F['fits']['lin']['coef'][0] + x['n'] * F['fits']['lin']['coef'][1] for x in e])
    yr = 2000 + (tt + 2457000 - 2451545.0) / 365.25
    ax.errorbar(yr + (0.03 if ph else 0), [x['oc_min'] for x in e], [x['sig_min'] for x in e], fmt='o', ms=2, color=c, ecolor=c, alpha=0.35, lw=0.6, label=f'{lab} (per event)')
for k, v in F['per_sector'].items():
    if not k.endswith('both'): continue
    tt = F['fits']['lin']['coef'][0] + v['mean_epoch'] * F['fits']['lin']['coef'][1]; yr = 2000 + (tt + 2457000 - 2451545.0) / 365.25
    ax.errorbar(yr, v['oc_min'], v['err_min'], fmt='s', ms=6, color=INK, ecolor=INK, lw=1.4)
xx = np.linspace(2018.4, 2026.3, 200)
for c_, ls in ((0.10, '--'), (-0.10, '--')):
    q_ = c_ * (xx - 2022.35) ** 2; q_ -= np.polyval(np.polyfit(xx, q_, 1), xx)
    ax.plot(xx, q_, color=C3, ls=ls, lw=1.2)
ax.text(2018.5, 4.2, 'green dashed: quadratic O-C with |c| = 0.10 min/yr^2 (recovered at delta chi2 = 8.8 in the injection test)', fontsize=7, color=INK2)
ax.set_ylim(-15, 15); ax.set_xlabel('year'); ax.set_ylabel('O - C (min)'); ax.legend(fontsize=7, loc='lower right')
fig.tight_layout(); fig.savefig(f'{V}/plots/etv.png', dpi=150); plt.close(fig)
# companion exclusion
import glob
fig, ax = plt.subplots(figsize=(6, 3.4))
for f, c, lab in ((f'{V}/results/companion_MEB0.85_M31.03_d9_e0.8.json', C1, 'delta chi2 > 9, e < 0.8'), (f'{V}/results/companion_MEB0.85_M31.03_d4_e0.8.json', C4, 'delta chi2 > 4, e < 0.8'),
                  (f'{V}/results/companion_MEB0.85_M31.03_d9_e0.5.json', C3, 'delta chi2 > 9, e < 0.5')):
    try: g = json.load(open(f))['grid']
    except Exception: continue
    ax.semilogx([x['a_au'] for x in g], [x['frac_excluded'] for x in g], 'o-', ms=3, color=c, lw=1.6, label=lab)
ax.set_xlabel('outer semi-major axis a (au), M3 = 1.03 Msun'); ax.set_ylabel('fraction of orbits excluded'); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(f'{V}/plots/companion.png', dpi=150); plt.close(fig)
print('plots written')
