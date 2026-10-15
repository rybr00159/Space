"""Diagnostic figures (matplotlib PNG) for the v6 light-curve model and timing analysis."""
import numpy as np, sys, json, pickle, os
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6/lib')
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
best = min(list(pk.values()) + list(pl.values()), key=lambda r: r['chi2'])
# ---------------- 1. folded eclipses with models along the degeneracy
models = [('best geometric (L3=%.3f, k=%.2f, b0=%.2f)' % (best['p']['L3'], best['p']['k'], best['p']['b0']), best['p'], C1)]
for L3v, c in ((0.975, C3), (0.5, C5), (0.0, C4)):
    if L3v in pl: models.append(('L3=%.3f, k=%.2f, b0=%.2f (dchi2=%.1f)' % (L3v, pl[L3v]['p']['k'], pl[L3v]['p']['b0'], pl[L3v]['chi2'] - best['chi2']), pl[L3v]['p'], c))
if 0.6 in pk: models.append(('k=0.60 (dchi2=%.0f)' % (pk[0.6]['chi2'] - best['chi2']), pk[0.6]['p'], C2))
chi, coefs, m = M.chi2(best['p'], return_coefs=True)
base = np.empty(d['t'].size)
for e, (a, b) in enumerate(zip(d['starts'], d['ends'])):
    dt = d['dt'][a:b]; base[a:b] = coefs[e, 0] + coefs[e, 1] * dt + coefs[e, 2] * dt ** 2
y = d['y'] / base
scale = np.ones(d['t'].size)
for k in sk: scale[d['key'] == k] = best['p'][f's_{k}']
y = 1 - (1 - y) / scale
P = best['p']['P']; T0 = best['p']['T0']
fig, axs = plt.subplots(2, 2, figsize=(10.5, 6.6), gridspec_kw=dict(height_ratios=[3, 1.3]), sharex='col')
for col, (ph, lab) in enumerate(((0.0, 'primary'), (0.5, 'secondary'))):
    sel = (d['ph'] == ph) & (d['texp'] < 0.003)   # 200-s data only (S68, S95, S105, S106)
    x = (((d['t'][sel] - T0) / P - ph + 0.5) % 1 - 0.5) * P * 24
    bins = np.arange(-4.8, 4.81, 0.1); idx = np.digitize(x, bins)
    ok = [i for i in range(1, bins.size) if (idx == i).sum() > 5]
    xb = np.array([x[idx == i].mean() for i in ok]); yb = np.array([y[sel][idx == i].mean() for i in ok])
    eb = np.array([y[sel][idx == i].std() / np.sqrt((idx == i).sum()) for i in ok])
    ax = axs[0, col]; ax.errorbar(xb, (yb - 1) * 1e3, eb * 1e3, fmt='o', ms=2.5, color=INK2, ecolor=GRID, lw=0.8, label='200-s sectors (S68, S95, S105, S106), 6-min bins')
    tt = T0 + ph * P + np.linspace(-0.2, 0.2, 2000)
    mref = None
    for name, p, c in models:
        f = M.eb(p, t=tt, texp=np.full(tt.size, 200 / 86400)); mm = 1 - (1 - p['L3']) * (1 - f)
        ax.plot((tt - T0 - ph * P) * 24, (mm - 1) * 1e3, color=c, lw=1.5, label=name)
        if mref is None: mref = mm
    ax.set_title(f'{lab} eclipse (SPOC68 depth scale)', loc='left', fontsize=10); ax.set_ylabel('flux - 1 (ppt)')
    ax.set_xlim(-2.0, 2.0)
    if col == 0: ax.legend(fontsize=6.5, loc='lower left')
    ax = axs[1, col]
    for name, p, c in models:
        f = M.eb(p, t=tt, texp=np.full(tt.size, 200 / 86400)); mm = 1 - (1 - p['L3']) * (1 - f)
        ax.plot((tt - T0 - ph * P) * 24, (mm - mref) * 1e6, color=c, lw=1.4)
    ax.errorbar(xb, (yb - np.interp(xb, (tt - T0 - ph * P) * 24, mref)) * 1e6, eb * 1e6, fmt='o', ms=2.5, color=INK2, ecolor=GRID, lw=0.8)
    ax.set_ylabel('minus best (ppm)'); ax.set_xlabel('hours from conjunction'); ax.set_ylim(-700, 700)
fig.tight_layout(); fig.savefig(f'{V}/plots/fold_models.png', dpi=150); plt.close(fig)
# ---------------- 2. profiles (real + injections)
fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
xs = sorted(pl); axs[0].plot(xs, [pl[x]['chi2'] - best['chi2'] for x in xs], 'o-', color=INK, ms=4, lw=2, label='real data')
for tag, c in (('injA', C3), ('injB', C1), ('injC', C5)):
    fn = f'{V}/results/profile_L3_{tag}.json'
    if not os.path.exists(fn): continue
    pr = json.load(open(fn))['profile']; tr = json.load(open(f'{V}/results/truth_{tag}.json'))['p']
    c0 = min(r['chi2'] for r in pr)
    axs[0].plot([r['val'] for r in pr], [r['chi2'] - c0 for r in pr], 's--', color=c, ms=3, lw=1.2, label=f'injection {tag[-1]} (true L3 = {tr["L3"]:.3f})')
    axs[0].axvline(tr['L3'], color=c, lw=1, ls=':')
xs = sorted(pk); axs[1].plot(xs, [pk[x]['chi2'] - best['chi2'] for x in xs], 'o-', color=INK, ms=4, lw=2)
for ax, nm in ((axs[0], 'L3 (fraction of TESS flux not from the EB)'), (axs[1], 'k = R2/R1')):
    ax.axhline(9, color=C2, lw=1, ls='--'); ax.axhline(4, color=C4, lw=1, ls=':'); ax.set_ylim(-1, 40); ax.set_xlabel(nm); ax.set_ylabel('profile delta chi2')
axs[0].legend(fontsize=7, loc='upper left'); axs[1].set_xscale('log')
fig.tight_layout(); fig.savefig(f'{V}/plots/profiles.png', dpi=150); plt.close(fig)
# ---------------- 3. O-C
F = json.load(open(f'{V}/results/etvfit_default.json'))
ev = F['events_used']; lin = F['fits']['lin']['coef']
def yr_of(n): return 2000 + (lin[0] + n * lin[1] + 2457000 - 2451545.0) / 365.25
fig, ax = plt.subplots(figsize=(9.5, 3.8))
for ph, c, lab in ((0.0, C1, 'primary'), (0.5, C2, 'secondary')):
    e = [x for x in ev if x['ph'] == ph]
    ax.errorbar(yr_of(np.array([x['n'] for x in e])) + (0.04 if ph else 0), [x['oc_min'] for x in e], [x['sig_min'] for x in e], fmt='o', ms=2, color=c, ecolor=c, alpha=0.35, lw=0.6, label=f'{lab} (per event)')
for k, v in F['per_sector'].items():
    if not k.endswith('both'): continue
    ax.errorbar(yr_of(v['mean_epoch']), v['oc_min'], v['err_min'], fmt='s', ms=6, color=INK, ecolor=INK, lw=1.4)
q = F['quadratic']; cul = abs(q['c_min_per_yr2']) + 1.96 * q['c_err_min_per_yr2']
xx = np.linspace(2018.4, 2026.4, 200); ym = np.mean([yr_of(x['n']) for x in ev])
for sgn in (1, -1):
    q_ = sgn * cul * (xx - ym) ** 2; q_ -= np.polyval(np.polyfit(xx, q_, 1), xx)
    ax.plot(xx, q_, color=C3, ls='--', lw=1.2)
ax.text(2018.5, 7.5, f'black squares: per-sector means; green dashed: quadratic O-C at the 95% limit |c| = {cul:.2f} min/yr$^2$ (linear part removed)', fontsize=7, color=INK2)
ax.set_ylim(-15, 10); ax.set_xlabel('year'); ax.set_ylabel('O - C (min)'); ax.legend(fontsize=7, loc='lower right')
fig.tight_layout(); fig.savefig(f'{V}/plots/etv.png', dpi=150); plt.close(fig)
# ---------------- 4. companion exclusion
import glob
fig, ax = plt.subplots(figsize=(6.4, 3.6))
fn = sorted(glob.glob(f'{V}/results/companion_MEB*_M31.03_d9_e0.8.json'))
if fn:
    g = json.load(open(fn[0]))['grid']; a = [x['a_au'] for x in g]
    ax.semilogx(a, [x['frac_excluded'] for x in g], 'o-', ms=3, color=C1, lw=1.6, label='eclipse timings (delta chi2 > 9)')
    ax.semilogx(a, [x['frac_excluded_gaiaRV'] for x in g], 's-', ms=3, color=C4, lw=1.2, label='Gaia DR3 RV constancy')
    ax.semilogx(a, [x['frac_excluded_ETV_or_gaiaRV'] for x in g], '-', color=INK, lw=2, label='either')
    ax.axvspan(297 * 1.0, 297 * 1.8, color=C3, alpha=0.15, lw=0)
    ax.text(310, 0.85, 'projected 1.0-1.8"\nat 297 pc', fontsize=7, color=INK2)
ax.set_xlabel('outer semi-major axis a (au); M3 = 1.03 Msun, e < 0.8, random orientation'); ax.set_ylabel('fraction of orbits excluded'); ax.legend(fontsize=7, loc='lower left')
fig.tight_layout(); fig.savefig(f'{V}/plots/companion.png', dpi=150); plt.close(fig)
# ---------------- 5. out-of-eclipse modulation folded at P_EB per sector
import importlib.util
spec = importlib.util.spec_from_file_location('ooe', f'{V}/scripts/e01_ooe.py'); sys.argv = [sys.argv[0], '3.0']
O = importlib.util.module_from_spec(spec); spec.loader.exec_module(O)
fig, axs = plt.subplots(2, 3, figsize=(10, 5.2), sharey=True)
for ax, (lab, lc) in zip(axs.ravel(), [('SPOC PDCSAP S1', D.load_spoc(1)), ('SPOC PDCSAP S28', D.load_spoc(28)), ('SPOC PDCSAP S68', D.load_spoc(68)),
                                      ('TESScut S95', D.load_tc(95)), ('TESScut S105', D.load_tc(105)), ('TESScut S106', D.load_tc(106))]):
    t, f = O.prep(D.clip(lc)); B = O.spline_basis(t, 3.0)
    cc, *_ = np.linalg.lstsq(B, f, rcond=None); r = f / (B @ cc) - 1
    ph = ((t - D.T00) / D.P0) % 1
    bins = np.linspace(0, 1, 26); idx = np.digitize(ph, bins)
    xb = [ph[idx == i].mean() for i in range(1, 26) if (idx == i).sum() > 3]; yb = [r[idx == i].mean() * 1e6 for i in range(1, 26) if (idx == i).sum() > 3]
    ax.plot(ph, r * 1e6, ',', color=GRID); ax.plot(xb, yb, 'o', ms=3, color=C1)
    ax.axvspan(0, 0.03, color=C2, alpha=0.15, lw=0); ax.axvspan(0.97, 1, color=C2, alpha=0.15, lw=0); ax.axvspan(0.47, 0.53, color=C2, alpha=0.15, lw=0)
    ax.set_title(lab, fontsize=9, loc='left'); ax.set_ylim(-1500, 1500)
for ax in axs[1]: ax.set_xlabel('phase (P_EB = 1.30367 d)')
for ax in axs[:, 0]: ax.set_ylabel('out-of-eclipse flux (ppm)')
fig.tight_layout(); fig.savefig(f'{V}/plots/ooe_fold.png', dpi=150); plt.close(fig)
print('plots written')
