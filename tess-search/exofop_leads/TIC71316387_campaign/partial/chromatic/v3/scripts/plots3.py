"""Figures for the v3 chromatic study."""
import sys, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/v3/scripts')
import lib3 as L

BLUE, ORANGE, AQUA, GREY, INK, MUTED = '#2a78d6', '#eb6834', '#1baf7a', '#9a9893', '#0b0b0b', '#52514e'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED, 'ytick.color': MUTED,
                     'axes.spines.top': False, 'axes.spines.right': False, 'grid.color': '#e4e3df', 'grid.linewidth': 0.6})
R = lambda f: json.load(open(L.V3 + '/results/' + f))
T3, P3, W3, WD, X3, C3, N3 = R('tess3.json'), R('predict3.json'), R('wasp3.json'), R('wasp3_diag.json'), R('wise3.json'), R('combine3.json'), R('null3.json')
E3 = R('ephem3.json')
shp = L.survey_shapes(T3['shapes'], E3)
AT = np.array([T3['adopted_deep']['mean'], T3['adopted_shallow']['mean']])
D = C3['measurements']['WASP']['D']
kb = C3['joint']['best_bound']
pb = P3['pairs'][kb]

# ---------------------------------------------------------------- 1. WASP folded
z = np.load(L.V3 + '/data/wasp3_default_clean.npz')
t, f, w, sg = z['t'], z['f'], z['w'], z['seg']
sd, ss = L.tmpl(t, shp, 'd', 0, 30 / 86400), L.tmpl(t, shp, 's', 0, 30 / 86400)
a, e, r, pr = L.fit_amps(t, f, w, sg, [sd, ss], 1)
base = f - r + a[0] * sd + a[1] * ss
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
for k, (which, nm) in enumerate((('d', 'deep eclipse (phase 0)'), ('s', 'shallow eclipse (phase 0.5)'))):
    j = shp[which]
    x = (((t - j['t0']) / j['per'] - j['phc'] + 0.5) % 1 - 0.5) * j['per'] * 24
    m = np.abs(x) < 4
    y = (f / base - 1)[m]; ww = w[m]
    bins = np.arange(-4, 4.001, 0.4)
    idx = np.digitize(x[m], bins)
    bc, bm, be = [], [], []
    for i in range(1, len(bins)):
        s = idx == i
        if s.sum() > 20:
            mu = np.sum(ww[s] * y[s]) / np.sum(ww[s])
            chi = L.mad_std((y[s] - mu) * np.sqrt(ww[s]))
            bc.append(0.5 * (bins[i - 1] + bins[i])); bm.append(mu)
            be.append(max(chi, 1.0) / np.sqrt(np.sum(ww[s])))
    ax[k].errorbar(bc, np.array(bm) * 1e3, np.array(be) * 1e3, fmt='o', ms=4, color=BLUE, ecolor=BLUE, elinewidth=1, label='SuperWASP 2006-07 (24-min bins)')
    xx = np.linspace(-4, 4, 800)
    prof = L.template(xx / 24.0, 0.0, j['T'], j['p'], j['q'], t0=0.0, per=L.P)
    i = 0 if which == 'd' else 1
    ax[k].plot(xx, -a[i] * prof * 1e3, color=BLUE, lw=1.5, label='template fit: %.0f ppm' % (a[i] * 1e6))
    rb = pb['ratio_deep' if which == 'd' else 'ratio_shallow']['WASP']
    ax[k].plot(xx, -rb * AT[i] / (1 + D) * prof * 1e3, color=ORANGE, lw=2, ls='--', label='bound %s prediction' % kb.replace('V', ''))
    rg = P3['background']['5900']['mid']['WASP']
    ax[k].plot(xx, -rg * AT[i] / (1 + D) * prof * 1e3, color=AQUA, lw=2, ls=':', label='G-type eclipsed star (5900 K)')
    ax[k].plot(xx, -AT[i] / (1 + D) * prof * 1e3, color=GREY, lw=1, label='TESS depth (undiluted by WASP)')
    ax[k].axhline(0, color='#cfcdc8', lw=0.8)
    ax[k].set_title(nm, fontsize=9, color=INK)
    ax[k].set_xlabel('hours from predicted eclipse centre')
ax[0].set_ylabel('relative flux - 1 (ppt)')
h, l = ax[0].get_legend_handles_labels()
fig.legend(h, l, fontsize=7.5, frameon=False, loc='lower center', ncol=5)
fig.suptitle('TIC 71316387: SuperWASP (2006-07) eclipses vs. predictions; times and shape fixed from TESS', fontsize=9)
plt.tight_layout(rect=(0, 0.08, 1, 1)); plt.savefig(L.V3 + '/plots/wasp_fold.png', dpi=130); plt.close()

# ---------------------------------------------------------------- 2. ratio vs wavelength
lam = {'WASP': 0.55, 'V': 0.55, 'TESS': 0.80, 'W1': 3.35, 'W2': 4.60, 'B': 0.44, 'R': 0.65, 'I': 0.80, 'Ks': 2.16, 'g': 0.47}
fig, ax = plt.subplots(1, 2, figsize=(10, 4), sharey=False)
bands = ['B', 'g', 'WASP', 'R', 'I', 'TESS', 'Ks', 'W1', 'W2']
for k, key in enumerate(('ratio_deep', 'ratio_shallow')):
    first = True
    for nm, p in P3['pairs'].items():
        if not p.get('feasible') or not C3['pairs'][nm]['cons']:
            continue
        xs = [lam[b] for b in bands]; ys = [p[key][b] if b != 'TESS' else 1.0 for b in bands]
        o = np.argsort(xs)
        ax[k].plot(np.array(xs)[o], np.array(ys)[o], color=ORANGE, lw=1.2 if nm != kb else 2.2, alpha=0.5 if nm != kb else 1,
                   label=('grazing bound pairs consistent with TESS s/d (K5-M0 + K7-M3)' if first else None))
        first = False
    EXT = json.load(open(L.V3 + '/results/ext_pairs3.json'))
    firstx = True
    for nm in ('M0V+M5V', 'M1V+M4.5V', 'M0.5V+M4.5V'):
        e_ = EXT[nm]
        bb = ['WASP', 'I', 'TESS', 'Ks', 'W1', 'W2']
        xs_ = [lam[b] for b in bb]; ys_ = [e_[key][b] for b in bb]
        o = np.argsort(xs_)
        ax[k].plot(np.array(xs_)[o], np.array(ys_)[o], color='#b3401a', lw=1.6, ls='--',
                   label=('bound pairs preferred by the TESS shape (M0-M1 + M4.5-M5, ~300 pc)' if firstx else None))
        firstx = False
    # background G/F envelope (eclipsed star 5600-7000 K for the deep eclipse; its TESS-s/d-implied secondary for the shallow)
    TEg = np.array(sorted(int(t) for t in P3['background']))
    def gen(te, b, g):
        y = np.array([P3['background'][str(t)][g][b] for t in TEg])
        return float(np.exp(np.interp(te, TEg, np.log(y))))
    t2map = {m['Teff1']: m['Teff2'] for m in C3['bg_models'].values()}
    xs = sorted(lam[b] for b in bands)
    order_b = [b for b in sorted(bands, key=lambda b: lam[b])]
    lo, hi = [], []
    for b in order_b:
        vals = []
        for t1 in (5600, 5900, 6200, 6600, 7000):
            te = t1 if k == 0 else t2map[t1]
            for g in ('grazing', 'mid', 'deep'):
                vals.append(1.0 if b == 'TESS' else gen(te, b, g))
        lo.append(min(vals)); hi.append(max(vals))
    ax[k].fill_between([lam[b] for b in order_b], lo, hi, color=AQUA, alpha=0.25, lw=0, label='background EB, G/F eclipsed star (5600-7000 K)')
    for nm, mk, col in (('WASP', 'o', BLUE), ('NEOWISE_W1', 's', BLUE), ('NEOWISE_W2', 'D', BLUE), ('ASASSN_V', '^', GREY)):
        mm = C3['measurements'][nm]
        x = lam[mm['band']] * (1.04 if nm == 'ASASSN_V' else 1)
        ax[k].errorbar([x], [mm['ratio'][k]], [mm['ratio_err'][k]], fmt=mk, ms=6, color=col, ecolor=col, elinewidth=1.2, capsize=0,
                       label=nm.replace('_', ' ') if k == 0 else None)
        ax[k].annotate(nm.split('_')[-1], (x, mm['ratio'][k]), xytext=(6, 4), textcoords='offset points', fontsize=7, color=MUTED)
    ax[k].plot([0.80], [1.0], 'k*', ms=8)
    ax[k].set_xscale('log')
    from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter, NullLocator
    ax[k].xaxis.set_minor_locator(NullLocator()); ax[k].xaxis.set_minor_formatter(NullFormatter())
    ax[k].xaxis.set_major_locator(FixedLocator([0.45, 0.55, 0.8, 1, 2, 3.35, 4.6]))
    ax[k].xaxis.set_major_formatter(FixedFormatter(['0.45', '0.55', '0.8', '1', '2', '3.4', '4.6']))
    ax[k].set_xlabel('wavelength (micron)'); ax[k].set_ylabel('depth(band) / depth(TESS)')
    ax[k].set_title(['deep eclipse', 'shallow eclipse'][k], fontsize=9)
    ax[k].axhline(0, color='#cfcdc8', lw=0.8); ax[k].grid(True, axis='y')
    ax[k].set_ylim(-1.5, 7.5)
ax[0].legend(fontsize=7, frameon=False, loc='upper left')
plt.tight_layout(); plt.savefig(L.V3 + '/plots/ratio_vs_band.png', dpi=130); plt.close()

# ---------------------------------------------------------------- 3. expected measurement distributions under each hypothesis
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
for k, (nm, key, mkey) in enumerate((('SuperWASP deep eclipse', 'WASP', 'WASP'), ('NEOWISE W1 deep eclipse', 'W1', 'NEOWISE_W1'))):
    mm = C3['measurements'][mkey]
    fac = (1 + mm['D']) / AT[0]
    tot = np.sqrt(np.array(mm['cov'])[0, 0]) / fac          # total error in amplitude units
    smp = np.array(N3[key]['samples'])[:, 0]
    zs = (smp - smp.mean()) / smp.std()                         # standardised fake-period null
    band = mm['band']
    hyps = [('bound %s' % kb.replace('V', ''), pb['ratio_deep'][band] / fac, ORANGE),
            ('G/F most favourable (5600 K, grazing)', C3['bg_models']['5600_grazing']['pred'][band][0] / fac, AQUA),
            ('G0 eclipsed star (5900 K)', P3['background']['5900']['mid'][band] / fac, GREY)]
    bins = np.linspace(-3 * tot, max(h[1] for h in hyps) + 3.5 * tot, 70) * 1e6
    for lab, mu, col in hyps:
        ax[k].hist((mu + zs * tot) * 1e6, bins=bins, color=col, alpha=0.35, density=True, label=lab)
    ax[k].axvline(mm['A'][0] * 1e6, color=BLUE, lw=2.2, label='observed')
    ax[k].set_xlabel('deep-eclipse template amplitude (ppm)'); ax[k].set_yticks([]); ax[k].set_title(nm, fontsize=9)
    ax[k].legend(fontsize=7, frameon=False)
fig.suptitle('Expected amplitude under each hypothesis = prediction + empirical (fake-period) null scaled to the total error budget', fontsize=9)
plt.tight_layout(); plt.savefig(L.V3 + '/plots/hypothesis_amplitudes.png', dpi=130); plt.close()

# ---------------------------------------------------------------- 4. Teff of the eclipsed star
TEg = np.array(sorted(int(t) for t in P3['background']))
tg = np.arange(3300, 7001, 25.0)
fig, ax = plt.subplots(figsize=(6.5, 3.6))
comb = np.ones_like(tg)
for mkey, col, ls in (('WASP', BLUE, '-'), ('NEOWISE_W1', ORANGE, '--'), ('NEOWISE_W2', GREY, ':')):
    mm = C3['measurements'][mkey]; b = mm['band']; r = mm['ratio'][0]; e = mm['ratio_err'][0]
    lk = np.zeros_like(tg)
    for g in ('grazing', 'mid', 'deep'):
        y = np.array([P3['background'][str(t)][g][b] for t in TEg])
        pr = np.exp(np.interp(tg, TEg, np.log(y)))
        lk += np.exp(-0.5 * ((r - pr) / e) ** 2) / 3
    comb *= lk
    ax.plot(tg, lk / lk.max(), color=col, ls=ls, lw=1.6, label=mkey.replace('_', ' '))
ax.plot(tg, comb / comb.max(), color=INK, lw=2.4, label='combined')
ax.axvspan(5300, 7000, color=AQUA, alpha=0.15, lw=0, label='G/F eclipsed star')
ax.axvline(4100, color=ORANGE, lw=0.8, alpha=0.6); ax.text(4120, 0.95, 'K7V', fontsize=7, color=MUTED)
ax.axvline(3850, color=ORANGE, lw=0.8, alpha=0.6); ax.text(3700, 0.95, 'M0V', fontsize=7, color=MUTED)
ax.set_xlabel('Teff of the star eclipsed at the deep eclipse (K)'); ax.set_ylabel('relative likelihood'); ax.set_yticks([0, 0.5, 1])
ax.legend(fontsize=7, frameon=False, loc='center right')
plt.tight_layout(); plt.savefig(L.V3 + '/plots/teff_eclipsed_star.png', dpi=130); plt.close()
print('plots written')
