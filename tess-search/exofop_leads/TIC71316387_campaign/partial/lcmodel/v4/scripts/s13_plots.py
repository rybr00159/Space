"""Figures: (1) folded primary/secondary eclipses (baseline-normalised per event, depth scale of each sector applied, binned
7.5 min; 200-s and 600-s sectors only so cadence smearing is small) with three models evaluated on the same cadence mix:
best total-occultation solution, best L3=0 grazing solution, and the L3=0.85 total solution; residuals below.
(2) profile likelihoods: Delta chi2 vs L3 (light curve only; + ellipsoidal limit), and vs r1 (density axis on top)."""
import os
os.environ['NUMBA_NUM_THREADS'] = '2'
import numpy as np, sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4/lib')
import fit3, fitlib as F, lclib as L
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
FD = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4'
res = json.load(open(f'{FD}/results/ml_default_recent.json'))
d = fit3.load_data(f'{FD}/results/data_default.pkl')
def pick(kind, branch, L3=None):
    rr = [x for x in res if x['kind'] == kind and x['branch'] == branch and (L3 is None or abs(x['fixed'].get('L3', -1) - L3) < 1e-9)]
    return min(rr, key=lambda x: x['tot'])
mods = [('best: total occultation, L3=%.3f k=%.2f b1=%.2f', pick('branch', 'total'), 'C3'),
        ('grazing, L3=0 (k=%.2f b1=%.2f)', pick('branch', 'partial'), 'C0'),
        ('total, L3=%.2f fixed (k=%.2f b1=%.2f)', pick('L3prof', 'total', 0.85), 'C2')]
gb = mods[0][1]['tot']
sel = np.isin(d['sector'], [28, 68, 95, 105, 106])
fig = plt.figure(figsize=(12, 8))
gs = fig.add_gridspec(2, 2, height_ratios=[3, 1.3])
for col, ph in enumerate((0.0, 0.5)):
    a0 = fig.add_subplot(gs[0, col]); a1 = fig.add_subplot(gs[1, col], sharex=a0)
    bins = np.arange(-0.2, 0.2001, 7.5 / 1440)
    for mi, (lab, p, c) in enumerate(mods):
        m = fit3.model_all(d, p)
        _, coefs = F.profile_chi2(d['y'], m, d['dt'], 1 / d['sig'] ** 2, d['starts'], d['ends'], 2)
        base = np.empty(d['y'].size)
        for e in range(len(d['starts'])):
            a, b = d['starts'][e], d['ends'][e]; dt = d['dt'][a:b]
            base[a:b] = coefs[e, 0] + coefs[e, 1] * dt + coefs[e, 2] * dt * dt
        scl = np.array([p.get(f's{s}', 1.0) for s in d['sector']])
        yn = 1 - (1 - d['y'] / base) / scl; mn = 1 - (1 - m) / scl
        s2 = sel & (d['ph'] == ph); x = d['dt'][s2]
        idx = np.digitize(x, bins)
        bx = np.array([x[idx == j].mean() for j in range(1, len(bins)) if (idx == j).sum() > 5])
        by = np.array([yn[s2][idx == j].mean() for j in range(1, len(bins)) if (idx == j).sum() > 5])
        be = np.array([yn[s2][idx == j].std() / np.sqrt((idx == j).sum()) for j in range(1, len(bins)) if (idx == j).sum() > 5])
        bm = np.array([mn[s2][idx == j].mean() for j in range(1, len(bins)) if (idx == j).sum() > 5])
        if mi == 0:
            a0.errorbar(bx * 24, (by - 1) * 1e6, be * 1e6, fmt='ko', ms=3, label='S28/68/95/105/106 binned 7.5 min')
        if mi == 0: txt = lab % (p['L3'], p['k'], p['b1'])
        elif mi == 1: txt = lab % (p['k'], p['b1'])
        else: txt = lab % (p['L3'], p['k'], p['b1'])
        txt += '  [dchi2 = %+.1f]' % (p['tot'] - gb)
        a0.plot(bx * 24, (bm - 1) * 1e6, '-', color=c, lw=1.5, label=txt)
        a1.plot(bx * 24, (by - bm) * 1e6, 'o-', color=c, ms=2.5, lw=0.8)
    a0.set_title(['primary eclipse', 'secondary eclipse'][col]); a0.legend(fontsize=7, loc='lower left' if col == 0 else 'lower right')
    a1.axhline(0, color='0.5'); a1.set_xlabel('hours from conjunction'); a1.set_ylim(-500, 500)
    a0.set_xlim(-2.5, 2.5)
    if col == 0: a0.set_ylabel('flux - 1 [ppm] (SPOC depth scale)'); a1.set_ylabel('data - model [ppm]')
plt.tight_layout(); plt.savefig(f'{FD}/plots/fold_models.png', dpi=110)
# profiles
pdv = json.load(open(f'{FD}/results/profile_derived_default_recent.json'))
fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
L3v = np.array([x['L3'] for x in pdv]); dl = np.array([x['dchi2_lc'] for x in pdv]); de = np.array([x['dchi2_ell_sectors'] for x in pdv])
dej = np.array([x['dchi2_ell_joint'] for x in pdv])
ax[0].plot(L3v, dl, 'ko-', label='light curve only (best of both branches)')
ax[0].plot(L3v, dl + de, 's--', color='C1', label='+ ellipsoidal limit (sector-scatter error, q from MS k-q)')
ax[0].plot(L3v, dl + dej, '^:', color='C3', label='+ ellipsoidal limit (fake-period error)')
ax[0].set_ylim(0, 40); ax[0].set_xlabel('L3 (fraction of TESS light not from the EB)'); ax[0].set_ylabel('Delta chi2 (profile)'); ax[0].legend(fontsize=7)
ax[0].axhline(9, color='0.6', ls=':')
for b, c in (('total', 'C3'), ('partial', 'C0')):
    rr = [x for x in res if x['kind'] == 'r1prof' and x['branch'] == b]
    gv = sorted(set(x['fixed']['r1'] for x in rr))
    yv = [min([x['tot'] for x in rr if x['fixed']['r1'] == v]) - gb for v in gv]
    ax[1].plot(gv, yv, 'o-', color=c, label=f'{b} branch')
ax[1].set_ylim(0, 40); ax[1].set_xlabel('r1 = R1/a'); ax[1].set_ylabel('Delta chi2 (profile)'); ax[1].legend(fontsize=8); ax[1].axhline(9, color='0.6', ls=':')
sec = ax[1].secondary_xaxis('top', functions=(lambda r: 3 * np.pi / (6.674e-8 * (L.P0 * 86400) ** 2) / np.maximum(r, 1e-3) ** 3 / 1.4 / 1.408,
                                              lambda rho: (3 * np.pi / (6.674e-8 * (L.P0 * 86400) ** 2) / np.maximum(rho, 1e-3) / 1.4 / 1.408) ** (1 / 3)))
sec.set_xlabel('rho1 [rho_sun] for q = 0.4')
plt.tight_layout(); plt.savefig(f'{FD}/plots/profiles.png', dpi=110)
print('ok')
