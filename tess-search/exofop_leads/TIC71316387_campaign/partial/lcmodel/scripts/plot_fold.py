"""Folded, baseline-normalised primary and secondary eclipses (all sectors, default products; each event divided by its
profiled quadratic baseline from the given model) with several model solutions overplotted, plus binned residuals."""
import os
os.environ['NUMBA_NUM_THREADS'] = '1'
import numpy as np, sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/scripts')
import fit3, lclib as L, fitlib as F
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d = fit3.load_data(L.HERE + '/results/data_default.pkl')
sols = json.load(open(sys.argv[1]))   # dict name -> param dict
out = sys.argv[2]
names = list(sols)
def detrended(p):
    m = fit3.model_all(d, p)
    c, coefs = F.profile_chi2(d['y'], m, d['dt'], 1 / d['sig'] ** 2, d['starts'], d['ends'], 2)
    base = np.empty_like(m)
    for e in range(len(d['starts'])):
        a, b = d['starts'][e], d['ends'][e]; dt = d['dt'][a:b]
        base[a:b] = coefs[e, 0] + coefs[e, 1] * dt + coefs[e, 2] * dt ** 2
    yy = d['y'] / base
    # rescale each sector's eclipse signal to the s=1 (SPOC) depth scale so that sectors can be co-added
    for sec, idx, texp, nsub in d['groups']:
        sc = p.get(f's{sec}', 1.0); yy[idx] = 1 - (1 - yy[idx]) / sc
    return yy, m, c
ref = sols[names[0]]
y, m0, c0 = detrended(ref)
ph = ((d['t'] - ref['T0']) / ref['P']) % 1
fig, ax = plt.subplots(2, 2, figsize=(11, 7), gridspec_kw=dict(height_ratios=[3, 1.3]), sharex='col')
cols = ['C3', 'C0', 'C2', 'C1', 'C4']
for j, cen in enumerate([0.0, 0.5]):
    x = ((ph - cen + 0.5) % 1 - 0.5) * ref['P'] * 24
    sel = np.abs(x) < 3.5
    # restrict to 200-s and 600-s sectors for the binned display (S1 30-min shown only through the model fit)
    hi = sel & (d['sector'] != 1)
    bins = np.arange(-3.5, 3.51, 0.125); idx = np.digitize(x[hi], bins)
    bx = np.array([x[hi][idx == k].mean() for k in range(1, len(bins)) if (idx == k).sum() > 5])
    by = np.array([y[hi][idx == k].mean() for k in range(1, len(bins)) if (idx == k).sum() > 5])
    be = np.array([y[hi][idx == k].std() / np.sqrt((idx == k).sum()) for k in range(1, len(bins)) if (idx == k).sum() > 5])
    ax[0, j].errorbar(bx, (by - 1) * 1e6, be * 1e6, fmt='o', color='k', ms=3, label='S28/68/95/105/106 binned (7.5 min), depth-scaled to s=1')
    tt = np.linspace(-3.5, 3.5, 600) / 24 + ref['T0'] + cen * ref['P']
    for k, n in enumerate(names):
        p = sols[n]
        u = fit3.q2u
        from ebmodel import eb_flux
        u1a, u2a = u(p['q1a'], p['q2a']); u1b, u2b = u(p['q1b'], p['q2b'])
        cosi = p['bn'] * p['rsum']
        mm = eb_flux(tt, p['T0'], p['P'], p['rsum'], p['k'], cosi, p['J'], p['L3'], p['ecosw'], p['esinw'], u1a, u2a, u1b, u2b, 200 / 86400, 2, fit3.XG, fit3.WG)
        ax[0, j].plot((tt - ref['T0'] - cen * ref['P']) * 24, (mm - 1) * 1e6, color=cols[k], lw=1.3, label=n)
        mi = np.interp(bx / 24 + ref['T0'] + cen * ref['P'], tt, mm)
        ax[1, j].plot(bx, (by - mi) * 1e6, 'o-', color=cols[k], ms=2.5, lw=0.8)
    ax[1, j].axhline(0, color='0.5'); ax[1, j].set_xlabel('hours from conjunction'); ax[0, j].set_title(['primary eclipse', 'secondary eclipse'][j])
    ax[1, j].set_ylim(-600, 600)
ax[0, 0].set_ylabel('flux - 1 (ppm)'); ax[1, 0].set_ylabel('data - model (ppm)'); ax[0, 0].legend(fontsize=7, loc='lower left')
plt.tight_layout(); plt.savefig(out, dpi=110)
print('saved', out)
