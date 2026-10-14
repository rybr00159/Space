"""TESS reference for the chromatic test (re-derivation, v2).

* Data: SPOC/TESS-SPOC PDCSAP S1 (1800 s), S28 (600 s), S68 (200 s): crowding-corrected by SPOC
  (CROWDSAP), so depths are relative to the flux of the target (+ anything unresolved within it).
  Cross-check only: TESScut aperture LCs of S95/S105/S106 from the lcmodel sub-campaign (their PRF
  contamination estimate is removed here); not used in the adopted reference.
* Windows of +-0.30 d around every deep and shallow eclipse, each with its own linear baseline,
  fitted simultaneously with the template amplitudes.
* Shape: V template (1-|x|^p)^q, fitted jointly to S28+S68 for each eclipse (T, p, q, phase offset).
* Errors: fake-phase null (windows around fake eclipse pairs, REAL eclipses masked +-0.045 d),
  plus beta scaling; adopted = max(null sd, formal*beta).
"""
import sys, json, glob
import numpy as np
from astropy.io import fits
from scipy.optimize import minimize
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/v2/scripts')
import lib2 as L

SPOC = L.BASE + '/exofop/agents/ctoi_vet/recent/data/mastDownload/HLSP/'
TC = L.BASE + '/tic713/lcmodel/data/'
HALF = 0.30
data = {}
for sec in (1, 28, 68):
    fn = glob.glob(SPOC + f'*0071316387-s{sec:04d}*/*lc.fits')[0]
    with fits.open(fn) as h:
        assert int(h[0].header['TICID']) == 71316387
        d = h[1].data
        t = d['TIME'] + 2457000.0
        f = d['PDCSAP_FLUX']
        ok = np.isfinite(t) & np.isfinite(f) & (d['QUALITY'] == 0)
        data[sec] = dict(t=t[ok], f=f[ok] / np.nanmedian(f[ok]), exptime={1: 1800., 28: 600., 68: 200.}[sec] / 86400,
                         src='SPOC PDCSAP', crowdsap=float(h[1].header['CROWDSAP']), flfrcsap=float(h[1].header['FLFRCSAP']))
for sec in (95, 105, 106):
    z = np.load(TC + f'lc_tc_s{sec}.npz', allow_pickle=True)
    t = z['t'] + 2457000.0; f = z['f'].astype(float); c = float(z['contam'])
    ok = np.isfinite(t) & np.isfinite(f)
    f = f[ok] / np.median(f[ok]); f = (f - c) / (1 - c)
    data[sec] = dict(t=t[ok], f=f, exptime=200 / 86400, src='TESScut aperture (lcmodel), contam removed', contam=c)


def windows(t, centres, half=HALF, mask_real=False):
    """keep points within +-half of the given phase centres; segment id = (centre idx, cycle)."""
    keep = np.zeros(len(t), bool); seg = np.full(len(t), -1)
    for j, c in enumerate(centres):
        ph = (t - L.T0) / L.P - c
        n = np.round(ph)
        dt = np.abs(ph - n) * L.P
        m = dt < half
        keep |= m
        seg[m] = (n[m].astype(int) * 2 + j)
    if mask_real:
        keep &= (L.dist_to_phase(t, 0.0) > 0.045) & (L.dist_to_phase(t, 0.5) > 0.045)
    _, sg = np.unique(seg[keep], return_inverse=True)
    return keep, sg


def shp_from(par):
    return {'d': dict(T=par[0], p=par[1], q=par[2], phc=par[3]), 's': dict(T=par[4], p=par[5], q=par[6], phc=0.5 + par[7])}


prep = {}
for sec in data:
    d = data[sec]
    k, sg = windows(d['t'], (0.0, 0.5))
    prep[sec] = (d['t'][k], d['f'][k], sg, d['exptime'])


SIG = {}
for s_ in prep:
    t_, f_, sg_, ex_ = prep[s_]
    oot = (L.dist_to_phase(t_, 0) > 0.06) & (L.dist_to_phase(t_, 0.5) > 0.06)
    SIG[s_] = L.mad_std(np.diff(f_[oot])) / np.sqrt(2)   # fixed per-sector white-noise sigma
print('per-sector sigma (ppm):', {k: round(v * 1e6) for k, v in SIG.items()})


def chi2(par, secs=(28, 68, 95, 105, 106)):
    T1, p1, q1, d1, T2, p2, q2, d2 = par
    if not (0.03 < T1 < 0.12 and 0.03 < T2 < 0.12 and 0.3 < p1 < 12 and 0.3 < p2 < 12 and 0.2 < q1 < 12 and 0.2 < q2 < 12
            and abs(d1) < 0.01 and abs(d2) < 0.01):
        return 1e30
    shp = shp_from(par)
    tot = 0
    for s in secs:
        t, f, sg, ex = prep[s]
        a, ae, r = L.measure(t, f, np.ones_like(t), sg, shp, ex, 1)
        tot += np.sum(r ** 2) / SIG[s] ** 2
    return tot


x0 = np.array([0.064, 1.2, 1.5, 0.0, 0.064, 1.2, 1.5, 0.0])
best = None
for start in (np.array([0.070, 3.0, 3.0, 0, 0.072, 3.0, 2.5, 0]), np.array([0.064, 1.3, 1.6, 0, 0.066, 1.3, 1.5, 0]), np.array([0.075, 4.0, 4.0, 0, 0.075, 4.0, 4.0, 0])):
    r = minimize(chi2, start, method='Nelder-Mead', options=dict(xatol=1e-6, fatol=1e-3, maxiter=8000, maxfev=8000))
    r = minimize(chi2, r.x, method='Nelder-Mead', options=dict(xatol=1e-6, fatol=1e-3, maxiter=8000, maxfev=8000))
    print('start', start[:3], '-> chi2 %.2f' % r.fun, np.round(r.x, 5), flush=True)
    if best is None or r.fun < best.fun:
        best = r
shp = shp_from(best.x)
for k in ('d', 's'):
    print(k, 'T=%.4f d (%.3f h) p=%.3f q=%.3f phase centre %.5f (%.2f min from nominal)' % (
        shp[k]['T'], shp[k]['T'] * 24, shp[k]['p'], shp[k]['q'], shp[k]['phc'], (shp[k]['phc'] - (0 if k == 'd' else 0.5)) * L.P * 1440))
out = dict(shapes=shp, shape_fit_chi2=best.fun, sectors={})

rng = np.random.default_rng(7)
for sec in sorted(data):
    t, f, sg, ex = prep[sec]
    a, ae, r = L.measure(t, f, np.ones_like(t), sg, shp, ex, 1)
    s = L.mad_std(r)
    aef = ae * s
    d = data[sec]
    fk, fke = [], []
    for dphi in L.fake_phases(150, rng):
        k, sgf = windows(d['t'], (dphi, 0.5 + dphi), mask_real=True)
        tt, ff = d['t'][k], d['f'][k]
        tem = [L.tmpl(tt, shp, 'd', dphi, ex), L.tmpl(tt, shp, 's', dphi, ex)]
        if not (np.sum(tem[0] > 0.3) > 5 and np.sum(tem[1] > 0.3) > 5):
            continue
        aa, ee, rr, _ = L.lsq(tt, ff, np.ones_like(tt), sgf, tem, 1)
        fk.append(aa); fke.append(ee * L.mad_std(rr))
    fk, fke = np.array(fk), np.array(fke)
    beta = np.std(fk / fke, axis=0, ddof=1)
    row = dict(src=d['src'], exptime_s=ex * 86400, n=int(len(t)), rms=float(s), deep=float(a[0]), shallow=float(a[1]),
               deep_err_formal=float(aef[0]), shallow_err_formal=float(aef[1]),
               null_mean=fk.mean(0).tolist(), null_sd=fk.std(0, ddof=1).tolist(), beta=beta.tolist(), n_fake=len(fk),
               crowdsap=d.get('crowdsap'), contam=d.get('contam'))
    row['deep_err'] = float(max(row['null_sd'][0], beta[0] * aef[0]))
    row['shallow_err'] = float(max(row['null_sd'][1], beta[1] * aef[1]))
    row['ratio_s_d'] = float(a[1] / a[0])
    out['sectors'][str(sec)] = row
    print('S%-3d %-40s deep %5.0f +- %3.0f (formal %3.0f, null sd %3.0f, null mean %4.0f, beta %.2f)  shallow %5.0f +- %3.0f (null mean %4.0f)  s/d %.3f' % (
        sec, d['src'][:40], a[0] * 1e6, row['deep_err'] * 1e6, aef[0] * 1e6, row['null_sd'][0] * 1e6, row['null_mean'][0] * 1e6, beta[0],
        a[1] * 1e6, row['shallow_err'] * 1e6, row['null_mean'][1] * 1e6, a[1] / a[0]), flush=True)


def combine(secs, key):
    v = np.array([out['sectors'][str(s)][key] for s in secs])
    e = np.array([out['sectors'][str(s)][key + '_err'] for s in secs])
    w = 1 / e ** 2
    m = np.sum(w * v) / np.sum(w)
    chi = np.sum(w * (v - m) ** 2)
    es = np.sqrt(1 / np.sum(w))
    err = es * np.sqrt(max(1.0, chi / (len(v) - 1)))
    return dict(mean=float(m), err=float(err), err_stat=float(es), chi2=float(chi), dof=len(v) - 1,
                sd=float(np.std(v, ddof=1)), unweighted=float(np.mean(v)), values=v.tolist(), errs=e.tolist())


for nm, secs in (('spoc', (1, 28, 68)), ('all', (1, 28, 68, 95, 105, 106)), ('tesscut', (95, 105, 106))):
    for key in ('deep', 'shallow'):
        c = combine(secs, key)
        out[f'{nm}_{key}'] = c
        print('%-8s %-8s mean %5.0f +- %3.0f ppm (chi2 %.1f/%d; sector sd %.0f)' % (nm, key, c['mean'] * 1e6, c['err'] * 1e6, c['chi2'], c['dof'], c['sd'] * 1e6))
# adopted reference: SPOC sectors; systematic floor = half the SPOC-vs-TESScut difference
for key in ('deep', 'shallow'):
    sp, tc = out['spoc_' + key], out['tesscut_' + key]
    sysf = 0.5 * abs(sp['mean'] - tc['mean'])
    out['adopted_' + key] = dict(mean=sp['mean'], err=float(np.hypot(sp['err'], sysf)), err_stat=sp['err'], sys_spoc_vs_tesscut=float(sysf))
    print('ADOPTED %s: %.0f +- %.0f ppm' % (key, sp['mean'] * 1e6, np.hypot(sp['err'], sysf) * 1e6))
json.dump(out, open(L.V2 + '/results/tess_ref2.json', 'w'), indent=1)

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
for sec in (28, 68, 95, 105, 106):
    t, f, sg, ex = prep[sec]
    tem = [L.tmpl(t, shp, 'd', 0, ex), L.tmpl(t, shp, 's', 0, ex)]
    a, ae, r, m = L.lsq(t, f, np.ones_like(t), sg, tem, 1)
    base = m + a[0] * tem[0] + a[1] * tem[1]
    fn = f / base
    for k, c in enumerate((0.0, 0.5)):
        dp = ((((t - L.T0) / L.P - c) + 0.5) % 1 - 0.5) * L.P * 24
        mm = np.abs(dp) < 3
        bins = np.arange(-3, 3.01, 0.1)
        idx = np.digitize(dp[mm], bins)
        bc = [0.5 * (bins[i - 1] + bins[i]) for i in range(1, len(bins)) if (idx == i).sum() > 3]
        bm = [np.mean(fn[mm][idx == i] - 1) for i in range(1, len(bins)) if (idx == i).sum() > 3]
        ax[k].plot(bc, np.array(bm) * 1e3, '.-', ms=3, lw=0.5, label=f'S{sec}')
for k, w in enumerate(('d', 's')):
    x = np.linspace(-3, 3, 600)
    A = out['spoc_' + ('deep' if w == 'd' else 'shallow')]['mean']
    ph0 = (shp[w]['phc'] - (0 if w == 'd' else 0.5)) * L.P * 24
    ax[k].plot(x, -A * 1e3 * L.C.vshape((x - ph0) / 24 / (shp[w]['T'] / 2), shp[w]['p'], shp[w]['q']), 'k-', lw=2, label='template (SPOC mean)')
    ax[k].set_xlabel('hours from nominal eclipse centre'); ax[k].set_title(['deep (phase 0)', 'shallow (phase 0.5)'][k])
ax[0].set_ylabel('relative flux - 1 (ppt), 6-min bins'); ax[0].legend(fontsize=7)
plt.tight_layout(); plt.savefig(L.V2 + '/plots/tess_ref2.png', dpi=110)
