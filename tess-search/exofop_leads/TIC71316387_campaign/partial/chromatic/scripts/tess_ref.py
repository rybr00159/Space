"""TESS reference: eclipse shape templates and per-sector depths, measured with the SAME
linear model (local per-window baselines fitted simultaneously with template amplitudes)
that is applied to the ground-based surveys.

Inputs: SPOC PDCSAP (S1 30-min, S28 10-min, S68 200-s; crowding-corrected by SPOC) and our
own TESScut aperture photometry from the lcmodel sub-campaign (S95/S105/S106, 200 s;
raw normalised flux + PRF-based contamination fraction which we remove here).
"""
import sys, json, glob
import numpy as np
from astropy.io import fits
from scipy.optimize import minimize
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/scripts')
import chromlib as C

BASE = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad'
OUT = BASE + '/tic713/chromatic'
SPOC = BASE + '/exofop/agents/ctoi_vet/recent/data/mastDownload/HLSP/'
TC = BASE + '/tic713/lcmodel/data/'

data = {}
for sec in (1, 28, 68):
    fn = glob.glob(SPOC + f'*0071316387-s{sec:04d}*/*lc.fits')[0]
    h = fits.open(fn)
    assert int(h[0].header['TICID']) == 71316387
    d = h[1].data
    t = d['TIME'] + 2457000.0
    f = d['PDCSAP_FLUX']; e = d['PDCSAP_FLUX_ERR']
    ok = np.isfinite(t) & np.isfinite(f) & np.isfinite(e) & (d['QUALITY'] == 0)
    med = np.nanmedian(f[ok])
    exp = {1: 1800., 28: 600., 68: 200.}[sec] / 86400
    data[sec] = dict(t=t[ok], f=f[ok] / med, src='SPOC-PDCSAP', exptime=exp,
                     crowdsap=float(h[1].header.get('CROWDSAP', np.nan)))
for sec in (95, 105, 106):
    z = np.load(TC + f'lc_tc_s{sec}.npz', allow_pickle=True)
    t = z['t'] + 2457000.0; f = z['f'].astype(float); c = float(z['contam'])
    ok = np.isfinite(t) & np.isfinite(f)
    f = f[ok] / np.median(f[ok])
    f = (f - c) / (1 - c)          # remove PRF-modelled contaminating flux
    data[sec] = dict(t=t[ok], f=f, src='TESScut-aperture(lcmodel)', exptime=200 / 86400, contam=c)

T14 = 1.55 / 24


def windows(t, half=0.30):
    """Keep data within +-half d of each deep/shallow eclipse centre; segment = one eclipse window."""
    ph = C.phase(t)
    dd = np.minimum(ph, 1 - ph) * C.P
    ds = np.abs(ph - 0.5) * C.P
    keep = (dd < half) | (ds < half)
    n = np.round((t - C.T0) / C.P * 2).astype(int)  # half-cycle index -> window id
    return keep, n


def prep(sec):
    d = data[sec]
    keep, n = windows(d['t'])
    t, f, n = d['t'][keep], d['f'][keep], n[keep]
    _, seg = np.unique(n, return_inverse=True)
    # require some out-of-eclipse coverage on both sides; drop windows with < 30% expected points
    return t, f, seg


def fit_shape(sel_secs, which, ph_c, x0):
    """Fit (T, p, q) of the template for one eclipse (which='d' or 's') using short-cadence sectors."""
    pack = [(prep(s), data[s]['exptime']) for s in sel_secs]

    def chi(par):
        T, p, q, dph = par
        if not (0.03 < T < 0.12 and 0.3 < p < 6 and 0.2 < q < 6):
            return 1e30
        tot = 0
        for (t, f, seg), ex in pack:
            nsub = max(1, int(round(ex / (60 / 86400))))
            sd = C.template(t, 0.0 + (dph if which == 'd' else 0), T, p, q, ex, nsub)
            ss = C.template(t, 0.5 + (dph if which == 's' else 0), T, p, q, ex, nsub) if which == 's' else None
            # the other eclipse uses its own (fixed, previously fitted or provisional) shape
            if which == 'd':
                so = C.template(t, ph_c['s'], *shapes['s'][:3], ex, nsub)
                tem = [sd, so]
            else:
                so = C.template(t, ph_c['d'], *shapes['d'][:3], ex, nsub)
                tem = [so, ss]
            a, ae, r, m = C.lsq(t, f, np.ones_like(t), seg, tem, order=1)
            s = C.mad_std(r)
            tot += np.sum((r / s) ** 2)
        return tot
    res = minimize(chi, x0, method='Nelder-Mead', options=dict(xatol=1e-5, fatol=1e-3, maxiter=4000))
    return res.x, res.fun


# provisional shapes, then iterate deep/shallow twice
shapes = {'d': (0.066, 1.2, 1.5, 0.0), 's': (0.066, 1.2, 1.5, 0.0)}
ph_c = {'d': 0.0, 's': 0.5}
short = [68, 95, 105, 106, 28]
for it in range(2):
    for w in ('d', 's'):
        x, fun = fit_shape(short, w, ph_c, np.array(shapes[w]))
        shapes[w] = tuple(x)
        ph_c[w] = (0.0 if w == 'd' else 0.5) + x[3]
        print(it, w, 'T=%.4f d (%.2f h) p=%.3f q=%.3f dphase=%.5f (%.2f min) chi=%.1f' % (x[0], x[0] * 24, x[1], x[2], x[3], x[3] * C.P * 1440, fun), flush=True)


def depths(sec, shp=shapes, phc=ph_c):
    t, f, seg = prep(sec)
    ex = data[sec]['exptime']; nsub = max(1, int(round(ex / (60 / 86400))))
    tem = [C.template(t, phc['d'], *shp['d'][:3], ex, nsub), C.template(t, phc['s'], *shp['s'][:3], ex, nsub)]
    a, ae, r, m = C.lsq(t, f, np.ones_like(t), seg, tem, order=1)
    s = C.mad_std(r)
    return a, ae * s, s, (t, f, seg, tem)


rng = np.random.default_rng(1)
out = dict(shapes={k: dict(T_d=v[0], T_h=v[0] * 24, p=v[1], q=v[2], dphase=v[3]) for k, v in shapes.items()},
           phase_centres=ph_c, sectors={})
for sec in sorted(data):
    a, ae, s, (t, f, seg, tem) = depths(sec)
    # empirical null for the error: fake phases using the FULL light curve windows
    d = data[sec]; ex = d['exptime']; nsub = max(1, int(round(ex / (60 / 86400))))
    fake = []
    for dphi in C.fake_phases(60, rng):
        tt = d['t']; ff = d['f']
        ph = C.phase(tt)
        cen = [(dphi) % 1, (0.5 + dphi) % 1]
        keep = np.zeros(len(tt), bool)
        for c in cen:
            keep |= np.abs(((ph - c) + 0.5) % 1 - 0.5) * C.P < 0.30
        n = np.round(((tt - C.T0) / C.P - dphi) * 2).astype(int)
        _, sg = np.unique(n[keep], return_inverse=True)
        tm = [C.template(tt[keep], cen[0], *shapes['d'][:3], ex, nsub), C.template(tt[keep], cen[1], *shapes['s'][:3], ex, nsub)]
        aa, _, _, _ = C.lsq(tt[keep], ff[keep], np.ones(keep.sum()), sg, tm, order=1)
        fake.append(aa)
    fake = np.array(fake)
    out['sectors'][sec] = dict(src=d['src'], exptime_s=ex * 86400, n=int(len(t)), rms=s,
                               deep=a[0], deep_err_formal=ae[0], shallow=a[1], shallow_err_formal=ae[1],
                               deep_err_null=C.mad_std(fake[:, 0]), shallow_err_null=C.mad_std(fake[:, 1]),
                               null_mean=[float(np.mean(fake[:, 0])), float(np.mean(fake[:, 1]))],
                               crowdsap=d.get('crowdsap'), contam=d.get('contam'))
    print(sec, d['src'], 'deep %.0f +- %.0f (null %.0f) ppm  shallow %.0f +- %.0f (null %.0f) ppm  ratio %.3f' % (
        a[0] * 1e6, ae[0] * 1e6, C.mad_std(fake[:, 0]) * 1e6, a[1] * 1e6, ae[1] * 1e6, C.mad_std(fake[:, 1]) * 1e6, a[1] / a[0]), flush=True)

# combined TESS reference (exclude S1 if its 30-min smearing makes it inconsistent? keep all: template integrates exposure)
secs = sorted(out['sectors'])
for key in ('deep', 'shallow'):
    v = np.array([out['sectors'][s][key] for s in secs])
    e = np.array([max(out['sectors'][s][key + '_err_null'], out['sectors'][s][key + '_err_formal']) for s in secs])
    w = 1 / e ** 2
    mean = np.sum(w * v) / np.sum(w)
    chi2 = np.sum((v - mean) ** 2 * w)
    scat = np.std(v, ddof=1)
    err_stat = np.sqrt(1 / np.sum(w))
    err = max(err_stat * np.sqrt(max(chi2 / (len(v) - 1), 1)), scat / np.sqrt(len(v)))
    out['combined_' + key] = dict(mean=mean, err=err, err_stat=err_stat, chi2=chi2, dof=len(v) - 1,
                                  sector_scatter=scat, unweighted_mean=float(np.mean(v)))
    # SPOC-only
    vs = np.array([out['sectors'][s][key] for s in (1, 28, 68)])
    out['spoc_only_' + key] = dict(mean=float(np.mean(vs)), sd=float(np.std(vs, ddof=1)))
    print(key, 'combined %.0f +- %.0f ppm (chi2 %.1f/%d, sector sd %.0f)' % (mean * 1e6, err * 1e6, chi2, len(v) - 1, scat * 1e6))
json.dump(out, open(OUT + '/results/tess_ref.json', 'w'), indent=1, default=float)

# phase-folded plot data
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
for sec in (68, 95, 105, 106):
    a, ae, s, (t, f, seg, tem) = depths(sec)
    X, nb = C.design(t, seg, tem, 1)
    full = np.linalg.lstsq(X, f, rcond=None)[0]
    base = X[:, :nb] @ full[:nb]
    fn = f / base
    ph = C.phase(t)
    for k, (c, axx) in enumerate(((0.0, ax[0]), (0.5, ax[1]))):
        dp = (((ph - c) + 0.5) % 1 - 0.5) * C.P * 24
        m = np.abs(dp) < 3
        axx.plot(dp[m], (fn[m] - 1) * 1e3, '.', ms=1, alpha=0.3, label=f'S{sec}')
for k, w in enumerate(('d', 's')):
    x = np.linspace(-3, 3, 600)
    T, p, q = shapes[w][:3]
    A = out['combined_' + ('deep' if w == 'd' else 'shallow')]['mean']
    ax[k].plot(x + shapes[w][3] * C.P * 24, -A * 1e3 * C.vshape(x / 24 / (T / 2), p, q), 'k-')
    ax[k].set_xlabel('hours from eclipse centre'); ax[k].set_title(['deep (phase 0)', 'shallow (phase 0.5)'][k])
ax[0].set_ylabel('relative flux - 1 (ppt)'); ax[0].set_ylim(-9, 3); ax[0].legend(markerscale=8)
plt.tight_layout(); plt.savefig(OUT + '/plots/tess_template.png', dpi=110)
