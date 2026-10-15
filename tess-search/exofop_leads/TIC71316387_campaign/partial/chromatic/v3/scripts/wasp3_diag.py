"""SuperWASP diagnostics (v3):
  1. airmass-based per-night baseline (c0 + c1*X [+ c2*X^2]) instead of polynomial in time;
  2. END-TO-END injection-recovery: a physically shaped eclipse pair (K7V+M2V model light curve in the
     WASP band from predict3, scaled) is injected into the RAW per-point flux at fake phases, then the
     full pipeline (normalisation, iterative clipping, weights, fit) is re-run -> tests clipping,
     template mismatch and normalisation (the cheap paired injection in run_survey only tests linearity);
  3. per-event amplitudes (each eclipse with its own night baseline) -> chi2 / outlier events;
  4. time-of-night / airmass distribution of in-eclipse points for deep vs shallow;
  5. phase-folded binned plot.
"""
import sys, json
import numpy as np
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/v3/scripts')
import lib3 as L
import wasp3 as W

T = json.load(open(L.V3 + '/results/tess3.json'))
E3 = json.load(open(L.V3 + '/results/ephem3.json'))
shp = L.survey_shapes(T['shapes'], E3)
sod = T['s_over_d']
EXPT = W.EXPT
out = {}


def airmass(t_bjd):
    from astropy.time import Time
    from astropy.coordinates import EarthLocation, AltAz
    import astropy.units as u
    loc = EarthLocation(lat=-32.3797 * u.deg, lon=20.8107 * u.deg, height=1798 * u.m)
    tt = Time(t_bjd - 70 / 86400.0, format='jd', scale='utc')     # ~UTC (sub-minute accuracy is plenty)
    aa = L._sky().transform_to(AltAz(obstime=tt, location=loc))
    return np.asarray(aa.secz), np.asarray(aa.alt.deg)


t, f, seg, cam, hjd = W.load(W.TARGET)
X, alt = airmass(t)
print('airmass range %.2f..%.2f; altitude min %.1f' % (X.min(), X.max(), alt.min()))

# ---- 1. airmass baseline: run the survey machinery with the regressor = airmass (Projector builds the
#         polynomial in whatever variable it is handed); templates still use time.
def run_x(t, f, seg, xvar, order, nfake=500, seed=3):
    rng = np.random.default_rng(seed)
    keep, sg, w = L.clean_and_weight(t, f, seg, shp, EXPT, 1)
    t, f, x = t[keep], f[keep], xvar[keep]
    pr = L.Projector(x, w, sg, order)
    a, e = L.measure_at(t, f, w, pr, shp, EXPT, 0.0)
    fk = []
    for d in L.fake_offsets(nfake, rng):
        r = L.measure_at(t, f, w, pr, shp, EXPT, d)
        if r is not None:
            fk.append(r[0])
    fk = np.array(fk)
    return dict(deep=float(a[0]), shallow=float(a[1]), deep_err=float(fk[:, 0].std(ddof=1)), shallow_err=float(fk[:, 1].std(ddof=1)),
                deep_null_mean=float(fk[:, 0].mean()), shallow_null_mean=float(fk[:, 1].mean()), n=int(len(t)))


for order in (1, 2):
    r = run_x(t, f, seg, X, order)
    out[f'airmass_o{order}'] = r
    print('TAM airmass-baseline order %d: deep %.0f +- %.0f  shallow %.0f +- %.0f  (null means %.0f, %.0f)' % (
        order, r['deep'] * 1e6, r['deep_err'] * 1e6, r['shallow'] * 1e6, r['shallow_err'] * 1e6, r['deep_null_mean'] * 1e6, r['shallow_null_mean'] * 1e6), flush=True)

# ---- 4. time-of-night / airmass of in-eclipse points
sd = L.tmpl(t, shp, 'd', 0, EXPT); ss = L.tmpl(t, shp, 's', 0, EXPT)
for nm, s in (('deep', sd), ('shallow', ss)):
    m = s > 0.5
    pos = np.zeros(len(t))
    for k in np.unique(seg):
        mm = seg == k
        pos[mm] = (t[mm] - t[mm].min()) / max(t[mm].max() - t[mm].min(), 1e-3)
    out[f'{nm}_inecl_airmass_median'] = float(np.median(X[m]))
    out[f'{nm}_inecl_nightpos_hist'] = np.histogram(pos[m], bins=5, range=(0, 1))[0].tolist()
    print('%-7s in-eclipse points: median airmass %.2f; fraction of night position histogram (5 bins) %s' % (
        nm, np.median(X[m]), np.histogram(pos[m], bins=5, range=(0, 1))[0].tolist()))
print('all points: median airmass %.2f' % np.median(X))

# ---- 3. per-event amplitudes (cleaned default data)
z = np.load(L.V3 + '/data/wasp3_default_clean.npz')
tc, fc, wc, sgc = z['t'], z['f'], z['w'], z['seg']
ev = {}
for which in 'ds':
    s = L.tmpl(tc, shp, which, 0, EXPT)
    j = shp[which]
    cyc = np.round((tc - j['t0']) / j['per'] - j['phc']).astype(int)
    rows = []
    for c in np.unique(cyc[s > 0.5]):
        segs = np.unique(sgc[(cyc == c) & (s > 0.05)])
        for k in segs:
            m = sgc == k
            if (s[m] > 0.5).sum() < 5:
                continue
            a, e, r, _ = L.fit_amps(tc[m], fc[m], wc[m], np.zeros(m.sum(), int), [s[m]], 1)
            # empirical error: formal x robust residual scale
            sc = L.mad_std(r * np.sqrt(wc[m]))
            rows.append((int(c), int(k), float(a[0]), float(e[0] * max(sc, 1)), int((s[m] > 0.5).sum()),
                         float((s[m] * (s[m] > 0.05)).sum())))
    rows = np.array(rows)
    wmean = np.sum(rows[:, 2] / rows[:, 3] ** 2) / np.sum(1 / rows[:, 3] ** 2)
    chi = np.sum(((rows[:, 2] - wmean) / rows[:, 3]) ** 2)
    ev[which] = dict(events=rows.tolist(), wmean=float(wmean), werr=float(1 / np.sqrt(np.sum(1 / rows[:, 3] ** 2))), chi2=float(chi), n=len(rows),
                     median=float(np.median(rows[:, 2])))
    # jackknife over events
    jk = []
    for i in range(len(rows)):
        mm = np.arange(len(rows)) != i
        jk.append(np.sum(rows[mm, 2] / rows[mm, 3] ** 2) / np.sum(1 / rows[mm, 3] ** 2))
    jk = np.array(jk)
    ev[which]['jackknife_sd'] = float(np.sqrt((len(jk) - 1) / len(jk) * np.sum((jk - jk.mean()) ** 2)))
    ev[which]['most_influential'] = rows[np.argsort(np.abs(jk - wmean))[-3:], :3].tolist()
    print('%s per-event (single-night) amplitudes: N=%d weighted mean %.0f +- %.0f (jackknife %.0f) chi2 %.1f/%d median %.0f' % (
        which, len(rows), wmean * 1e6, ev[which]['werr'] * 1e6, ev[which]['jackknife_sd'] * 1e6, chi, len(rows) - 1, np.median(rows[:, 2]) * 1e6))
out['events'] = ev

# ---- 2. end-to-end injection-recovery with a physical eclipse shape
P3 = json.load(open(L.V3 + '/results/predict3.json'))
pr = P3['pairs']['K7V+M2V']
import predict3 as PR
import sed3 as S
tab = S.mamajek()
tg = np.linspace(-0.12, 0.12, 2401)
r1, r2 = tab['K7V'], tab['M2V']
ld1 = S.ld_quadratic(r1['Teff'], 4.65, 'WASP'); ld2 = S.ld_quadratic(r2['Teff'], 4.8, 'WASP')
k1, k2 = PR.lightcurves(pr['R1'], pr['R2'], pr['a_Rsun'], pr['inc'], ld1, ld2, tg)
prof_d = k1 / k1.max(); prof_s = k2 / k2.max()
# the template amplitude of these unit-peak profiles
amp_d = PR.template_amp(tg, prof_d, shp, 'd'); amp_s = PR.template_amp(tg, prof_s, shp, 's')
print('physical K7+M2 WASP-band profiles: template amplitude / true peak = %.3f (deep) %.3f (shallow)' % (amp_d, amp_s))
out['phys_template_amp_over_peak'] = [float(amp_d), float(amp_s)]
rng = np.random.default_rng(21)
raw_t, raw_f, raw_seg, _, _ = W.load(W.TARGET)


def inject(tt, ff, dphi, Dd, Ds):
    j = shp['d']
    dtd = ((tt - j['t0']) / j['per'] - j['phc'] - dphi + 0.5) % 1 - 0.5
    dts = ((tt - j['t0']) / j['per'] - shp['s']['phc'] - dphi + 0.5) % 1 - 0.5
    pd_ = np.interp(dtd * j['per'], tg, prof_d, left=0, right=0)
    ps_ = np.interp(dts * j['per'], tg, prof_s, left=0, right=0)
    return ff * (1 - Dd * pd_ - Ds * ps_)


keep0, sg0, w0 = L.clean_and_weight(raw_t, raw_f, raw_seg, shp, EXPT, 1)
tt0, ff0 = raw_t[keep0], raw_f[keep0]
prj0 = L.Projector(tt0, w0, sg0, 1)
for D in (0.0025, 0.0060):
    rec = []
    for dphi in L.fake_offsets(60, rng):
        fi = inject(raw_t, raw_f, dphi, D, D * 0.30)
        keep, sg, w = L.clean_and_weight(raw_t, fi, raw_seg, shp, EXPT, 1)
        tt, ff = raw_t[keep], fi[keep]
        prj = L.Projector(tt, w, sg, 1)
        a1, _ = L.measure_at(tt, ff, w, prj, shp, EXPT, dphi)
        # same phase without injection (paired) using the same cleaning of uninjected data
        a0, _ = L.measure_at(tt0, ff0, w0, prj0, shp, EXPT, dphi)
        rec.append((a1[0], a1[1], a0[0], a0[1]))
    rec = np.array(rec)
    gd = (rec[:, 0] - rec[:, 2]) / (D * amp_d)
    gs = (rec[:, 1] - rec[:, 3]) / (0.30 * D * amp_s)
    out[f'e2e_inj_{D}'] = dict(peak_depth=D, expected_template_amp=D * amp_d, rec_mean=float(rec[:, 0].mean()), rec_sd=float(rec[:, 0].std(ddof=1)),
                               gain_deep=float(gd.mean()), gain_deep_sem=float(gd.std(ddof=1) / np.sqrt(len(gd))),
                               gain_shallow=float(gs.mean()), gain_shallow_sem=float(gs.std(ddof=1) / np.sqrt(len(gs))),
                               frac_recovered_3sig=float(np.mean(rec[:, 0] > 3 * 571e-6)))
    print('E2E injection peak %.4f (template amp %.0f ppm): recovered %.0f +- %.0f ppm; gain deep %.3f +- %.3f, shallow %.3f +- %.3f; frac > 3 sigma %.2f' % (
        D, D * amp_d * 1e6, rec[:, 0].mean() * 1e6, rec[:, 0].std(ddof=1) * 1e6, gd.mean(), gd.std(ddof=1) / np.sqrt(len(gd)),
        gs.mean(), gs.std(ddof=1) / np.sqrt(len(gs)), np.mean(rec[:, 0] > 3 * 571e-6)), flush=True)

# ---- 5. folded binned light curve (default cleaned data)
s_d = L.tmpl(tc, shp, 'd', 0, EXPT); s_s = L.tmpl(tc, shp, 's', 0, EXPT)
a, e, r, prj = L.fit_amps(tc, fc, wc, sgc, [s_d, s_s], 1)
base = fc - r + a[0] * (-s_d) * -1 + a[1] * (-s_s) * -1     # model = base - a s  =>  base = (f - r) + a s
fold = {}
for nm, which in (('deep', 'd'), ('shallow', 's')):
    j = shp[which]
    x = (((tc - j['t0']) / j['per'] - j['phc'] + 0.5) % 1 - 0.5) * j['per'] * 24
    m = np.abs(x) < 4
    y = (fc / base - 1)[m]; ww = wc[m]
    bins = np.arange(-4, 4.001, 0.5)
    idx = np.digitize(x[m], bins)
    rows = []
    for i in range(1, len(bins)):
        s = idx == i
        if s.sum() > 10:
            mu = np.sum(ww[s] * y[s]) / np.sum(ww[s])
            rows.append((0.5 * (bins[i - 1] + bins[i]), mu, 1 / np.sqrt(np.sum(ww[s])) * max(1, L.mad_std(y[s] * np.sqrt(ww[s]))), int(s.sum())))
    fold[nm] = rows
out['fold'] = fold
json.dump(out, open(L.V3 + '/results/wasp3_diag.json', 'w'), indent=1)
print('saved')
