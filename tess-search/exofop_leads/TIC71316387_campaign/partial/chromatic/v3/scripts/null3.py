"""Fake-PERIOD null for SuperWASP and NEOWISE W1/W2 (complements the fake-phase null).

The fake-phase null keeps the real period, so nearby offsets give correlated amplitudes: only
~1/(eclipse width in phase) ~ 30-60 of them are independent, which limits the tail probabilities it can
resolve.  Here each trial uses a random period P' in [0.6, 3.5] d (avoiding +-1.5% around P, P/2, 2P,
3P/2, 2P/3 and 1 d aliases of these) and a random epoch; the deep+shallow template pair is fitted at P'
while the REAL eclipses are modelled at the true ephemeris.  Trials are then (nearly) independent draws
of the noise in an eclipse-shaped matched filter with the survey's own cadence and systematics.
Output: null samples and empirical tail probabilities for the observed amplitudes and for the
amplitudes predicted by each hypothesis.
"""
import sys, json
import numpy as np
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/v3/scripts')
import lib3 as L

T = json.load(open(L.V3 + '/results/tess3.json'))
E3 = json.load(open(L.V3 + '/results/ephem3.json'))
shp = L.survey_shapes(T['shapes'], E3)
rng = np.random.default_rng(77)


def bad_period(p):
    P = L.P
    for m in (1.0, 0.5, 2.0, 1.5, 2 / 3., 3.0, 1 / 3.):
        if abs(p / (P * m) - 1) < 0.015:
            return True
    # 1-day aliases of P
    for k in (-2, -1, 1, 2):
        f = 1 / P + k
        if f > 0 and abs(p * f - 1) < 0.015:
            return True
    return False


def fake_period_null(t, f, w, seg, order, exptime, ntrial=3000):
    pr = L.Projector(t, w, seg, order)
    real = [L.tmpl(t, shp, 'd', 0, exptime), L.tmpl(t, shp, 's', 0, exptime)]
    out = []
    while len(out) < ntrial:
        p = rng.uniform(0.6, 3.5)
        if bad_period(p):
            continue
        t0 = rng.uniform(0, p) + L.T0
        sh = {'d': dict(shp['d'], per=p, t0=t0, phc=0.0), 's': dict(shp['s'], per=p, t0=t0, phc=0.5)}
        cd, cs = L.tmpl(t, sh, 'd', 0, exptime), L.tmpl(t, sh, 's', 0, exptime)
        if (cd > 0.5).sum() < 3 or (cs > 0.5).sum() < 3:
            continue
        a, e, r, _ = L.fit_amps(t, f, w, None, [cd, cs] + real, proj=pr)
        out.append((a[0], a[1], e[0], e[1], p))
    return np.array(out)


res = {}
z = np.load(L.V3 + '/data/wasp3_default_clean.npz')
W = json.load(open(L.V3 + '/results/wasp3.json'))['default']
nl = fake_period_null(z['t'], z['f'], z['w'], z['seg'], 1, 30 / 86400, 3000)
res['WASP'] = dict(samples=nl[:, :2].tolist(), obs=[W['deep'], W['shallow']], formal_real=[W['deep_formal'], W['shallow_formal']])
print('WASP fake-period null: deep mean %.0f sd %.0f (mad %.0f) kurt %.2f; shallow sd %.0f; corr %.2f; formal-err-scaled sd (beta) %.2f/%.2f' % (
    nl[:, 0].mean() * 1e6, nl[:, 0].std() * 1e6, L.mad_std(nl[:, 0]) * 1e6,
    np.mean(((nl[:, 0] - nl[:, 0].mean()) / nl[:, 0].std()) ** 4) - 3, nl[:, 1].std() * 1e6, np.corrcoef(nl[:, 0], nl[:, 1])[0, 1],
    np.std(nl[:, 0] / nl[:, 2]), np.std(nl[:, 1] / nl[:, 3])))
# null sd scaled to the real-phase formal error (coverage-aware) via beta
beta_d, beta_s = np.std(nl[:, 0] / nl[:, 2]), np.std(nl[:, 1] / nl[:, 3])
res['WASP']['beta'] = [float(beta_d), float(beta_s)]
res['WASP']['err_scaled'] = [float(beta_d * W['deep_formal']), float(beta_s * W['shallow_formal'])]
res['WASP']['corr'] = float(np.corrcoef(nl[:, 0], nl[:, 1])[0, 1])
print('   -> beta-scaled errors at the real phase: deep %.0f shallow %.0f ppm' % (beta_d * W['deep_formal'] * 1e6, beta_s * W['shallow_formal'] * 1e6))
zs = (nl[:, 0] / nl[:, 2]); zo = W['deep'] / W['deep_formal']
print('   observed deep S/N(formal) %.2f; fraction of fake-period trials with S/N >= that: %d/%d' % (zo, np.sum(zs >= zo), len(zs)))
res['WASP']['frac_snr_ge_obs'] = float((np.sum(zs >= zo) + 1) / (len(zs) + 1))

# NEOWISE W1 / W2 (strict + AllWISE) -- rebuild the cleaned arrays exactly as wise3.py
import wise3_load as WL
for band in ('w1', 'w2'):
    t, f, w, seg, r = WL.cleaned(band)
    nl = fake_period_null(t, f, w, seg, 0, 7.7 / 86400, 3000)
    beta_d, beta_s = np.std(nl[:, 0] / nl[:, 2]), np.std(nl[:, 1] / nl[:, 3])
    res[band.upper()] = dict(samples=nl[:, :2].tolist(), obs=[r['deep'], r['shallow']], beta=[float(beta_d), float(beta_s)],
                             err_scaled=[float(beta_d * r['deep_formal']), float(beta_s * r['shallow_formal'])],
                             corr=float(np.corrcoef(nl[:, 0], nl[:, 1])[0, 1]))
    zs = nl[:, 0] / nl[:, 2]; zo = r['deep'] / r['deep_formal']
    res[band.upper()]['frac_snr_ge_obs'] = float((np.sum(zs >= zo) + 1) / (len(zs) + 1))
    print('%s fake-period null: deep sd %.0f kurt %.2f; shallow sd %.0f; beta-scaled errors deep %.0f shallow %.0f; obs deep %.0f S/N %.2f, frac S/N>=obs %d/%d' % (
        band.upper(), nl[:, 0].std() * 1e6, np.mean(((nl[:, 0] - nl[:, 0].mean()) / nl[:, 0].std()) ** 4) - 3, nl[:, 1].std() * 1e6,
        beta_d * r['deep_formal'] * 1e6, beta_s * r['shallow_formal'] * 1e6, r['deep'] * 1e6, zo, np.sum(zs >= zo), len(zs)))
json.dump(res, open(L.V3 + '/results/null3.json', 'w'))
print('saved')
