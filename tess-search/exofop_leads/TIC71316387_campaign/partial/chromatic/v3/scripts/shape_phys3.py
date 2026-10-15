"""Does a main-sequence K/M pair reproduce the TESS eclipse SHAPE?  For a given pair (Mamajek radii, masses,
TESS LD), the deep-eclipse profile shape depends only on the inclination; its amplitude scales with the
eclipsed star's light fraction f1.  Fit inclination + amplitude (+ time offset) to the 200-s TESS sectors
(chunk baselines, fixed white-noise weights as in shape_grid3), compare chi2 with the empirical template,
and convert the fitted amplitude into the light fraction f1 and hence the distance at which an MS pair
of that type would sit (f1 scales as d^-2 relative to its value at 297 pc)."""
import sys, json
import numpy as np
from scipy.optimize import minimize_scalar
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/v3/scripts')
import lib3 as L
import tess3 as TS
import predict3 as PR
import sed3 as S

data = TS.load_all()
SECS = (68, 95, 105, 106)
T3 = json.load(open(L.V3 + '/results/tess3.json'))
P3 = json.load(open(L.V3 + '/results/predict3.json'))
shp = T3['shapes']
pre = {}
for s in SECS:
    d = data[s]
    oot = (np.abs(L.dphase(d['t'], 0)) > 0.08) & (np.abs(L.dphase(d['t'], 0.5)) > 0.08)
    sig = L.mad_std(np.diff(d['f'][oot])) / np.sqrt(2)
    a, e, r = TS.fit_sector(d, shp)
    k = np.abs(r) < 5 * L.mad_std(r)
    t, f = d['t'][k], d['f'][k]
    w = np.ones(k.sum()) / sig ** 2
    pr = L.Projector(t, w, TS.chunks(t), 1)
    pre[s] = dict(t=t, w=w, pr=pr, fp=pr(f), ex=d['ex'])
tab = S.mamajek()
tg = np.linspace(-0.15, 0.15, 1201)


def cols_from_profile(t, ex, prof_d, prof_s, dt_d, dt_s):
    out = []
    for prof, phc, dt in ((prof_d, shp['d']['phc'], dt_d), (prof_s, shp['s']['phc'], dt_s)):
        acc = np.zeros(len(t)); nsub = 4
        for o in ((np.arange(nsub) + 0.5) / nsub - 0.5) * ex:
            x = L.dphase(t + o, phc) - dt
            acc += np.interp(x, tg, prof, left=0, right=0)
        out.append(acc / nsub)
    return out


def chi2_cols(colfun):
    tot = 0; amps = []
    for s in SECS:
        P = pre[s]
        cd, cs = colfun(P['t'], P['ex'])
        Sd, Ss = P['pr'](-cd), P['pr'](-cs)
        G = np.array([[np.sum(P['w'] * Sd * Sd), np.sum(P['w'] * Sd * Ss)], [np.sum(P['w'] * Sd * Ss), np.sum(P['w'] * Ss * Ss)]])
        b = np.array([np.sum(P['w'] * Sd * P['fp']), np.sum(P['w'] * Ss * P['fp'])])
        a = np.linalg.solve(G, b)
        tot += np.sum(P['w'] * P['fp'] ** 2) - b @ a
        amps.append(a)
    return tot, np.mean(amps, axis=0)


# reference: empirical template
c_t, a_t = chi2_cols(lambda t, ex: (L.tmpl(t, shp, 'd', 0, ex), L.tmpl(t, shp, 's', 0, ex)))
print('empirical template: chi2 %.1f, amplitudes %.0f / %.0f ppm' % (c_t, a_t[0] * 1e6, a_t[1] * 1e6))
out = dict(template_chi2=c_t, pairs={})
for pair in ('K7V+M2V', 'M0V+M3V', 'K6V+M1V', 'K5V+M1V'):
    p = P3['pairs'][pair]
    r1, r2 = tab[p['spt1']], tab[p['spt2']]
    ld1 = S.ld_quadratic(r1['Teff'], 4.65, 'TESS'); ld2 = S.ld_quadratic(r2['Teff'], 4.8, 'TESS')
    res = []
    imin = np.degrees(np.arccos((p['R1'] + p['R2']) / p['a_Rsun'])) + 0.2
    for inc in np.arange(imin, 90.01, 0.5):
        k1, k2 = PR.lightcurves(p['R1'], p['R2'], p['a_Rsun'], inc, ld1, ld2, tg)
        if k1.max() <= 0:
            continue
        pd_, ps_ = k1 / k1.max(), k2 / k2.max()
        c, a = chi2_cols(lambda t, ex: cols_from_profile(t, ex, pd_, ps_, 0.0, 0.0))
        res.append((inc, c, a[0], a[1], k1.max(), k2.max()))
    res = np.array(res)
    i = np.argmin(res[:, 1])
    inc, c, Ad, As, k1m, k2m = res[i]
    f1 = Ad / k1m; f2 = As / k2m
    fMS = p['frac_TESS'][0]
    dist = 297.4 * np.sqrt(fMS / f1)
    # inclination giving the observed depth at the MS light fraction (the predict3 solution), and its chi2
    j = np.argmin(np.abs(res[:, 2] / res[:, 4] * 0 + res[:, 4] * fMS - a_t[0]))
    out['pairs'][pair] = dict(grid=res.tolist(), best_inc=float(inc), best_chi2=float(c), dchi2_vs_template=float(c - c_t), f1_fit=float(f1),
                              f2_fit=float(f2), f1_MS_297pc=float(fMS), implied_distance_pc=float(dist),
                              inc_depth_matched=float(res[j, 0]), chi2_depth_matched=float(res[j, 1]), dchi2_depth_matched=float(res[j, 1] - c))
    print('%-8s best inc %.1f: chi2 %.1f (template %+.1f); fitted f1 = %.4f vs MS@297pc %.4f -> MS distance %.0f pc; '
          'depth-matched (MS@297pc) inc %.1f has chi2 worse by %.1f' % (pair, inc, c, c - c_t, f1, fMS, dist, res[j, 0], res[j, 1] - c), flush=True)
json.dump(out, open(L.V3 + '/results/shape_phys3.json', 'w'), indent=1)


# ---------------------------------------------------------------- radius-scaling scan with the bound (297 pc) depth constraint
def chi2_fixed_deep(colfun, Ad):
    tot = 0
    for s in SECS:
        P = pre[s]
        cd, cs = colfun(P['t'], P['ex'])
        Sd, Ss = P['pr'](-cd), P['pr'](-cs)
        Gdd, Gds, Gss = np.sum(P['w'] * Sd * Sd), np.sum(P['w'] * Sd * Ss), np.sum(P['w'] * Ss * Ss)
        bd, bs = np.sum(P['w'] * Sd * P['fp']), np.sum(P['w'] * Ss * P['fp'])
        As = (bs - Ad * Gds) / Gss
        tot += np.sum(P['w'] * P['fp'] ** 2) - 2 * Ad * bd - 2 * As * bs + Ad ** 2 * Gdd + 2 * Ad * As * Gds + As ** 2 * Gss
    return tot


scan = {}
for pair in ('K7V+M2V', 'M0V+M3V'):
    p = P3['pairs'][pair]
    r1, r2 = tab[p['spt1']], tab[p['spt2']]
    ld1 = S.ld_quadratic(r1['Teff'], 4.65, 'TESS'); ld2 = S.ld_quadratic(r2['Teff'], 4.8, 'TESS')
    fMS = p['frac_TESS'][0]
    rows = []
    for sc in (0.8, 0.9, 1.0, 1.1, 1.2, 1.35, 1.5, 1.75, 2.0):
        R1, R2 = p['R1'] * sc, p['R2'] * sc
        imin = np.degrees(np.arccos(min((R1 + R2) / p['a_Rsun'], 1))) + 0.1
        best_free, best_bound = (np.inf, None), (np.inf, None)
        for inc in np.arange(imin, 90.01, 0.25):
            k1, k2 = PR.lightcurves(R1, R2, p['a_Rsun'], inc, ld1, ld2, tg)
            if k1.max() <= 0:
                continue
            cf = lambda t, ex: cols_from_profile(t, ex, k1, k2, 0.0, 0.0)
            # free amplitude (in units of the occulted fraction): chi2 and implied f1
            c, a = chi2_cols(lambda t, ex: cols_from_profile(t, ex, k1 / k1.max(), k2 / k2.max(), 0.0, 0.0))
            if c < best_free[0]:
                best_free = (c, inc, a[0] / k1.max())
            # bound at 297 pc with luminosity scaling as R^2 at fixed Teff (colour): f1 = fMS * sc^2
            cb = chi2_fixed_deep(cf, fMS * sc ** 2)
            if cb < best_bound[0]:
                best_bound = (cb, inc)
        rows.append((sc, best_free[0] - c_t, best_free[1], best_free[2], best_bound[0] - c_t, best_bound[1]))
        print('%-8s radii x%.2f: free-amplitude best inc %.2f dchi2(vs template) %+.1f, implied f1 %.4f (bound-MS value %.4f); '
              'bound-at-297pc constrained best inc %.2f dchi2 %+.1f' % (pair, sc, best_free[1], best_free[0] - c_t, best_free[2], fMS * sc ** 2,
                                                                     best_bound[1], best_bound[0] - c_t), flush=True)
    scan[pair] = rows
out['radius_scan'] = scan
json.dump(out, open(L.V3 + '/results/shape_phys3.json', 'w'), indent=1)
