"""Template-shape search (generalised V, (1-|x|^p)^q; deep and shallow independent) on the 200-s sectors
(S68 SPOC, S95/S105/S106 TESScut) with FIXED per-sector white-noise sigma, so chi2 is comparable between
shapes.  Projectors (quarter-orbit chunk baselines) are precomputed once; for each trial shape only the
template columns are recomputed.  Multi-start Nelder-Mead over (T, p, q, dphase) per eclipse, plus a
coarse (p, q) grid with T and dphase optimised, to map the chi2 surface."""
import sys, json
import numpy as np
from scipy.optimize import minimize
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/v3/scripts')
import lib3 as L
import tess3 as TS

data = TS.load_all()
SECS = (68, 95, 105, 106)
pre = {}
for s in SECS:
    d = data[s]
    oot = (np.abs(L.dphase(d['t'], 0)) > 0.08) & (np.abs(L.dphase(d['t'], 0.5)) > 0.08)
    sig = L.mad_std(np.diff(d['f'][oot])) / np.sqrt(2)
    sh0 = {'d': dict(T=0.065, p=1.5, q=1.0, phc=0.0), 's': dict(T=0.065, p=1.5, q=1.0, phc=0.5)}
    a, e, r = TS.fit_sector(d, sh0)
    k = np.abs(r) < 5 * L.mad_std(r)
    t, f = d['t'][k], d['f'][k]
    w = np.ones(k.sum()) / sig ** 2
    pr = L.Projector(t, w, TS.chunks(t), 1)
    fp = pr(f)
    pre[s] = dict(t=t, w=w, pr=pr, fp=fp, ex=d['ex'], c0=float(np.sum(w * fp * fp)))
    print('S%d sigma %.0f ppm n %d' % (s, sig * 1e6, k.sum()), flush=True)


def col(t, ex, T, p, q, phc):
    return L.template(t, phc, T, p, q, ex)


def chi2_pair(sd, ss):
    tot = 0.0
    for s in SECS:
        P = pre[s]
        S = [P['pr'](-col(P['t'], P['ex'], *sd)), P['pr'](-col(P['t'], P['ex'], *ss))]
        G = np.array([[np.sum(P['w'] * S[i] * S[j]) for j in range(2)] for i in range(2)])
        b = np.array([np.sum(P['w'] * S[i] * P['fp']) for i in range(2)])
        tot += P['c0'] - b @ np.linalg.solve(G, b)
    return tot


cur = {'d': (0.060, 1.7, 1.0, 0.0), 's': (0.065, 1.7, 1.0, 0.5)}


def fit_one(which, start):
    def f(x):
        T, p, q, dph = x
        if not (0.03 < T < 0.25 and 0.4 < p < 12 and 0.2 < q < 60 and abs(dph) < 0.01):
            return 1e30
        sh = dict(cur); sh[which] = (T, p, q, (0.0 if which == 'd' else 0.5) + dph)
        return chi2_pair(sh['d'], sh['s'])
    r = minimize(f, start, method='Nelder-Mead', options=dict(xatol=1e-6, fatol=1e-3, maxiter=3000, maxfev=3000))
    r = minimize(f, r.x, method='Nelder-Mead', options=dict(xatol=1e-6, fatol=1e-3, maxiter=3000, maxfev=3000))
    return r


out = {}
for it in range(2):
    for which in 'ds':
        best = None
        for st in ([0.060, 1.7, 1.0, 0.0], [0.07, 2.0, 1.0, 0.0], [0.09, 3.4, 10.0, 0.0], [0.055, 1.2, 0.8, 0.0], [0.08, 1.5, 2.5, 0.0]):
            r = fit_one(which, st)
            print('  it%d %s start %s -> chi2 %.2f  T=%.4f p=%.3f q=%.3f dph=%.5f' % (it, which, st[:3], r.fun, *r.x), flush=True)
            if best is None or r.fun < best.fun:
                best = r
        T, p, q, dph = best.x
        cur[which] = (T, p, q, (0.0 if which == 'd' else 0.5) + dph)
        out[f'{which}_it{it}'] = dict(chi2=float(best.fun), x=best.x.tolist())
        print('iter %d %s BEST chi2 %.2f T=%.4f (%.2f h) p=%.3f q=%.3f dphase %.5f (%.2f min)' % (it, which, best.fun, T, T * 24, p, q, dph, dph * L.P * 1440), flush=True)
T3 = json.load(open(L.V3 + '/results/tess3.json'))['shapes']
c3 = chi2_pair(tuple(T3['d'][k] for k in ('T', 'p', 'q', 'phc')), tuple(T3['s'][k] for k in ('T', 'p', 'q', 'phc')))
cb = chi2_pair(cur['d'], cur['s'])
npts = sum(len(pre[s]['t']) for s in SECS)
print('best chi2 %.2f, tess3 shape chi2 %.2f (delta %.2f), N=%d' % (cb, c3, c3 - cb, npts))
# coarse map of the deep-eclipse (p, q) surface with T, dphase optimised
grid = []
for p in (1.0, 1.5, 2.0, 3.0):
    for q in (0.7, 1.0, 1.6, 3.0, 10.0):
        def f(x):
            T, dph = x
            if not (0.03 < T < 0.25 and abs(dph) < 0.01):
                return 1e30
            return chi2_pair((T, p, q, dph), cur['s'])
        r = minimize(f, [cur['d'][0], cur['d'][3]], method='Nelder-Mead', options=dict(xatol=1e-6, fatol=1e-3))
        grid.append((p, q, r.x[0], r.fun - cb))
        print('   deep grid p=%.1f q=%4.1f T=%.4f dchi2=%.1f' % (p, q, r.x[0], r.fun - cb), flush=True)
best_shapes = {'d': dict(T=cur['d'][0], p=cur['d'][1], q=cur['d'][2], phc=cur['d'][3]),
               's': dict(T=cur['s'][0], p=cur['s'][1], q=cur['s'][2], phc=cur['s'][3])}
json.dump(dict(best=best_shapes, chi2_best=cb, chi2_tess3=c3, npts=npts, fits=out, deep_grid=grid), open(L.V3 + '/results/shape_grid3.json', 'w'), indent=1)
