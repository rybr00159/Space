"""Transit shape in the raw-pixel aperture light curves (s03 'mask'), own code.
Exposure-integrated trapezoid: depth, T14, flat fraction f = T23/T14, mid-time offset. Per-sector white errors from
the out-of-transit robust scatter; a red-noise factor beta from 300 fake epochs (single-epoch depth scatter vs
formal error) per sector multiplies the errors.  Profiles chi2 in f and T14; V shape = f 0.
Also: odd/even depths and the phase-0.5 depth at fixed adopted shape.
Excludes S104 E62 (gap edge; rejected by both pipelines).
Output: out/s05_shape_{ap}_k{kn}.json
"""
import sys, numpy as np
sys.path.insert(0, 'tic214311908/analysis/panel_systematics/scripts')
from ps_lib import ROOT, jdump, T0, P, trap
from scipy.optimize import minimize

ap = sys.argv[1] if len(sys.argv) > 1 else 'box'
kn = sys.argv[2] if len(sys.argv) > 2 else '0.75'
d = np.load(f'{ROOT}/out/s03_lc_mask_k{kn}.npz')
t, f, s = d['t'], d[ap], d['sec']
ok = np.isfinite(f)
t, f, s = t[ok], f[ok], s[ok]
E = np.round((t - T0) / P).astype(int)
ok = ~((s == 104) & (E == 62))
t, f, s, E = t[ok], f[ok], s[ok], E[ok]
dtc = {12: 30 / 1440, 39: 10 / 1440}
expo = np.array([dtc.get(int(x), 200 / 86400) for x in s])
ph = t - (T0 + E * P)
win = np.abs(ph) < 0.30
sig = np.zeros(len(t)); beta = {}
rng = np.random.default_rng(3)
for sec in np.unique(s):
    m = s == sec
    oot = m & (np.abs(ph) > 0.1)
    r = 1.4826 * np.median(np.abs(f[oot] - np.median(f[oot])))
    sig[m] = r
    # red noise: fake epochs (phase 0.15-0.85), single-window box depth (adopted shape) vs formal
    tm, fm = t[m], f[m]
    dd, ee = [], []
    for k in range(300):
        tc = tm.min() + rng.uniform(0, tm.max() - tm.min())
        phs = ((tc - T0) / P) % 1
        if phs < 0.12 or phs > 0.88:
            continue
        w = np.abs(tm - tc) < 0.30
        if w.sum() < (12 if sec == 12 else 20):
            continue
        mod = trap(tm[w], tc, dtc.get(int(sec), 200 / 86400))
        if (mod > 0.5).sum() < (2 if sec == 12 else 5):
            continue
        X = np.column_stack([np.ones(w.sum()), tm[w] - tc, -mod])
        b, *_ = np.linalg.lstsq(X, fm[w], rcond=None)
        cov = np.linalg.inv(X.T @ X) * r ** 2
        dd.append(b[2]); ee.append(np.sqrt(cov[2, 2]))
    dd, ee = np.array(dd), np.array(ee)
    beta[int(sec)] = float(1.4826 * np.median(np.abs(dd - np.median(dd))) / np.median(ee))
    print(sec, "rms", round(r, 2), "beta", beta[int(sec)], len(dd))
    sig[m] *= max(beta[int(sec)], 1.0) if np.isfinite(beta[int(sec)]) else 1.5
tw, fw, sw, Ew, xw = t[win], f[win], sig[win], E[win], expo[win]
ue = np.unique(Ew)


def model(p, tt=tw, EE=Ew, xx=xw):
    dep, t14, fl, dtm = p
    m = np.zeros(len(tt))
    for xe in np.unique(xx):
        k = xx == xe
        m[k] = trap(tt[k], T0 + EE[k] * P + dtm, xe, t14=t14, flat=fl)
    return -dep * m


def chi2(p, fix=None):
    dep, t14, fl, dtm = p
    if not (0 < t14 < 0.2 and 0 <= fl < 1):
        return 1e12
    m = model(p)
    # per-window linear baseline profiled
    r = fw - m
    c = 0.0
    for e in ue:
        k = Ew == e
        X = np.column_stack([np.ones(k.sum()), tw[k] - tw[k].mean()]) / sw[k][:, None]
        y = r[k] / sw[k]
        b, *_ = np.linalg.lstsq(X, y, rcond=None)
        c += np.sum((y - X @ b) ** 2)
    return c


best = None
for t14 in (1.4 / 24, 1.7 / 24, 2.0 / 24):
    for fl in (0.3, 0.7):
        o = minimize(chi2, [15, t14, fl, 0.0], method='Nelder-Mead', options=dict(xatol=1e-5, fatol=1e-3, maxiter=4000))
        if best is None or o.fun < best.fun:
            best = o
out = dict(ap=ap, knots=kn, beta=beta, n_windows=int(len(ue)), best=dict(depth_ppt=best.x[0], T14_h=best.x[1] * 24, flat=best.x[2],
                                                                             dt_min=best.x[3] * 1440, chi2=best.fun, npts=int(len(tw))))
prof = {}
for fl in (0.0, 0.15, 0.3, 0.45, 0.6, 0.72, 0.85):
    o = min((minimize(lambda q: chi2([q[0], q[1], fl, q[2]]), [15, t0, 0.0], method='Nelder-Mead',
                      options=dict(xatol=1e-5, fatol=1e-3, maxiter=3000)) for t0 in (1.5 / 24, 1.8 / 24, 2.2 / 24)),
            key=lambda o: o.fun)
    prof[fl] = dict(dchi2=o.fun - best.fun, depth=o.x[0], T14_h=o.x[1] * 24)
out['profile_flat'] = prof
pt = {}
for t14h in (1.3, 1.5, 1.7, 1.9, 2.1):
    o = min((minimize(lambda q: chi2([q[0], t14h / 24, q[1], q[2]]), [15, f0, 0.0], method='Nelder-Mead',
                      options=dict(xatol=1e-5, fatol=1e-3, maxiter=3000)) for f0 in (0.3, 0.7)), key=lambda o: o.fun)
    pt[t14h] = dict(dchi2=o.fun - best.fun, depth=o.x[0], flat=o.x[1])
out['profile_T14'] = pt
# adopted shape: depth, odd/even
def depth_fixed(sel):
    X = []
    m = model([1.0, 1.70 / 24, 0.72, 0.0])
    cols, y, wts = [], [], []
    A = np.zeros((sel.sum(), 1 + 2 * len(np.unique(Ew[sel]))))
    for j, e in enumerate(np.unique(Ew[sel])):
        k = Ew[sel] == e
        A[k, 1 + 2 * j] = 1; A[k, 2 + 2 * j] = tw[sel][k] - tw[sel][k].mean()
    A[:, 0] = m[sel]
    Aw = A / sw[sel][:, None]; yw = fw[sel] / sw[sel]
    b, *_ = np.linalg.lstsq(Aw, yw, rcond=None)
    cov = np.linalg.inv(Aw.T @ Aw)
    return float(-b[0]), float(np.sqrt(cov[0, 0]))
out['depth_adopted_shape'] = depth_fixed(np.ones(len(tw), bool))
od = depth_fixed(Ew % 2 == 1); ev = depth_fixed(Ew % 2 == 0)
out['odd'] = od; out['even'] = ev; out['odd_minus_even'] = [od[0] - ev[0], float(np.hypot(od[1], ev[1]))]
print(ap, kn, 'beta', {k: round(v, 2) for k, v in beta.items()})
print('best', {k: round(v, 3) for k, v in out['best'].items()})
print('profile flat', {k: round(v['dchi2'], 1) for k, v in prof.items()})
print('profile T14', {k: round(v['dchi2'], 1) for k, v in pt.items()})
print('depth adopted shape', out['depth_adopted_shape'], 'odd-even', out['odd_minus_even'])
jdump(out, f'{ROOT}/out/s05_shape_{ap}_k{kn}.json')
