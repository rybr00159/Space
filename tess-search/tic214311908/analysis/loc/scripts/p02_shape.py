"""Transit-shape cross-check (input to the physical-plausibility test, item 4).  Matched-filter light curve at the
M3V's apparent position: L(t) = sum_pix w p (y - b) / sum_pix w p^2 over the fit region (p = PRF at the M3V,
w = 1/sigma^2 from the DIFF noise map, b = per-cadence background plane already removed), i.e. the flux (e-/s) of a
point source at the M3V.  Per transit window (+-0.45 d) a quadratic baseline; common trapezoid model (depth,
T14, flat fraction f = T23/T14, global time shift), exposure-integrated.  Profile chi2 over f and over T14.
Noise: per-window OOT scatter, inflated by the binned-residual (red-noise) factor at 1.6 h.
Output: out/shape_profile.json, figs/shape_fold.png"""
import sys
sys.path.insert(0, 'tic214311908/analysis/loc/scripts')
from l2fit import *
from scipy.optimize import least_squares
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

segs = []
for sec in SECTORS:
    M = SectorModel(sec, 'DIFF'); Dd = M.Dd
    i3 = int(np.where(M.hyp.name == 'M3V')[0][0]); p = M.H[i3]; w = M.w
    L = (Dd.Y[:, M.mask] @ (w * p)) / (w * p * p).sum()
    for tc in real_tcs(Dd.t):
        m = np.abs(Dd.t - tc) < 0.45
        if m.sum() < 20: continue
        oot = m & (np.abs(Dd.t - tc) > 0.07)
        if (Dd.t[oot] < tc).sum() < 5 or (Dd.t[oot] > tc).sum() < 5: continue
        if (np.abs(Dd.t[m] - tc) < 0.03).sum() < 1: continue
        segs.append(dict(sec=sec, tc=tc, t=Dd.t[m], L=L[m], exp=Dd.exptime))
print(len(segs), 'transit windows')


def model_seg(sg, depth, t14, f, dt):
    tpl = template(sg['t'], [sg['tc'] + dt], sg['exp'], t14=t14, flat=f)
    return -depth * tpl


def solve(params, segs, ret=False):
    depth, t14, f, dt = params
    res = []; fits = []
    for sg in segs:
        mdl = model_seg(sg, depth, t14, f, dt)
        x = sg['t'] - sg['tc']
        X = np.column_stack([np.ones_like(x), x, x * x])
        y = sg['L'] - mdl
        c, *_ = np.linalg.lstsq(X, y, rcond=None)
        r = (y - X @ c) / sg['sig']
        res.append(r); fits.append(X @ c)
    return (np.concatenate(res), fits) if ret else np.concatenate(res)


# per-window noise from a first fit, with red-noise inflation from 1.6-h binned residuals
for sg in segs: sg['sig'] = np.std(sg['L'])
p0 = [15.0, T14, 0.7, 0.0]
r0 = least_squares(lambda q: solve(q, segs), p0, bounds=([0, 0.03, 0.0, -0.02], [60, 0.12, 0.98, 0.02]))
_, fits = solve(r0.x, segs, ret=True)
for sg, fb in zip(segs, fits):
    rr = sg['L'] - fb - model_seg(sg, *r0.x)
    s = 1.4826 * np.median(np.abs(rr - np.median(rr)))
    nb = max(1, int(round(T14 / sg['exp'])))
    k = len(rr) // nb
    binned = rr[:k * nb].reshape(k, nb).mean(1) if k >= 3 else rr
    beta = max(1.0, np.std(binned) / (s / np.sqrt(nb))) if k >= 3 else 1.0
    sg['sig'] = s * beta; sg['beta'] = beta
best = least_squares(lambda q: solve(q, segs), r0.x, bounds=([0, 0.03, 0.0, -0.02], [60, 0.12, 0.98, 0.02]))
chi_best = float((best.fun ** 2).sum())
prof_f = []
for f in np.arange(0.0, 0.96, 0.05):
    rf = least_squares(lambda q: solve([q[0], q[1], f, q[2]], segs), [best.x[0], best.x[1], best.x[3]],
                       bounds=([0, 0.03, -0.02], [60, 0.12, 0.02]))
    prof_f.append([float(f), float((rf.fun ** 2).sum() - chi_best), float(rf.x[0]), float(rf.x[1] * 24)])
prof_t = []
for t14 in np.arange(0.045, 0.10, 0.0025):
    rt = least_squares(lambda q: solve([q[0], t14, q[1], q[2]], segs), [best.x[0], best.x[2], best.x[3]],
                       bounds=([0, 0.0, -0.02], [60, 0.98, 0.02]))
    prof_t.append([float(t14 * 24), float((rt.fun ** 2).sum() - chi_best), float(rt.x[1])])
J = best.jac; cov = np.linalg.pinv(J.T @ J)
out = dict(n_windows=len(segs), depth=float(best.x[0]), depth_err=float(np.sqrt(cov[0, 0])), T14_h=float(best.x[1] * 24),
           T14_err_h=float(np.sqrt(cov[1, 1]) * 24), flat=float(best.x[2]), flat_err=float(np.sqrt(cov[2, 2])), dt_min=float(best.x[3] * 1440),
           chi2=chi_best, npts=int(len(best.fun)), profile_flat=prof_f, profile_T14=prof_t,
           beta_median=float(np.median([sg['beta'] for sg in segs])))
to_json(out, f'{BASE}/out/shape_profile.json')
print(f"depth {best.x[0]:.1f}+-{np.sqrt(cov[0,0]):.1f} e/s (point source at the M3V); T14 {best.x[1]*24:.2f}+-{np.sqrt(cov[1,1])*24:.2f} h; "
      f"T23/T14 {best.x[2]:.2f}+-{np.sqrt(cov[2,2]):.2f}; time shift {best.x[3]*1440:+.1f} min; chi2 {chi_best:.0f}/{len(best.fun)}")
print('profile over flat fraction (f, dchi2, depth, T14[h]):')
for r in prof_f: print('   %.2f %6.2f %5.1f %.2f' % tuple(r))
print('profile over T14 (h, dchi2, f):')
for r in prof_t[::2]: print('   %.2f %6.2f %.2f' % tuple(r))
# fold figure
fig, ax = plt.subplots(figsize=(7, 4))
_, fits = solve(best.x, segs, ret=True)
ph = np.concatenate([(sg['t'] - sg['tc']) * 24 for sg in segs]); yy = np.concatenate([sg['L'] - fb for sg, fb in zip(segs, fits)])
o = np.argsort(ph); ph, yy = ph[o], yy[o]
bins = np.arange(-10, 10.01, 0.25); idx = np.digitize(ph, bins)
bx = [ph[idx == i].mean() for i in range(1, len(bins)) if (idx == i).sum() > 3]; by = [yy[idx == i].mean() for i in range(1, len(bins)) if (idx == i).sum() > 3]
ax.plot(ph, yy, ',', color='0.7'); ax.plot(bx, by, 'ko', ms=3)
tt = np.linspace(-10, 10, 2000) / 24
ax.plot(tt * 24, -best.x[0] * template(tt, [best.x[3]], 0.0001, t14=best.x[1], flat=best.x[2]), 'r-')
ax.plot(tt * 24, -best.x[0] * template(tt, [best.x[3]], 0.0001, t14=best.x[1], flat=0.0), 'b--', lw=0.8)
ax.set_xlim(-6, 6); ax.set_ylim(-40, 25); ax.set_xlabel('hours from mid-transit'); ax.set_ylabel('matched-filter flux at M3V (e/s)')
fig.tight_layout(); fig.savefig(f'{BASE}/figs/shape_fold.png', dpi=90)
