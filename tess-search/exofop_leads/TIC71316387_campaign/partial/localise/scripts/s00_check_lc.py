"""Sanity check: aperture photometry from the TESScut cubes; fold on the EB ephemeris; fit a triangular
(V-shaped) eclipse with free depth, width and centre for the deep (phase 0) and shallow (phase 0.5) eclipses."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from loclib import *
from scipy.optimize import least_squares
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

out = {}
fig, axs = plt.subplots(2, 6, figsize=(24, 7))
for k, sec in enumerate(SECTORS):
    S = Sector(sec)
    g = S.good
    xt, yt = S.xt_wcs, S.yt_wcs
    ap = circle_mask((S.ny, S.nx), xt, yt, 1.6)
    ann = circle_mask((S.ny, S.nx), xt, yt, 8) & ~circle_mask((S.ny, S.nx), xt, yt, 4)
    bk = np.median(S.flux[:, ann], axis=1)
    f = (S.flux[:, ap] - bk[:, None]).sum(1)
    t = S.t[g]; f = f[g]
    # local detrend: running median 0.5 d excluding eclipses
    ph = ((t - T0) / P) % 1.0
    ine = (np.abs(((ph + 0.5) % 1) - 0.5) * P < 0.6 * T14) | (np.abs(ph - 0.5) * P < 0.6 * T14)
    trend = np.zeros_like(f)
    for i in range(len(t)):
        m = (np.abs(t - t[i]) < 0.3) & ~ine
        trend[i] = np.median(f[m]) if m.sum() > 3 else np.nan
    fn = f / trend
    ok = np.isfinite(fn)
    t, fn, ph = t[ok], fn[ok], ph[ok]
    res = dict(sector=sec, camera=S.camera, ccd=S.ccd, exptime_min=S.exptime * 1440, n_good=int(g.sum()),
               xt_wcs=xt, yt_wcs=yt, median_ap_flux=float(np.median(f)))
    for name, pc in [('deep', 0.0), ('shallow', 0.5)]:
        dph = ((ph - pc + 0.5) % 1) - 0.5
        dtt = dph * P
        sel = np.abs(dtt) < 3 * T14
        def model(p):
            d, w, c, z = p
            return z * (1 - d * eclipse_shape_w(dtt[sel], c, w, S.exptime))
        def eclipse_shape_w(x, c, w, ex):
            sub = np.linspace(-0.5, 0.5, 11) * ex
            a = np.abs(x[:, None] + sub[None, :] - c)
            return np.clip(1 - a / (w / 2), 0, 1).mean(1)
        r = least_squares(lambda p: model(p) - fn[sel], [0.006, T14, 0.0, 1.0],
                          bounds=([0, 0.02, -0.02, 0.9], [0.05, 0.2, 0.02, 1.1]))
        J = r.jac; cov = np.linalg.pinv(J.T @ J) * (r.fun**2).sum() / (sel.sum() - 4)
        e = np.sqrt(np.diag(cov))
        res[name] = dict(depth_ppm=r.x[0] * 1e6, depth_err=e[0] * 1e6, T14_h=r.x[1] * 24, T14_err_h=e[1] * 24,
                         dt_min=r.x[2] * 1440, dt_err_min=e[2] * 1440, rms_ppm=float(np.std(r.fun) * 1e6))
        ax = axs[0 if name == 'deep' else 1, k]
        ax.plot(dtt[sel] * 24, fn[sel], '.', ms=2, alpha=0.4)
        xs = np.linspace(-3 * T14, 3 * T14, 400)
        ax.plot(xs * 24, r.x[3] * (1 - r.x[0] * eclipse_shape_w(xs, r.x[2], r.x[1], S.exptime)), 'r-')
        ax.set_title(f'S{sec} {name}: {r.x[0]*1e6:.0f}+-{e[0]*1e6:.0f} ppm, T14={r.x[1]*24:.2f} h')
        ax.set_ylim(0.988, 1.006)
    print(res, flush=True)
    out[sec] = res
plt.tight_layout(); plt.savefig(f'{L}/plots/s00_aperture_fold.png', dpi=80)
to_json(out, f'{L}/results/s00_check_lc.json')
