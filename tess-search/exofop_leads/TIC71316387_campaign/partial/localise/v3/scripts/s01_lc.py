"""Sanity check of each TESScut cube: aperture light curve, background, per-sector deep/shallow depth and timing.
Depth is fitted per window as linear baseline + A * V-shape (T14 = 1.55 h), using the same window geometry as the
difference images (in-eclipse cadences excluded from the baseline).  Timing: grid over tc shifts of +-15 min."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v3lib import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

out = {}
fig, axs = plt.subplots(2, 6, figsize=(24, 7))
for k, sec in enumerate(SECTORS):
    C = Cube(sec)
    yy, xx = np.mgrid[0:C.ny, 0:C.nx]
    rr = np.hypot(xx - C.xt, yy - C.yt)
    ap = rr <= 2.0
    # background: median of pixels > 6 px from any T<17 star
    g = C.g
    far = np.ones((C.ny, C.nx), bool)
    for i in np.where(g['T'] < 17.5)[0]:
        far &= np.hypot(xx - g.x.values[i], yy - g.y.values[i]) > 5
    bk = np.median(C.f[:, far], axis=1)
    lc = (C.f[:, ap] - bk[:, None]).sum(1)
    gd = C.good
    med = np.median(lc[gd])
    fl = lc / med
    res = {}
    for kind, ph in (('deep', 0.0), ('shallow', 0.5)):
        tcs, ns = C.epochs(ph)
        rows = []
        for shift in np.arange(-15, 15.1, 1.0) / 1440:
            A_all, s_all = [], []
            for tc in tcs:
                tcs_ = tc + shift
                m = gd & (np.abs(C.t - tcs_) < OOT_B)
                if m.sum() < 8:
                    continue
                tt = C.t[m] - tcs_
                sh = vshape(C.t[m], tcs_, C.dt)
                if (sh > 0.3).sum() < 1:
                    continue
                X = np.column_stack([np.ones(m.sum()), tt, -sh])
                c, *_ = np.linalg.lstsq(X, fl[m], rcond=None)
                r = fl[m] - X @ c
                cv = np.linalg.inv(X.T @ X) * np.var(r)
                A_all.append(c[2]); s_all.append(np.sqrt(cv[2, 2]))
            A_all = np.array(A_all); s_all = np.array(s_all)
            w = 1 / s_all ** 2
            rows.append((shift * 1440, (A_all * w).sum() / w.sum(), 1 / np.sqrt(w.sum()), len(A_all)))
        rows = np.array(rows)
        j = np.argmax(rows[:, 1])
        res[kind] = dict(depth_ap_ppm=float(rows[15, 1] * 1e6), err_ppm=float(rows[15, 2] * 1e6), n=int(rows[15, 3]),
                         best_shift_min=float(rows[j, 0]), depth_at_best_ppm=float(rows[j, 1] * 1e6))
    # fold plot
    ph = ((C.t - T0) / P + 0.25) % 1 - 0.25
    ax = axs[0, k]; ax.plot(ph[gd], fl[gd], ',', alpha=0.4)
    bins = np.linspace(-0.25, 0.75, 400); bc = 0.5 * (bins[1:] + bins[:-1])
    bm = [np.median(fl[gd][(ph[gd] >= a) & (ph[gd] < b)]) if ((ph[gd] >= a) & (ph[gd] < b)).sum() > 2 else np.nan for a, b in zip(bins[:-1], bins[1:])]
    ax.plot(bc, bm, 'r-'); ax.set_ylim(0.988, 1.006); ax.set_title(f'S{sec} cam{C.cam} ccd{C.ccd}')
    ax = axs[1, k]; ax.plot(C.t[gd], fl[gd], ',', alpha=0.5); ax.set_ylim(0.985, 1.01)
    ax.set_title(f'cut {C.ncut} cadences')
    # contamination of the 2-px aperture by Gaia neighbours (TESS-band flux ratio, PRF-free estimate via 2-px radius)
    out[sec] = dict(cam=C.cam, ccd=C.ccd, cadence_min=C.dt * 1440, n_good=int(gd.sum()), n_centroid_cut=C.ncut,
                    xt=C.xt, yt=C.yt, epoch=C.epoch, bkg_med=float(np.median(bk[gd])), ap_flux_med=float(med),
                    jac=C.jac.tolist(), **res)
    print(sec, out[sec]['cam'], out[sec]['ccd'], f"{out[sec]['cadence_min']:.1f}min good={out[sec]['n_good']} cut={C.ncut} bkg={out[sec]['bkg_med']:.0f}",
          {k: (round(v['depth_ap_ppm']), round(v['err_ppm']), v['n'], v['best_shift_min']) for k, v in res.items()}, flush=True)
plt.tight_layout(); plt.savefig(f'{V3}/plots/s01_lc.png', dpi=60)
to_json(out, f'{V3}/results/s01_lc.json')
