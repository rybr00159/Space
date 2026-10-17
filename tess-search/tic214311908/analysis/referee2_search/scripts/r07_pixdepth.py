"""r07: model-free whole-cutout check that no distant pixel carries the 5.69-d dip (independent of the phot pipeline,
which can couple distant variables into the M3V light curve, as the RRab shows).

For each sector, from the raw TESScut cube (FLUX, QUALITY == 0, finite), per pixel and per transit: depth =
(local linear baseline fitted at 0.06 < |dt| < 0.30 d) - mean flux at |dt| < 0.022 d (flat part of a 1.6-h transit),
averaged over the sector's transits; the same at 80 fake epochs (>= 0.4 d from any transit) gives a per-pixel
empirical error. Output per sector: S/N image, the pixels with |S/N| > 4, and the distance of each from the M3V;
summed depth in the 3x3 box around the M3V; the strongest positive-S/N pixel beyond 2.5 px of the M3V.
A per-cadence common mode (median of all pixels > 6 px from the Mira) is subtracted first.
Uses the phot ephemeris (T0, P) and the TESScut TIME (FFI-barycentred; target offset < 31 s, negligible).
Output: out/r07_pixdepth.json, figs/r07_pixdepth.png
"""
import sys, glob
sys.path.insert(0, "tic214311908/analysis/referee2_search/scripts")
from rlib2 import *
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

CUBES = "/home/user/work/tesscut5/214311906"
M3V = (260.969674, -51.323821)       # Gaia DR3 2016.0 (PM ~0.14"/yr: < 1.5" by 2026, irrelevant at the px level)
MIRA = (260.9776003, -51.3237201)
PAIR = None
rng = np.random.default_rng(707)


def epoch_depth(t, F, tc):
    k = np.abs(t - tc) < 0.30
    inn = np.abs(t - tc) < 0.022
    base = k & (np.abs(t - tc) > 0.06)
    if inn.sum() < 1 or base.sum() < 8: return None
    X = np.column_stack([np.ones(base.sum()), t[base] - tc])
    c, *_ = np.linalg.lstsq(X, F[base], rcond=None)
    return c[0] - F[inn].mean(0)            # positive = dimming (per pixel)


def main():
    out = {}
    fig, axs = plt.subplots(2, 3, figsize=(14, 8)); axs = axs.ravel()
    for i, sec in enumerate(SECTORS):
        fn = glob.glob(f"{CUBES}/tess-s{sec:04d}-*_astrocut.fits")[0]
        with fits.open(fn, memmap=True) as h:
            tab = h[1].data; w = WCS(h[2].header)
            t = np.array(tab["TIME"], float); F = np.array(tab["FLUX"], np.float64); q = np.array(tab["QUALITY"], int)
        ok = (q == 0) & np.isfinite(t) & np.all(np.isfinite(F.reshape(len(t), -1)), axis=1)
        t, F = t[ok], F[ok].reshape(ok.sum(), -1)
        ny, nx = 21, 21
        mx, my = w.all_world2pix([M3V], 0)[0]; Mx, My = w.all_world2pix([MIRA], 0)[0]
        # common-mode removal: per-cadence median of all pixels > 6 px from the Mira (background / scattered-light steps)
        yy0, xx0 = np.mgrid[0:ny, 0:nx]
        cm_pix = (np.hypot(xx0 - Mx, yy0 - My) > 6).ravel()
        F = F - np.median(F[:, cm_pix], axis=1)[:, None]
        es = [e for e in range(-450, 100) if t.min() + 0.3 < T0 + e * P < t.max() - 0.3]
        D = [epoch_depth(t, F, T0 + e * P) for e in es]
        D = [x for x in D if x is not None]
        if len(D) == 0: continue
        Dm = np.mean(D, 0)
        fk = []
        while len(fk) < 80:
            tc = rng.uniform(t.min() + 0.3, t.max() - 0.3)
            if abs(phase_dist(tc)) < 0.4: continue
            # fake 'sector' of the same number of events: average len(D) random epochs
            vals = []
            for _ in range(len(D)):
                tc2 = rng.uniform(t.min() + 0.3, t.max() - 0.3)
                if abs(phase_dist(tc2)) < 0.4: continue
                v = epoch_depth(t, F, tc2)
                if v is not None: vals.append(v)
            if len(vals) == len(D): fk.append(np.mean(vals, 0))
        fk = np.array(fk)
        E = 1.4826 * np.median(np.abs(fk - np.median(fk, 0)), 0)
        SN = (Dm / E).reshape(ny, nx)
        yy, xx = np.mgrid[0:ny, 0:nx]
        dist_m3v = np.hypot(xx - mx, yy - my) * 21.0
        dist_mira = np.hypot(xx - Mx, yy - My) * 21.0
        far = dist_m3v.ravel() > 2.5 * 21
        sn = SN.ravel()
        hi = np.where(np.abs(sn) > 4)[0]
        # empirical tail: how often does a fake-epoch S/N image have a far pixel above the observed far max?
        fkSN = fk / E
        far_max_fake = np.max(fkSN[:, far], 1)
        jfar = int(np.argmax(np.where(far, sn, -np.inf)))
        box = (np.abs(xx - mx) <= 1.5) & (np.abs(yy - my) <= 1.5)
        res = dict(n_transits=len(D), m3v_xy=[float(mx), float(my)], mira_xy=[float(Mx), float(My)],
                   sum_depth_3x3_e_s=float(Dm.reshape(ny, nx)[box].sum()),
                   sum_depth_3x3_err=float(np.std(fk.reshape(len(fk), ny, nx)[:, box].sum(1))),
                   max_SN_pixel_dist_m3v_arcsec=float(dist_m3v.ravel()[np.argmax(sn)]), max_SN=float(sn.max()),
                   far_max_SN=float(sn[jfar]), far_max_dist_m3v_arcsec=float(dist_m3v.ravel()[jfar]),
                   far_max_dist_mira_arcsec=float(dist_mira.ravel()[jfar]),
                   far_max_fake_p50=float(np.median(far_max_fake)), far_max_fake_p95=float(np.percentile(far_max_fake, 95)),
                   p_far=float(np.mean(far_max_fake >= sn[jfar])),
                   err_far_pixels_median_e_s=float(np.median(E[far])), err_far_pixels_p90_e_s=float(np.percentile(E[far], 90)),
                   err_m3v_pixel_e_s=float(E.reshape(ny, nx)[int(round(my)), int(round(mx))]),
                   pixels_absSN_gt4=[dict(x=int(j % nx), y=int(j // nx), SN=float(sn[j]), d_m3v=float(dist_m3v.ravel()[j]),
                                          d_mira=float(dist_mira.ravel()[j])) for j in hi])
        out[sec] = res
        print(sec, {k: v for k, v in res.items() if k != "pixels_absSN_gt4"}, "n|SN|>4:", len(hi),
              [(p["x"], p["y"], round(p["SN"], 1), round(p["d_m3v"])) for p in res["pixels_absSN_gt4"]], flush=True)
        im = axs[i].imshow(SN, origin="lower", cmap="RdBu_r", vmin=-6, vmax=6)
        axs[i].plot(mx, my, "k+", ms=12); axs[i].plot(Mx, My, "gx", ms=10)
        axs[i].set_title(f"S{sec} dip S/N (+ = dimming)")
        del F
    fig.colorbar(im, ax=axs.tolist(), shrink=0.6)
    fig.savefig(f"{ME}/figs/r07_pixdepth.png", dpi=80)
    jdump(out, f"{ME}/out/r07_pixdepth.json")


if __name__ == "__main__":
    main()
