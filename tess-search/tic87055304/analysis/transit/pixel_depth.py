"""Per-pixel transit flux drop (Delta F, e-/s) in the TESScut FFI cutouts, every sector; background-independent
check of the sector-to-sector depth consistency and of the intrinsic depth.

For each sector and each pixel of the central 9x9: all transits with coverage >= 0.7 and >= 0.5 h of data on both
sides inside +-0.25 d windows (restricted to one momentum-dump segment) are fitted jointly with
  F_pix(t) = c0_i + c1_i (t - Tc_i) - dF_pix * s(t),
s = (1 - template) / template_depth (unit-depth, exposure-integrated free-density MCMC shape), per-transit linear
baselines c0_i, c1_i. dF_pix is the mid-transit flux drop in e-/s. Errors: white residual scatter per pixel x the
time-averaging beta of the summed 5x5 residuals.
Because dF is a flux difference, it does not depend on the background level or on contaminating (constant) flux.
Intrinsic depth on the target = dF_sum / (F_T10 * 10^(-0.4 (T - 10)) * PRF fraction), F_T10 = 15000 e-/s
(TESS instrument handbook; ~10 % uncertainty), T = 11.7249 (TIC 8.2).
PRF fits of the dF map: scale-only fits of the target PRF and of each Gaia neighbour's PRF (positions fixed),
chi^2 compared (the separation of the 9" neighbour, 0.43 px, is too small to be resolved this way).
Writes pixel_depth.json and pixel_depth.png.
"""
import os, json, glob, warnings
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import tmodel as tm
from transit_fits import Template
from vshape import beta_time_avg
from aperture_depth import gaia, TARGET_ID, NB_ID, RA_T, DEC_T, raw

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = "/home/user/work/tic87/data"
mc = json.load(open(os.path.join(HERE, "mcmc_free.json")))
th = np.array([mc["params"][n]["med"] for n in ["T0", "P", "k", "b", "rho", "u1", "u2"]])
tmpl = Template(th); T0, P = th[0], th[1]
F_T10 = 15000.0; TMAG = 11.7249
F_TARGET = F_T10 * 10 ** (-0.4 * (TMAG - 10))
W, HALF = 0.25, 4


def main():
    import PRF
    res = dict(F_target_eps=F_TARGET, sectors={})
    files = sorted(glob.glob(f"{DATA}/tesscut/*.fits"))
    fig, axes = plt.subplots(2, len(files), figsize=(3.3 * len(files), 6.4))
    for j, fpath in enumerate(files):
        s = int(os.path.basename(fpath).split("-")[1][1:])
        h = fits.open(fpath); d = h[1].data; w = WCS(h[2].header)
        cam, ccd = h[0].header["CAMERA"], h[0].header["CCD"]
        col0, row0 = h[1].header["1CRV4P"] - h[1].header["1CRPX4"], h[1].header["2CRV4P"] - h[1].header["2CRPX4"]
        t = d["TIME"]; flux = d["FLUX"].astype(float); q = d["QUALITY"]
        dumps = np.unique(np.round(np.asarray(raw[s]["TIME"])[(np.asarray(raw[s]["QUALITY"]) & 32) > 0], 3))
        ok = (q == 0) & np.isfinite(t) & np.all(np.isfinite(flux.reshape(len(t), -1)), axis=1)
        ok &= ~((t > 4178.80) & (t < 4179.20))
        if len(dumps):
            ok &= np.min(np.abs(t[:, None] - dumps[None, :]), axis=1) >= 1 / 24
        t, flux = t[ok], flux[ok]
        tx, ty = w.world_to_pixel_values(RA_T, DEC_T)
        med = np.median(flux, axis=0)
        yy, xx = np.mgrid[0:med.shape[0], 0:med.shape[1]]
        far = np.hypot(xx - tx, yy - ty) >= 5
        sky = far & (med < np.percentile(med[far], 35))
        flux = flux - np.median(flux[:, sky], axis=1)[:, None, None]
        E = tm.EXPT[s]
        # segments
        lab = np.zeros(len(t), int); k = 0
        for i in range(1, len(t)):
            if t[i] - t[i - 1] > 0.5 or np.any((dumps > t[i - 1]) & (dumps <= t[i])):
                k += 1
            lab[i] = k
        ep_all = np.arange(np.ceil((t.min() - T0) / P), np.floor((t.max() - T0) / P) + 1).astype(int)
        blocks = []
        for e in ep_all:
            Tc = T0 + e * P
            m = np.abs(t - Tc) < W
            if m.sum() < 10: continue
            inn = m & (np.abs(t - Tc) < tmpl.T14 / 2)
            if inn.sum() == 0: continue
            vals, cnt = np.unique(lab[inn], return_counts=True)
            m &= lab == vals[np.argmax(cnt)]
            cov = min(1, np.sum(m & (np.abs(t - Tc) < tmpl.T14 / 2)) * E / 86400 / tmpl.T14)
            pre = np.sum(m & (t < Tc - tmpl.T14 / 2)) * E / 86400; post = np.sum(m & (t > Tc + tmpl.T14 / 2)) * E / 86400
            if cov >= 0.7 and pre >= 0.5 / 24 and post >= 0.5 / 24:
                blocks.append((e, np.where(m)[0], Tc))
        idx = np.concatenate([b[1] for b in blocks])
        nb = len(blocks)
        X = np.zeros((len(idx), 2 * nb + 1)); r0 = 0
        for i, (e, ii, Tc) in enumerate(blocks):
            n = len(ii)
            X[r0:r0 + n, 2 * i] = 1; X[r0:r0 + n, 2 * i + 1] = t[ii] - Tc
            X[r0:r0 + n, -1] = -(1 - tmpl(t[ii] - Tc, E)) / tmpl.depth
            r0 += n
        ix, iy = int(np.round(tx)), int(np.round(ty))
        ys, xs = slice(iy - HALF, iy + HALF + 1), slice(ix - HALF, ix + HALF + 1)
        cube = flux[idx][:, ys, xs]
        Y = cube.reshape(len(idx), -1)
        coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
        R = Y - X @ coef
        dof = len(idx) - X.shape[1]
        sig = np.sqrt(np.sum(R ** 2, axis=0) / dof)
        XtXi = np.linalg.inv(X.T @ X)
        dF = coef[-1].reshape(2 * HALF + 1, 2 * HALF + 1)
        edF = (sig * np.sqrt(XtXi[-1, -1])).reshape(dF.shape)
        # red noise on the summed 5x5
        c5 = slice(HALF - 2, HALF + 3)
        Rsum = R.reshape(len(idx), 2 * HALF + 1, 2 * HALF + 1)[:, c5, c5].sum(axis=(1, 2))
        beta, _ = beta_time_avg(Rsum, t[idx], E)
        # summed-aperture fits (proper error from the summed light curve)
        sums = {}
        for name, half in (("3x3", 1), ("5x5", 2), ("7x7", 3), ("9x9", 4)):
            cc = slice(HALF - half, HALF + half + 1)
            ysum = cube[:, cc, cc].sum(axis=(1, 2))
            cs, *_ = np.linalg.lstsq(X, ysum, rcond=None)
            rs = ysum - X @ cs
            ss = np.sqrt(np.sum(rs ** 2) / dof)
            bs, _ = beta_time_avg(rs, t[idx], E)
            sums[name] = dict(dF=float(cs[-1]), err=float(ss * np.sqrt(XtXi[-1, -1]) * bs), beta=float(bs))
        # PRF predictions
        prf = PRF.TESS_PRF(cam, ccd, s, col0 + tx, row0 + ty)
        def prf_img(ra, dec):
            x, y = w.world_to_pixel_values(ra, dec)
            img = prf.locate(float(x), float(y), med.shape)
            return img[ys, xs]
        tgt = gaia[gaia.source_id == TARGET_ID].iloc[0]
        Pt = prf_img(tgt.ra, tgt.dec)
        fits_ = {}
        wgt = 1 / (edF * beta) ** 2
        for r in gaia[(gaia.sep < 90) & (gaia.Tmag < 17.5)].itertuples():
            Pn = prf_img(r.ra, r.dec)
            if Pn.sum() < 0.05: continue
            a = np.sum(dF * Pn * wgt) / np.sum(Pn ** 2 * wgt)
            chi2 = float(np.sum((dF - a * Pn) ** 2 * wgt))
            Fn = F_T10 * 10 ** (-0.4 * (r.Tmag - 10))
            fits_[int(r.source_id)] = dict(sep=float(r.sep), Tmag=float(r.Tmag), chi2=chi2, dF_total=float(a), depth_if_host_ppm=float(a / Fn * 1e6))
        ct = fits_[int(TARGET_ID)]
        sres = dict(n_transits=nb, epochs=[int(b[0]) for b in blocks], beta=float(beta), sums=sums,
                    dF_sum_9x9_err=float(np.sqrt(np.sum((edF * beta) ** 2))),
                    target_prf_fraction_3x3=float(Pt[HALF - 1:HALF + 2, HALF - 1:HALF + 2].sum()), target_prf_fraction_9x9=float(Pt.sum()),
                    intrinsic_depth_target_ppm=float(ct["dF_total"] / F_TARGET * 1e6),
                    intrinsic_depth_target_err_ppm=float(sums["5x5"]["err"] / Pt[HALF - 2:HALF + 3, HALF - 2:HALF + 3].sum() / F_TARGET * 1e6),
                    chi2_target=ct["chi2"], npix=int(dF.size),
                    best_hosts=sorted(fits_.items(), key=lambda kv: kv[1]["chi2"])[:6],
                    dF_map=dF.tolist(), dF_err_map=edF.tolist())
        res["sectors"][s] = sres
        print(s, nb, "beta %.2f" % beta, {k: (round(v["dF"], 2), round(v["err"], 2)) for k, v in sums.items()},
              "intrinsic(target) %.0f +- %.0f ppm" % (sres["intrinsic_depth_target_ppm"], sres["intrinsic_depth_target_err_ppm"]),
              "chi2_t %.1f / %d" % (ct["chi2"], dF.size - 1), "best:", [(k, round(v["sep"], 1), round(v["chi2"], 1)) for k, v in sres["best_hosts"][:3]], flush=True)
        ax = axes[0, j]
        vmax = np.nanmax(np.abs(dF))
        ax.imshow(dF, origin="lower", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.plot(tx - ix + HALF, ty - iy + HALF, "k+", ms=10)
        ax.set_title(f"S{s} dF map (e/s), n={nb}", fontsize=9)
        ax = axes[1, j]
        ax.imshow(dF - ct["dF_total"] * Pt, origin="lower", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_title(f"residual after target PRF; chi2={ct['chi2']:.0f}/{dF.size - 1}", fontsize=8)
    # consistency of the background-independent dF across sectors
    v = np.array([res["sectors"][s]["sums"]["5x5"]["dF"] for s in res["sectors"]]); e = np.array([res["sectors"][s]["sums"]["5x5"]["err"] for s in res["sectors"]])
    wm = np.sum(v / e ** 2) / np.sum(1 / e ** 2)
    from scipy import stats
    chi2 = float(np.sum(((v - wm) / e) ** 2))
    res["dF5x5_consistency"] = dict(mean=float(wm), err=float(1 / np.sqrt(np.sum(1 / e ** 2))), chi2=chi2, dof=len(v) - 1, p=float(stats.chi2.sf(chi2, len(v) - 1)))
    v3 = np.array([res["sectors"][s]["sums"]["3x3"]["dF"] for s in res["sectors"]]); e3 = np.array([res["sectors"][s]["sums"]["3x3"]["err"] for s in res["sectors"]])
    wm3 = np.sum(v3 / e3 ** 2) / np.sum(1 / e3 ** 2); chi3 = float(np.sum(((v3 - wm3) / e3) ** 2))
    res["dF3x3_consistency"] = dict(mean=float(wm3), err=float(1 / np.sqrt(np.sum(1 / e3 ** 2))), chi2=chi3, dof=len(v3) - 1, p=float(stats.chi2.sf(chi3, len(v3) - 1)))
    print(res["dF5x5_consistency"], res["dF3x3_consistency"])
    json.dump(res, open(os.path.join(HERE, "pixel_depth.json"), "w"), indent=1)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "pixel_depth.png"), dpi=120)


if __name__ == "__main__":
    main()
