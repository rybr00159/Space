"""s11: where does the ~8.5-9.2-d periodicity of the TARS light curve of TIC 214311908 come from?

TARS uses a 1-pixel aperture (header APERTURE = 1). We rebuild 1-pixel light curves from the TESScut cubes
(read-only, /home/user/work/tesscut5) for every pixel of the 21x21 cutout in S12/S39/S66/S93 (the TARS sectors),
background-subtracted with the phot faint-pixel background (phot/data/jitter_s*.npz), QUALITY == 0, 1-h bins,
each normalised by its median (as a normalised 1-px light curve). For each pixel: Lomb-Scargle power and semi-amplitude
(in e-/s and relative) at the TARS best period of that sector (from s06) and the best period in 7-11 d.
(1) The pixel whose light curve correlates best with the TARS light curve is identified (it should hold the M3V).
(2) The e-/s semi-amplitude image at the TARS period is fitted with the Mira's empirical image (phot b-image, s04b),
    the M3V PRF, the pair PRF and a constant: which star carries the 8.7-d signal?
(3) The same LS on the deblended light curves (raw model-B M3V, model-A M3V and pair, Mira) regressed on background and
    pointing (as phot p12), amplitudes in e-/s at the TARS period.
Output: out/tars_pixels.json, figs/tars_pixels.png
"""
import sys; sys.path.insert(0, "tic214311908/analysis/search/scripts")
sys.path.insert(1, "tic214311908/analysis/phot/scripts")
from slib import *
from astropy.io import fits
from astropy.timeseries import LombScargle
import glob
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt


def bin1h(t, Y):
    k = np.floor((t - t[0]) * 24).astype(int); u, inv = np.unique(k, return_inverse=True); n = np.bincount(inv)
    tb = np.bincount(inv, t) / n
    if Y.ndim == 1: return tb, np.bincount(inv, Y) / n, n
    Yb = np.zeros((len(u),) + Y.shape[1:]); np.add.at(Yb, inv, Y); Yb /= n.reshape((-1,) + (1,) * (Y.ndim - 1))
    return tb, Yb, n


def sinamp(t, y, P):
    X = np.column_stack([np.ones(len(t)), np.sin(2 * np.pi * t / P), np.cos(2 * np.pi * t / P), (t - t.mean()) / 10])
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    return float(np.hypot(c[1], c[2]))


def main():
    import photlib
    pp = json.load(open(f"{ROOT}/out/public_periodograms.json"))
    out = {}
    fig, axs = plt.subplots(2, 4, figsize=(17, 8))
    for col, s in enumerate((12, 39, 66, 93)):
        C = photlib.Cutout(s)
        J = np.load(f"{PHOT}/data/jitter_s{s}.npz")
        ok = C.fin & (C.q == 0) & np.isfinite(J["bkg"])
        t = C.t[ok]; cube = C.flux[ok] - J["bkg"][ok][:, None, None]
        tb, cb, nb_ = bin1h(t, cube)
        med = np.nanmedian(cb, 0)
        # TARS light curve
        f = glob.glob(f"{BASE}/ctx/lc/hlsp_tars_tess_ffi_s{s:04d}-0000000214311908_tess_v01_lc.fits")[0]
        d = fits.open(f)[1].data; tt = np.array(d["TIME"], float); ft = np.array(d["FLUX"], float)
        m = np.isfinite(tt) & np.isfinite(ft); tt, ft = tt[m], ft[m] / np.nanmedian(ft[m])
        Ptars = pp[f"TARS_M3V_S{s}"]["top"][0][0]
        # correlation of each pixel's normalised 1-h light curve with TARS (interpolated to the pixel bins)
        ftb = np.interp(tb, tt, ft, left=np.nan, right=np.nan)
        k = np.isfinite(ftb)
        corr = np.full(med.shape, np.nan); amp_e = np.full(med.shape, np.nan); amp_rel = np.full(med.shape, np.nan)
        pbest = np.full(med.shape, np.nan); pow_t = np.full(med.shape, np.nan)
        for iy in range(med.shape[0]):
            for ix in range(med.shape[1]):
                y = cb[:, iy, ix]
                if not np.all(np.isfinite(y)) or med[iy, ix] <= 0: continue
                yn = y / med[iy, ix]
                corr[iy, ix] = np.corrcoef(yn[k], ftb[k])[0, 1]
                amp_e[iy, ix] = sinamp(tb, y, Ptars); amp_rel[iy, ix] = amp_e[iy, ix] / med[iy, ix]
                fr = np.linspace(1 / 11, 1 / 7, 300); L = LombScargle(tb, yn); pw = L.power(fr)
                pbest[iy, ix] = 1 / fr[np.argmax(pw)]; pow_t[iy, ix] = L.power(1 / Ptars)
        z = np.load(f"{ROOT}/data/pixfull_s{s}.npz", allow_pickle=True)
        m3 = z["m3v_xy"]; mi = z["mira_xy"]; pa = z["pair_xy"]
        iy3, ix3 = int(round(m3[1])), int(round(m3[0]))
        jbest = np.unravel_index(np.nanargmax(corr), corr.shape)
        # fit of the e-/s amplitude image (signed: project on the phase of the M3V pixel) with templates
        def signed_amp_image():
            X = np.column_stack([np.ones(len(tb)), np.sin(2 * np.pi * tb / Ptars), np.cos(2 * np.pi * tb / Ptars), (tb - tb.mean()) / 10])
            c, *_ = np.linalg.lstsq(X, cb.reshape(len(tb), -1), rcond=None)
            a, b = c[1].reshape(med.shape), c[2].reshape(med.shape)
            ph = np.arctan2(b[iy3, ix3], a[iy3, ix3])
            return a * np.cos(ph) + b * np.sin(ph)
        S = signed_amp_image()
        from s05b_pixfull_periodic import PRFs
        prfs = PRFs(s)
        bimg = z["b"] / np.nansum(z["b"])
        yy, xx = np.mgrid[0:med.shape[0], 0:med.shape[1]]
        fm = (np.hypot(xx - mi[0], yy - mi[1]) <= 4) & np.isfinite(S)
        # noise per pixel ~ amplitude scatter at fake periods
        fakes = np.random.default_rng(s).uniform(3.0, 6.0, 40)
        SF = []
        for pf in fakes:
            X = np.column_stack([np.ones(len(tb)), np.sin(2 * np.pi * tb / pf), np.cos(2 * np.pi * tb / pf), (tb - tb.mean()) / 10])
            c, *_ = np.linalg.lstsq(X, cb.reshape(len(tb), -1), rcond=None)
            SF.append(np.hypot(c[1], c[2]).reshape(med.shape))
        E = np.std(SF, 0) + 1e-3
        fits_ = {}
        for nm, T in (("Mira_b", bimg), ("M3V", prfs.img(*m3)), ("pair", prfs.img(*pa))):
            w = 1 / E[fm]
            X = np.column_stack([T[fm], np.ones(fm.sum())]) * w[:, None]
            c, *_ = np.linalg.lstsq(X, S[fm] * w, rcond=None)
            fits_[nm] = dict(chi2=float(np.sum((X @ c - S[fm] * w) ** 2)), amp=float(c[0]))
        # joint Mira + M3V + pair
        X = np.column_stack([bimg[fm], prfs.img(*m3)[fm], prfs.img(*pa)[fm], np.ones(fm.sum())]) / E[fm][:, None]
        c, *_ = np.linalg.lstsq(X, S[fm] / E[fm], rcond=None)
        cov = np.linalg.pinv(X.T @ X)
        fits_["joint"] = dict(amp_Mira=float(c[0]), amp_M3V=float(c[1]), amp_pair=float(c[2]),
                              err=[float(np.sqrt(cov[i, i])) for i in range(3)], chi2=float(np.sum((X @ c - S[fm] / E[fm]) ** 2)), npix=int(fm.sum()))
        out[s] = dict(P_tars=Ptars, tars_pixel_best_corr=[int(jbest[1]), int(jbest[0])], corr_best=float(np.nanmax(corr)),
                      m3v_pixel=[ix3, iy3], corr_m3v_pixel=float(corr[iy3, ix3]),
                      amp_rel_m3v_pixel=float(amp_rel[iy3, ix3]), amp_e_m3v_pixel=float(amp_e[iy3, ix3]),
                      median_amp_rel_within_4px_of_Mira=float(np.nanmedian(amp_rel[fm])),
                      frac_pixels_within_4px_with_Pbest_7p5_10=float(np.nanmean((pbest[fm] > 7.5) & (pbest[fm] < 10))),
                      template_fits=fits_)
        print(f"S{s}: P_tars {Ptars:.2f}; best-corr pixel {jbest[::-1]} (r {np.nanmax(corr):.2f}), M3V pixel ({ix3},{iy3}) r {corr[iy3, ix3]:.2f}; "
              f"amp at M3V pixel {amp_e[iy3, ix3]:.1f} e/s ({amp_rel[iy3, ix3]*1e3:.1f} ppt); joint fit amps Mira {c[0]:.0f}+-{np.sqrt(cov[0,0]):.0f} "
              f"M3V {c[1]:.1f}+-{np.sqrt(cov[1,1]):.1f} pair {c[2]:.1f}+-{np.sqrt(cov[2,2]):.1f} e/s; single chi2 " +
              ", ".join(f"{k} {v['chi2']:.0f}" for k, v in fits_.items() if k != "joint") + f" / {fm.sum()} px; frac pix P 7.5-10 d {out[s]['frac_pixels_within_4px_with_Pbest_7p5_10']:.2f}", flush=True)
        ax = axs[0, col]; im = ax.imshow(S / E, origin="lower", cmap="RdBu_r", vmin=-8, vmax=8)
        ax.plot(*m3, "g+", ms=10); ax.plot(*mi, "kx", ms=10); ax.plot(*pa, "m+", ms=8); ax.set_title(f"S{s}: signed amp/err at P_TARS {Ptars:.2f} d")
        ax = axs[1, col]; im2 = ax.imshow(pbest, origin="lower", vmin=7, vmax=11); ax.plot(*m3, "r+", ms=10); ax.plot(*mi, "kx", ms=10)
        ax.set_title("best LS period 7-11 d per pixel"); plt.colorbar(im2, ax=ax, fraction=0.04)
    fig.tight_layout(); fig.savefig(f"{ROOT}/figs/tars_pixels.png", dpi=75); plt.close(fig)
    save_json({str(k): v for k, v in out.items()}, f"{ROOT}/out/tars_pixels.json")


if __name__ == "__main__":
    main()
