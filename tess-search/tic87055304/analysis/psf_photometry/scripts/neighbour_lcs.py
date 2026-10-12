"""Task 4: the neighbour TIC 87055315's own S13 QLP and eleanor-lite light curves, compared with the target's.

For each file: quality == 0, fold on the ephemeris, depth with depthfit (per-transit quadratic baseline, control-phase
errors).  Then test how the neighbour's light curve relates to the target's:
  - identical inputs? (eleanor-lite: same RAW_FLUX array?)
  - linear map nb - 1 = k (tgt - 1) + c over common cadences (k = ratio of normalisations if both come from nearly the
    same aperture flux; a dilution/contamination correction F' = (F - C)/(F_med - C) multiplies fractional variations
    by F_med/(F_med - C));
  - ratio of the QLP flux errors.
Prediction for the depth in an aperture centred on the neighbour (QLP-like circular aperture, r = 2.5 px, at S13
positions), using the PRF scene: under H_target / H_neighbour, with and without a contamination correction that
removes the flux of all other stars (TIC-style: nb's own flux fraction only).
"""
import os, sys, json, glob
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from astropy.io import fits
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import *
from depthfit import summarize, detrended, WIN
from psfphot import load_scene

MD = f"{DATA}/mast_dl/mastDownload/HLSP"
OUT = os.path.join(ROOT, "neighbour_lcs"); os.makedirs(OUT, exist_ok=True)


def load(kind, tic):
    if kind == "qlp":
        f = glob.glob(f"{MD}/hlsp_qlp_tess_ffi_s0013-{tic:016d}_tess_v01_llc/*.fits")[0]
        h = fits.open(f); d = h[1].data
        ok = (d["QUALITY"] == 0) & np.isfinite(d["SAP_FLUX"])
        return dict(t=np.array(d["TIME"][ok], float), sap=np.array(d["SAP_FLUX"][ok], float), ksp=np.array(d["KSPSAP_FLUX"][ok], float),
                    err=np.array(d["KSPSAP_FLUX_ERR"][ok], float), bkg=np.array(d["SAP_BKG"][ok], float),
                    x=np.array(d["SAP_X"][ok], float), y=np.array(d["SAP_Y"][ok], float), cad=np.array(d["CADENCENO"][ok]),
                    tmag=h[0].header["TESSMAG"], file=f)
    f = glob.glob(f"{MD}/hlsp_gsfc-eleanor-lite_tess_ffi_s0013-{tic:016d}_tess_v1.0_lc/*.fits")[0]
    h = fits.open(f); d = h[1].data
    ok = (d["QUALITY"] == 0) & np.isfinite(d["CORR_FLUX"]) & np.isfinite(d["RAW_FLUX"])
    return dict(t=np.array(d["TIME"][ok], float), raw=np.array(d["RAW_FLUX"][ok], float), corr=np.array(d["CORR_FLUX"][ok], float),
                pca=np.array(d["PCA_FLUX"][ok], float), bkg=np.array(d["FLUX_BKG"][ok], float), cad=np.array(d["FFIINDEX"][ok]),
                xc=np.array(d["X_CENTROID"][ok], float), yc=np.array(d["Y_CENTROID"][ok], float), aperture=h[0].header.get("APERTURE"),
                ap_mask=np.array(h[2].data[h[0].header.get("APERTURE")]), tpf_shape=None, file=f)


def main():
    res = {}
    lcs = {}
    for kind in ("qlp", "eleanor"):
        for name, tic in (("target", 87055304), ("neighbour", 87055315)):
            lcs[(kind, name)] = load(kind, tic)
    ex = 1800.0
    # identical-input checks
    et, en = lcs[("eleanor", "target")], lcs[("eleanor", "neighbour")]
    common = np.intersect1d(et["cad"], en["cad"])
    it = np.isin(et["cad"], common); inn = np.isin(en["cad"], common)
    res["eleanor_identical_raw_flux"] = bool(np.allclose(et["raw"][it], en["raw"][inn]))
    res["eleanor_identical_corr_flux"] = bool(np.allclose(et["corr"][it], en["corr"][inn], rtol=1e-6))
    res["eleanor_identical_ap_mask"] = bool(np.array_equal(et["ap_mask"], en["ap_mask"]))
    res["eleanor_centroid_diff_px"] = [float(np.median(en["xc"][inn] - et["xc"][it])), float(np.median(en["yc"][inn] - et["yc"][it]))]
    qt, qn = lcs[("qlp", "target")], lcs[("qlp", "neighbour")]
    common = np.intersect1d(qt["cad"], qn["cad"])
    it = np.isin(qt["cad"], common); inn = np.isin(qn["cad"], common)
    a = qt["sap"][it] - 1; b = qn["sap"][inn] - 1
    good = (np.abs(a) < 0.02) & (np.abs(b) < 0.1)
    k, c = np.polyfit(a[good], b[good], 1)
    res["qlp_linear_map_nb_vs_target"] = dict(k=float(k), c=float(c), corr=float(np.corrcoef(a[good], b[good])[0, 1]),
                                              rms_resid=float(np.std(b[good] - (k * a[good] + c))), rms_nb=float(np.std(b[good])))
    res["qlp_err_ratio_nb_over_target"] = float(np.median(qn["err"]) / np.median(qt["err"]))
    res["qlp_aperture_centre_diff_px"] = [float(np.median(qn["x"][inn] - qt["x"][it])), float(np.median(qn["y"][inn] - qt["y"][it]))]
    res["qlp_bkg_ratio"] = float(np.median(qn["bkg"][inn]) / np.median(qt["bkg"][it]))
    # depths
    series = {("qlp", "target"): ["sap", "ksp"], ("qlp", "neighbour"): ["sap", "ksp"],
              ("eleanor", "target"): ["raw", "corr", "pca"], ("eleanor", "neighbour"): ["raw", "corr", "pca"]}
    fig, axs = plt.subplots(2, 2, figsize=(11, 7))
    for j, (key, cols) in enumerate(series.items()):
        L = lcs[key]
        ys = [L[c] / np.nanmedian(L[c]) for c in cols]
        r = summarize(L["t"], ys, ex, cols)
        res[f"{key[0]}_{key[1]}"] = {c: dict(depth=r[c]["A"], err=r[c]["err"], n_tr=r[c]["n_tr"]) for c in cols}
        ax = axs.ravel()[j]
        for c, y, col in zip(cols, ys, ["k", "C1", "C2"]):
            dd = detrended(L["t"], y, ex)
            if dd is None: continue
            dtk, yd, A, keep = dd
            ax.plot(dtk * 24, yd * 1e3, ".", color=col, ms=3, alpha=0.6, label=f"{c}: depth {r[c]['A']*1e3:.2f}+-{r[c]['err']*1e3:.2f} ppt")
        ax.set_title(f"S13 {key[0]} light curve of the {key[1]} (TIC {87055304 if key[1]=='target' else 87055315})", fontsize=9)
        ax.set_xlabel("hours from mid-transit"); ax.set_ylabel("ppt"); ax.legend(fontsize=7)
        ax.set_ylim(-40, 15) if key == ("qlp", "neighbour") else ax.set_ylim(-12, 6)
    plt.tight_layout(); plt.savefig(os.path.join(ROOT, "neighbour_own_lcs_s13.png"), dpi=90); plt.close()

    # predictions from the PRF scene for a QLP-like circular aperture centred on each star (S13)
    sec = load_sector(13); sc, js = load_scene(13)
    prf = PRFModel(sec, sig_px=js["sig_jit"], xt=js["target_xy"][0], yt=js["target_xy"][1])
    X = sc.x.values + js["dx0"]; Y = sc.y.values + js["dy0"]; F = sc.flux.values
    src = sc.Source.values; itg = np.where(src == TGT)[0][0]; inb = np.where(src == NB)[0][0]
    imgs = {}
    near = np.where(np.hypot(X - X[itg], Y - Y[itg]) < 9)[0]
    tot = np.zeros((21, 21))
    for kk in near:
        im = F[kk] * prf.img(X[kk], Y[kk]); tot += im; imgs[kk] = im
    yy, xx = np.mgrid[0:21, 0:21]
    # observed pair depth (S13 PSF pair) and the PSF-split amplitudes
    dep = json.load(open(os.path.join(ROOT, "lc", "nominal", "depths.json")))["sectors"]["13"]
    dF_pair = dep["dF_t"] + dep["dF_n"]
    pred = {}
    for apname, (cx, cy) in (("centred_on_target", (X[itg], Y[itg])), ("centred_on_nb", (X[inb], Y[inb]))):
        for rad in (1.25, 2.0, 2.5):
            ap = np.hypot(xx - cx, yy - cy) <= rad
            Ftot = tot[ap].sum(); Ft_ = imgs[itg][ap].sum(); Fn_ = imgs[inb][ap].sum()
            # TIC-style contamination correction: keep only the named star's flux fraction
            own = Fn_ if apname == "centred_on_nb" else Ft_
            pred[f"{apname}_r{rad}"] = dict(
                frac_target=float(Ft_ / Ftot), frac_nb=float(Fn_ / Ftot),
                depth_H_target_raw=float(dF_pair * Ft_ / F[itg] / Ftot), depth_H_nb_raw=float(dF_pair * Fn_ / F[inb] / Ftot),
                depth_H_target_decontaminated=float(dF_pair * Ft_ / F[itg] / own), depth_H_nb_decontaminated=float(dF_pair * Fn_ / F[inb] / own),
                dilution_factor=float(Ftot / own))
    res["aperture_predictions_s13"] = pred
    res["note"] = ("QLP/eleanor flux fractions use the S13 PRF scene (TESS_PRF at Gaia positions, field registration);"
                   " 'decontaminated' = depth after removing all flux except the named star's own, i.e. what a"
                   " contamination-corrected light curve of that star would show.")
    json.dump(res, open(os.path.join(OUT, "neighbour_lcs.json"), "w"), indent=1, default=float)
    print(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
