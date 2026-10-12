"""Task 5: secondary eclipse and odd/even depths in the separated PSF light curves (nominal photometry).

Quantities are flux decrements dF (e-/s) from depthfit (template shape and duration of the primary, per-transit quadratic
baseline), converted to depths relative to the Gaia-predicted flux of each star (F_pair_PSF * r/(1+r), r = Gaia/TIC flux
ratio) -- the 'catalogue normalisation' that is insensitive to the PSF flux-split bias of the static image.
Errors: control-phase / epoch-scatter / white maximum (depthfit.summarize).
 - secondary at phase 0.5 (circular orbit);
 - eccentric secondary: scan of phases 0.30-0.70 (step 0.005) with the primary's duration, per star, combined over sectors;
 - odd/even: depth on the 2P ephemeris for odd and even epochs.
"""
import os, sys, json
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import *
from depthfit import summarize, depth_fit
from analyze_lc import masks


def combine(vals, errs):
    v = np.array(vals); e = np.array(errs); w = 1 / e ** 2
    m = np.sum(w * v) / np.sum(w)
    return float(m), float(1 / np.sqrt(np.sum(w))), float(np.sum(w * (v - m) ** 2))


def main(tag="nominal"):
    od = os.path.join(ROOT, "lc", tag)
    res = dict(per_sector={}, combined={})
    acc = {k: ([], []) for k in ["sec_t", "sec_n", "odd_t", "odd_n", "even_t", "even_n", "pri_t", "pri_n"]}
    scan_ph = np.arange(0.30, 0.7001, 0.005)
    scan = {"t": [], "n": [], "te": [], "ne": []}
    for s in SECTORS:
        d = np.load(os.path.join(od, f"psf_s{s}.npz"))
        js = json.load(open(os.path.join(ROOT, "static", f"static_s{s}.json")))
        m = masks(d)
        t = d["t"][m]; ex = float(d["exptime"])
        ft = d["flux"][m, int(d["jt"])]; fn = d["flux"][m, int(d["jn"])]
        msum = np.median(ft + fn); rg = js["ratio_gaia"]
        Fn_cat = msum * rg / (1 + rg); Ft_cat = msum / (1 + rg)
        out = {}
        # primary (for reference) and secondary
        for name, t0x, Px in [("pri", T0, P), ("sec", T0 + 0.5 * P, P), ("odd", T0 + P, 2 * P), ("even", T0, 2 * P)]:
            r = summarize(t, [ft, fn], ex, ["t", "n"], P_=Px, t0=t0x)
            out[name] = dict(dFt=r["t"]["A"], dFt_err=r["t"]["err"], dFn=r["n"]["A"], dFn_err=r["n"]["err"],
                             depth_t_cat=r["t"]["A"] / Ft_cat, depth_t_cat_err=r["t"]["err"] / Ft_cat,
                             depth_n_cat=r["n"]["A"] / Fn_cat, depth_n_cat_err=r["n"]["err"] / Fn_cat,
                             n_tr=r["t"]["n_tr"])
            acc[f"{name}_t"][0].append(out[name]["depth_t_cat"]); acc[f"{name}_t"][1].append(out[name]["depth_t_cat_err"])
            acc[f"{name}_n"][0].append(out[name]["depth_n_cat"]); acc[f"{name}_n"][1].append(out[name]["depth_n_cat_err"])
        # phase scan (white-noise errors scaled by the sector's control/white ratio at phase 0.5)
        ksc_t = out["sec"]["dFt_err"] / max(depth_fit(t, ft, ex, P, T0 + 0.5 * P)["err"], 1e-9)
        ksc_n = out["sec"]["dFn_err"] / max(depth_fit(t, fn, ex, P, T0 + 0.5 * P)["err"], 1e-9)
        st, sn, ste, sne = [], [], [], []
        for ph in scan_ph:
            a = depth_fit(t, ft, ex, P, T0 + ph * P); b = depth_fit(t, fn, ex, P, T0 + ph * P)
            st.append(a["A"] / Ft_cat); sn.append(b["A"] / Fn_cat)
            ste.append(a["err"] * ksc_t / Ft_cat); sne.append(b["err"] * ksc_n / Fn_cat)
        scan["t"].append(st); scan["n"].append(sn); scan["te"].append(ste); scan["ne"].append(sne)
        res["per_sector"][str(s)] = out
        print(f"S{s}: sec dFn {out['sec']['dFn']:.2f}+-{out['sec']['dFn_err']:.2f} (nb cat depth {out['sec']['depth_n_cat']*1e3:.2f}+-{out['sec']['depth_n_cat_err']*1e3:.2f} ppt) "
              f"| odd n {out['odd']['depth_n_cat']*1e3:.1f}+-{out['odd']['depth_n_cat_err']*1e3:.1f} even n {out['even']['depth_n_cat']*1e3:.1f}+-{out['even']['depth_n_cat_err']*1e3:.1f} ppt", flush=True)
    for k, (v, e) in acc.items():
        if all(np.isfinite(v)):
            mm, ee, c2 = combine(v, e); res["combined"][k] = dict(depth=mm, err=ee, chi2=c2, dof=len(v) - 1)
    c = res["combined"]
    for st_ in ("t", "n"):
        dif = c[f"odd_{st_}"]["depth"] - c[f"even_{st_}"]["depth"]
        de = np.hypot(c[f"odd_{st_}"]["err"], c[f"even_{st_}"]["err"])
        c[f"odd_minus_even_{st_}"] = dict(diff=dif, err=float(de), sigma=float(dif / de))
        c[f"secondary_over_primary_{st_}"] = dict(ratio=c[f"sec_{st_}"]["depth"] / c[f"pri_{st_}"]["depth"],
                                                   upper2sig_depth=c[f"sec_{st_}"]["depth"] + 2 * c[f"sec_{st_}"]["err"])
    # combined scan
    for st_ in ("t", "n"):
        V = np.array(scan[st_]); E = np.array(scan[st_ + "e"]); W = 1 / E ** 2
        mv = np.sum(W * V, 0) / np.sum(W, 0); me = 1 / np.sqrt(np.sum(W, 0))
        snr = mv / me; k = np.argmax(np.abs(snr))
        c[f"scan_{st_}"] = dict(max_abs_snr=float(snr[k]), phase=float(scan_ph[k]), depth=float(mv[k]), err=float(me[k]),
                                 n_phases=len(scan_ph), max_depth_2sig=float(np.max(mv + 2 * me)))
        res[f"scan_{st_}"] = dict(phase=scan_ph.tolist(), depth=mv.tolist(), err=me.tolist())
    print(json.dumps(res["combined"], indent=1))
    json.dump(res, open(os.path.join(ROOT, "lc", tag, "secondary_oddeven.json"), "w"), indent=1, default=float)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    for a, st_, nm in zip(ax, ("t", "n"), ("target", "neighbour")):
        sc_ = res[f"scan_{st_}"]
        a.errorbar(np.array(sc_["phase"]), np.array(sc_["depth"]) * 1e3, np.array(sc_["err"]) * 1e3, fmt=".", ms=3, color="k")
        a.axhline(0, color="0.5"); a.axvline(0.5, color="C3", ls=":")
        a.axhline(c[f"pri_{st_}"]["depth"] * 1e3, color="C0", ls="--", label=f"primary depth {c[f'pri_{st_}']['depth']*1e3:.1f} ppt")
        a.set_xlabel("phase"); a.set_ylabel(f"eclipse depth rel. to Gaia-predicted {nm} flux (ppt)")
        a.set_title(f"{nm} PSF light curve: secondary-eclipse scan (all sectors)"); a.legend(fontsize=8)
    plt.tight_layout(); plt.savefig(os.path.join(ROOT, "secondary_scan_psf.png"), dpi=90); plt.close()


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "nominal")
