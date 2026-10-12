"""Summary figures: injection leakage, robustness of f_target across photometry variants, per-transit amplitudes,
depth versus aperture position."""
import os, json, glob
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import ROOT, SECTORS

C_T, C_N = "#2a6fb0", "#c0392b"


def fig_injection():
    R = json.load(open(os.path.join(ROOT, "injection", "injection_results.json")))
    cfgs = ["same", "smear0.4", "gauss", "prf_ccd300", "pos+0.05", "pos-0.05"]
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), sharey=True)
    real = json.load(open(os.path.join(ROOT, "lc", "nominal", "depths.json")))
    for ax, host in zip(axs, ("target", "neighbour")):
        for j, c in enumerate(cfgs):
            rr = [r for r in R if r["host"] == host and r["inj_prf"] == c]
            ft = np.array([r["f_target_resp"] for r in rr])
            xs = j + np.linspace(-0.25, 0.25, len(ft))
            ax.plot(xs, ft, "o", color=C_T if host == "target" else C_N, ms=4, alpha=0.8)
        ax.axhline(1, color=C_T, ls=":", lw=1); ax.axhline(0, color=C_N, ls=":", lw=1)
        fr = real["combined"]["f_target"]
        ax.axhspan(fr["f_lo"], fr["f_hi"], color="0.85", zorder=0, label=f"real data (combined): f_target = {fr['f_best']:.2f}")
        ax.set_xticks(range(len(cfgs))); ax.set_xticklabels(["same PRF", "PRF smeared 0.4px", "Gaussian PSF", "PRF of CCD (300,300)",
                                                             "star +0.05px to nb", "star -0.05px from nb"], rotation=20, fontsize=8)
        ax.set_title(f"injected into the {host} (pixel level; 5 sectors x 2 phases)")
        ax.legend(fontsize=8, loc="center right")
    axs[0].set_ylabel("recovered f_target = dF_t / (dF_t + dF_n)")
    plt.tight_layout(); plt.savefig(os.path.join(ROOT, "injection_leakage.png"), dpi=100); plt.close()


def fig_variants():
    tags = sorted(glob.glob(os.path.join(ROOT, "lc", "*", "depths.json"))) + [os.path.join(ROOT, "lc", "nominal", "depths_regs.json")]
    rows = []
    for f in tags:
        d = json.load(open(f)); c = d["combined"]["f_target"]
        rows.append((d["tag"] + ("+regs" if d.get("regs") else ""), c["f_best"], c["f_lo"], c["f_hi"], d["combined"]["delta_chi2_target_minus_neighbour"]))
    order = ["nominal"] + sorted([r[0] for r in rows if r[0] != "nominal"])
    rows = sorted(rows, key=lambda r: order.index(r[0]))
    fig, ax = plt.subplots(figsize=(8, 0.38 * len(rows) + 1.5))
    for i, (tg, f, lo, hi, dc) in enumerate(rows):
        ax.errorbar(f, i, xerr=[[f - lo], [hi - f]], fmt="o", color="k" if tg == "nominal" else "0.35", ms=5)
        ax.text(2.02, i, f"dchi2(T-N) = {dc:+.0f}", va="center", fontsize=7)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows], fontsize=8)
    ax.axvline(1, color=C_T, ls="--", label="target hosts (f=1)"); ax.axvline(0, color=C_N, ls="--", label="neighbour hosts (f=0)")
    ax.set_xlim(-0.6, 2.0); ax.invert_yaxis(); ax.set_xlabel("f_target (combined over 5 sectors, 1-sigma)")
    ax.legend(fontsize=8, loc="lower right"); ax.set_title("Fraction of the transit flux decrement on the target: photometry variants", fontsize=9)
    plt.tight_layout(); plt.savefig(os.path.join(ROOT, "variants_f_target.png"), dpi=100); plt.close()


def fig_epochs():
    d = json.load(open(os.path.join(ROOT, "lc", "nominal", "depths.json")))
    fig, ax = plt.subplots(figsize=(6, 5))
    for k, (s, v) in enumerate(d["sectors"].items()):
        at = np.array(v["A_epoch_t"]); an = np.array(v["A_epoch_n"])
        ax.plot(an, at, "o", ms=5, label=f"S{s}", color=f"C{k}")
    lim = 40
    S = np.mean([v["dF_t"] + v["dF_n"] for v in d["sectors"].values()])
    ax.plot([0], [S], "*", color=C_T, ms=16, mec="k", label=f"expected if target hosts (dF = {S:.1f})")
    ax.plot([S], [0], "*", color=C_N, ms=16, mec="k", label="expected if neighbour hosts")
    x = np.linspace(-15, lim, 10)
    ax.plot(x, S - x, color="0.6", ls=":", lw=1, label="dF_t + dF_n = const (flux-transfer noise)")
    ax.axhline(0, color="0.85", lw=0.8); ax.axvline(0, color="0.85", lw=0.8)
    ax.set_xlabel("single-transit dF on neighbour PSF flux (e-/s)"); ax.set_ylabel("single-transit dF on target PSF flux (e-/s)")
    ax.set_xlim(-15, lim); ax.set_ylim(-25, lim); ax.legend(fontsize=7)
    ax.set_title("Per-transit flux decrements (nominal PSF photometry)", fontsize=9)
    plt.tight_layout(); plt.savefig(os.path.join(ROOT, "per_transit_dF.png"), dpi=100); plt.close()


def fig_apertures():
    fs = [os.path.join(ROOT, "aperture", f"aperture_s{s}.json") for s in SECTORS]
    if not all(os.path.exists(f) for f in fs): return
    names = ["1px_opposite", "1px_target", "1px_neighbour", "2x2_opposite", "2x2_target", "2x2_neighbour", "3x3_target"]
    fig, axs = plt.subplots(1, 5, figsize=(20, 4.2), sharey=False)
    for ax, f in zip(axs, fs):
        a = json.load(open(f)); ap = dict(a["apertures"])
        for k_, v_ in a.get("aliases", {}).items(): ap[k_] = ap[v_]
        y = np.array([ap[n]["depth"] for n in names]); e = np.array([ap[n]["depth_err"] for n in names])
        ft = np.array([ap[n]["frac_target"] for n in names]); fn = np.array([ap[n]["frac_nb"] for n in names])
        Fap = np.array([ap[n]["F_ap"] for n in names])
        tn = a["tests"]["named"]["eps0.05"]
        x = np.arange(len(names))
        ax.errorbar(x, y * 1e3, e * 1e3, fmt="ko", ms=5, label="measured", zorder=3)
        ax.plot(x, tn["S_target"] * ft / Fap * 1e3, "s-", color=C_T, ms=5, label=f"target hosts (chi2 {tn['chi2_target']:.0f})")
        ax.plot(x, tn["S_nb"] * fn / Fap * 1e3, "^-", color=C_N, ms=5, label=f"neighbour hosts (chi2 {tn['chi2_nb']:.0f})")
        lab = [n.replace("_", "\n") + ("\n(=" + a["aliases"][n].split("_")[1] + ")" if n in a.get("aliases", {}) else "") for n in names]
        ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=7)
        ax.set_title(f"S{a['sector']}: depth vs aperture (dof {tn['dof']})", fontsize=9); ax.legend(fontsize=7)
        ax.set_ylabel("transit depth in aperture (ppt)")
    plt.tight_layout(); plt.savefig(os.path.join(ROOT, "depth_vs_aperture.png"), dpi=100); plt.close()
    # pixel-map localisation
    fig, axs = plt.subplots(1, 5, figsize=(20, 4.2))
    for ax, s in zip(axs, SECTORS):
        g = np.load(os.path.join(ROOT, "aperture", f"grid_s{s}.npz"))
        a = json.load(open(os.path.join(ROOT, "aperture", f"aperture_s{s}.json")))
        grid, gx, nb = g["grid"], g["gx"], g["nb"]
        dchi = grid - grid.min()
        h = ax.contourf(gx, gx, np.clip(dchi, 0, 60), levels=[0, 2.3, 6.2, 11.8, 20, 40, 60], cmap="Greys_r")
        ax.plot(0, 0, "+", color=C_T, ms=14, mew=2, label="target"); ax.plot(nb[0], nb[1], "x", color=C_N, ms=11, mew=2, label="neighbour")
        fp = a["free_position"]; ax.plot(fp["dx_px"], fp["dy_px"], "o", mfc="none", color="orange", ms=8, label="best fit")
        ax.set_aspect("equal"); ax.set_xlabel("column offset from target (px)"); ax.set_ylabel("row offset (px)")
        ax.set_title(f"S{s} pixel-map source position (dchi2 contours 2.3/6.2/11.8)", fontsize=8); ax.legend(fontsize=7)
    plt.tight_layout(); plt.savefig(os.path.join(ROOT, "pixelmap_localisation.png"), dpi=100); plt.close()


if __name__ == "__main__":
    for f in (fig_injection, fig_variants, fig_epochs, fig_apertures):
        try:
            f(); print("ok", f.__name__)
        except Exception as e:
            print("failed", f.__name__, e)


def fig_pixmap_scatter(sub="aperture"):
    fig, axs = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
    for k, s in enumerate(SECTORS):
        p = np.load(os.path.join(ROOT, sub, f"pixmap_s{s}.npz"))
        a = json.load(open(os.path.join(ROOT, sub, f"aperture_s{s}.json")))
        names = list(p["names"]); sel = np.array([i for i, n in enumerate(names) if n.startswith("pix_")])
        tt = a["tests"]["pixelmap"]["eps0.05"]
        y = p["y"][sel] / p["Fap"][sel] * 1e3; e = p["err"][sel] / p["Fap"][sel] * 1e3
        pt = tt["S_target"] * p["fr_t"][sel] / p["Fap"][sel] * 1e3; pn = tt["S_nb"] * p["fr_n"][sel] / p["Fap"][sel] * 1e3
        bright = p["Fap"][sel] > 0.08 * p["Fap"][sel].max()
        for ax, pr in zip(axs, (pt, pn)):
            ax.errorbar(pr[bright], y[bright], e[bright], fmt="o", ms=4, color=f"C{k}", label=f"S{s}", alpha=0.85)
    tot_t = sum(json.load(open(os.path.join(ROOT, sub, f"aperture_s{s}.json")))["tests"]["pixelmap"]["eps0.05"]["chi2_target"] for s in SECTORS)
    tot_n = sum(json.load(open(os.path.join(ROOT, sub, f"aperture_s{s}.json")))["tests"]["pixelmap"]["eps0.05"]["chi2_nb"] for s in SECTORS)
    for ax, nm, c2 in zip(axs, ("target", "neighbour"), (tot_t, tot_n)):
        lim = [0, 6.5]; ax.plot(lim, lim, "k--", lw=1); ax.set_xlim(lim); ax.set_ylim(-1, 7)
        ax.set_xlabel(f"predicted pixel depth if the {nm} hosts (ppt)")
        ax.set_title(f"single pixels within 2.5 px (flux > 8% of peak); total chi2 = {c2:.0f}", fontsize=9)
    axs[0].set_ylabel("measured pixel depth (ppt)"); axs[1].legend(fontsize=8)
    plt.tight_layout(); plt.savefig(os.path.join(ROOT, "pixel_depth_vs_prediction.png"), dpi=100); plt.close()


if __name__ == "__main__":
    fig_pixmap_scatter()
    print("ok fig_pixmap_scatter")
