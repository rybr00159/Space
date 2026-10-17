"""c13: figures (figs/*.png). INV = blue (#2a78d6), REF = orange (#eb6834), uncleaned = gray; reference palette slots 1-2."""
import sys
sys.path.insert(0, "tic214311908/analysis/clean/scripts")
from cleanlib import *
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tfit import Model, mid_depth

C_INV, C_REF, C_GRAY, C_INK, C_MUTED = "#2a78d6", "#eb6834", "#8a8984", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": C_MUTED, "axes.labelcolor": C_INK, "xtick.color": C_MUTED, "ytick.color": C_MUTED,
                     "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False})
EPH = json.load(open(f"{CL}/out/ephemeris_final.json"))["combined"]
P, T0 = EPH["P_d"], EPH["T0_epoch0_BTJD_TDB"]


def fold_bins(t, y, Pp, Tm, nb=25):
    ph = ((t - Tm) / Pp) % 1; k = np.floor(ph * nb).astype(int)
    return (np.arange(nb) + 0.5) / nb, np.array([np.nanmean(y[k == j]) if np.any(k == j) else np.nan for j in range(nb)]), \
        np.array([np.nanstd(y[k == j]) / np.sqrt(max(np.sum(k == j), 1)) for j in range(nb)])


def fig_rrab():
    fig, axs = plt.subplots(2, 3, figsize=(11, 5.5), sharex=True); axs = axs.ravel()
    L = json.load(open(f"{CL}/out/leak_inv.json"))
    for i, s in enumerate(SECTORS):
        ax = axs[i]
        for mode, col, lab in (("none", C_GRAY, "phase-1 (uncleaned)"), ("all", C_INV, "cleaned")):
            z = np.load(f"{CL}/data/lcinv_{mode}_s{s}.npz")
            ok = z["quality_ok"].astype(bool) & ~z["intransit_window"].astype(bool) & np.isfinite(z["det_M3V_B"])
            x, y, e = fold_bins(z["t_tdb"][ok], (z["det_M3V_B"][ok] - 1) * 1e3, P_RR, TMAX_RR)
            ax.errorbar(x, y, e, fmt="o-", ms=3, lw=1.2, color=col, label=f"{lab}: A1 {L[str(s)][mode]['M3V_B']['RR']['amp1']:.2f} ppt")
        ax.axhline(0, color=C_MUTED, lw=0.6)
        ax.set_title(f"S{s}", fontsize=9, color=C_INK); ax.legend(fontsize=7, frameon=False, loc="lower left")
        if i % 3 == 0: ax.set_ylabel("M3V flux - 1 (ppt)")
        if i >= 3: ax.set_xlabel("RRab phase (Gaia P = 0.6485275 d, 0 = max light)")
    fig.suptitle("INV M3V light curve folded on the RRab period: before / after excluding variables from the field estimators", fontsize=10)
    fig.tight_layout(); fig.savefig(f"{CL}/figs/rrab_leak.png", dpi=130); plt.close(fig)


def fig_fold():
    fig, axs = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, pipe, fn, col in ((axs[0], "INV", "win_inv_raw.npz", C_INV), (axs[1], "REF", "win_ref.npz", C_REF)):
        z = np.load(f"{CL}/data/{fn}")
        g = json.load(open(f"{CL}/out/global_{pipe.lower()}_single.json"))
        m = g["map"]
        # remove per-window linear baseline (fit out of transit)
        f = z["f"].copy(); t = z["t"]; w = z["w"]
        for i in np.unique(w):
            ii = w == i; tc = z["tc0"][i]; oot = ii & (np.abs(t - tc) > 0.07)
            c = np.polyfit(t[oot] - tc, f[oot], 1); f[ii] = f[ii] / np.polyval(c, t[ii] - tc)
        ph = (((t - m["T0"]) / m["P"] + 0.5) % 1 - 0.5) * m["P"] * 24
        ax.plot(ph, (f - 1) * 1e3, ".", ms=1.5, color=C_GRAY, alpha=0.35, rasterized=True)
        bins = np.arange(-7.2, 7.21, 1 / 6)
        k = np.digitize(ph, bins)
        bx = 0.5 * (bins[1:] + bins[:-1]); by = np.array([np.mean((f[k == j] - 1) * 1e3) if np.sum(k == j) > 3 else np.nan for j in range(1, len(bins))])
        be = np.array([np.std((f[k == j] - 1) * 1e3) / np.sqrt(np.sum(k == j)) if np.sum(k == j) > 3 else np.nan for j in range(1, len(bins))])
        ax.errorbar(bx, by, be, fmt="o", ms=3.5, color=col, lw=1, label="10-min bins")
        tt = np.linspace(-0.3, 0.3, 2000) + m["T0"]
        mod = Model(tt, 120 / 86400).flux(m["T0"], m["P"], m["k"], m["b"], m["aR"], m["u1"], m["u2"])
        ax.plot((tt - m["T0"]) * 24, (mod - 1) * 1e3, color=C_INK, lw=1.5, label=f"batman (rho prior): depth {g['depth_mid_ppt']['med']:.1f} ppt")
        ax.set_xlim(-7.2, 7.2); ax.set_ylim(-35, 20); ax.set_xlabel("hours from mid-transit"); ax.set_title(f"{pipe} cleaned ({len(np.unique(w))} transits)", fontsize=10)
        ax.legend(fontsize=8, frameon=False, loc="lower left")
    axs[0].set_ylabel("relative flux - 1 (ppt of 979 e-/s)")
    fig.tight_layout(); fig.savefig(f"{CL}/figs/fold_clean.png", dpi=130); plt.close(fig)


def fig_depth_epoch():
    D = json.load(open(f"{CL}/out/depth_epoch.json"))
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    ax = axs[0]
    for k_, (pipe, col, off) in enumerate((("inv", C_INV, -0.12), ("ref", C_REF, 0.12))):
        d = D[pipe]["per_sector"]
        secs = sorted(int(s) for s in d)
        x = np.arange(len(secs)) + off
        ax.errorbar(x, [d[str(s)][0] for s in secs], [d[str(s)][1] for s in secs], fmt="o", ms=6, color=col, lw=1.5, capsize=0,
                    label=f"{pipe.upper()} cleaned (chi2 {D[pipe]['chi2']:.1f}/5)")
        ax.axhline(D[pipe]["all"][0], color=col, lw=0.8, ls="--")
    ax.set_xticks(range(6)); ax.set_xticklabels([f"S{s}" for s in SECTORS]); ax.set_ylabel("mid-transit depth (ppt)")
    ax.set_title("Per-sector depth, fixed shape, fake-epoch errors", fontsize=10); ax.legend(fontsize=8, frameon=False)
    ax = axs[1]
    I = pd.read_csv(f"{CL}/out/inject_summary.csv")
    for pipe, col, mk in (("inv", C_INV, "o"), ("ref", C_REF, "s")):
        for var, ls in (("invprf" if pipe == "inv" else "refprf", "-"), ("sx+", ":"), ("sx-", ":"), ("blur", "--")):
            x = I[I.variant == var].sort_values("sector")
            ax.plot(range(6), x[f"{pipe}_resp_mean"].values, marker=mk, ms=4, ls=ls, lw=1, color=col,
                    label=f"{pipe.upper()} {var}" if var in ("invprf", "refprf", "blur") or (var == "sx+") else None)
    obs = D["inv"]["per_sector"]; mu = D["inv"]["all"][0]
    ax.errorbar(np.arange(6) + 0.15, [obs[str(s)][0] / 15.7 for s in SECTORS], [obs[str(s)][1] / 15.7 for s in SECTORS], fmt="D", ms=4,
                color=C_GRAY, lw=1, label="observed INV / 15.7 ppt")
    ax.set_xticks(range(6)); ax.set_xticklabels([f"S{s}" for s in SECTORS]); ax.set_ylabel("recovered / injected depth")
    ax.set_title("Pixel-level injection of a constant 15.7-ppt transit", fontsize=10); ax.legend(fontsize=7, frameon=False, ncol=2)
    fig.tight_layout(); fig.savefig(f"{CL}/figs/depth_epoch.png", dpi=130); plt.close(fig)


def fig_oc():
    D = pd.read_csv(f"{CL}/out/per_transit_clean.csv")
    fig, ax = plt.subplots(figsize=(9, 3.6))
    for pipe, col, off in (("inv", C_INV, -1.5), ("ref", C_REF, 1.5)):
        d = D[D.pipeline == pipe]
        oc = (d.tc - (T0 + P * d.epoch)) * 1440
        ax.errorbar(d.epoch + off, oc, d.tc_err_cal * 1440, fmt="o", ms=4, color=col, lw=1, label=f"{pipe.upper()} (injection-calibrated errors)")
    ax.axhline(0, color=C_MUTED, lw=0.8); ax.set_xlabel("epoch"); ax.set_ylabel("O - C (min), combined ephemeris"); ax.legend(fontsize=8, frameon=False)
    fig.tight_layout(); fig.savefig(f"{CL}/figs/oc_clean.png", dpi=130); plt.close(fig)


def fig_scan():
    fig, axs = plt.subplots(2, 1, figsize=(10, 5), sharex=True)
    for ax, pipe, col in ((axs[0], "inv", C_INV), (axs[1], "ref", C_REF)):
        for suf, ls, a in (("", "-", 1.0), ("_noS12", ":", 0.8)):
            d = pd.read_csv(f"{CL}/out/phase_scan_{pipe}{suf}.csv")
            ax.plot(d.phase, d.depth_ppt / (d.err_ppt * (d.depth_ppt / d.err_ppt)[np.abs(((d.phase + 0.5) % 1) - 0.5) > 0.03].pipe(robust_std)),
                    ls=ls, color=col, lw=1, alpha=a, label=f"{pipe.upper()}{' without S12' if suf else ''}")
        for x_ in (0.4925, 0.576): ax.axvline(x_, color=C_MUTED, lw=0.7, ls="--")
        ax.axhline(0, color=C_MUTED, lw=0.6); ax.set_ylabel("dip S/N (calibrated)"); ax.legend(fontsize=8, frameon=False, loc="upper right")
    axs[1].set_xlabel("orbital phase")
    fig.suptitle("Phase scans of the cleaned light curves (dashed: phases 0.4925 and 0.576 flagged in phase 2)", fontsize=10)
    fig.tight_layout(); fig.savefig(f"{CL}/figs/phase_scan.png", dpi=130); plt.close(fig)


if __name__ == "__main__":
    for f in (fig_rrab, fig_fold, fig_depth_epoch, fig_oc, fig_scan):
        try:
            f(); print("ok", f.__name__)
        except Exception as ex:
            print("FAILED", f.__name__, repr(ex))
