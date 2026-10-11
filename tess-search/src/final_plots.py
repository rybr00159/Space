"""Standardised figures for the final EB table: every system folded at its adopted period and epoch."""
import os, sys, warnings
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
import stage2
from eb_analysis import fold_profile

W = "/home/user/work"


def figure(row, out, ax_row=None):
    secs = stage2.load_all(int(row.tic))
    t = np.concatenate([s["t"] for s in secs]); f = np.concatenate([s["f"] for s in secs])
    P, T0 = float(row.P), float(row.T0_BTJD)
    ph = ((t - T0) / P + 0.25) % 1.0 - 0.25
    own = ax_row is None
    if own:
        fig, ax_row = plt.subplots(1, 3, figsize=(13, 3.4))
    pb, prof, _ = fold_profile(t, f, P, T0 - 0.25 * P, nb=400)
    ax_row[0].plot(ph, f, ",", color="0.6"); ax_row[0].plot(pb - 0.25, prof, "-", lw=0.8, color="C0")
    ax_row[0].set_title(f"TIC {row.tic}  P = {P:.6f} d", fontsize=9); ax_row[0].set_xlabel("phase", fontsize=8)
    dur = 0.08 * P if P < 2 else 0.25
    for a, c0, lab in ((ax_row[1], 0.0, "primary"), (ax_row[2], float(row.phase_sec), "secondary")):
        x = (((t - T0) / P - c0 + 0.5) % 1.0 - 0.5) * P * 24
        m = np.abs(x) < 3 * dur * 24
        a.plot(x[m], f[m], ".", ms=1, color="0.6")
        hb = np.linspace(-3 * dur * 24, 3 * dur * 24, 61); idx = np.digitize(x[m], hb)
        bc = [(hb[q - 1] + hb[q]) / 2 for q in range(1, len(hb)) if (idx == q).sum() > 3]
        bv = [f[m][idx == q].mean() for q in range(1, len(hb)) if (idx == q).sum() > 3]
        a.plot(bc, bv, "o-", ms=2, lw=0.8, color="C0"); a.set_xlabel("hours", fontsize=8)
        a.set_title(f"{lab} (phase {c0:.3f})", fontsize=9)
    if own:
        fig.tight_layout(); fig.savefig(out, dpi=100); plt.close(fig)


if __name__ == "__main__":
    d = pd.read_csv(f"{W}/vet/new_ebs_final.csv")
    d = d.sort_values(["kind", "primary_pct"], ascending=[True, False])
    os.makedirs(f"{W}/vet/final_figs", exist_ok=True)
    for r in d.itertuples():
        figure(r, f"{W}/vet/final_figs/EB_TIC{r.tic}.png")
    # contact sheets: full-phase panel only, 15 per page
    for page in range(0, len(d), 15):
        sub = d.iloc[page:page + 15]; nr = (len(sub) + 2) // 3
        fig, axes = plt.subplots(nr, 3, figsize=(15, 2.6 * nr))
        for a in axes.ravel(): a.axis("off")
        for a, r in zip(axes.ravel(), sub.itertuples()):
            a.axis("on")
            secs = stage2.load_all(int(r.tic))
            t = np.concatenate([s["t"] for s in secs]); f = np.concatenate([s["f"] for s in secs])
            ph = ((t - r.T0_BTJD) / r.P + 0.25) % 1.0 - 0.25
            pb, prof, _ = fold_profile(t, f, r.P, r.T0_BTJD - 0.25 * r.P, nb=300)
            a.plot(ph, f, ",", color="0.6"); a.plot(pb - 0.25, prof, "-", lw=0.8)
            a.set_title(f"TIC {r.tic}  P={r.P:.5f} d  {r.kind}", fontsize=8); a.tick_params(labelsize=7)
        fig.tight_layout(); fig.savefig(f"{W}/vet/final_figs/contact_sheet_{page // 15 + 1}.png", dpi=85); plt.close(fig)
    print("done")
