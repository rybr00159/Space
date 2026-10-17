import sys; sys.path.insert(0, ".")
from photlib import *
from wotan import flatten
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, axs = plt.subplots(6, 2, figsize=(16, 20), gridspec_kw=dict(width_ratios=[3, 1]))
for row, s in enumerate(SECTORS):
    z = np.load(f"{ROOT}/data/phot_s{s}.npz"); names = list(z["names"])
    t = z["t"]; ok = z["qok"] & (z["bkg"] < 3 * np.median(z["bkg"][z["q0"]]))
    ph = ((t - T0_EPH) / P_EPH + 0.5) % 1 - 0.5
    intr = np.abs(ph * P_EPH) < 0.6 * T14_EPH
    for k, nm in enumerate(("M3V", "pair", "Mira")):
        f = z["F"][:, names.index(nm)]
        m = ok & np.isfinite(f) & ~(np.abs(ph * P_EPH) < 0.1)
        tr = flatten(t[m], f[m], method="biweight", window_length=0.75, return_trend=True, break_tolerance=0.3, edge_cutoff=0)[1]
        kk = np.isfinite(tr); r = f - np.interp(t, t[m][kk], tr[kk])
        s_ = robust_std(r[m])
        # box depth
        sel = ok & intr & np.isfinite(r); d = -np.mean(r[sel]); e = s_ / np.sqrt(sel.sum())
        print(f"S{s} {nm}: rstd {s_:.1f} e/s per cadence; in-transit mean dip {d:.1f} +- {e:.1f} e/s (n={sel.sum()})")
        if nm in ("M3V", "pair"):
            tb = np.arange(t[ok].min(), t[ok].max(), 1 / 24)
            idx = np.digitize(t[ok], tb); rb = np.array([np.nanmean(r[ok][idx == i]) if np.sum(idx == i) > 2 else np.nan for i in range(1, len(tb))])
            axs[row, 0].plot(tb[:-1], rb + (0 if nm == "M3V" else 60), ".", ms=3, label=nm)
            pb = np.linspace(-0.5, 0.5, 201) * 0.6
            x = ph[ok] * P_EPH; xb = np.digitize(x, pb)
            yb = np.array([np.nanmean(r[ok][xb == i]) if np.sum(xb == i) > 3 else np.nan for i in range(1, len(pb))])
            axs[row, 1].plot(pb[:-1] * 24, yb + (0 if nm == "M3V" else 40), ".-", ms=3, label=nm)
    for tc in T0_EPH + P_EPH * np.arange(-300, 100):
        if t[ok].min() < tc < t[ok].max(): axs[row, 0].axvline(tc, color="r", lw=0.5)
    axs[row, 0].set_title(f"S{s}"); axs[row, 0].set_ylim(-80, 120); axs[row, 1].set_ylim(-60, 80); axs[row, 1].axvline(-0.8, color='r', lw=.5); axs[row, 1].axvline(0.8, color='r', lw=.5)
axs[0, 0].legend(); fig.tight_layout(); fig.savefig(f"{ROOT}/explore/e03_quick.png", dpi=60)
