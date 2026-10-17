import sys; sys.path.insert(0, ".")
from photlib import *
from wotan import flatten
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
model = sys.argv[1] if len(sys.argv) > 1 else "B"
fig, axs = plt.subplots(6, 2, figsize=(16, 22), gridspec_kw=dict(width_ratios=[3, 1]))
def tmpl(t):
    ph = ((t - T0_EPH) / P_EPH + 0.5) % 1 - 0.5; a = np.abs(ph * P_EPH); T14 = 0.067; T23 = 0.78 * T14
    return np.where(a <= T23 / 2, 1.0, np.where(a >= T14 / 2, 0.0, (T14 / 2 - a) / ((T14 - T23) / 2)))
for row, s in enumerate(SECTORS):
    z = np.load(f"{ROOT}/data/phot_s{s}.npz"); J = json.load(open(f"{ROOT}/out/phot_s{s}.json"))
    t = z["t"]; ok = z["qok"] & (z["bkg"] < J["bkg_lim"]) & z["jit_ok"]
    ph = ((t - T0_EPH) / P_EPH + 0.5) % 1 - 0.5
    y = tmpl(t)
    lines = []
    for k, (mk, nm) in enumerate((("B", "M3V"), ("Bp", "pair"), ("B", "Mira"), ("A", "M3V"), ("A", "pair"), ("A", "Mira"))):
        if model != "all" and mk not in ("B", "Bp"): continue
        names = list(z[f"names_{mk}"]); f = z[f"F_{mk}"][:, names.index(nm)]
        m = ok & np.isfinite(f) & ~(np.abs(ph * P_EPH) < 0.1)
        tr = flatten(t[m], f[m], method="biweight", window_length=0.75, return_trend=True, break_tolerance=0.5, edge_cutoff=0)[1]
        kk = np.isfinite(tr); r = f - np.interp(t, t[m][kk], tr[kk])
        s_ = robust_std(r[m])
        sel = ok & np.isfinite(r)
        A = np.column_stack([np.ones(sel.sum()), -y[sel]]); c, *_ = np.linalg.lstsq(A, r[sel], rcond=None)
        e = s_ / np.sqrt(np.sum((y[sel] - y[sel].mean()) ** 2))
        lines.append(f"{mk}:{nm} rstd {s_:.1f} dip {c[1]:.1f}+-{e:.1f} e/s")
        if mk in ("B", "Bp") and nm in ("M3V", "pair"):
            off = 0 if nm == "M3V" else 80
            w = 1 / 24 if s != 12 else 1 / 12
            tb = np.arange(t[ok].min(), t[ok].max() + w, w); idx = np.digitize(t[ok], tb)
            rb = np.array([np.nanmean(r[ok][idx == i]) if np.sum(idx == i) >= 2 else np.nan for i in range(1, len(tb))])
            axs[row, 0].plot(tb[:-1] + w / 2, rb + off, ".", ms=3, label=nm)
            pb = np.arange(-0.3, 0.3001, 0.01); x = ph[ok] * P_EPH; xb = np.digitize(x, pb)
            yb = np.array([np.nanmean(r[ok][xb == i]) if np.sum(xb == i) >= 2 else np.nan for i in range(1, len(pb))])
            axs[row, 1].plot((pb[:-1] + 0.005) * 24, yb + off * 0.6, ".-", ms=4, label=nm)
    print(f"S{s}: " + "; ".join(lines))
    for tc in T0_EPH + P_EPH * np.arange(-300, 100):
        if t[ok].min() < tc < t[ok].max(): axs[row, 0].axvline(tc, color="r", lw=0.5)
    axs[row, 0].set_title(f"S{s} (model B: M3V; Bp: pair +80)"); axs[row, 0].set_ylim(-60, 140); axs[row, 1].set_ylim(-50, 90)
    for v in (-0.8, 0.8): axs[row, 1].axvline(v, color="r", lw=0.5)
axs[0, 0].legend(); fig.tight_layout(); fig.savefig(f"{ROOT}/explore/e06_quick2_{model}.png", dpi=60)
