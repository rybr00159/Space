"""Summary figures for the five vetted signals (static PNGs for the repo README)."""
import sys, json, warnings, os
sys.path.insert(0, "/home/user/work/src"); warnings.filterwarnings("ignore")
sys.argv = [sys.argv[0], "90086266"]
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from wotan import flatten
import stage2
from lcan import read, bin10, in_mask
import loc5
from loc5 import Cut, component_images, localise

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF, "text.color": INK,
                     "axes.labelcolor": INK, "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
                     "legend.frameon": False, "lines.linewidth": 1.8})


def lc(tic, sector, prov, col, P, t0, dur, win=1.0):
    p = stage2.products(tic); r = p[(p.sector == sector) & (p.prov == prov)].sort_values("exptime").iloc[0]
    t, q, cols, bkg, aper, crowd = read(stage2.fetch(r.uri), prov)
    ok = (q == 0) & np.isfinite(t) & np.isfinite(cols[col]) & (cols[col] > 0)
    tt, ff = bin10(t[ok], cols[col][ok] / np.nanmedian(cols[col][ok]), r.exptime)
    m = ~in_mask(tt, P, t0, dur, 0.8)
    tr = flatten(tt[m], ff[m], method="biweight", window_length=win, return_trend=True, break_tolerance=0.4)[1]
    g = np.isfinite(tr)
    return tt, ff / np.interp(tt, tt[m][g], tr[g])


def fold_binned(t, f, P, t0, half_h, nb=40):
    h = (((t - t0) / P + 0.5) % 1 - 0.5) * P * 24
    m = np.abs(h) < half_h; e = np.linspace(-half_h, half_h, nb + 1); i = np.digitize(h[m], e)
    x = np.array([(e[k - 1] + e[k]) / 2 for k in range(1, nb + 1) if np.sum(i == k) > 2])
    y = np.array([np.mean(f[m][i == k]) for k in range(1, nb + 1) if np.sum(i == k) > 2])
    return x, (y - 1) * 1e6


def snr_panel(ax, tic, sector, cfg, mark=None, title=""):
    C = Cut(tic, sector, cfg)
    comps, O = component_images(C, cfg)
    q = np.sqrt(sum((A / S) ** 2 for A, S in comps)) if cfg["mode"] == "fourier_free" else comps[0][0] / comps[0][1]
    v = np.nanmax(np.abs(q))
    im = ax.imshow(q, origin="lower", cmap="Greys" if cfg["mode"] == "fourier_free" else "RdBu_r",
                   vmin=0 if cfg["mode"] == "fourier_free" else -v, vmax=v)
    g = C.g; m = (g.x > -0.5) & (g.x < C.nx - 0.5) & (g.y > -0.5) & (g.y < C.ny - 0.5) & (g["T"] < g["T"][C.it] + 5)
    ax.scatter(g.x[m], g.y[m], s=10, facecolors="none", edgecolors=INK2, lw=0.6)
    ax.plot(C.xt, C.yt, "+", color=BLUE, ms=14, mew=2.2)
    ax.annotate("target", (C.xt, C.yt), xytext=(6, 6), textcoords="offset points", color=INK, fontsize=8)
    if mark is not None:
        r = g[g.source_id == mark].iloc[0]
        ax.plot(r.x, r.y, "o", ms=13, mfc="none", mec=ORANGE, mew=2)
        ax.annotate("Gaia DR3 …0577280\n(RS CVn variable)", (r.x, r.y), xytext=(8, -22), textcoords="offset points", color=INK, fontsize=8)
    ax.set_title(title); ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03, label="signal S/N per pixel")


os.makedirs("figs", exist_ok=True)
cfgs = json.load(open("loccfg.json"))

# ---- TIC 90086266
P, T0, T14 = 2.339469, 4180.4740, 0.128
fig, ax = plt.subplots(1, 3, figsize=(14, 4.2), gridspec_kw={"width_ratios": [1.15, 1.15, 1]})
for (sec, prov, col, c, lab) in ((13, "TESS-SPOC", "SAP_FLUX", BLUE, "2019 (S13)"), (104, "QLP", "SAP_FLUX", ORANGE, "2026 (S104)")):
    t, f = lc(90086266, sec, prov, col, P, T0, T14)
    x, y = fold_binned(t, f, P, T0, 4.5, 54)
    ax[0].plot(x, y, "-o", color=c, ms=3.5, lw=1.5, label=lab)
ax[0].set_xlabel("hours from mid-eclipse"); ax[0].set_ylabel("relative flux (ppm)"); ax[0].legend(loc="lower left")
ax[0].set_title("A 3% eclipse every 2.33947 d")
oe = json.load(open("out/90086266/oddeven.json"))
for k, (key, lab) in enumerate((("S13 TESS-SPOC SAP_FLUX", "2019 (S13)"), ("S104 QLP SAP_FLUX", "2026 (S104)"))):
    ev = np.array(oe[key]["events"]); x0 = k * 12 + (ev[:, 0] - ev[:, 0].min())
    for par, c, mk, nm in ((0, BLUE, "o", "even"), (1, ORANGE, "s", "odd")):
        s = ev[:, 0] % 2 == par
        ax[1].errorbar(x0[s], ev[s, 1] * 1e6, ev[s, 2] * 1e6, fmt=mk, color=c, ms=6, lw=1.2, label=nm if k == 0 else None)
    ax[1].text(k * 12 + 4.5, 0.96, lab, transform=ax[1].get_xaxis_transform(), ha="center", color=INK, fontsize=9)
    if k == 1: ax[1].axvline(10.5, color=INK2, lw=0.8, ls="--")
ax[1].set_xticks([]); ax[1].set_xlabel("successive eclipses"); ax[1].set_ylabel("eclipse depth (ppm, SAP flux)")
ax[1].legend(loc="lower left", ncol=2); ax[1].set_title("Alternate eclipses: equal in 2019, not in 2026")
snr_panel(ax[2], 90086266, 104, cfgs["90086266"], title="Dimming is centred on the target (S104)")
fig.tight_layout(); fig.savefig("figs/tic90086266.png", dpi=110); plt.close(fig)

if os.environ.get("ONLY90"): sys.exit(0)
# ---- TIC 61332759
fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.2), gridspec_kw={"width_ratios": [1, 1.1]})
snr_panel(ax[0], 61332759, 67, cfgs["61332759"], mark=6711181657530577280, title="Variability source (S67): the neighbour")
P2 = 1.480598
t, f = lc(61332759, 67, "QLP", "SAP_FLUX", P2, 3148.7857, 0.0001, win=4.0)
ph = ((t - 3148.7857) / P2) % 1; e = np.linspace(0, 1, 51); i = np.digitize(ph, e)
y = np.array([np.mean(f[i == k]) for k in range(1, 51)])
ax[1].plot((e[:-1] + e[1:]) / 2, (y - 1) * 1e6, "-o", color=BLUE, ms=3.5)
ax[1].set_xlabel("phase (P = 1.4806 d)"); ax[1].set_ylabel("relative flux in the QLP aperture (ppm)")
ax[1].set_title("Two unequal dips per 1.48-day cycle")
fig.tight_layout(); fig.savefig("figs/tic61332759.png", dpi=110); plt.close(fig)

# ---- TIC 61141412
ev = [e for e in json.load(open("out/61141412/etv_events.json")) if e["snr"] > 4]
k = np.array([e["k"] for e in ev]); tc = np.array([e["tc"] for e in ev]); s = np.array([e["err"] for e in ev])
A = np.column_stack([np.ones(len(k)), k]); W = 1 / s; c1, *_ = np.linalg.lstsq(A * W[:, None], tc * W, rcond=None)
A2 = np.column_stack([np.ones(len(k)), k, k ** 2]); c2, *_ = np.linalg.lstsq(A2 * W[:, None], tc * W, rcond=None)
oc = (tc - A @ c1) * 1440
fig, ax = plt.subplots(1, 2, figsize=(11, 4.0))
ax[0].errorbar(tc, oc, s * 1440, fmt="o", color=BLUE, ms=6, lw=1.2)
kk = np.linspace(k.min(), k.max(), 200)
ax[0].plot(c1[0] + c1[1] * kk, (np.column_stack([np.ones(200), kk, kk ** 2]) @ c2 - (c1[0] + c1[1] * kk)) * 1440, color=ORANGE, lw=1.8)
ax[0].annotate("quadratic fit", (c1[0] + c1[1] * kk[120], ((c2[0] + c2[1] * kk[120] + c2[2] * kk[120] ** 2) - (c1[0] + c1[1] * kk[120])) * 1440),
               xytext=(10, 12), textcoords="offset points", color=INK, fontsize=8)
ax[0].axhline(0, color=INK2, lw=0.8)
ax[0].set_xlabel("BTJD (2019 → 2026)"); ax[0].set_ylabel("observed − linear prediction (min)")
ax[0].set_title("Dip times drift: a third body")
for sec, prov, col, cc, lab in ((67, "TESS-SPOC", "PDCSAP_FLUX", BLUE, "S67 (2023)"), (27, "TESS-SPOC", "PDCSAP_FLUX", ORANGE, "S27 (2020)")):
    t, f = lc(61141412, sec, prov, col, 6.528604, 3143.003, 0.05, win=0.6)
    # fold on the quadratic ephemeris
    kk_ = np.round((t - c2[0]) / c2[1]); tpred = c2[0] + c2[1] * kk_ + c2[2] * kk_ ** 2
    h = (t - tpred) * 24; m = np.abs(h) < 4; e = np.linspace(-4, 4, 33); i = np.digitize(h[m], e)
    x = [(e[q - 1] + e[q]) / 2 for q in range(1, 33) if np.sum(i == q) > 2]
    y = [np.mean(f[m][i == q] - 1) * 1e6 for q in range(1, 33) if np.sum(i == q) > 2]
    ax[1].plot(x, y, "-o", color=cc, ms=3.5, label=lab)
ax[1].legend(loc="lower left"); ax[1].set_xlabel("hours from predicted dip"); ax[1].set_ylabel("relative flux (ppm)")
ax[1].set_title("~1-hour dips every 6.5285 d")
fig.tight_layout(); fig.savefig("figs/tic61141412.png", dpi=110); plt.close(fig)

# ---- giants
fig, ax = plt.subplots(1, 2, figsize=(11, 4.0))
for a, (tic, P, t0, dur, secs, title) in zip(ax, ((214311906, 5.69149, 3840.913, 0.067, (39, 66, 93, 104), "TIC 214311906 (Mira): flat-bottomed, 1.6 h"),
                                                   (55829034, 0.616146, 3837.8667, 0.085, (66, 93, 104), "TIC 55829034 (bulge giant): V-shaped, ~2 h"))):
    T, F = [], []
    for sec in secs:
        t, f = lc(tic, sec, "QLP", "SAP_FLUX", P, t0, dur, win=1.0); T.append(t); F.append(f)
    x, y = fold_binned(np.concatenate(T), np.concatenate(F), P, t0, 5 * dur * 24, 40)
    a.plot(x, y, "-o", color=BLUE, ms=3.5)
    a.set_xlabel("hours from mid-eclipse"); a.set_ylabel("relative flux (ppm)"); a.set_title(title)
fig.tight_layout(); fig.savefig("figs/giants.png", dpi=110); plt.close(fig)
print("done")
