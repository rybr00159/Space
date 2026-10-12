"""Summary figure of the pixel-level localisation: per sector, the fractional flux drop per pixel (dF / O) in the
central 3x3 pixels, observed vs. predicted (TESS PRF + Gaia) for the signal on the target and on the 9" neighbour,
plus the fitted source positions. Writes localization_summary.png."""
import os, json, glob
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import prf_localize as pl
from prf_localize import gaia, TARGET_ID, NB_ID

HERE = os.path.dirname(os.path.abspath(__file__))
tgt = gaia[gaia.source_id == TARGET_ID].iloc[0]; nbr = gaia[gaia.source_id == NB_ID].iloc[0]
H = pl.HALF
files = sorted(glob.glob("/home/user/work/tic87/data/tesscut/*.fits"))
fig, axes = plt.subplots(3, len(files), figsize=(3.0 * len(files), 9.0))
rows = []
for j, f in enumerate(files):
    S = pl.load_sector(f); M = pl.PRFModel(S)
    dF, edF, O = pl.dF_map(S, S["blocks"])
    xt, yt = map(float, S["w"].world_to_pixel_values(tgt.ra, tgt.dec)); xn, yn = map(float, S["w"].world_to_pixel_values(nbr.ra, nbr.dec))
    sl = slice(H - 1, H + 2)
    R = (dF / O)[sl, sl]; Re = (edF / O)[sl, sl]
    Pt = (M.img(xt, yt) / O)[sl, sl]; Pn = (M.img(xn, yn) / O)[sl, sl]
    # scale predictions to the observed (weighted) ratio level
    w = 1 / Re ** 2
    at = np.sum(R * Pt * w) / np.sum(Pt * Pt * w); an = np.sum(R * Pn * w) / np.sum(Pn * Pn * w)
    chi_t = float(np.sum((R - at * Pt) ** 2 * w)); chi_n = float(np.sum((R - an * Pn) ** 2 * w))
    rows.append(dict(sector=S["s"], chi2_ratio_target=chi_t, chi2_ratio_neighbour=chi_n, dof=8))
    vmax = np.max(R) * 1e3
    ext = [-1.5, 1.5, -1.5, 1.5]
    for i, (img, lab) in enumerate(((R, f"S{S['s']} observed dF/O"), (at * Pt, f"if on target (chi2={chi_t:.0f}/8)"), (an * Pn, f"if on neighbour (chi2={chi_n:.0f}/8)"))):
        ax = axes[i, j]
        ax.imshow(img * 1e3, origin="lower", cmap="magma", vmin=0, vmax=vmax, extent=ext)
        ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.5, 1.5)
        for (yy, xx), v in np.ndenumerate(img * 1e3):
            ax.text(xx - 1, yy - 1, f"{v:.1f}", ha="center", va="center", fontsize=7, color="w" if v < 0.6 * vmax else "k")
        ax.plot(xt - S["ix"], yt - S["iy"], "c+", ms=12, mew=2)
        ax.plot(xn - S["ix"], yn - S["iy"], "gx", ms=10, mew=2)
        # North / East arrows
        e = np.array(pl.to_EN(S["w"], xt + 1, yt, xt, yt)); n = np.array(pl.to_EN(S["w"], xt, yt + 1, xt, yt))
        J = np.array([[e[0], n[0]], [e[1], n[1]]]) / 21.0  # (E,N) per pixel
        Ji = np.linalg.inv(J)
        for vec, c, l in ((np.array([0, 1]), "w", "N"), (np.array([1, 0]), "y", "E")):
            d = Ji @ vec
            d = d / np.hypot(*d) * 0.45
            x0_, y0_ = -0.95, -0.95
            ax.annotate("", xy=(x0_ + d[0], y0_ + d[1]), xytext=(x0_, y0_), arrowprops=dict(arrowstyle="->", color=c, lw=1.5))
            ax.text(x0_ + d[0] * 1.35, y0_ + d[1] * 1.35, l, color=c, fontsize=8, ha="center", va="center", fontweight="bold")
        ax.set_title(lab, fontsize=8); ax.set_xticks([]); ax.set_yticks([])
fig.suptitle("Fractional flux drop per pixel (ppt) in the central 3x3 TESScut pixels; cyan + = TIC 87055304, green x = 9\" neighbour TIC 87055315", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "localization_summary.png"), dpi=120)
json.dump(rows, open(os.path.join(HERE, "localization_ratio_maps.json"), "w"), indent=1)
print(rows)
