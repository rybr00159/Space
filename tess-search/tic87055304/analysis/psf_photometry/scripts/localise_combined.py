"""Combine the per-sector pixel-map localisation grids (aperture.py: chi2 of a single PRF at a free position fitted to the
per-pixel transit flux decrements, control-phase covariance) in sky coordinates.

Each sector's grid (pixel offsets from the target) is mapped to (dE, dN) arcsec with the cutout WCS (local Jacobian at the
target), interpolated onto a common sky grid and summed.  Reports the best position, its distance from the target and the
neighbour, and chi2 at the two stars' positions.
"""
import os, sys, json
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import *
from psfphot import load_scene

SUB = sys.argv[1] if len(sys.argv) > 1 else "aperture"


def main():
    E = np.arange(-18, 18.01, 0.2); N = np.arange(-18, 18.01, 0.2)
    EE, NN = np.meshgrid(E, N)
    tot = np.zeros_like(EE); per = {}
    nb_sky = None
    for s in SECTORS:
        g = np.load(os.path.join(ROOT, SUB, f"grid_s{s}.npz"))
        grid, gx = g["grid"], g["gx"]          # grid[b, a]: oy = gx[b], ox = gx[a]
        sec = load_sector(s); sc, js = load_scene(s)
        w = sec["wcs"]
        xt, yt = js["target_xy"]
        r0 = np.array(w.pixel_to_world_values(xt, yt)); h = 0.1
        rx = np.array(w.pixel_to_world_values(xt + h, yt)); ry = np.array(w.pixel_to_world_values(xt, yt + h))
        cosd = np.cos(np.radians(r0[1]))
        J = np.array([[(rx[0] - r0[0]) * cosd, (ry[0] - r0[0]) * cosd], [rx[1] - r0[1], ry[1] - r0[1]]]) * 3600 / h  # arcsec per px
        Ji = np.linalg.inv(J)
        # pixel offsets for each sky grid point
        pix = Ji @ np.vstack([EE.ravel(), NN.ravel()])
        f = RegularGridInterpolator((gx, gx), grid, bounds_error=False, fill_value=np.nan)
        val = f(np.vstack([pix[1], pix[0]]).T).reshape(EE.shape)
        val = val - np.nanmin(val)
        per[s] = val
        tot += np.nan_to_num(val, nan=1e3)
        nbp = np.array(js["nb_xy"]) - np.array(js["target_xy"])
        nb_sky = J @ nbp
    k = np.unravel_index(np.argmin(tot), tot.shape)
    best = (EE[k], NN[k])
    interp = RegularGridInterpolator((N, E), tot)
    chi_t = float(interp([0.0, 0.0])[0] - tot.min()); chi_n = float(interp([nb_sky[1], nb_sky[0]])[0] - tot.min())
    m1 = tot - tot.min() < 2.30
    res = dict(best_offset_arcsec_EN=[float(best[0]), float(best[1])], neighbour_offset_arcsec_EN=nb_sky.tolist(),
               dist_best_to_nb_arcsec=float(np.hypot(best[0] - nb_sky[0], best[1] - nb_sky[1])),
               dist_best_to_target_arcsec=float(np.hypot(*best)),
               dchi2_at_target=chi_t, dchi2_at_neighbour=chi_n,
               one_sigma_E_range=[float(EE[m1].min()), float(EE[m1].max())], one_sigma_N_range=[float(NN[m1].min()), float(NN[m1].max())])
    print(json.dumps(res, indent=1))
    json.dump(res, open(os.path.join(ROOT, SUB, "localisation_combined.json"), "w"), indent=1)
    fig, ax = plt.subplots(figsize=(6, 5.5))
    d = tot - tot.min()
    cs = ax.contourf(E, N, np.clip(d, 0, 100), levels=[0, 2.3, 6.2, 11.8, 25, 50, 100], cmap="Greys_r")
    plt.colorbar(cs, ax=ax, label="delta chi2 (sum of 5 sectors)")
    ax.plot(0, 0, "+", color="#2a6fb0", ms=16, mew=2.5, label="target (Gaia DR3)")
    ax.plot(nb_sky[0], nb_sky[1], "x", color="#c0392b", ms=12, mew=2.5, label="neighbour TIC 87055315")
    ax.plot(*best, "o", mfc="none", mec="orange", ms=10, mew=2, label="best-fit source")
    ax.set_xlim(18, -18); ax.set_xlabel("dE (arcsec, east left)"); ax.set_ylabel("dN (arcsec)"); ax.set_aspect("equal")
    ax.set_title("Source of the transit flux decrement (pixel maps, 5 sectors)\ncontours: delta chi2 = 2.3, 6.2, 11.8, 25, 50", fontsize=9)
    ax.legend(fontsize=8, loc="lower left")
    plt.tight_layout(); plt.savefig(os.path.join(ROOT, f"localisation_combined{'_regs' if 'regs' in SUB else ''}.png"), dpi=100); plt.close()


if __name__ == "__main__":
    main()
