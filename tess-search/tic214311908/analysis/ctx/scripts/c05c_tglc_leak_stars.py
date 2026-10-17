"""Follow-up of c05: TGLC S39 leak model evaluated at every Gaia star within 30" of the M3V, plus a noise simulation
(signal injected at the M3V position in the leak model + real noise vectors from 400 random phases) to calibrate the
free-position fit. Also robustness variants (drop the Mira's own TGLC light curve; drop the reference stars).
Output: out/tglc_leak_stars.json, fig/tglc_leak_map_s39.png
"""
import os, sys, json, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import c05_tglc_leak as L

sector = 39; out = {}
pos, prf, J, w = L.positions(sector); sxy = L.tglc_star_xy(sector); Ji = np.linalg.inv(J)
xm, ym = pos["Mira"]
g = pd.read_csv(L.GAIA); ep = 2021.42; dt = ep - 2016.0
g["ra_ep"] = g.ra + np.nan_to_num(g.pmra) * dt / 3.6e6 / np.cos(np.radians(g.dec)); g["dec_ep"] = g.dec + np.nan_to_num(g.pmdec) * dt / 3.6e6
xy = w.all_world2pix(np.column_stack([g.ra_ep, g.dec_ep]), 0); g["x"], g["y"] = xy[:, 0], xy[:, 1]
EN = (J @ np.vstack([g.x - xm, g.y - ym])).T; g["E"], g["N"] = EN[:, 0], EN[:, 1]
m3 = g[g.source_id == L.STARS["M3V"]].iloc[0]; g["dM3V"] = np.hypot(g.E - m3.E, g.N - m3.N)
for col in ("psf_flux", "aperture_flux"):
    obs, rnd, lvl = L.observed(sector, col)
    for variant, names in (("all6", list(L.STARS)), ("no_Mira", ["M3V", "pairA", "pairB", "ref1", "ref2"]),
                           ("no_refs", ["M3V", "pairA", "pairB", "Mira"])):
        ov = np.array([obs[k] for k in names]); R = np.array([rnd[k] for k in names]); cov = np.cov(R)
        def model(x, y):
            cp, ca = L.leak_coeffs(pos, prf, sxy, (x, y)); c = cp if col == "psf_flux" else ca
            return np.array([c[k] for k in names])
        rows = []
        for _, s in g[g.dM3V < 30].sort_values("dM3V").iterrows():
            dF, e, chi2 = L.gls(ov, cov, model(s.x, s.y))
            rows.append(dict(source_id=int(s.source_id), G=float(s.G), T=float(s["T"]), E=round(float(s.E), 1), N=round(float(s.N), 1),
                             dM3V=round(float(s.dM3V), 1), dF=round(dF, 2), dF_err=round(e, 2), chi2=round(chi2, 2),
                             need_depth=round(dF / (15000 * 10 ** (-0.4 * (s["T"] - 10))), 3)))
        df = pd.DataFrame(rows); chimin_star = df.chi2.min()
        # free-position grid
        grid = np.arange(-45, 45.1, 1.0); chi = np.full((len(grid), len(grid)), np.nan)
        for i, dn in enumerate(grid):
            for j, de in enumerate(grid):
                dx, dy = Ji @ np.array([de, dn]); chi[i, j] = L.gls(ov, cov, model(xm + dx, ym + dy))[2]
        i0, j0 = np.unravel_index(np.nanargmin(chi), chi.shape)
        res = dict(stars=df.head(15).to_dict("records"), grid_best_EN=[float(grid[j0]), float(grid[i0])], grid_chi2_min=float(np.nanmin(chi)),
                   dof=len(names) - 3)
        # noise simulation: inject at the M3V with the fitted dF, real noise vectors (random phases, mean removed)
        if variant == "all6":
            mv = model(m3.x, m3.y); dF0 = L.gls(ov, cov, mv)[0]; Rc = R - R.mean(axis=1, keepdims=True)
            sub = np.arange(-30, 0.1, 1.5); best = []; chi_m3 = []
            # coarse grid near the stars of interest for speed
            gx = np.arange(-40, 10.1, 1.5); gy = np.arange(-20, 25.1, 1.5)
            Mgrid = {(de, dn): model(xm + (Ji @ np.array([de, dn]))[0], ym + (Ji @ np.array([de, dn]))[1]) for dn in gy for de in gx}
            rng = np.random.default_rng(5)
            for k in rng.choice(R.shape[1], 300, replace=False):
                sim = dF0 * mv + Rc[:, k]
                c2 = {key: L.gls(sim, cov, mm)[2] for key, mm in Mgrid.items()}
                kb = min(c2, key=c2.get); best.append(kb); chi_m3.append(L.gls(sim, cov, mv)[2] - c2[kb])
            best = np.array(best); chi_m3 = np.array(chi_m3)
            d_off = np.hypot(best[:, 0] - m3.E, best[:, 1] - m3.N)
            res["sim"] = dict(n=len(best), dF_injected=dF0, median_offset_arcsec=float(np.median(d_off)),
                              p68_offset=float(np.percentile(d_off, 68)), p95_offset=float(np.percentile(d_off, 95)),
                              frac_dchi2_at_truth_ge_observed=float(np.mean(chi_m3 >= (L.gls(ov, cov, mv)[2] - np.nanmin(chi)))),
                              dchi2_truth_p50_p95=[float(np.percentile(chi_m3, 50)), float(np.percentile(chi_m3, 95))],
                              best_EN_sample=best[:20].tolist())
            if col == "psf_flux":
                fig, ax = plt.subplots(figsize=(6.5, 6))
                im = ax.imshow(chi - np.nanmin(chi), origin="lower", extent=[grid[0], grid[-1], grid[0], grid[-1]], vmin=0, vmax=60, cmap="viridis_r")
                ax.contour(grid, grid, chi - np.nanmin(chi), levels=[2.3, 6.2, 11.8], colors="w", linewidths=0.8)
                near = g[(np.abs(g.E) < 45) & (np.abs(g.N) < 45)]
                ax.scatter(near.E, near.N, s=np.clip(30 - 2 * (near["T"] - 11) ** 1.5, 2, 60), facecolors="none", edgecolors="r", lw=0.7)
                for k in ("M3V", "pairA", "Mira"):
                    s = g[g.source_id == L.STARS[k]].iloc[0]; ax.annotate(k, (s.E, s.N), color="k", fontsize=8)
                ax.scatter(best[:, 0], best[:, 1], s=4, c="orange", alpha=0.4, label="sim: injected on M3V")
                ax.plot(grid[j0], grid[i0], "w+", ms=12); ax.set_xlabel("E offset from Mira (arcsec)"); ax.set_ylabel("N offset (arcsec)")
                ax.set_title("TGLC S39 psf_flux leak model: chi2 - min"); plt.colorbar(im); ax.legend(loc="lower left", fontsize=7)
                plt.tight_layout(); plt.savefig("fig/tglc_leak_map_s39.png", dpi=90)
        out[f"{col}_{variant}"] = res
        print(col, variant, "best grid", res["grid_best_EN"], "chi2min", round(res["grid_chi2_min"], 2), "dof", res["dof"])
        print(df.head(8)[["source_id", "G", "E", "N", "dM3V", "dF", "chi2", "need_depth"]].to_string(index=False))
        if "sim" in res: print("   sim:", {k: v for k, v in res["sim"].items() if k != "best_EN_sample"})
json.dump(out, open("out/tglc_leak_stars.json", "w"), indent=1, default=float)
