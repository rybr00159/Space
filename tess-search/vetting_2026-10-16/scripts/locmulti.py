"""Localise with component images (phase-free Fourier or eclipse template) in each sector, plus a joint sky chi2 map."""
import sys, json, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from loc5 import Cut, component_images, fit_psf_multi, chi2_map, direct_offset, HERE
tic = int(sys.argv[1]); sectors = [int(s) for s in sys.argv[2].split(",")]
cfg_all = json.load(open(f"{HERE}/loccfg.json")); cfg0 = cfg_all[str(tic)]
R = cfg0.get("map_radius", 90); step = cfg0.get("map_step", 3)
ge, gn = np.meshgrid(np.arange(-R, R + step, step), np.arange(-R, R + step, step))
grid = np.column_stack([ge.ravel(), gn.ravel()])
tot = np.zeros(len(grid)); per = []; star_tot = None
for s in sectors:
    cfg = dict(cfg0); cfg.update(cfg0.get("per_sector", {}).get(str(s), {}))
    C = Cut(tic, s, cfg); psf = C.prf()
    comps, O = component_images(C, cfg)
    woff = direct_offset(C, psf, O)
    snr = max(float(np.nanmax(np.abs(A / S))) for A, S in comps)
    # free fit started at the strongest pixel of the summed component S/N image
    q = sum((A / S) ** 2 for A, S in comps); yy, xx = np.unravel_index(np.nanargmax(q), q.shape)
    F = fit_psf_multi(comps, psf, float(xx), float(yy))
    pos = np.array([F["x"], F["y"]]) - woff; dEN = C.J @ (pos - np.array([C.xt, C.yt])); cov = C.J @ F["cov"] @ C.J.T
    m = chi2_map(comps, psf, C, woff, grid); tot += m - m.min()
    g = C.g[np.hypot(C.g.dE, C.g.dN) < R].copy()
    sm = chi2_map(comps, psf, C, woff, g[["dE", "dN"]].values)
    g[f"c2_s{s}"] = sm - m.min()
    star_tot = g[["source_id", "dE", "dN", "T", f"c2_s{s}"]] if star_tot is None else star_tot.merge(g[["source_id", f"c2_s{s}"]], on="source_id", how="outer")
    per.append(dict(sector=s, snr=snr, dE=float(dEN[0]), dN=float(dEN[1]), eE=float(np.sqrt(cov[0, 0])), eN=float(np.sqrt(cov[1, 1]))))
    print(f"S{s}: peak S/N {snr:.1f}; free fit offset ({dEN[0]:+.1f}, {dEN[1]:+.1f})\" ± ({np.sqrt(cov[0,0]):.1f}, {np.sqrt(cov[1,1]):.1f}); map minimum at {grid[np.argmin(m)]}", flush=True)
i = int(np.argmin(tot))
cols = [c for c in star_tot.columns if c.startswith("c2_s")]
star_tot["c2_sum"] = star_tot[cols].sum(1)
star_tot["sep"] = np.hypot(star_tot.dE, star_tot.dN)
star_tot = star_tot.sort_values("c2_sum")
tgt = star_tot[star_tot.sep < 1.0]
print(f"joint map minimum at offset {grid[i]} arcsec; target delta-chi2 {float(tgt.c2_sum.iloc[0]) if len(tgt) else float('nan'):.1f}")
print(star_tot.head(10)[["source_id", "sep", "T", "c2_sum"] + cols].round(1).to_string(index=False))
star_tot.to_csv(f"{HERE}/out/{tic}/loc/joint_stars.csv", index=False)
np.save(f"{HERE}/out/{tic}/loc/joint_map.npy", np.column_stack([grid, tot]))
json.dump(per, open(f"{HERE}/out/{tic}/loc/per_sector_multi.json", "w"), indent=1)
