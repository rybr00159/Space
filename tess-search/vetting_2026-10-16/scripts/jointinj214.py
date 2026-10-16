"""Does the joint chi2 map resolve sources near TIC 214311906? Inject the signal (other phase) at known stars, build the
same joint map as for the real signal, and report the map minimum and the delta-chi2 at each candidate star."""
import sys, json, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from loc5 import Cut, component_images, chi2_map, direct_offset, template
cfg = json.load(open("loccfg.json"))["214311906"]; cfg = dict(cfg, P=5.691471, t0=3840.9101)
stars = {"Mira": 5925373268528684288, "pair": 5925373272850782208, "M3V": 5925373272850784768, "bright41": 5925373204131295872}
R, step = 60, 3
ge, gn = np.meshgrid(np.arange(-R, R + step, step), np.arange(-R, R + step, step)); grid = np.column_stack([ge.ravel(), gn.ravel()])
secs = (39, 66, 93, 104)
cache = {}
for s in secs:
    C = Cut(214311906, s, cfg); psf = C.prf(); comps, O = component_images(C, cfg); woff = direct_offset(C, psf, O)
    amp = np.nansum(np.clip(comps[0][0], 0, None)[np.hypot(*np.meshgrid(np.arange(C.nx) - C.xt, np.arange(C.ny) - C.yt)) < 3])
    cache[s] = (C, psf, woff, amp)
results = {}
for inj in ("pair", "M3V", "bright41", None):
    tot = np.zeros(len(grid)); star_c2 = {k: 0.0 for k in stars}
    for s in secs:
        C, psf, woff, amp = cache[s]
        cfg2 = dict(cfg)
        fl0 = C.flux.copy()
        if inj is not None:
            cfg2["t0"] = cfg["t0"] + 0.37 * cfg["P"]
            r = C.g[C.g.source_id == stars[inj]].iloc[0]
            img = psf(r.x + woff[0], r.y + woff[1]); img = img / img.sum()
            C.flux = fl0 - amp * template(C.t, cfg2)[:, None, None] * img[None]
        comps, _ = component_images(C, cfg2); C.flux = fl0
        m = chi2_map(comps, psf, C, woff, grid); tot += m - m.min()
        for k, gid in stars.items():
            r = C.g[C.g.source_id == gid].iloc[0]
            star_c2[k] += float(chi2_map(comps, psf, C, woff, np.array([[r.dE, r.dN]]))[0] - m.min())
    i = int(np.argmin(tot)); best = min(star_c2.values())
    label = f"injected at {inj}" if inj else "real signal"
    results[label] = dict(map_min=grid[i].tolist(), star_dchi2={k: round(v - best, 1) for k, v in star_c2.items()})
    print(f"{label:20s}: joint map minimum at {grid[i]}; delta-chi2 at stars " + ", ".join(f"{k} {v - best:.1f}" for k, v in star_c2.items()), flush=True)
json.dump(results, open("out/214311906/joint_injection_test.json", "w"), indent=1)
