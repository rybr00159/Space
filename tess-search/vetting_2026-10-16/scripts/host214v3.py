"""TIC 214311906 field with the corrected ephemeris: (1) absolute dips in r = 0.75 px apertures at the M3V, pair, Mira
and the mirror point; (2) joint chi2 map (S39, S66, S93, S103, S104) calibrated with injections at each star."""
import sys, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
from wotan import flatten
from loc5 import Cut, circle_mask, component_images, chi2_map, direct_offset, template
E = json.load(open("out/214311906/ephemeris_v2.json")); P, T0 = E["P"], E["T0"]
cfg = json.load(open("loccfg.json"))["214311906"]; cfg = dict(cfg, P=P, t0=T0)
IDS = {"M3V": 5925373272850784768, "pair": 5925373272850782208, "Mira": 5925373268528684288, "bright41": 5925373204131295872}
SECS = (39, 66, 93, 103, 104)
def trap(t):
    tc = T0 + np.round((t - T0) / P) * P; a = np.abs(t - tc); T14, T23 = 0.067, 0.067 * 0.78
    return np.where(a <= T23 / 2, 1.0, np.where(a >= T14 / 2, 0.0, (T14 / 2 - a) / ((T14 - T23) / 2)))
tot = {}; cache = {}
for sec in SECS + (12,):
    C = Cut(214311906, sec, cfg); g = C.good; t = C.t[g]; cube = C.flux[g]
    med = np.median(cube, 0); bk = np.median(cube[:, med < np.percentile(med, 25)], 1)
    pos = {k: (float(C.g[C.g.source_id == v].x.iloc[0]), float(C.g[C.g.source_id == v].y.iloc[0])) for k, v in IDS.items() if k != "bright41"}
    pos["mirror"] = (2 * pos["Mira"][0] - pos["M3V"][0], 2 * pos["Mira"][1] - pos["M3V"][1])
    row = []
    for name, (x, y) in pos.items():
        ap = circle_mask((C.ny, C.nx), x, y, 0.75); f = cube[:, ap].sum(1) - bk * ap.sum()
        m = trap(t) == 0
        tr = flatten(t[m], f[m], method="biweight", window_length=1.0, return_trend=True, break_tolerance=0.4)[1]
        k_ = np.isfinite(tr); r_ = f - np.interp(t, t[m][k_], tr[k_]); y_ = trap(t)
        A = np.column_stack([np.ones(len(t)), y_]); c, *_ = np.linalg.lstsq(A, r_, rcond=None)
        s = 1.4826 * np.median(np.abs(r_[m] - np.median(r_[m]))); e = s / np.sqrt(np.sum(y_ ** 2))
        row.append(f"{name} {-c[1]:5.1f}±{e:3.1f} [{int(ap.sum())}px]"); tot.setdefault(name, []).append((-c[1], e))
    print(f"S{sec:3d}: " + " | ".join(row), flush=True)
for name, v in tot.items():
    v = np.array(v); w = 1 / v[:, 1] ** 2; print(f"  weighted mean {name}: {np.sum(w*v[:,0])/w.sum():.2f} ± {1/np.sqrt(w.sum()):.2f} e/s")
# joint map with injections
R, step = 60, 3
ge, gn = np.meshgrid(np.arange(-R, R + step, step), np.arange(-R, R + step, step)); grid = np.column_stack([ge.ravel(), gn.ravel()])
for s in SECS:
    C = Cut(214311906, s, cfg); psf = C.prf(); comps, O = component_images(C, cfg); woff = direct_offset(C, psf, O)
    cache[s] = (C, psf, woff)
res = {}
for inj in (None, "M3V", "pair", "Mira"):
    totm = np.zeros(len(grid)); sc = {k: 0.0 for k in IDS}
    for s in SECS:
        C, psf, woff = cache[s]; cfg2 = dict(cfg); fl0 = C.flux.copy()
        if inj is not None:
            cfg2["t0"] = T0 + 0.37 * P
            r = C.g[C.g.source_id == IDS[inj]].iloc[0]; img = psf(r.x + woff[0], r.y + woff[1]); img = img / img.sum()
            C.flux = fl0 - 16.5 * template(C.t, cfg2)[:, None, None] * img[None]     # 16.5 e/s total, as measured by the referee
        comps, _ = component_images(C, cfg2); C.flux = fl0
        m = chi2_map(comps, psf, C, woff, grid); totm += m - m.min()
        for k, gid in IDS.items():
            r = C.g[C.g.source_id == gid].iloc[0]
            sc[k] += float(chi2_map(comps, psf, C, woff, np.array([[r.dE, r.dN]]))[0] - m.min())
    best = min(sc.values()); i = int(np.argmin(totm))
    lab = "real signal" if inj is None else f"injected at {inj}"
    res[lab] = dict(map_min=grid[i].tolist(), dchi2={k: round(v - best, 1) for k, v in sc.items()})
    print(f"{lab:18s}: map minimum {grid[i]}; delta-chi2 " + ", ".join(f"{k} {v - best:.1f}" for k, v in sc.items()), flush=True)
json.dump(dict(abs_dip={k: v for k, v in tot.items()}, joint=res), open("out/214311906/host_v3.json", "w"), indent=1, default=float)
