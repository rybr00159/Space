"""TIC 214311906 field, redone with barycentred TESScut times:
(1) per-event times from small-aperture photometry at the M3V and Mira positions -> linear ephemeris;
(2) absolute dip (e/s) in r = 0.75 px apertures at the Mira, the pair (TIC 214311910), the M3V (TIC 214311908) and the
    mirror point east of the Mira, per sector."""
import sys, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
from scipy.optimize import least_squares
from wotan import flatten
from loc5 import Cut, circle_mask
cfg = json.load(open("loccfg.json"))["214311906"]
IDS = {"Mira": 5925373268528684288, "pair": 5925373272850782208, "M3V": 5925373272850784768}
P0, T00, T14, FL = 5.69149, 3840.913, 0.067, 0.78
def trap(t, tc, T14=T14, fl=FL):
    a = np.abs(t - tc); T23 = T14 * fl
    return np.where(a <= T23 / 2, 1.0, np.where(a >= T14 / 2, 0.0, (T14 / 2 - a) / ((T14 - T23) / 2)))
apflux = {}; events = []
for sec in (12, 39, 66, 93, 103, 104):
    C = Cut(214311906, sec, cfg); g = C.good; t = C.t[g]; cube = C.flux[g]
    med = np.median(cube, 0); bk = np.median(cube[:, med < np.percentile(med, 25)], 1)
    pos = {k: (float(C.g[C.g.source_id == v].x.iloc[0]), float(C.g[C.g.source_id == v].y.iloc[0])) for k, v in IDS.items()}
    mx, my = 2 * pos["Mira"][0] - pos["M3V"][0], 2 * pos["Mira"][1] - pos["M3V"][1]; pos["mirror"] = (mx, my)
    for name, (x, y) in pos.items():
        ap = circle_mask((C.ny, C.nx), x, y, 0.75)
        f = cube[:, ap].sum(1) - bk * ap.sum()
        apflux[(sec, name)] = (t, f, float(np.median(f)))
    # per-event times on the M3V aperture (r = 1 px) and the Mira aperture
    for name in ("M3V", "Mira"):
        x, y = pos[name]; ap = circle_mask((C.ny, C.nx), x, y, 1.0); f = cube[:, ap].sum(1) - bk * ap.sum(); f = f / np.median(f)
        ph = np.abs(((t - T00) / P0 + 0.5) % 1 - 0.5) * P0; m = ph > 0.12
        tr = flatten(t[m], f[m], method="biweight", window_length=1.0, return_trend=True, break_tolerance=0.4)[1]
        k_ = np.isfinite(tr); fd = f / np.interp(t, t[m][k_], tr[k_]); s = 1.4826 * np.median(np.abs(fd[m] - np.median(fd[m])))
        for e in np.unique(np.round((t - T00) / P0)):
            tc = T00 + e * P0; w = np.abs(t - tc) < 0.3
            if np.sum(np.abs(t[w] - tc) < T14) < 4 or w.sum() < 15: continue
            grid = np.arange(-0.15, 0.1501, 0.002); c2 = []
            for dt in grid:
                A = np.column_stack([np.ones(w.sum()), trap(t[w], tc + dt)]); c, *_ = np.linalg.lstsq(A, fd[w], rcond=None)
                c2.append(np.sum(((fd[w] - A @ c) / s) ** 2) if c[1] < 0 else np.inf)
            j = int(np.argmin(c2))
            if not np.isfinite(c2[j]): continue
            r = least_squares(lambda x_: (fd[w] - x_[1] + x_[2] * trap(t[w], tc + x_[0])) / s, [grid[j], 1.0, 0.002], bounds=([-0.16, 0.9, 0], [0.16, 1.1, 0.1]))
            J = r.jac; cov = np.linalg.pinv(J.T @ J) * max(1, 2 * r.cost / (w.sum() - 3))
            snr = r.x[2] / np.sqrt(cov[2, 2])
            if snr > 3: events.append(dict(sector=sec, ap=name, k=int(e), tc=float(tc + r.x[0]), err=float(np.sqrt(cov[0, 0])), depth=float(r.x[2]), snr=float(snr)))
# ephemeris from the M3V-aperture events
out = {}
for name in ("M3V", "Mira"):
    ev = [e for e in events if e["ap"] == name]
    k = np.array([e["k"] for e in ev]); tt = np.array([e["tc"] for e in ev]); s = np.array([e["err"] for e in ev])
    A = np.column_stack([np.ones(len(k)), k]); W = 1 / s; c, *_ = np.linalg.lstsq(A * W[:, None], tt * W, rcond=None)
    cov = np.linalg.pinv((A * W[:, None]).T @ (A * W[:, None])); res = tt - A @ c; chi2 = float(np.sum((res * W) ** 2)); dof = len(tt) - 2
    cov *= max(1, chi2 / dof)
    out[name] = dict(P=float(c[1]), P_err=float(np.sqrt(cov[1, 1])), T0=float(c[0]), T0_err=float(np.sqrt(cov[0, 0])), cov=cov.tolist(), chi2=chi2, dof=dof, n=len(tt))
    print(f"{name} aperture: {len(tt)} events; P = {c[1]:.6f} ± {np.sqrt(cov[1,1]):.6f} d, T0 = {c[0]:.4f} ± {np.sqrt(cov[0,0]):.4f}; chi2 {chi2:.1f}/{dof}; "
          f"by sector {dict(zip(*np.unique([e['sector'] for e in ev], return_counts=True)))}; mean S/N {np.mean([e['snr'] for e in ev]):.1f}")
P, T0 = out["M3V"]["P"], out["M3V"]["T0"]
# absolute dips with the new ephemeris
print("absolute dip (e/s) in r=0.75 px apertures [aperture median flux]:")
tot = {}
for sec in (12, 39, 66, 93, 103, 104):
    row = []
    for name in ("M3V", "pair", "Mira", "mirror"):
        t, f, lvl = apflux[(sec, name)]
        ph = np.abs(((t - T0) / P + 0.5) % 1 - 0.5) * P; m = ph > 0.12
        tr = flatten(t[m], f[m], method="biweight", window_length=1.0, return_trend=True, break_tolerance=0.4)[1]
        k_ = np.isfinite(tr); r_ = f - np.interp(t, t[m][k_], tr[k_])
        y = trap(t, T0 + np.round((t - T0) / P) * P); A = np.column_stack([np.ones(len(t)), y])
        c, *_ = np.linalg.lstsq(A, r_, rcond=None); s = 1.4826 * np.median(np.abs(r_[m] - np.median(r_[m]))); e = s / np.sqrt(np.sum(y ** 2))
        row.append(f"{name} {-c[1]:5.1f}±{e:3.1f} [{lvl:6.0f}]"); tot.setdefault(name, []).append((-c[1], e))
    print(f"  S{sec:3d}: " + " | ".join(row))
for name, v in tot.items():
    v = np.array(v); w = 1 / v[:, 1] ** 2; print(f"  weighted mean {name}: {np.sum(w*v[:,0])/w.sum():.2f} ± {1/np.sqrt(w.sum()):.2f} e/s")
json.dump(dict(ephemeris=out, events=events, abs_dip={k: v for k, v in tot.items()}), open("out/214311906/host_v2.json", "w"), indent=1, default=float)
