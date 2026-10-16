"""Per-event times for the 6.5286-d dips of TIC 61141412 (shape fixed, time and depth free), then linear and
quadratic ephemerides."""
import sys, json, warnings
sys.path.insert(0, "/home/user/work/src"); warnings.filterwarnings("ignore")
sys.argv = [sys.argv[0], "61141412"]
import numpy as np
from scipy.optimize import least_squares
from wotan import flatten
import stage2
from lcan import read, bin10, trap
P0, T0, T14, FL = 6.528604, 3143.003, 0.042, 0.3
best = {13: ("SPOC", "PDCSAP_FLUX"), 27: ("TESS-SPOC", "PDCSAP_FLUX"), 67: ("TESS-SPOC", "PDCSAP_FLUX"), 94: ("QLP", "SAP_FLUX"), 104: ("QLP", "SAP_FLUX")}
p = stage2.products(61141412); ev = []
for sec, (prov, col) in best.items():
    r = p[(p.sector == sec) & (p.prov == prov)].sort_values("exptime").iloc[0]
    t, q, cols, bkg, aper, crowd = read(stage2.fetch(r.uri), prov)
    ok = (q == 0) & np.isfinite(t) & np.isfinite(cols[col]) & (cols[col] > 0)
    tt, ff = bin10(t[ok], cols[col][ok] / np.nanmedian(cols[col][ok]), r.exptime)
    fl = flatten(tt, ff, method="biweight", window_length=0.6, break_tolerance=0.4); m = np.isfinite(fl); tt, fl = tt[m], fl[m]
    s = 1.4826 * np.median(np.abs(fl - np.median(fl)))
    for k in np.unique(np.round((tt - T0) / P0)):
        tc = T0 + k * P0; w = np.abs(tt - tc) < 0.35
        if np.sum(np.abs(tt[w] - tc) < 0.1) < 4 or w.sum() < 20: continue
        # scan time offsets +-0.12 d, then refine
        grid = np.arange(-0.12, 0.1201, 0.004); c2 = []
        for dt in grid:
            y = trap(tt[w], tc + dt, 1.0, T14, T14 * FL)
            A = np.column_stack([np.ones(w.sum()), 1 - y]); c, *_ = np.linalg.lstsq(A, fl[w], rcond=None)
            c2.append(np.sum(((fl[w] - A @ c) / s) ** 2))
        c2 = np.array(c2); j = int(np.argmin(c2))
        fit = least_squares(lambda x: (fl[w] - x[1] - (-x[2]) * (1 - trap(tt[w], tc + x[0], 1.0, T14, T14 * FL))) / s,
                            [grid[j], 1.0, -0.0012], bounds=([-0.13, 0.99, -0.01], [0.13, 1.01, 0.01]))
        J = fit.jac; cov = np.linalg.pinv(J.T @ J) * max(1, 2 * fit.cost / (w.sum() - 3))
        d = fit.x[2]; e_d = np.sqrt(cov[2, 2]); e_t = np.sqrt(cov[0, 0])
        dchi = c2.max() - c2.min()
        ev.append(dict(sector=sec, k=int(k), tc=tc + fit.x[0], err=e_t, depth=d, depth_err=e_d, snr=d / e_d))
for e in ev:
    print(f"S{e['sector']:3d} k={e['k']:4d} t={e['tc']:.4f}±{e['err']*1440:4.1f} min  depth {e['depth']*1e6:5.0f}±{e['depth_err']*1e6:4.0f} (S/N {e['snr']:4.1f})")
good = [e for e in ev if e["snr"] > 4 and e["err"] < 0.02]
k = np.array([e["k"] for e in good]); t = np.array([e["tc"] for e in good]); w = 1 / np.array([e["err"] for e in good])
for deg in (1, 2):
    A = np.vander(k, deg + 1, increasing=True)
    c, *_ = np.linalg.lstsq(A * w[:, None], t * w, rcond=None); res = t - A @ c
    chi2 = np.sum((res * w) ** 2)
    print(f"degree {deg}: coeffs {np.round(c, 7).tolist()}  chi2 {chi2:.1f} for {len(t) - deg - 1} dof; O-C (min) by sector: " +
          ", ".join(f"S{s}: {np.mean(res[[g['sector'] == s for g in good]]) * 1440:+.1f}" for s in sorted(set(e['sector'] for e in good))))
json.dump(ev, open("out/61141412/etv_events.json", "w"), indent=1, default=float)
