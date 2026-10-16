"""Alternate-eclipse test for TIC 90086266 per sector and product: per-event trapezoid fits (shape fixed per sector,
depth and time free), then even vs odd epochs on the linked ephemeris P = 2.339469 d."""
import sys, json, warnings
sys.path.insert(0, "/home/user/work/src"); warnings.filterwarnings("ignore")
sys.argv = [sys.argv[0], "90086266"]
import numpy as np
from scipy.optimize import least_squares
from wotan import flatten
import stage2
from lcan import read, bin10, trap, in_mask
P, T0, T14, FL = 2.339469, 4180.4740, 0.128, 0.12
p = stage2.products(90086266); res = {}
for _, r in p.iterrows():
    t, q, cols, bkg, aper, crowd = read(stage2.fetch(r.uri), r.prov)
    for c in ("PDCSAP_FLUX", "SAP_FLUX"):
        if c not in cols: continue
        ok = (q == 0) & np.isfinite(t) & np.isfinite(cols[c]) & (cols[c] > 0)
        tt, ff = bin10(t[ok], cols[c][ok] / np.nanmedian(cols[c][ok]), r.exptime)
        m = ~in_mask(tt, P, T0, T14, 0.8)
        tr = flatten(tt[m], ff[m], method="biweight", window_length=1.0, return_trend=True, break_tolerance=0.4)[1]
        g = np.isfinite(tr); f = ff / np.interp(tt, tt[m][g], tr[g])
        s = 1.4826 * np.median(np.abs(f[m] - np.median(f[m])))
        ev = []
        for k in np.unique(np.round((tt - T0) / P)):
            tc = T0 + k * P; w = np.abs(tt - tc) < 0.25
            if np.sum(np.abs(tt[w] - tc) < T14 / 2) < 0.5 * T14 / (max(600, r.exptime) / 86400) or w.sum() < 15: continue
            fit = least_squares(lambda x: (f[w] - trap(tt[w], tc + x[0], x[1], T14 * x[2], T14 * x[2] * FL)) / s,
                                [0, 0.03, 1.0], bounds=([-0.02, 0, 0.7], [0.02, 0.1, 1.3]))
            J = fit.jac; cov = np.linalg.pinv(J.T @ J) * max(1, 2 * fit.cost / (w.sum() - 3))
            ev.append((int(k), fit.x[1], np.sqrt(cov[1, 1]), fit.x[2] * T14 * 24, fit.x[0] * 1440))
        ev = np.array(ev)
        if len(ev) < 4: continue
        even = ev[ev[:, 0] % 2 == 0]; odd = ev[ev[:, 0] % 2 == 1]
        def wm(a):
            w = 1 / a[:, 2] ** 2; mu = np.sum(w * a[:, 1]) / w.sum()
            # error inflated by the scatter of the events about their mean
            chi = np.sum(((a[:, 1] - mu) / a[:, 2]) ** 2) / max(1, len(a) - 1)
            return mu, np.sqrt(1 / w.sum() * max(1, chi)), np.mean(a[:, 3]), np.std(a[:, 3]) / np.sqrt(len(a))
        me, mo = wm(even), wm(odd)
        key = f"S{int(r.sector)} {r.prov} {c}"
        res[key] = dict(even=me, odd=mo, diff_sigma=(me[0] - mo[0]) / np.hypot(me[1], mo[1]), n_even=len(even), n_odd=len(odd),
                        events=ev.tolist())
        print(f"{key:30s} even {me[0]*1e6:7.0f}±{me[1]*1e6:4.0f} ({len(even)})  odd {mo[0]*1e6:7.0f}±{mo[1]*1e6:4.0f} ({len(odd)})  "
              f"diff {(me[0]-mo[0])*1e6:+6.0f} ppm = {res[key]['diff_sigma']:+.1f} sigma ({(me[0]-mo[0])/mo[0]*100:+.1f}%);  T14 even {me[2]:.2f}±{me[3]:.2f} h odd {mo[2]:.2f}±{mo[3]:.2f} h")
json.dump(res, open("out/90086266/oddeven.json", "w"), indent=1, default=float)
