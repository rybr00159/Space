"""Independent aperture photometry of TIC 90086266 from TESScut cubes; odd/even per sector and aperture."""
import sys, json, warnings
sys.path.insert(0, "/home/user/work/src"); warnings.filterwarnings("ignore")
sys.argv = [sys.argv[0], "90086266"]
import numpy as np
from scipy.optimize import least_squares
from wotan import flatten
from lcan import trap, in_mask
from loc5 import Cut, circle_mask
cfg = json.load(open("loccfg.json"))["90086266"]
P, T0, T14, FL = 2.339469, 4180.4740, 0.128, 0.12
out = {}
for sec in (13, 104):
    C = Cut(90086266, sec, cfg); g = C.good; t = C.t[g]; cube = C.flux[g]
    med = np.median(cube, 0); bkgpix = med < np.percentile(med, 30)
    bkg = np.median(cube[:, bkgpix], 1)
    for rad in (1.0, 1.5, 2.5):
        ap = circle_mask((C.ny, C.nx), C.xt, C.yt, rad)
        f = cube[:, ap].sum(1) - bkg * ap.sum()
        if C.exptime < 500 / 86400:
            k = np.floor((t - t[0]) / (600 / 86400)).astype(int); _, inv, cnt = np.unique(k, return_inverse=True, return_counts=True)
            tb = np.bincount(inv, t) / cnt; fb = np.bincount(inv, f) / cnt; tb, fb = tb[cnt >= 2], fb[cnt >= 2]
        else:
            tb, fb = t, f
        fb = fb / np.median(fb)
        m = ~in_mask(tb, P, T0, T14, 0.8)
        tr = flatten(tb[m], fb[m], method="biweight", window_length=1.0, return_trend=True, break_tolerance=0.4)[1]
        ok = np.isfinite(tr); fd = fb / np.interp(tb, tb[m][ok], tr[ok])
        s = 1.4826 * np.median(np.abs(fd[m] - np.median(fd[m])))
        ev = []
        for kk in np.unique(np.round((tb - T0) / P)):
            tc = T0 + kk * P; w = np.abs(tb - tc) < 0.25
            if np.sum(np.abs(tb[w] - tc) < T14 / 2) < 6 or w.sum() < 15: continue
            fit = least_squares(lambda x: (fd[w] - trap(tb[w], tc + x[0], x[1], T14, T14 * FL)) / s, [0, 0.03], bounds=([-0.02, 0], [0.02, 0.1]))
            J = fit.jac; cov = np.linalg.pinv(J.T @ J) * max(1, 2 * fit.cost / (w.sum() - 2))
            ev.append((int(kk), fit.x[1], np.sqrt(cov[1, 1])))
        ev = np.array(ev)
        def wm(a):
            w = 1 / a[:, 2] ** 2; mu = np.sum(w * a[:, 1]) / w.sum(); chi = np.sum(((a[:, 1] - mu) / a[:, 2]) ** 2) / max(1, len(a) - 1)
            return mu, np.sqrt(max(1, chi) / w.sum())
        e, o = wm(ev[ev[:, 0] % 2 == 0]), wm(ev[ev[:, 0] % 2 == 1])
        d = (e[0] - o[0]) / np.hypot(e[1], o[1])
        out[f"S{sec}_r{rad}"] = dict(even=e, odd=o, sigma=d, events=ev.tolist())
        print(f"S{sec} aperture r={rad} px: depth even {e[0]*1e6:.0f}±{e[1]*1e6:.0f} odd {o[0]*1e6:.0f}±{o[1]*1e6:.0f} -> {d:+.1f} sigma ({(e[0]-o[0])/o[0]*100:+.1f}%)  "
              f"[{len(ev)} events]")
json.dump(out, open("out/90086266/oddeven_tesscut.json", "w"), indent=1, default=float)
