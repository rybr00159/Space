"""Are QLP S103/S104 times early? Fit the same TIC 90086266 eclipses in QLP SAP and in TESScut aperture photometry."""
import sys, json, warnings
sys.path.insert(0, "/home/user/work/src"); warnings.filterwarnings("ignore")
sys.argv = [sys.argv[0], "90086266"]
import numpy as np
from scipy.optimize import least_squares
from wotan import flatten
import stage2
from lcan import read, bin10, trap, in_mask
from loc5 import Cut, circle_mask
P, T0, T14, FL = 2.339469, 4180.4740, 0.128, 0.12
cfg = json.load(open("loccfg.json"))["90086266"]

def times(t, f, ex):
    m = ~in_mask(t, P, T0, T14, 0.8)
    tr = flatten(t[m], f[m], method="biweight", window_length=1.0, return_trend=True, break_tolerance=0.4)[1]
    k = np.isfinite(tr); f = f / np.interp(t, t[m][k], tr[k]); s = 1.4826 * np.median(np.abs(f[m] - np.median(f[m])))
    out = {}
    for e in np.unique(np.round((t - T0) / P)):
        tc = T0 + e * P; w = np.abs(t - tc) < 0.25
        if np.sum(np.abs(t[w] - tc) < T14 / 2) < 5 or w.sum() < 15: continue
        r = least_squares(lambda x: (f[w] - trap(t[w], tc + x[0], x[1], T14, T14 * FL)) / s, [0, 0.03], bounds=([-0.02, 0], [0.02, 0.1]))
        J = r.jac; cov = np.linalg.pinv(J.T @ J) * max(1, 2 * r.cost / (w.sum() - 2))
        out[int(e)] = (tc + r.x[0], np.sqrt(cov[0, 0]))
    return out

p = stage2.products(90086266)
for sec in (13, 104):
    r = p[(p.sector == sec) & (p.prov == "QLP")].iloc[0]
    t, q, cols, bkg, aper, crowd = read(stage2.fetch(r.uri), "QLP")
    ok = (q == 0) & np.isfinite(t) & np.isfinite(cols["SAP_FLUX"]) & (cols["SAP_FLUX"] > 0)
    tq, fq = bin10(t[ok], cols["SAP_FLUX"][ok] / np.nanmedian(cols["SAP_FLUX"][ok]), r.exptime)
    A = times(tq, fq, r.exptime)
    C = Cut(90086266, sec, cfg); g = C.good; tt = C.t[g]; cube = C.flux[g]
    ap = circle_mask((C.ny, C.nx), C.xt, C.yt, 1.5); med = np.median(cube, 0); bk = np.median(cube[:, med < np.percentile(med, 30)], 1)
    fc = cube[:, ap].sum(1) - bk * ap.sum()
    tb, fb = bin10(tt, fc / np.median(fc), C.exptime * 86400)
    B = times(tb, fb, C.exptime * 86400)
    common = sorted(set(A) & set(B))
    d = np.array([(A[e][0] - B[e][0]) * 1440 for e in common]); s = np.array([np.hypot(A[e][1], B[e][1]) * 1440 for e in common])
    w = 1 / s ** 2
    print(f"S{sec}: {len(common)} eclipses; QLP minus TESScut eclipse time = {np.sum(w*d)/w.sum():+.2f} ± {1/np.sqrt(w.sum()):.2f} min (per-eclipse {np.round(d,1).tolist()})")
