"""Diagnostic: single-window noise of delta(p) in fake vs real windows (robust vs plain std)."""
import numpy as np
from rcommon import *
from r3_pixel_ratio import collect, PHASES
for s in SECTORS:
    c = Cube(s); yy, xx = np.mgrid[0:c.ny, 0:c.nx]; far = np.hypot(xx-c.xt, yy-c.yt) > 5
    tp = (int(round(c.xt)), int(round(c.yt))); a = (int(round(c.xt + c.u[0])), int(round(c.yt + c.u[1])))
    real = collect(c); fake = [w for ph in PHASES for w in collect(c, ph)]
    def vals(ws):
        out = []
        for w in ws:
            D = w['D'] - np.median(w['D'][far]); O = w['O'] - np.percentile(w['O'], 10)
            out.append((D[tp[1], tp[0]] / O[tp[1], tp[0]] * 1e3, D[a[1], a[0]] / O[a[1], a[0]] * 1e3, w['sD'][tp[1], tp[0]], np.median(w['D'][far])))
        return np.array(out)
    R, F = vals(real), vals(fake)
    mad = lambda v: 1.4826 * np.median(np.abs(v - np.median(v)))
    print(f"S{s} real n={len(R)} delta_t {np.round(R[:,0],2)} delta_a {np.round(R[:,1],2)} | fake n={len(F)}: std(dt)={F[:,0].std():.2f} mad={mad(F[:,0]):.2f}  std(da)={F[:,1].std():.2f} mad={mad(F[:,1]):.2f} ppt; sD_t real med {np.median(R[:,2]):.2f} fake med {np.median(F[:,2]):.2f}; worst fake |da| {np.sort(np.abs(F[:,1]))[-3:].round(1)}; floor max {np.abs(F[:,3]).max():.1f}")
