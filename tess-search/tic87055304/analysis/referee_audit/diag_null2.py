"""Diagnostic: single-window delta(a)-delta(t) in clean fake windows: tails and their causes."""
import numpy as np, os
os.environ.setdefault('R3_CLEAN', 'mad10')
from rcommon import *
from r3_pixel_ratio import collect, clean, PHASES
for s in SECTORS:
    c = Cube(s); tp = (int(round(c.xt)), int(round(c.yt))); a = (int(round(c.xt + c.u[0])), int(round(c.yt + c.u[1])))
    real = collect(c); fake = [w for ph in PHASES for w in collect(c, ph)]
    real, fake, dr = clean(real, fake)
    def v(w):
        D = w['D'] - w['floor']; O = w['O'] - np.percentile(w['O'], 10)
        return 1e3 * (D[a[1], a[0]] / O[a[1], a[0]] - D[tp[1], tp[0]] / O[tp[1], tp[0]]), 1e3 * D[tp[1], tp[0]] / O[tp[1], tp[0]], w['tc'], w['floor_val'], w['sD'][tp[1], tp[0]]
    F = np.array([v(w) for w in fake]); R = np.array([v(w) for w in real])
    mad = 1.4826 * np.median(np.abs(F[:, 0] - np.median(F[:, 0])))
    o = np.argsort(np.abs(F[:, 0]))[-4:]
    print(f"S{s}: n_fake={len(F)} dropped real {dr}; L1 single std {F[:,0].std():.2f} mad {mad:.2f}; real L1 {np.round(R[:,0],2)}; tails:",
          [(round(F[k,0],1), round(F[k,2],2), round(F[k,3],2), round(F[k,4],2)) for k in o])
