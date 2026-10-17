import sys, time
sys.path.insert(0, "tic214311908/analysis/phot/scripts")
from photlib import *
from p05_phot import setup
for s in SECTORS:
    S = setup(s); C, W, J, g, prf, rBA = S
    mira = g[g.name == "Mira"].iloc[0]; m3v = g[g.name == "M3V"].iloc[0]
    cand = []
    for r in g.itertuples():
        if not (11.5 <= r.T <= 15.0): continue
        if not (2.5 < r.x < C.nx - 3.5 and 2.5 < r.y < C.ny - 3.5): continue
        dm = np.hypot(r.x - mira.x, r.y - mira.y)
        if dm < 3.0: continue
        d = np.hypot(g.x - r.x, g.y - r.y)
        if np.any((d > 0) & (d < 1.5) & (g["T"] < r.T + 1.5)): continue
        cand.append((r.source_id, round(r.x,1), round(r.y,1), round(r.T,2), round(dm,1)))
    print(s, len(cand), cand)
