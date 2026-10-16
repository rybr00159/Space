"""Depth of the 5.6915 d dips in the QLP light curves of the Mira (TIC 214311906), the M-dwarf pair TIC 214311910
(6.6") and the M3V TIC 214311908 (18"), per sector; plus the stars' relative positions."""
import sys, json, warnings
sys.path.insert(0, "/home/user/work/src"); warnings.filterwarnings("ignore")
sys.argv = [sys.argv[0], "214311906"]
import numpy as np, pandas as pd
from wotan import flatten
import stage2
from lcan import read, bin10, in_mask, trap
P, T0, T14, FL = 5.69149, 3840.913, 0.067, 0.78
g = pd.read_csv("out/214311906/gaia_150arcsec_all.csv")
ids = {214311906: 5925373268528684288, 214311910: 5925373272850782208, 214311908: 5925373272850784768}
t0r = g[g.source_id == ids[214311906]].iloc[0]
for tic, gid in ids.items():
    r = g[g.source_id == gid].iloc[0]
    dE = (r.ra - t0r.ra) * np.cos(np.radians(t0r.dec)) * 3600; dN = (r.dec - t0r.dec) * 3600
    print(f"TIC {tic}: offset from Mira (dE, dN) = ({dE:+.1f}, {dN:+.1f})\"  G {r.G:.2f} T {r['T']:.2f}")
rows = []
for tic in ids:
    p = stage2.products(tic)
    for _, r in p[p.prov == "QLP"].sort_values("sector").iterrows():
        t, q, cols, bkg, aper, crowd = read(stage2.fetch(r.uri), "QLP")
        for col in ("SAP_FLUX", "DET_FLUX"):
            if col not in cols: continue
            ok = (q == 0) & np.isfinite(t) & np.isfinite(cols[col]) & (cols[col] > 0)
            lvl = np.nanmedian(cols[col][ok])
            tt, ff = bin10(t[ok], cols[col][ok] / lvl, r.exptime)
            m = ~in_mask(tt, P, T0, T14, 1.0)
            tr = flatten(tt[m], ff[m], method="biweight", window_length=1.0, return_trend=True, break_tolerance=0.4)[1]
            k = np.isfinite(tr); f = ff / np.interp(tt, tt[m][k], tr[k])
            s = 1.4826 * np.median(np.abs(f[m] - np.median(f[m])))
            y = 1 - trap(tt, T0 + np.round((tt - T0) / P) * P, 1.0, T14, T14 * FL)
            A = np.column_stack([np.ones(len(tt)), y]); c, *_ = np.linalg.lstsq(A, f, rcond=None)
            e = s / np.sqrt(max(1, np.sum(y ** 2)))
            rows.append(dict(tic=tic, sector=int(r.sector), col=col, level=lvl, depth_ppm=-c[1] * 1e6, err_ppm=e * 1e6))
d = pd.DataFrame(rows)
pv = d[d.col == "SAP_FLUX"].pivot(index="sector", columns="tic", values=["depth_ppm", "err_ppm", "level"]).round(0)
print(pv.to_string())
print(d[d.col == "DET_FLUX"].pivot(index="sector", columns="tic", values="depth_ppm").round(0).to_string())
d.to_csv("out/214311906/trio_depths.csv", index=False)
