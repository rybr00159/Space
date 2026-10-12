"""For each final EB: Gaia DR3 stars within 42" (2 TESS px) bright enough to produce the primary eclipse if fully eclipsed."""
import json, sys, time, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from vet import gaia_neighbours

IN = sys.argv[1] if len(sys.argv) > 1 else "/home/user/work/vet/new_ebs_final.csv"
OUT = sys.argv[2] if len(sys.argv) > 2 else "/home/user/work/vet/eb_neighbours.json"
d = pd.read_csv(IN)
out = {}
for r in d.itertuples():
    for a in range(4):
        try:
            g = gaia_neighbours(r.ra, r.dec, r.Tmag, r.primary_pct / 100, radius=42.0); break
        except Exception as e:
            g = None; time.sleep(10 * (a + 1))
    if g is None: out[int(r.tic)] = "query failed"; continue
    g = g[g.sep > 2.0]  # exclude the target itself
    m = g[g.can_mimic]
    out[int(r.tic)] = [dict(sep=round(float(x.sep), 1), G=round(float(x.phot_g_mean_mag), 2), Test=round(float(x.tmag_est), 2)) for x in m.itertuples()]
    print(r.tic, f"depth {r.primary_pct:.2f}%:", out[int(r.tic)] if out[int(r.tic)] else "no Gaia neighbour within 42\" can produce it", flush=True)
json.dump(out, open(OUT, "w"), indent=1)
