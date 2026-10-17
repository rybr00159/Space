import sys; sys.path.insert(0, "tic214311908/analysis/search/scripts")
from slib import *
import json
P = 0.64853; 
for star, model in [("M3V","B"),("M3V","A"),("pair","A"),("Mira","B")]:
    d = load_star(star, model)
    m = ~d.intr.values & np.isfinite(d.det.values)
    t = d.t.values[m]; f = d.det.values[m]-1; s = d.sector.values[m]
    for PP in (P, 2*P):
        ph = ((t - 3840.0)/PP) % 1
        nb = 20
        k = (ph*nb).astype(int)
        row=[]
        for ss in (39,66,93,103,104,'all'):
            sel = (s==ss) if ss!='all' else np.ones(len(s),bool)
            mb = np.array([np.median(f[sel & (k==i)]) for i in range(nb)])*1e3
            row.append(f"S{ss}: " + " ".join(f"{v:+.1f}" for v in mb))
        print(star, model, "P=%.5f"%PP); print("\n".join("   "+r for r in row))
