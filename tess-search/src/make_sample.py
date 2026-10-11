"""Select search targets from a TIC table: dwarfs/subgiants T<=12 plus cool dwarfs T<=13.5."""
import sys, pandas as pd, numpy as np
d = pd.read_parquet(sys.argv[1])
for c in ["Tmag", "Teff", "rad", "logg", "d"]: d[c] = pd.to_numeric(d[c], errors="coerce")
ok = (d.objType == "STAR") & ~d.disposition.isin(["DUPLICATE", "SPLIT", "ARTIFACT"])
main = ok & (d.Tmag <= 12.0) & ((d.rad <= 2.5) | d.rad.isna())
cool = ok & (d.Tmag <= 13.5) & (d.Teff < 4000) & (d.rad < 0.7)
s = d[main | cool].copy()
s["ID"] = s.ID.astype(np.int64)
s = s.sort_values("Tmag")
np.savetxt(sys.argv[2], s.ID.values, fmt="%d")
print(len(s), "targets;", int(cool.sum()), "cool dwarfs")
