"""Details of the non-trivial VizieR catalogue entries near the M3V found by c08: El-Badry+2021 wide-binary catalogue
(J/MNRAS/506/2269), eROSITA eRASS1 main (J/A+A/682/A34), and the descriptions of the other catalogues listing the M3V.
Output: out/vizier_extra.json"""
import json, warnings
warnings.filterwarnings("ignore")
import numpy as np
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import astropy.units as u
c = SkyCoord(260.9696737, -51.3238212, unit="deg")
v = Vizier(columns=["**", "+_r"], row_limit=10)
out = {}
for cat, rad in (("J/MNRAS/506/2269/catalog", 5), ("J/A+A/682/A34/erass1-m", 20), ("J/AJ/161/234/table2", 5), ("J/AJ/161/234/table3", 5),
                 ("J/AJ/161/234/table4", 5), ("J/AJ/165/267/table1", 5), ("J/A+A/704/A237/80pcbpok", 5), ("J/AJ/158/93/table2", 5),
                 ("J/AJ/159/84/table2", 5), ("J/ApJ/891/58/table6", 5), ("J/A+A/653/A98/aspic1_1", 5), ("J/A+A/657/A7/tablea3", 5),
                 ("J/A+A/684/A121/like_ctp", 5), ("J/A+A/692/A150/allskytg", 5), ("J/A+A/669/A104/catalog", 5), ("II/374/table2", 5)):
    try:
        r = v.query_region(c, radius=rad * u.arcsec, catalog=cat)
        if len(r) == 0:
            out[cat] = None; continue
        t = r[0]; desc = t.meta.get("description", "")
        rows = [{k: (str(t[k][i])) for k in t.colnames} for i in range(len(t))]
        out[cat] = dict(description=desc, rows=rows)
        print("==", cat, "|", desc, "| n =", len(t))
        print("  ", {k: rows[0][k] for k in list(rows[0])[:40]})
    except Exception as e:
        out[cat] = f"error {e}"; print(cat, "error", e)
json.dump(out, open("out/vizier_extra.json", "w"), indent=1)
