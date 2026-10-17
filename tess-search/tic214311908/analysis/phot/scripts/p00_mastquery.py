"""List all MAST time-series products for the three TICs (QLP, TESS-SPOC, SPOC, TGLC, GSFC-ELEANOR-LITE, CDIPS, PATHOS, ...)."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from astroquery.mast import Observations
rows = []
for tic, (ra, dec) in {214311908: (260.970485, -51.324202), 214311906: (260.977648, -51.323688), 214311910: (260.974900, -51.324069)}.items():
    obs = Observations.query_criteria(target_name=[str(tic), f"TIC {tic}", f"{tic:016d}"], dataproduct_type="timeseries")
    obs2 = Observations.query_region(f"{ra} {dec}", radius="0.003 deg")
    obs2 = obs2[obs2["dataproduct_type"] == "timeseries"]
    for o in list(obs) + list(obs2):
        rows.append(dict(tic=tic, obs_collection=o["obs_collection"], provenance_name=o["provenance_name"], target_name=o["target_name"],
                         sequence_number=o["sequence_number"], t_exptime=o["t_exptime"], obsid=o["obsid"], dataURL=o["dataURL"]))
d = pd.DataFrame(rows).drop_duplicates(subset=["obsid"])
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 80)
print(d.sort_values(["provenance_name", "sequence_number"]).to_string())
d.to_csv("tic214311908/analysis/phot/data/mast/mast_timeseries.csv", index=False)
