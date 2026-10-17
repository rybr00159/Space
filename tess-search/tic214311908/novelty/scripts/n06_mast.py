"""n06: MAST inventory of time-series products within 5" of the M3V (and by TIC name), to see whether any product or
pipeline output (e.g. SPOC/TESS-SPOC DV reports, a new HLSP) has appeared since the ctx/ check. Read-only.
Output: out/mast_inventory.csv, out/mast_status.json
Usage: python -I n06_mast.py <novelty_dir>
"""
import sys, os, json, datetime
import pandas as pd
from astroquery.mast import Observations
import astropy.units as u
from astropy.coordinates import SkyCoord

D = sys.argv[1]
OUT = os.path.join(D, "out")
NOW = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
c = SkyCoord(260.970484934795, -51.3242024657333, unit="deg")
obs = Observations.query_region(c, radius=5 * u.arcsec).to_pandas()
keep = [k for k in ("obs_collection", "provenance_name", "project", "target_name", "sequence_number", "dataproduct_type",
                    "t_min", "t_max", "s_ra", "s_dec", "obs_id", "instrument_name") if k in obs]
o = obs[keep].copy()
o = o[~o.obs_collection.isin(["TESS"]) | o.dataproduct_type.isin(["timeseries"])]   # drop FFI image rows
o.to_csv(os.path.join(OUT, "mast_inventory.csv"), index=False)
summ = (o.groupby(["obs_collection", "provenance_name", "target_name"], dropna=False)
        .agg(n=("obs_id", "size"), sectors=("sequence_number", lambda s: sorted(set(int(x) for x in s if pd.notna(x)))))
        .reset_index())
print(summ.to_string())
spoc = obs[(obs.obs_collection == "TESS") & (obs.dataproduct_type == "timeseries")]
res = {"queried_utc": NOW, "n_rows": int(len(obs)), "summary": summ.astype(str).to_dict("records"),
       "tess_mission_timeseries_rows": spoc[keep].astype(str).to_dict("records")}
json.dump(res, open(os.path.join(OUT, "mast_status.json"), "w"), indent=1)
print("TESS mission time-series rows (SPOC 2-min/20-s):", len(spoc))
