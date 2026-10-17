"""Query MAST for all time-series products (TESS SPOC/TESS-SPOC/QLP/TGLC/ELEANOR/CDIPS/PATHOS/etc.) within 60" of the Mira,
and list which TIC IDs they belong to."""
import warnings, json, sys
warnings.filterwarnings("ignore")
from astroquery.mast import Observations, Catalogs
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
c = SkyCoord(260.977650, -51.323690, unit="deg")
obs = Observations.query_region(c, radius=60*u.arcsec)
print(len(obs))
obs = obs[obs["dataproduct_type"] == "timeseries"]
cols = ["obs_collection", "provenance_name", "target_name", "sequence_number", "t_exptime", "s_ra", "s_dec", "obsid", "project", "dataURL"]
obs.sort(["provenance_name", "target_name", "sequence_number"])
for r in obs:
    sep = SkyCoord(r["s_ra"], r["s_dec"], unit="deg").separation(c).arcsec
    print(r["obs_collection"], r["provenance_name"], r["target_name"], r["sequence_number"], r["t_exptime"], f"{sep:.1f}", r["obsid"])
obs.write("out/mast_timeseries_60arcsec.ecsv", overwrite=True)
