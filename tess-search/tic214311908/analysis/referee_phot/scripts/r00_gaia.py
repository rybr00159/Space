"""Referee: own Gaia DR3 query (330" around the Mira) + TIC query for nearby stars (TIC Tmag used as catalogue flux).
Output: data/gaia330.csv, data/tic330.csv"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from astroquery.gaia import Gaia
from astroquery.mast import Catalogs
from astropy.coordinates import SkyCoord
import astropy.units as u
D = "tic214311908/analysis/referee_phot/data"
RA, DE = 260.977650, -51.323690
q = f"""SELECT source_id, ra, dec, parallax, parallax_error, pmra, pmdec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag,
 bp_rp, ruwe, has_xp_sampled, phot_bp_rp_excess_factor, ipd_frac_multi_peak
 FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {RA}, {DE}, {330/3600}))"""
g = Gaia.launch_job_async(q).get_results().to_pandas()
g.columns = [c.lower() for c in g.columns]
g.to_csv(f"{D}/gaia330.csv", index=False)
print(len(g), "Gaia sources")
t = Catalogs.query_region(SkyCoord(RA, DE, unit="deg"), radius=330 * u.arcsec, catalog="TIC").to_pandas()
keep = ["ID", "GAIA", "ra", "dec", "Tmag", "e_Tmag", "Teff", "rad", "mass", "Jmag", "Hmag", "Kmag", "GAIAmag", "gaiabp", "gaiarp", "objType", "disposition", "duplicate_id", "pmRA", "pmDEC", "TESSflag"]
t = t[[c for c in keep if c in t.columns]]
t.to_csv(f"{D}/tic330.csv", index=False)
print(len(t), "TIC sources")
print(t[t.ID.astype(str).isin(["214311906", "214311908", "214311910"])].to_string())
