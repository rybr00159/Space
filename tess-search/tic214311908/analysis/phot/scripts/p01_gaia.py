"""Gaia DR3 sources within 360" of the Mira (covers the full 21x21 px TESScut cutouts), with TESS magnitudes
(Stassun et al. 2019 colour relation; BP-RP clipped to [-0.2, 3.5]; G+0.5 fallback without colour)
and TIC IDs/Tmag for the three stars of interest from the TIC (MAST)."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from astroquery.gaia import Gaia
from astroquery.mast import Catalogs
D = "tic214311908/analysis/phot/data/gaia"
RA, DE = 260.977650, -51.323690
q = f"""SELECT source_id, ra, dec, ra_error, dec_error, parallax, parallax_error, pmra, pmdec, phot_g_mean_mag AS G,
 phot_bp_mean_mag AS BP, phot_rp_mean_mag AS RP, bp_rp, ruwe, ipd_frac_multi_peak, ipd_gof_harmonic_amplitude,
 phot_variable_flag, astrometric_params_solved, duplicated_source
 FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {RA}, {DE}, {360/3600}))"""
g = Gaia.launch_job_async(q).get_results().to_pandas()
g.columns = [c.lower() if c not in ("G", "BP", "RP") else c for c in g.columns]
if "g" in g.columns and "G" not in g.columns: g = g.rename(columns={"g": "G", "bp": "BP", "rp": "RP"})
c = np.clip(g.bp_rp.values, -0.2, 3.5)
T = g.G - 0.00522555 * c**3 + 0.0891337 * c**2 - 0.633923 * c + 0.0324473
T = np.where(np.isfinite(g.bp_rp), T, g.G - 0.5)
g["T"] = T
dE = (g.ra - RA) * np.cos(np.radians(DE)) * 3600; dN = (g.dec - DE) * 3600
g["dE"], g["dN"], g["sep"] = dE, dN, np.hypot(dE, dN)
g = g.sort_values("G").reset_index(drop=True)
g.to_csv(f"{D}/gaia_360arcsec.csv", index=False)
print(len(g), "sources;", (g["T"] < 16).sum(), "with T<16")
t = Catalogs.query_criteria(catalog="TIC", ID=[214311906, 214311908, 214311910]).to_pandas()
print(t[["ID", "GAIA", "Tmag", "e_Tmag", "Teff", "rad", "mass", "contratio", "ra", "dec", "pmRA", "pmDEC"]].to_string())
t.to_csv(f"{D}/tic_three.csv", index=False)
