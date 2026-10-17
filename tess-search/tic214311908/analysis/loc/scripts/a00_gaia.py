"""Gaia DR3 sources within 330 arcsec of the Mira (covers the 21x21 px cutouts), plus astrophysical parameters,
and TIC rows of the three named stars. Output: data/gaia_330.csv, data/tic_trio.csv"""
import numpy as np, pandas as pd
from astroquery.gaia import Gaia
from astroquery.mast import Catalogs
OUT = 'tic214311908/analysis/loc/data'
RA, DEC = 260.977650, -51.323690
q = f"""SELECT g.source_id, g.ra, g.dec, g.ra_error, g.dec_error, g.parallax, g.parallax_error, g.pmra, g.pmdec, g.pmra_error, g.pmdec_error,
 g.phot_g_mean_mag, g.phot_bp_mean_mag, g.phot_rp_mean_mag, g.bp_rp, g.ruwe, g.ipd_frac_multi_peak, g.ipd_gof_harmonic_amplitude,
 g.phot_bp_rp_excess_factor, g.astrometric_params_solved, g.phot_variable_flag, g.duplicated_source, g.radial_velocity,
 a.teff_gspphot, a.logg_gspphot, a.radius_gspphot, a.distance_gspphot, a.classprob_dsc_combmod_star
FROM gaiadr3.gaia_source g LEFT JOIN gaiadr3.astrophysical_parameters a ON g.source_id = a.source_id
WHERE 1=CONTAINS(POINT('ICRS', g.ra, g.dec), CIRCLE('ICRS', {RA}, {DEC}, {330/3600}))"""
j = Gaia.launch_job_async(q); t = j.get_results().to_pandas()
t.columns = [c.lower() for c in t.columns]
t['dE'] = (t.ra - RA) * np.cos(np.radians(DEC)) * 3600; t['dN'] = (t.dec - DEC) * 3600; t['sep'] = np.hypot(t.dE, t.dN)
t.sort_values('sep').to_csv(f'{OUT}/gaia_330.csv', index=False)
print(len(t), 'Gaia sources; G<16:', (t.phot_g_mean_mag < 16).sum())
c = Catalogs.query_criteria(catalog='TIC', ID=[214311906, 214311908, 214311910])
c.to_pandas().to_csv(f'{OUT}/tic_trio.csv', index=False)
print(c['ID', 'Tmag', 'e_Tmag', 'Teff', 'rad', 'mass', 'GAIA', 'Jmag', 'Kmag', 'contratio'])
