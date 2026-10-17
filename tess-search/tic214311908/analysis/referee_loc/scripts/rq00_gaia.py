"""Independent Gaia DR3 query (cone 400" around the Mira) -> data/gaia_400.csv"""
import warnings; warnings.filterwarnings('ignore')
from astroquery.gaia import Gaia
Gaia.MAIN_GAIA_TABLE = 'gaiadr3.gaia_source'
q = """SELECT source_id, ra, dec, ra_error, dec_error, parallax, parallax_error, pmra, pmdec, phot_g_mean_mag,
 phot_bp_mean_mag, phot_rp_mean_mag, bp_rp, ruwe, ipd_frac_multi_peak, ipd_gof_harmonic_amplitude,
 phot_bp_rp_excess_factor, astrometric_params_solved, phot_variable_flag, teff_gspphot
 FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec), CIRCLE('ICRS',260.977650,-51.323690,0.1111))"""
j = Gaia.launch_job_async(q)
t = j.get_results().to_pandas()
t.to_csv('tic214311908/analysis/referee_loc/data/gaia_400.csv', index=False)
print(len(t))
