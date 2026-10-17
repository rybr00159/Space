"""Gaia DR3 sources within 330" of the M3V (for plate calibration and star counts). Output: data/gaia_330.csv"""
from astroquery.gaia import Gaia
Gaia.ROW_LIMIT = -1
q = ("select source_id, ra, dec, parallax, parallax_error, pmra, pmdec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, bp_rp, ruwe, "
     "ipd_frac_multi_peak, astrometric_params_solved, DISTANCE(POINT('ICRS', ra, dec), POINT('ICRS', 260.969674, -51.323821))*3600 as sep_m3v "
     "from gaiadr3.gaia_source where 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', 260.969674, -51.323821, 0.09167))")
t = Gaia.launch_job_async(q).get_results().to_pandas()
t.to_csv('tic214311908/analysis/star/data/gaia_330.csv', index=False); print(len(t))
