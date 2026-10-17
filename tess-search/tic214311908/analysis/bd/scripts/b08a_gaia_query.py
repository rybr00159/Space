"""b08a: Gaia DR3 comparison sample for the RVS radial-velocity error test (read-only archive query).
M dwarfs like the target: parallax > 8 mas, 2.55 < bp_rp < 2.95, 12.2 < grvs_mag < 12.9, rv_method_used = 2 (combined-CCF
RV, as the target), with RV error, expected S/N, number of (deblended) transits, template, RUWE.
Output: data/gaia/rv_comparison.csv"""
import sys; sys.path.insert(0, 'tic214311908/analysis/bd/scripts')
from blib import *
from astroquery.gaia import Gaia
Gaia.ROW_LIMIT = -1
q = """SELECT source_id, parallax, bp_rp, phot_g_mean_mag, grvs_mag, radial_velocity, radial_velocity_error, rv_method_used,
rv_nb_transits, rv_nb_deblended_transits, rv_visibility_periods_used, rv_expected_sig_to_noise, rv_template_teff,
rv_template_logg, rv_template_fe_h, rv_atm_param_origin, ruwe, ipd_frac_multi_peak, non_single_star, l, b
FROM gaiadr3.gaia_source
WHERE parallax > 8 AND bp_rp BETWEEN 2.55 AND 2.95 AND grvs_mag BETWEEN 12.2 AND 12.9 AND rv_method_used = 2
AND radial_velocity_error IS NOT NULL"""
j = Gaia.launch_job_async(q)
t = j.get_results().to_pandas()
t.to_csv(f'{BD}/data/gaia/rv_comparison.csv', index=False)
print(len(t), t.radial_velocity_error.describe())
