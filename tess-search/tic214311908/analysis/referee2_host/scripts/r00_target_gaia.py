"""r00: fetch the M3V's Gaia DR3 row (+ 2MASS/AllWISE best neighbours) directly from the Gaia archive to verify the
inputs used by the host analysis. Output: out/r00_target_gaia.json"""
import json, numpy as np
from astroquery.gaia import Gaia
B = 'tic214311908/analysis/referee2_host'
sid = 5925373272850784768
q = f"select * from gaiadr3.gaia_source where source_id = {sid}"
r = Gaia.launch_job_async(q).get_results().to_pandas().iloc[0]
d = {k: (None if (isinstance(v, float) and not np.isfinite(v)) else (v.item() if hasattr(v, 'item') else (v.decode() if isinstance(v, bytes) else v))) for k, v in r.items()}
keys = ['source_id', 'ra', 'dec', 'parallax', 'parallax_error', 'pmra', 'pmdec', 'ruwe', 'phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag', 'bp_rp',
        'phot_bp_rp_excess_factor', 'phot_bp_n_blended_transits', 'phot_rp_n_blended_transits', 'phot_bp_n_obs', 'phot_rp_n_obs', 'ipd_frac_multi_peak',
        'ipd_gof_harmonic_amplitude', 'ipd_frac_odd_win', 'astrometric_excess_noise', 'astrometric_excess_noise_sig', 'non_single_star', 'visibility_periods_used',
        'radial_velocity', 'radial_velocity_error', 'rv_method_used', 'rv_nb_transits', 'rv_nb_deblended_transits', 'rv_visibility_periods_used',
        'rv_expected_sig_to_noise', 'rv_template_teff', 'rv_template_logg', 'rv_template_fe_h', 'rv_chisq_pvalue', 'rv_renormalised_gof', 'rv_amplitude_robust',
        'grvs_mag', 'grvs_mag_error', 'phot_variable_flag', 'has_xp_continuous', 'has_xp_sampled', 'j_m', 'j_msigcom', 'h_m', 'h_msigcom', 'ks_m', 'ks_msigcom', 'ph_qual',
        'phot_g_mean_flux_over_error', 'phot_bp_mean_flux_over_error', 'phot_rp_mean_flux_over_error', 'teff_gspphot', 'mh_gspphot', 'logg_gspphot', 'distance_gspphot']
out = {k: d.get(k) for k in keys}
json.dump(out, open(f'{B}/out/r00_target_gaia.json', 'w'), indent=1, default=str)
for k, v in out.items(): print(k, v)
