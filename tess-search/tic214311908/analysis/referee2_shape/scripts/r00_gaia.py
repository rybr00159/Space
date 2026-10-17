"""r00: Gaia DR3 astrometric details of the M3V (for the RUWE / white-dwarf test). Output: out/gaia_astrom.json"""
import json
from astroquery.gaia import Gaia
q = """SELECT source_id, phot_g_mean_mag, bp_rp, parallax, parallax_error, ruwe, astrometric_n_obs_al, astrometric_n_good_obs_al,
astrometric_n_bad_obs_al, astrometric_chi2_al, astrometric_excess_noise, astrometric_excess_noise_sig, astrometric_params_solved,
astrometric_matched_transits, visibility_periods_used, astrometric_sigma5d_max, scan_direction_strength_k1, scan_direction_strength_k2,
scan_direction_strength_k3, scan_direction_strength_k4, scan_direction_mean_k1, scan_direction_mean_k2, scan_direction_mean_k3,
scan_direction_mean_k4, ipd_gof_harmonic_amplitude, ipd_gof_harmonic_phase, ipd_frac_multi_peak, ipd_frac_odd_win, nu_eff_used_in_astrometry,
pmra_error, pmdec_error, ra_error, dec_error, astrometric_gof_al, ecl_lat, ecl_lon, radial_velocity, radial_velocity_error, rv_nb_transits,
rv_amplitude_robust, rv_chisq_pvalue, rv_renormalised_gof, rv_expected_sig_to_noise, grvs_mag, rv_template_teff
FROM gaiadr3.gaia_source WHERE source_id = 5925373272850784768"""
r = Gaia.launch_job(q).get_results().to_pandas().iloc[0].to_dict()
r = {k: (float(v) if hasattr(v, '__float__') else str(v)) for k, v in r.items()}
json.dump(r, open('out/gaia_astrom.json', 'w'), indent=1)
for k, v in r.items(): print(k, v)
