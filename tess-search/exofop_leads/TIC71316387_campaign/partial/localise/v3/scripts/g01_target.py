import json, sys
import pyvo
svc = pyvo.dal.TAPService('https://gea.esac.esa.int/tap-server/tap')
q = "SELECT * FROM gaiadr3.gaia_source WHERE source_id = 6571343840265900544"
r = svc.run_sync(q).to_table()
row = {c: (r[c][0].item() if hasattr(r[c][0], 'item') else str(r[c][0])) for c in r.colnames}
json.dump(row, open('target_gaia_source_v3.json', 'w'), indent=1, default=str)
keys = ['ra','dec','parallax','parallax_error','pmra','pmdec','ruwe','astrometric_excess_noise','astrometric_excess_noise_sig',
 'astrometric_n_obs_al','astrometric_n_good_obs_al','astrometric_n_bad_obs_al','astrometric_gof_al','astrometric_chi2_al',
 'astrometric_params_solved','visibility_periods_used','astrometric_sigma5d_max','nu_eff_used_in_astrometry','pseudocolour',
 'ipd_gof_harmonic_amplitude','ipd_gof_harmonic_phase','ipd_frac_multi_peak','ipd_frac_odd_win','scan_direction_strength_k1',
 'scan_direction_strength_k2','scan_direction_strength_k3','scan_direction_strength_k4','scan_direction_mean_k1','scan_direction_mean_k2',
 'scan_direction_mean_k3','scan_direction_mean_k4','duplicated_source','matched_transits','new_matched_transits','matched_transits_removed',
 'phot_g_n_obs','phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag','phot_bp_rp_excess_factor','phot_bp_n_blended_transits',
 'phot_rp_n_blended_transits','phot_bp_n_contaminated_transits','phot_rp_n_contaminated_transits','phot_bp_n_obs','phot_rp_n_obs','bp_rp',
 'phot_variable_flag','non_single_star','rv_nb_transits','radial_velocity','radial_velocity_error','rv_chisq_pvalue','rv_renormalised_gof',
 'rv_amplitude_robust','rv_expected_sig_to_noise','rv_template_teff','classprob_dsc_combmod_star','teff_gspphot','logg_gspphot','mh_gspphot',
 'ag_gspphot','distance_gspphot','phot_g_mean_flux_over_error','phot_g_mean_flux_error','phot_g_mean_flux']
for k in keys:
    print(f'{k:36s} {row.get(k)}')
