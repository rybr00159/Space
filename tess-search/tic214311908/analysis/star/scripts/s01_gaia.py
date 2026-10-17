"""Gaia DR3 data for TIC 214311908 (M3V) and its neighbourhood.
Queries (Gaia archive TAP via astroquery):
  - gaia_source within 120" of the M3V (all columns needed for census and diagnostics)
  - astrophysical_parameters for the three named stars + the 0.83" companion + the background star near the M3V
  - nss_two_body_orbit / nss_acceleration / nss_non_linear_spectro / vari_classifier / vari_summary for the same
  - dr2_neighbourhood -> DR2 astrometry for the same (DR2-DR3 PM difference)
  - external crossmatches: 2MASS, AllWISE best neighbours
Outputs: data/gaia_120.csv, data/gaia_ap.csv, data/gaia_nss_vari.json, data/gaia_dr2.csv, data/gaia_xm.csv
"""
import os, json, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from astroquery.gaia import Gaia
Gaia.ROW_LIMIT = -1
D = 'tic214311908/analysis/star/data'
M3V = 5925373272850784768; PAIRA = 5925373272850782208; PAIRB = 5925373272827253632; MIRA = 5925373268528684288
BKG = 5925373272827260544; REF = 5925373204131295872
ids = [M3V, PAIRA, PAIRB, MIRA, BKG, REF]
idl = ','.join(str(i) for i in ids)

def q(adql):
    j = Gaia.launch_job_async(adql)
    return j.get_results().to_pandas()

r = q(f"select ra, dec from gaiadr3.gaia_source where source_id={M3V}")
ra0, de0 = float(r.ra[0]), float(r.dec[0])
cols = ("source_id, ra, dec, ra_error, dec_error, parallax, parallax_error, pmra, pmra_error, pmdec, pmdec_error, ruwe, "
        "astrometric_params_solved, astrometric_excess_noise, astrometric_excess_noise_sig, astrometric_n_good_obs_al, "
        "astrometric_n_obs_al, astrometric_chi2_al, visibility_periods_used, astrometric_gof_al, nu_eff_used_in_astrometry, "
        "pseudocolour, ecl_lat, astrometric_sigma5d_max, "
        "ipd_frac_multi_peak, ipd_frac_odd_win, ipd_gof_harmonic_amplitude, ipd_gof_harmonic_phase, "
        "phot_g_mean_mag, phot_g_mean_flux, phot_g_mean_flux_error, phot_g_n_obs, phot_bp_mean_mag, phot_rp_mean_mag, "
        "phot_bp_mean_flux, phot_bp_mean_flux_error, phot_rp_mean_flux, phot_rp_mean_flux_error, "
        "phot_bp_n_obs, phot_rp_n_obs, phot_bp_n_blended_transits, phot_rp_n_blended_transits, phot_bp_n_contaminated_transits, "
        "phot_rp_n_contaminated_transits, phot_bp_rp_excess_factor, bp_rp, phot_variable_flag, "
        "radial_velocity, radial_velocity_error, rv_nb_transits, rv_amplitude_robust, rv_chisq_pvalue, rv_renormalised_gof, "
        "rv_method_used, rv_template_teff, grvs_mag, "
        "non_single_star, has_xp_continuous, has_xp_sampled, has_rvs, has_epoch_photometry, duplicated_source, "
        "teff_gspphot, logg_gspphot, mh_gspphot, distance_gspphot, l, b, "
        f"DISTANCE(POINT('ICRS', ra, dec), POINT('ICRS', {ra0}, {de0}))*3600 as sep_m3v")
g = q(f"select {cols} from gaiadr3.gaia_source where 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra0}, {de0}, {120/3600}))")
g.to_csv(f'{D}/gaia_120.csv', index=False)
print('gaia_source within 120":', len(g))
ap = q(f"select * from gaiadr3.astrophysical_parameters where source_id in ({idl})")
ap.to_csv(f'{D}/gaia_ap.csv', index=False)
out = {}
for tab in ['nss_two_body_orbit', 'nss_acceleration_astro', 'nss_non_linear_spectro', 'nss_vim_fl', 'vari_classifier_result',
            'vari_summary', 'vari_long_period_variable', 'vari_rotation_modulation', 'vari_eclipsing_binary']:
    try:
        t = q(f"select * from gaiadr3.{tab} where source_id in ({idl})")
        out[tab] = t.astype(object).where(pd.notnull(t), None).to_dict('records')
        print(tab, len(t))
    except Exception as e:
        out[tab] = f'error {e}'; print(tab, 'error', e)
json.dump(out, open(f'{D}/gaia_nss_vari.json', 'w'), indent=1, default=str)
dr2 = q(f"select n.dr3_source_id, n.dr2_source_id, n.angular_distance, n.magnitude_difference, d.ra, d.dec, d.parallax, "
        f"d.parallax_error, d.pmra, d.pmra_error, d.pmdec, d.pmdec_error, d.astrometric_excess_noise, d.phot_g_mean_mag, "
        f"d.bp_rp, d.radial_velocity, d.radial_velocity_error, d.astrometric_params_solved "
        f"from gaiadr3.dr2_neighbourhood as n join gaiadr2.gaia_source as d on n.dr2_source_id = d.source_id "
        f"where n.dr3_source_id in ({idl})")
dr2.to_csv(f'{D}/gaia_dr2.csv', index=False)
print(dr2.to_string())
xm = []
for t, nb in [('tmass_psc_xsc_best_neighbour', 'tmass_psc_xsc_join'), ('allwise_best_neighbour', None)]:
    try:
        x = q(f"select * from gaiadr3.{t} where source_id in ({idl})"); x['table'] = t; xm.append(x)
    except Exception as e:
        print(t, e)
pd.concat(xm).to_csv(f'{D}/gaia_xm.csv', index=False)
print(pd.concat(xm).to_string())
