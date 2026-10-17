"""Comparison sample of nearby M dwarfs with colours like TIC 214311908 (BP-RP 2.756).
Gaia DR3: parallax > 10 mas, parallax_over_error > 20, BP-RP 2.45-3.05, in three separate queries
(gaia_source; 2MASS photometry via tmass_psc_xsc_best_neighbour + tmass_original_valid; DR2 astrometry via
dr2_neighbourhood; DSC class probabilities), merged on source_id.
Used for: the empirical main sequence (over-luminosity / photometric binarity), the distribution of RUWE, IPD
diagnostics and C* of comparable stars, the calibration of DR2-DR3 proper-motion differences, and the DSC
'binary star' probability of normal M dwarfs.
Output: data/mdwarf_sample.csv"""
import warnings, time
warnings.filterwarnings('ignore')
import pandas as pd
from astroquery.gaia import Gaia
Gaia.ROW_LIMIT = -1
D = 'tic214311908/analysis/star/data'
W = "g.parallax > 10 and g.parallax_over_error > 20 and g.bp_rp between 2.45 and 3.05"
def run(q, lab):
    for k in range(3):
        try:
            t = Gaia.launch_job_async(q).get_results().to_pandas(); print(lab, len(t), flush=True); return t
        except Exception as e:
            print(lab, 'retry', e, flush=True); time.sleep(20)
    raise RuntimeError(lab)
g = run(f"""select g.source_id, g.ra, g.dec, g.l, g.b, g.parallax, g.parallax_error, g.pmra, g.pmra_error, g.pmdec, g.pmdec_error,
g.phot_g_mean_mag, g.bp_rp, g.phot_bp_mean_mag, g.phot_rp_mean_mag, g.phot_bp_rp_excess_factor, g.ruwe,
g.ipd_frac_multi_peak, g.ipd_gof_harmonic_amplitude, g.ipd_frac_odd_win, g.phot_bp_n_blended_transits, g.phot_bp_n_obs,
g.astrometric_n_good_obs_al, g.visibility_periods_used, g.astrometric_excess_noise, g.non_single_star,
g.radial_velocity, g.radial_velocity_error, g.rv_nb_transits, g.nu_eff_used_in_astrometry, g.ecl_lat, g.astrometric_params_solved, g.pseudocolour
from gaiadr3.gaia_source as g where {W}""", 'gaia')
t = run(f"""select g.source_id, t.ks_m, t.ks_msigcom, t.j_m, t.j_msigcom, t.h_m, t.h_msigcom, t.ph_qual
from gaiadr3.gaia_source as g join gaiadr3.tmass_psc_xsc_best_neighbour as x on g.source_id = x.source_id
join gaiadr1.tmass_original_valid as t on x.original_ext_source_id = t.designation where {W}""", '2mass')
d = run(f"""select g.source_id, n.angular_distance as dr2_dist, n.magnitude_difference as dr2_dmag, d.pmra as pmra_dr2,
d.pmdec as pmdec_dr2, d.pmra_error as pmra_error_dr2, d.pmdec_error as pmdec_error_dr2, d.parallax as plx_dr2, d.parallax_error as plx_error_dr2
from gaiadr3.gaia_source as g join gaiadr3.dr2_neighbourhood as n on g.source_id = n.dr3_source_id
join gaiadr2.gaia_source as d on n.dr2_source_id = d.source_id where {W}""", 'dr2')
d = d.sort_values('dr2_dist').drop_duplicates('source_id')
a = run(f"""select g.source_id, a.classprob_dsc_specmod_binarystar, a.classprob_dsc_combmod_binarystar
from gaiadr3.gaia_source as g join gaiadr3.astrophysical_parameters as a on g.source_id = a.source_id where {W}""", 'ap')
m = g.merge(t, on='source_id', how='left').merge(d, on='source_id', how='left').merge(a, on='source_id', how='left')
m.to_csv(f'{D}/mdwarf_sample.csv', index=False)
print(len(m))
