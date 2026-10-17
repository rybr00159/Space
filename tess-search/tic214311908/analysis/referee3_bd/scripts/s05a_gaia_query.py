"""s05a: Gaia DR3 archive queries (read-only): (1) the target and the wide companions' RV columns; (2) an independent RV-error
comparison sample: rv_method_used = 2, |grvs_mag - target| < 0.3, rv_template_teff <= 3700, parallax > 5 mas,
2.4 < bp_rp < 3.1, 15 <= rv_nb_transits <= 30.
Outputs: data/gaia_target.csv, data/gaia_rvcomp.csv"""
import sys; sys.path.insert(0, 'tic214311908/analysis/referee3_bd/scripts')
from rlib3 import *
from astroquery.gaia import Gaia
Gaia.ROW_LIMIT = -1
cols = """source_id, ra, dec, parallax, bp_rp, phot_g_mean_mag, grvs_mag, radial_velocity, radial_velocity_error, rv_method_used,
rv_nb_transits, rv_nb_deblended_transits, rv_visibility_periods_used, rv_expected_sig_to_noise, rv_renormalised_gof, rv_chisq_pvalue,
rv_amplitude_robust, rv_template_teff, rv_template_logg, rv_template_fe_h, rv_atm_param_origin, ruwe, ipd_frac_multi_peak, non_single_star"""
q1 = f"""SELECT {cols} FROM gaiadr3.gaia_source WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', 260.969674, -51.323821, 0.0062))"""
t1 = Gaia.launch_job(q1).get_results().to_pandas()
t1.to_csv(f'{R3}/data/gaia_target.csv', index=False)
print(t1[['source_id', 'phot_g_mean_mag', 'grvs_mag', 'radial_velocity', 'radial_velocity_error', 'rv_method_used', 'rv_nb_transits', 'rv_expected_sig_to_noise', 'rv_template_teff', 'ruwe', 'parallax']].to_string(), flush=True)
tg = t1[t1.source_id == 5925373272850784768].iloc[0]
q2 = f"""SELECT {cols} FROM gaiadr3.gaia_source WHERE rv_method_used = 2 AND radial_velocity_error IS NOT NULL
AND grvs_mag BETWEEN {tg.grvs_mag - 0.3} AND {tg.grvs_mag + 0.3} AND rv_template_teff <= 3700 AND parallax > 5
AND bp_rp BETWEEN 2.4 AND 3.1 AND rv_nb_transits BETWEEN 15 AND 30"""
t2 = Gaia.launch_job_async(q2).get_results().to_pandas()
t2.to_csv(f'{R3}/data/gaia_rvcomp.csv', index=False)
print(len(t2), t2.radial_velocity_error.describe(), flush=True)
