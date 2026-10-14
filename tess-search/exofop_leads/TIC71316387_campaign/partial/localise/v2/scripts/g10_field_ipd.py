"""Gaia DR3 field comparison: IPD / blending diagnostics for stars of similar G within 3 deg of the target,
plus their nearest Gaia neighbour distance (to separate 'isolated' from resolved pairs)."""
import pyvo, sys, time
import numpy as np
R = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v2'
svc = pyvo.dal.TAPService('https://gea.esac.esa.int/tap-server/tap')
cols = """g.source_id, g.ra, g.dec, g.phot_g_mean_mag, g.bp_rp, g.parallax, g.ruwe, g.astrometric_excess_noise, g.astrometric_excess_noise_sig,
g.ipd_frac_multi_peak, g.ipd_frac_odd_win, g.ipd_gof_harmonic_amplitude, g.ipd_gof_harmonic_phase, g.visibility_periods_used,
g.astrometric_n_good_obs_al, g.astrometric_matched_transits, g.phot_bp_n_obs, g.phot_rp_n_obs, g.phot_bp_n_blended_transits,
g.phot_rp_n_blended_transits, g.phot_bp_n_contaminated_transits, g.phot_rp_n_contaminated_transits, g.phot_bp_rp_excess_factor,
g.duplicated_source, g.non_single_star, g.astrometric_params_solved, g.nu_eff_used_in_astrometry, g.pseudocolour, g.ecl_lat"""
q = f"""SELECT {cols} FROM gaiadr3.gaia_source AS g
WHERE 1=CONTAINS(POINT(g.ra,g.dec), CIRCLE(330.887948,-42.835511, 3.0)) AND g.phot_g_mean_mag BETWEEN 10.75 AND 12.75"""
t0 = time.time()
job = svc.submit_job(q, maxrec=200000)
job.run()
job.wait(phases=['COMPLETED', 'ERROR', 'ABORTED'], timeout=600)
r = job.fetch_result().to_table()
print(len(r), 'stars', time.time() - t0)
r.write(f'{R}/gaia/field3deg_G10.75-12.75.csv', overwrite=True, format='csv')
