"""r06a: download Gaia DR3 NSS spectroscopic binaries (SB1/SB2) and eclipsing binaries among nearby M dwarfs
(parallax > 10, BP-RP 2.0-3.4) with 2MASS photometry, as a validation sample for the colour-residual twin discriminant.
Output: data/sb2/nss_mdwarfs.csv"""
import pandas as pd, time
from astroquery.gaia import Gaia
B = 'tic214311908/analysis/referee2_host'
Gaia.ROW_LIMIT = -1
q = """
select n.source_id, n.nss_solution_type, n.period, n.mass_ratio, n.mass_ratio_error, n.semi_amplitude_primary, n.semi_amplitude_secondary,
 n.eccentricity, n.goodness_of_fit, n.significance,
 g.ra, g.dec, g.parallax, g.parallax_error, g.pmra, g.pmdec, g.phot_g_mean_mag, g.phot_bp_mean_mag, g.phot_rp_mean_mag, g.bp_rp,
 g.phot_bp_rp_excess_factor, g.ruwe, g.ipd_frac_multi_peak, g.ipd_gof_harmonic_amplitude, g.radial_velocity, g.radial_velocity_error, g.grvs_mag,
 t.j_m, t.j_msigcom, t.h_m, t.h_msigcom, t.ks_m, t.ks_msigcom, t.ph_qual
from gaiadr3.nss_two_body_orbit as n
join gaiadr3.gaia_source as g on g.source_id = n.source_id
join gaiadr3.tmass_psc_xsc_best_neighbour as xt on xt.source_id = g.source_id
join gaiadr1.tmass_original_valid as t on t.designation = xt.original_ext_source_id
where g.parallax > 10 and g.bp_rp between 2.0 and 3.4
"""
for attempt in range(4):
    try:
        t = Gaia.launch_job_async(q).get_results().to_pandas(); break
    except Exception as e:
        print('retry', str(e)[:200], flush=True); time.sleep(20)
t.to_csv(f'{B}/data/sb2/nss_mdwarfs.csv', index=False)
print(len(t)); print(t.nss_solution_type.value_counts())
