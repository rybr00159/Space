"""h03: Field comparison sample of nearby M dwarfs from the Gaia DR3 archive (own query):
parallax > 10 mas, parallax_over_error > 20, 2.30 < BP-RP < 3.25, with 2MASS PSC and AllWISE best neighbours.
Columns: astrometry/photometry/quality (RUWE, IPD, C* inputs, blended transits), RV stats (method, error, nb transits,
grvs_mag, expected S/N), has_xp_continuous, 2MASS JHKs (+ph_qual), AllWISE W1-W3 (+cc_flags, ext_flag, ph_qual).
Output: data/field_mdwarfs.csv"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from astroquery.gaia import Gaia
B = 'tic214311908/analysis/host'
Gaia.ROW_LIMIT = -1
q = """
select g.source_id, g.ra, g.dec, g.l, g.b, g.parallax, g.parallax_error, g.pmra, g.pmdec, g.phot_g_mean_mag, g.phot_bp_mean_mag,
 g.phot_rp_mean_mag, g.bp_rp, g.phot_bp_rp_excess_factor, g.phot_bp_n_blended_transits, g.phot_rp_n_blended_transits, g.phot_bp_n_obs,
 g.ruwe, g.ipd_frac_multi_peak, g.ipd_gof_harmonic_amplitude, g.astrometric_excess_noise, g.non_single_star, g.visibility_periods_used,
 g.radial_velocity, g.radial_velocity_error, g.rv_method_used, g.rv_nb_transits, g.rv_nb_deblended_transits, g.rv_expected_sig_to_noise,
 g.rv_template_teff, g.rv_chisq_pvalue, g.rv_renormalised_gof, g.rv_amplitude_robust, g.grvs_mag, g.has_xp_continuous,
 t.j_m, t.j_msigcom, t.h_m, t.h_msigcom, t.ks_m, t.ks_msigcom, t.ph_qual, g.phot_rp_n_obs
from gaiadr3.gaia_source as g
join gaiadr3.tmass_psc_xsc_best_neighbour as xt on xt.source_id = g.source_id
join gaiadr1.tmass_original_valid as t on t.designation = xt.original_ext_source_id
where g.parallax > 10 and g.parallax_over_error > 20 and g.bp_rp between 2.30 and 3.25 and g.dec between {d0} and {d1}
"""
import time
T = []
for d0 in range(-90, 90, 30):
    for attempt in range(4):
        try:
            T.append(Gaia.launch_job_async(q.format(d0=d0, d1=d0 + 30)).get_results().to_pandas()); break
        except Exception as e:
            print('retry', d0, str(e)[:100], flush=True); time.sleep(20)
    print('dec band', d0, len(T[-1]), flush=True)
t = pd.concat(T, ignore_index=True).drop_duplicates('source_id')
# AllWISE best neighbours by source_id chunks
W = []
ids = t.source_id.values
for k in range(0, len(ids), 3000):
    sub = ','.join(str(int(s)) for s in ids[k:k + 3000])
    qq = f"""select xw.source_id, w.w1mpro, w.w1mpro_error, w.w2mpro, w.w2mpro_error, w.w3mpro, w.cc_flags as wise_cc, w.ext_flag as wise_ext,
 w.ph_qual as wise_ph_qual from gaiadr3.allwise_best_neighbour as xw join gaiadr1.allwise_original_valid as w
 on w.designation = xw.original_ext_source_id where xw.source_id in ({sub})"""
    for attempt in range(4):
        try:
            W.append(Gaia.launch_job_async(qq).get_results().to_pandas()); break
        except Exception as e:
            print('retry w', k, str(e)[:100], flush=True); time.sleep(20)
W = pd.concat(W, ignore_index=True).drop_duplicates('source_id')
t = t.merge(W, on='source_id', how='left')
print(len(t))
t.to_csv(f'{B}/data/field_mdwarfs.csv', index=False)
