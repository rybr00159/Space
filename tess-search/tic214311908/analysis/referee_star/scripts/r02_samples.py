"""Referee: comparison samples from the Gaia archive (independent of the investigator's s04).
(a) nearby M dwarfs: parallax > 12 mas, parallax_over_error > 20, BP-RP 2.55-2.95, with 2MASS (best neighbour) photometry.  -> data/mdw_cmd.csv
(b) local DR2-DR3 proper-motion calibration: all 5-parameter DR3 sources with G 13.5-15.5 within 1.0 deg of the M3V
    that have a DR2 counterpart (dr2_neighbourhood, angular_distance < 0.5"), plus the DR2 PM and errors. -> data/pm_local.csv"""
import warnings, time, sys; warnings.filterwarnings('ignore')
from astroquery.gaia import Gaia
Gaia.ROW_LIMIT = -1
D = 'tic214311908/analysis/referee_star/data'
def run(q, lab):
    for k in range(4):
        try:
            t = Gaia.launch_job_async(q).get_results().to_pandas(); print(lab, len(t), flush=True); return t
        except Exception as e: print(lab, 'retry', e, flush=True); time.sleep(20)
which = sys.argv[1]
if which == 'a':
    q = """select g.source_id, g.ra, g.dec, g.l, g.b, g.parallax, g.parallax_error, g.phot_g_mean_mag, g.bp_rp, g.phot_bp_mean_mag, g.phot_rp_mean_mag,
    g.phot_bp_rp_excess_factor, g.ruwe, g.ipd_frac_multi_peak, g.ipd_gof_harmonic_amplitude, g.radial_velocity, g.radial_velocity_error, g.grvs_mag,
    g.astrometric_params_solved, g.non_single_star,
    t.j_m, t.j_msigcom, t.h_m, t.h_msigcom, t.ks_m, t.ks_msigcom, t.ph_qual
    from gaiadr3.gaia_source as g
    join gaiadr3.tmass_psc_xsc_best_neighbour as x on g.source_id = x.source_id
    join gaiadr1.tmass_original_valid as t on x.original_ext_source_id = t.designation
    where g.parallax > 12 and g.parallax_over_error > 20 and g.bp_rp between 2.55 and 2.95"""
    run(q, 'mdw').to_csv(f'{D}/mdw_cmd.csv', index=False)
else:
    q = """select g.source_id, g.ra, g.dec, g.phot_g_mean_mag, g.bp_rp, g.parallax, g.pmra, g.pmdec, g.pmra_error, g.pmdec_error, g.ruwe,
    d.pmra as pmra2, d.pmdec as pmdec2, d.pmra_error as pmra2_err, d.pmdec_error as pmdec2_err, n.angular_distance
    from gaiadr3.gaia_source as g join gaiadr3.dr2_neighbourhood as n on g.source_id = n.dr3_source_id
    join gaiadr2.gaia_source as d on n.dr2_source_id = d.source_id
    where 1=contains(point(g.ra,g.dec), circle(260.969674,-51.323821,1.0)) and g.phot_g_mean_mag between 13.5 and 15.5
    and g.astrometric_params_solved > 3 and d.pmra is not null and n.angular_distance < 500"""
    run(q, 'pm').to_csv(f'{D}/pm_local.csv', index=False)
