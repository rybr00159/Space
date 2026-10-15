"""Query Gaia DR3 'twin' comparison stars (similar colour/brightness, high latitude near the target) with 2MASS and AllWISE photometry."""
import sys, warnings; warnings.filterwarnings('ignore')
from astroquery.gaia import Gaia
Gaia.ROW_LIMIT=-1
out=sys.argv[1]; test=len(sys.argv)>2
base="""SELECT g.source_id,g.ra,g.dec,g.l,g.b,g.parallax,g.parallax_error,g.pmra,g.pmdec,g.ruwe,g.astrometric_excess_noise,g.astrometric_excess_noise_sig,
g.ipd_frac_multi_peak,g.ipd_gof_harmonic_amplitude,g.phot_g_mean_mag,g.phot_bp_mean_mag,g.phot_rp_mean_mag,g.bp_rp,g.phot_bp_rp_excess_factor,
g.phot_bp_n_blended_transits,g.phot_bp_n_obs,g.radial_velocity,g.radial_velocity_error,g.rv_nb_transits,g.rv_chisq_pvalue,g.rv_amplitude_robust,
g.teff_gspphot,g.logg_gspphot,g.mh_gspphot,g.ebpminrp_gspphot,g.non_single_star,g.has_xp_sampled,g.astrometric_n_good_obs_al,g.visibility_periods_used,
tm.j_m,tm.j_msigcom,tm.h_m,tm.h_msigcom,tm.ks_m,tm.ks_msigcom,tm.ph_qual AS tm_qual,
aw.w1mpro,aw.w1mpro_error,aw.w2mpro,aw.w2mpro_error,aw.w3mpro,aw.w3mpro_error,aw.ph_qual AS aw_qual,aw.cc_flags,aw.ext_flag
FROM gaiadr3.gaia_source AS g
JOIN gaiadr3.tmass_psc_xsc_best_neighbour AS xm ON g.source_id=xm.source_id
JOIN gaiadr1.tmass_original_valid AS tm ON xm.original_ext_source_id=tm.designation
JOIN gaiadr3.allwise_best_neighbour AS axm ON g.source_id=axm.source_id
JOIN gaiadr1.allwise_original_valid AS aw ON axm.allwise_oid=aw.allwise_oid
WHERE """
if test:
    q=base+"g.source_id=6571343840265900544"
else:
    q=base+("1=CONTAINS(POINT(g.ra,g.dec),CIRCLE(330.887948,-42.835511,7)) "
            "AND g.bp_rp BETWEEN 0.66 AND 0.92 AND g.phot_g_mean_mag BETWEEN 10.5 AND 13.0 AND g.parallax_over_error>30 "
            "AND g.parallax BETWEEN 1.5 AND 8")
j=Gaia.launch_job(q); t=j.get_results(); t.write(out,overwrite=True); print(len(t))
if test: print(t)
