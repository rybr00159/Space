import sys, json, signal
from astroquery.gaia import Gaia
from astropy.table import Table
Gaia.ROW_LIMIT=-1
out=sys.argv[1]
sid=6571343840265900544
qs={
 'dr3_source':f"select * from gaiadr3.gaia_source where source_id={sid}",
 'dr3_ap':f"select * from gaiadr3.astrophysical_parameters where source_id={sid}",
 'dr3_nss':f"select * from gaiadr3.nss_two_body_orbit where source_id={sid}",
 'dr3_nss_acc':f"select * from gaiadr3.nss_acceleration_astro where source_id={sid}",
 'dr3_vari':f"select * from gaiadr3.vari_summary where source_id={sid}",
 'dr3_varclass':f"select * from gaiadr3.vari_classifier_result where source_id={sid}",
 'dr3_eb':f"select * from gaiadr3.vari_eclipsing_binary where source_id={sid}",
 'dr2_source':f"select * from gaiadr2.gaia_source where source_id={sid}",
 'dr3_neigh':"select source_id,ra,dec,phot_g_mean_mag,bp_rp,parallax,parallax_error,pmra,pmdec,ruwe,ipd_frac_multi_peak, DISTANCE(POINT(ra,dec),POINT(330.887948,-42.835511))*3600 as sep from gaiadr3.gaia_source where 1=CONTAINS(POINT(ra,dec),CIRCLE(330.887948,-42.835511,0.0333)) order by sep",
 'dr3_xp':f"select * from gaiadr3.xp_summary where source_id={sid}",
 'dr3_tmass':f"select * from gaiadr3.tmass_psc_xsc_best_neighbour where source_id={sid}",
}
for k,q in qs.items():
    try:
        j=Gaia.launch_job(q); t=j.get_results()
        t.write(f"{out}/{k}.ecsv",overwrite=True); print(k,len(t))
    except Exception as e: print(k,'ERR',e)
