"""q01: Gaia DR3 query (read-only ADQL): astrometry/photometry of all sources within 25" of the M3V, XP synthetic
photometry (gaiadr3.synthetic_photometry_gspc) where available, and the Mira's vari_long_period_variable row.
Output: data/gaia_cone25.csv, data/gaia_gspc.csv, data/gaia_mira_lpv.csv"""
import os
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'): os.environ[k]='1'
from astroquery.gaia import Gaia
Gaia.ROW_LIMIT=-1
D=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+'/data'
q1="""SELECT source_id, ra, dec, ref_epoch, parallax, parallax_error, pmra, pmdec, radial_velocity, radial_velocity_error, rv_method_used,
phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, bp_rp, ruwe, astrometric_params_solved, has_xp_continuous, phot_variable_flag,
DISTANCE(POINT(ra,dec), POINT(260.969674,-51.323821))*3600 AS sep
FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT(ra,dec), CIRCLE(260.969674,-51.323821, 25/3600.))"""
r=Gaia.launch_job(q1).get_results().to_pandas(); r.to_csv(f'{D}/gaia_cone25.csv',index=False)
print(r[['source_id','sep','phot_g_mean_mag','bp_rp','parallax','pmra','pmdec','ruwe','has_xp_continuous','radial_velocity']].sort_values('sep').to_string())
ids=','.join(str(x) for x in r.source_id)
q2=f"""SELECT * FROM gaiadr3.synthetic_photometry_gspc WHERE source_id IN ({ids})"""
s=Gaia.launch_job(q2).get_results().to_pandas(); s.to_csv(f'{D}/gaia_gspc.csv',index=False)
cols=[c for c in s.columns if c.endswith('_mag') and any(b in c for b in ('johnson_v','cousins_i','sdss_g','sdss_r','sdss_i','sdss_z','decam_z','ps1_z','cousins_r'))]
print(s[['source_id']+cols].to_string())
q3=f"""SELECT * FROM gaiadr3.vari_long_period_variable WHERE source_id IN ({ids})"""
l=Gaia.launch_job(q3).get_results().to_pandas(); l.to_csv(f'{D}/gaia_mira_lpv.csv',index=False); print(l.T)
