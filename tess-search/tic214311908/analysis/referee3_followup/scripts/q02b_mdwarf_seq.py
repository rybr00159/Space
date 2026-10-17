"""q02b: M-dwarf main sequence with XP synthetic photometry and 2MASS Ks across BP-RP 2.2-3.8 (d < 20 pc), read-only ADQL.
Two-step: gaia_source x synthetic_photometry_gspc, then the 2MASS best-neighbour + tmass_original_valid for those ids.
Output: data/mdw_seq.csv"""
import os
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'): os.environ[k]='1'
import pandas as pd, numpy as np
from astroquery.gaia import Gaia
Gaia.ROW_LIMIT=-1
D=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+'/data'
q="""SELECT g.source_id, g.parallax, g.phot_g_mean_mag, g.bp_rp, g.ruwe,
s.g_sdss_mag, s.r_sdss_mag, s.i_sdss_mag, s.z_sdss_mag, s.v_jkc_mag, s.i_jkc_mag, s.y_ps1_mag
FROM gaiadr3.gaia_source AS g JOIN gaiadr3.synthetic_photometry_gspc AS s ON s.source_id = g.source_id
WHERE g.parallax > 50 AND g.ruwe < 1.2 AND g.ipd_frac_multi_peak = 0 AND g.non_single_star = 0 AND g.bp_rp > 2.2 AND g.bp_rp < 3.8"""
a=Gaia.launch_job_async(q).get_results().to_pandas(); a.to_csv(f'{D}/mdw_seq_step1.csv',index=False); print('gaia+gspc', len(a), flush=True)
out=[]
ids=a.source_id.values
for i in range(0,len(ids),800):
    q2=("SELECT x.source_id, t.ks_m, t.ks_msigcom, t.j_m, t.ph_qual FROM gaiadr3.tmass_psc_xsc_best_neighbour AS x "
        "JOIN gaiadr1.tmass_original_valid AS t ON t.designation = x.original_ext_source_id WHERE x.source_id IN (%s)"%','.join(str(v) for v in ids[i:i+800]))
    out.append(Gaia.launch_job(q2).get_results().to_pandas())
t=pd.concat(out); r=a.merge(t,on='source_id'); r=r[r.ph_qual.astype(str).str[2]=='A']
r['MKs']=r.ks_m+5*np.log10(r.parallax/100.); r.to_csv(f'{D}/mdw_seq.csv',index=False); print('final', len(r))
