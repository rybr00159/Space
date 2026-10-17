"""q02: empirical colours of nearby M dwarfs. Sample = ../referee_star/data/mdw_cmd.csv (Gaia DR3 + 2MASS, written by
referee_star), cleaned (RUWE < 1.15, ipd_frac_multi_peak = 0, non_single_star = 0, 2MASS AAA); XP synthetic photometry
(gaiadr3.synthetic_photometry_gspc) fetched for these source_ids (read-only ADQL, chunks of 1000).
Output: data/mdw_gspc.csv"""
import os
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'): os.environ[k]='1'
import pandas as pd, numpy as np
from astroquery.gaia import Gaia
Gaia.ROW_LIMIT=-1
B='tic214311908/analysis'
D=f'{B}/referee3_followup/data'
m=pd.read_csv(f'{B}/referee_star/data/mdw_cmd.csv')
m=m[(m.ruwe<1.15)&(m.ipd_frac_multi_peak==0)&(m.non_single_star==0)&(m.ph_qual=='AAA')&(m.parallax/m.parallax_error>20)]
m['MKs']=m.ks_m+5*np.log10(m.parallax/100.)
m=m[(m.MKs>5.0)&(m.MKs<9.5)]
print('clean sample', len(m))
out=[]
ids=m.source_id.values
for i in range(0,len(ids),1000):
    q="SELECT source_id, g_sdss_mag, r_sdss_mag, i_sdss_mag, z_sdss_mag, v_jkc_mag, i_jkc_mag, r_jkc_mag, y_ps1_mag FROM gaiadr3.synthetic_photometry_gspc WHERE source_id IN (%s)"%','.join(str(x) for x in ids[i:i+1000])
    out.append(Gaia.launch_job(q).get_results().to_pandas())
s=pd.concat(out); r=m.merge(s,on='source_id'); r.to_csv(f'{D}/mdw_gspc.csv',index=False); print(len(r))
