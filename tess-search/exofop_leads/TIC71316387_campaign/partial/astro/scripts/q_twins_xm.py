import sys, warnings, time; warnings.filterwarnings('ignore')
from astroquery.gaia import Gaia
from astropy.table import Table, vstack, join
Gaia.ROW_LIMIT=-1
T=Table.read(sys.argv[1]); ids=list(T['source_id'])+[6571343840265900544]
ids=sorted(set(ids)); out=[]
for i in range(0,len(ids),400):
    ch=','.join(str(x) for x in ids[i:i+400])
    q=f"""SELECT xm.source_id, tm.j_m,tm.j_msigcom,tm.h_m,tm.h_msigcom,tm.ks_m,tm.ks_msigcom,tm.ph_qual AS tm_qual, xm.number_of_neighbours AS tm_nn,
    aw.w1mpro,aw.w1mpro_error,aw.w2mpro,aw.w2mpro_error,aw.w3mpro,aw.w3mpro_error,aw.w4mpro,aw.ph_qual AS aw_qual,aw.cc_flags,aw.ext_flag
    FROM gaiadr3.tmass_psc_xsc_best_neighbour AS xm
    JOIN gaiadr1.tmass_original_valid AS tm ON xm.original_ext_source_id=tm.designation
    JOIN gaiadr3.allwise_best_neighbour AS axm ON xm.source_id=axm.source_id
    JOIN gaiadr1.allwise_original_valid AS aw ON axm.allwise_oid=aw.allwise_oid
    WHERE xm.source_id IN ({ch})"""
    for a in range(4):
        try:
            t=Gaia.launch_job(q).get_results(); out.append(t); print(i,len(t),flush=True); break
        except Exception as e: print('retry',str(e)[:80],flush=True); time.sleep(5)
X=vstack(out); X.write(sys.argv[2],overwrite=True); print('total',len(X))
