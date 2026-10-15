import warnings; warnings.filterwarnings('ignore')
import pyvo
s=pyvo.dal.TAPService('http://archive.eso.org/tap_obs')
q="SELECT dp_id,instrument_name,dataproduct_type,obs_release_date,s_ra,s_dec,t_min,t_exptime,target_name FROM ivoa.ObsCore WHERE INTERSECTS(s_region,CIRCLE('ICRS',330.887948,-42.835511,10.0/3600))=1"
t=s.search(q).to_table(); print('ESO obscore phase3 products covering position:',len(t))
import collections; print(collections.Counter([(r['instrument_name'],r['dataproduct_type']) for r in t]))
sp=t[[str(x)=='spectrum' for x in t['dataproduct_type']]]; print(sp)
try:
    q2="SELECT dp_id,instrument,dp_cat,dp_type,ra,dec,target,exposure,date_obs FROM dbo.raw WHERE ra BETWEEN 330.885 AND 330.891 AND dec BETWEEN -42.838 AND -42.833"
    t2=s.search(q2).to_table(); print('ESO raw frames within ~10":',len(t2)); print(t2[:20])
except Exception as e: print('raw ERR',str(e)[:200])
