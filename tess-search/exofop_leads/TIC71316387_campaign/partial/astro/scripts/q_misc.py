import warnings; warnings.filterwarnings('ignore')
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord
c=SkyCoord(330.887948,-42.835511,unit='deg')
Vizier.TIMEOUT=60
v=Vizier(columns=['*','+_r'],row_limit=-1)
for cat in ['B/wds/wds','B/sb9/main','I/358/vclassre','I/358/veb','J/MNRAS/522/29/table2']:
    try:
        r=v.query_region(c,radius=30*u.arcsec,catalog=cat)
        print(cat,'->',[(k,len(r[k])) for k in r.keys()] if len(r) else 'NONE within 30"',flush=True)
    except Exception as e: print(cat,'ERR',str(e)[:100],flush=True)
import pyvo
try:
    s=pyvo.dal.TAPService('http://archive.eso.org/tap_obs')
    q="SELECT dp_id,instrument_name,dataproduct_type,obs_release_date,s_ra,s_dec,t_min,exptime FROM ivoa.ObsCore WHERE INTERSECTS(s_region,CIRCLE('ICRS',330.887948,-42.835511,10.0/3600))=1"
    t=s.search(q).to_table(); print('ESO obscore',len(t)); print(t[:30])
except Exception as e: print('ESO ERR',str(e)[:200])
