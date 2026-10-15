import warnings; warnings.filterwarnings('ignore')
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord
c=SkyCoord(330.887948,-42.835511,unit='deg')
Vizier.TIMEOUT=90
v=Vizier(columns=['*','+_r'],row_limit=-1)
for cat in ['J/other/PASA/42.51','J/MNRAS/478/4513','V/153','J/MNRAS/506/150']:
    try:
        r=v.query_region(c,radius=30*u.arcsec,catalog=cat)
        print(cat,'->',[(k,len(r[k])) for k in r.keys()] if len(r) else 'NONE within 30"')
        for k in r.keys(): print(r[k][:3])
    except Exception as e: print(cat,'ERR',str(e)[:150])
# list tables of GALAH DR4 to confirm the query hit the star table
try:
    cats=Vizier(row_limit=1).get_catalogs('J/other/PASA/42.51')
    for k in cats.keys(): print('table',k,cats[k].meta.get('description','')[:80],len(cats[k].colnames))
except Exception as e: print('meta ERR',str(e)[:100])
