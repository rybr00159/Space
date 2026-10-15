import warnings; warnings.filterwarnings('ignore')
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord
c=SkyCoord(330.887948,-42.835511,unit='deg')
Vizier.TIMEOUT=90
v=Vizier(columns=['*','+_r'],row_limit=-1)
for cat in ['III/283/ravedr6','III/279/rave_dr5','J/MNRAS/506/150/stars','III/286/catalog','III/286/allvis','J/A+A/666/A121','J/A+A/659/A95/sos','J/ApJS/258/16/tess-ebs','J/A+A/674/A16','I/358/vclassre','I/358/veb','J/AJ/161/154','J/MNRAS/522/29/table2','B/wds/wds','J/AJ/155/42','B/sb9/main','J/A+A/664/A4']:
    try:
        r=v.query_region(c,radius=30*u.arcsec,catalog=cat)
        print(cat,'->',[ (k,len(r[k])) for k in r.keys()] if len(r) else 'NONE within 30"')
        for k in r.keys(): print(r[k][:3])
    except Exception as e: print(cat,'ERR',str(e)[:100])
