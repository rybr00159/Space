import warnings; warnings.filterwarnings('ignore')
from astroquery.vizier import Vizier
import astropy.units as u, numpy as np
from astropy.coordinates import SkyCoord
c=SkyCoord(330.887948,-42.835511,unit='deg')
Vizier.TIMEOUT=90
v=Vizier(columns=['**','+_r'],row_limit=-1)
for cat in ['I/320/spm4','I/331/apop','I/319/xpm','I/259/tyc2','I/337/tgas','I/340/ucac5','I/339/hsoy','I/322A/out','I/317/sample','I/312/sample','I/350/gaiaedr3','I/345/gaia2']:
    try:
        r=v.query_region(c,radius=3*u.arcsec,catalog=cat)
        for t in r:
            row=t[0]
            keys=[k for k in t.colnames if any(s in k.lower() for s in ['pm','ep','e_ra','e_de','ra','de','plx','_r','mag'])]
            print('=====',cat,' | '.join(f'{k}={row[k]}' for k in keys if not np.ma.is_masked(row[k]))[:1500])
    except Exception as e: print(cat,'ERR',str(e)[:100])
