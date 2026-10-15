import warnings, sys; warnings.filterwarnings('ignore')
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord
c=SkyCoord(330.887948,-42.835511,unit='deg')
Vizier.TIMEOUT=200
v=Vizier(columns=['*','+_r'],row_limit=5)
res=v.query_region(c,radius=3*u.arcsec)
out=open(sys.argv[1],'w')
for k in res.keys():
    t=res[k]
    out.write(f"{k} | {t.meta.get('description','')[:100]} | n={len(t)} | cols={','.join(t.colnames[:40])}\n")
print(len(res.keys()))
