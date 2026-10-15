import warnings,sys; warnings.filterwarnings('ignore')
from astroquery.vizier import Vizier
Vizier.TIMEOUT=60
for kw in sys.argv[1:]:
    r=Vizier.find_catalogs(kw)
    for k,v in list(r.items())[:12]: print(kw,'|',k,'|',v.description[:110])
