import warnings,sys; warnings.filterwarnings('ignore')
from astroquery.vizier import Vizier
Vizier.TIMEOUT=90
for c in sys.argv[1:]:
    try:
        v=Vizier(row_limit=2,columns=['**'])
        r=v.get_catalogs(c)
        for k in r.keys(): print(c,'|',k,'|',len(r[k]),'|',r[k].meta.get('description','')[:90],'|',','.join(r[k].colnames[:30]))
    except Exception as e: print(c,'ERR',str(e)[:150])
