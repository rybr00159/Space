from astroquery.vizier import Vizier
import warnings; warnings.filterwarnings('ignore')
import sys
for cid in sys.argv[1:]:
    try:
        v = Vizier(row_limit=3, columns=['**'])
        res = v.get_catalogs(cid)
    except Exception as e:
        print(cid, 'ERR', e); continue
    print('=====', cid, len(res))
    for t in res:
        print('--', t.meta.get('name'), len(t), '|', t.meta.get('description', '')[:100])
        print('   cols:', ', '.join(t.colnames)[:600])
