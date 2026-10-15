import warnings, sys; warnings.filterwarnings('ignore')
from astroquery.vizier import Vizier
import astropy.units as u, numpy as np
from astropy.coordinates import SkyCoord
c=SkyCoord(330.887948,-42.835511,unit='deg')
Vizier.TIMEOUT=120
v=Vizier(columns=['**','+_r'],row_limit=-1)
cats=['J/A+A/659/A95/sos','J/A+A/684/A29/main','J/A+A/684/A29/extend','J/A+A/685/A59/catalog','J/ApJ/984/58/table1','J/MNRAS/524/1855/xpparams',
'J/A+A/692/A115/pristine_gaia','J/A+A/695/A75/catalog','J/A+A/711/A54/table1','J/MNRAS/522/29/table2','J/ApJS/279/50/table2','J/MNRAS/520/5635/sources',
'J/MNRAS/471/770/table1','J/AJ/154/259/table6','V/160/shboost','I/354/starhorse2021','I/355/paramp','I/355/paramsup','I/363/gaia_rv_dr3','J/MNRAS/529/1802/targets',
'J/A+A/677/A137/catalog','J/A+A/692/A150/allskytg','II/374/table2','J/ApJ/925/164/catalog','V/136/tycall','J/A+A/653/A98/aspic1_1','IV/39/tic82']
out=open(sys.argv[1],'w')
for cat in cats:
    try:
        res=v.query_region(c,radius=3*u.arcsec,catalog=cat)
        for t in res:
            r=t[0]
            s=' | '.join(f'{k}={r[k]}' for k in t.colnames if not np.ma.is_masked(r[k]))
            out.write(f"===== {cat} {t.meta.get('description','')[:120]}\n{s}\n")
    except Exception as e: out.write(f"===== {cat} ERR {str(e)[:150]}\n")
out.close()
