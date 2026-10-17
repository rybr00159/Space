"""More VizieR look-ups (correct table names): El-Badry+21 wide binaries, Zhang+23 XP parameters, eRASS1 main/soft,
Andrae+23 XP metallicity (if on VizieR), plus the full eRASS1 record. Positions as in s02_phot.py.
Outputs: data/vizier2_<name>.csv"""
import warnings, time
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord
D = 'tic214311908/analysis/star/data'
c0 = SkyCoord(260.969674, -51.323821, unit='deg')
v = Vizier(columns=['**', '_RAJ2000', '_DEJ2000', '_r'], row_limit=-1, timeout=120)
for name, cid, rad in [('elbadry21', 'J/MNRAS/506/2269/catalog', 30), ('zhang23', 'J/MNRAS/524/1855/xpparams', 20),
                       ('erass1_main', 'J/A+A/682/A34/erass1-m', 120), ('erass1_soft', 'J/A+A/682/A34/erass1-s', 120),
                       ('erass1_hard', 'J/A+A/682/A34/erass1-h', 120)]:
    try:
        r = v.query_region(c0, radius=rad * u.arcsec, catalog=cid)
        t = r[0].to_pandas() if len(r) else pd.DataFrame()
    except Exception as e:
        print(name, 'error', e); continue
    t.to_csv(f'{D}/vizier2_{name}.csv', index=False)
    print('=====', name, cid, len(t))
    if len(t):
        print(t.T.to_string()[:6000])
