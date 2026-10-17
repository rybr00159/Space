"""h02: list every VizieR catalogue with an entry within 3" of the M3V (Gaia DR3 2016 position and the J2000 position
propagated from it), to find literature metallicities / XP-based parameters. Output: data/vizier_allcats.txt,
data/vizier_allcats_rows.csv (all rows, long format)"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import astropy.units as u
B = 'tic214311908/analysis/host'
ra, dec = 260.969674, -51.323821
pmra, pmdec = -113.671152, 85.54355
cosd = np.cos(np.radians(dec))
c16 = SkyCoord(ra, dec, unit='deg')
c00 = SkyCoord(ra - pmra * 16 / 3.6e6 / cosd, dec - pmdec * 16 / 3.6e6, unit='deg')
v = Vizier(row_limit=50, columns=['*', '_r'], timeout=600)
lines, rows = [], []
seen = set()
for tag, c in [('2016', c16), ('J2000', c00)]:
    res = v.query_region(c, radius=3 * u.arcsec)
    for k in res.keys():
        t = res[k]
        if k in seen:
            continue
        seen.add(k)
        lines.append(f'{k}\t{len(t)}\t{t.meta.get("description", "")}')
        df = t.to_pandas()
        for i, r in df.iterrows():
            for col, val in r.items():
                rows.append(dict(cat=k, row=i, col=col, val=str(val)))
open(f'{B}/data/vizier_allcats.txt', 'w').write('\n'.join(lines))
pd.DataFrame(rows).to_csv(f'{B}/data/vizier_allcats_rows.csv', index=False)
print('\n'.join(lines))
