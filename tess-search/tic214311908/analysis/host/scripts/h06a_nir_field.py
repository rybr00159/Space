"""h06a: local photometric calibration checks around the M3V (radius 8 arcmin):
DENIS (B/denis) vs 2MASS PSC (II/246) for the same stars (zero point of the DENIS frame that measured the M3V,
image 1724791 / strip 9861, and of the whole neighbourhood); unWISE (II/363) vs AllWISE (II/328); VHS DR5 (II/367)
vs 2MASS for unsaturated stars; Gaia DR3 for colours. Output: data/local_{denis,2mass,unwise,allwise,vhs,gaia}.csv"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import astropy.units as u
B = 'tic214311908/analysis/host'
c = SkyCoord(260.969674, -51.323821, unit='deg')
v = Vizier(row_limit=-1, columns=['**', '_RAJ2000', '_DEJ2000', '_r'], timeout=600)
for lab, cat in [('denis', 'B/denis/denis'), ('2mass', 'II/246/out'), ('unwise', 'II/363/unwise'), ('allwise', 'II/328/allwise'),
                 ('vhs', 'II/367/vhs_dr5'), ('gaia', 'I/355/gaiadr3')]:
    r = v.query_region(c, radius=8 * u.arcmin, catalog=cat)
    t = r[0].to_pandas()
    t.to_csv(f'{B}/data/local_{lab}.csv', index=False)
    print(lab, len(t))
