"""TIC-8 entries of the candidate stars (Tmag, 2MASS, Teff, radius, mass, contamination ratio)."""
import sys; sys.path.insert(0, 'tic214311908/analysis/referee_loc/scripts')
from rl_lib import *
from astroquery.mast import Catalogs
t = Catalogs.query_region('260.97765 -51.32369', radius=0.012, catalog='TIC').to_pandas()
cols = ['ID', 'GAIA', 'ra', 'dec', 'Tmag', 'e_Tmag', 'Jmag', 'Hmag', 'Kmag', 'GAIAmag', 'gaiabp', 'gaiarp', 'Teff', 'rad', 'mass',
        'contratio', 'dstArcSec', 'objType', 'lumclass', 'd', 'TESSflag', 'disposition', 'duplicate_id']
t = t[[c for c in cols if c in t.columns]].sort_values('dstArcSec')
t.to_csv(f'{ROOT}/data/tic_near.csv', index=False)
print(t[t.Tmag < 18].head(25).to_string())
