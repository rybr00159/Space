"""Gaia DR3 sources within 330" of the target (covers the 21x21 TESScut cutouts), G < 19, via VizieR (I/355/gaiadr3)."""
import numpy as np, pandas as pd, sys
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord
Vizier.ROW_LIMIT = -1
v = Vizier(columns=["Source", "RA_ICRS", "DE_ICRS", "pmRA", "pmDE", "Plx", "Gmag", "BPmag", "RPmag", "RUWE", "Epoch"],
           column_filters={"Gmag": "<19"}, timeout=300, row_limit=-1)
c = SkyCoord(272.60647748, -43.99546646, unit="deg")
for a in range(4):
    try:
        r = v.query_region(c, radius=330 * u.arcsec, catalog="I/355/gaiadr3"); break
    except Exception as e:
        print("retry", e)
t = r[0].to_pandas()
t.to_csv(sys.argv[1], index=False)
print(len(t), t.Gmag.min())
