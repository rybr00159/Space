"""Own VizieR query of Gaia DR3 (I/355/gaiadr3), G < 18, within 330" of the TIC position (covers the 21x21 stamp)."""
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord
v = Vizier(columns=['Source', 'RA_ICRS', 'DE_ICRS', 'pmRA', 'pmDE', 'Gmag', 'BP-RP', 'Plx', 'RUWE', 'Epoch'],
           column_filters={'Gmag': '<18'}, row_limit=-1, timeout=300)
r = v.query_region(SkyCoord(272.606557, -43.995371, unit='deg'), radius=330 * u.arcsec, catalog='I/355/gaiadr3')
t = r[0]
t.write('/home/user/work/tic87/referee_audit/gaia/gaia_dr3_330as_G18_vizier.csv', format='csv', overwrite=True)
print(len(t))
