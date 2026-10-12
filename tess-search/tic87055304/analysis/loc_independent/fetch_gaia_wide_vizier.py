"""Fallback: Gaia DR3 (I/355/gaiadr3) G<17.5 within 330" from VizieR."""
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
v = Vizier(columns=['Source', 'RA_ICRS', 'DE_ICRS', 'e_RA_ICRS', 'e_DE_ICRS', 'pmRA', 'pmDE', 'Plx', 'Gmag', 'BPmag', 'RPmag', 'BP-RP', 'RUWE'],
           column_filters={'Gmag': '<17.5'}, row_limit=-1, timeout=300)
r = v.query_region(SkyCoord(272.606557, -43.995371, unit='deg'), radius=330 * u.arcsec, catalog='I/355/gaiadr3')[0]
import pandas as pd
d = pd.DataFrame(dict(source_id=np.array(r['Source']), ra=np.array(r['RA_ICRS']), dec=np.array(r['DE_ICRS']),
                      pmra=np.array(r['pmRA'].filled(np.nan)), pmdec=np.array(r['pmDE'].filled(np.nan)),
                      parallax=np.array(r['Plx'].filled(np.nan)), ref_epoch=2016.0,
                      phot_g_mean_mag=np.array(r['Gmag'].filled(np.nan)), phot_bp_mean_mag=np.array(r['BPmag'].filled(np.nan)),
                      phot_rp_mean_mag=np.array(r['RPmag'].filled(np.nan)), bp_rp=np.array(r['BP-RP'].filled(np.nan)),
                      ruwe=np.array(r['RUWE'].filled(np.nan))))
c = SkyCoord(d.ra.values, d.dec.values, unit='deg')
d['sep'] = c.separation(SkyCoord(272.606557, -43.995371, unit='deg')).arcsec
d.to_csv('/home/user/work/tic87/loc_independent/gaia_wide/gaia_dr3_330arcsec_G17p5_vizier.csv', index=False)
print(len(d))
