"""v4: base rate of the target's Gaia diagnostics among field stars of similar brightness around the target
(G 11.3-12.3, within 6 deg) from VizieR I/355/gaiadr3; neighbours within 3" (any G) via CDS XMatch."""
import numpy as np, time
import astropy.units as u
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
from astroquery.xmatch import XMatch
O = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v4/gaia'
v = Vizier(columns=['Source', 'RA_ICRS', 'DE_ICRS', 'Gmag', 'BP-RP', 'Plx', 'RUWE', 'IPDfmp', 'IPDgofha', 'IPDfow', 'o_BPmag', 'NBPblend',
                    'o_RPmag', 'NRPblend', 'E(BP/RP)', 'epsi', 'sepsi', 'Dup'], column_filters={'Gmag': '11.3..12.3'}, row_limit=-1, timeout=600)
c = SkyCoord(330.887948, -42.835511, unit='deg')
r1 = v.query_region(c, radius=6 * u.deg, catalog='I/355/gaiadr3')[0]
print('stars', len(r1), flush=True)
sub = Table([r1['Source'], r1['RA_ICRS'], r1['DE_ICRS']], names=['sid', 'ra', 'dec'])
outs = []
for i0 in range(0, len(sub), 4000):
    for k in range(4):
        try:
            x = XMatch.query(cat1=sub[i0:i0 + 4000], cat2='vizier:I/355/gaiadr3', max_distance=3 * u.arcsec, colRA1='ra', colDec1='dec')
            break
        except Exception as e:
            print('retry', e); time.sleep(10)
    outs.append(x)
X = vstack(outs)
n3 = {}
for s, src in zip(X['sid'], X['Source']):
    if int(src) != int(s):
        n3[int(s)] = n3.get(int(s), 0) + 1
r1['n3'] = [n3.get(int(s), 0) for s in r1['Source']]
r1.write(f'{O}/field_G11.3-12.3_6deg.ecsv', overwrite=True, format='ascii.ecsv')
print('done', len(r1), 'with neighbours <3":', int((np.array(r1['n3']) > 0).sum()))
