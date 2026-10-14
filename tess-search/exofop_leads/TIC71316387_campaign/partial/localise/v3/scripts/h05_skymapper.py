"""SkyMapper DR2 Shallow Survey (5-20 s, unsaturated for the target) companion search with the h01 machinery:
ePSF from unsaturated isolated Gaia stars in the 6' cutout, scan 0.6-4", injection-recovery (dm 1-4.5, sep 0.8-3")."""
import sys, os, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import h01_vhs as h
from v3lib import V3, to_json
h.SATLIM = 60000.0; h.GAIN = 1.0; h.RN = 8.0; h.EP = 2015.0
from astropy.io import fits
from astropy.wcs import WCS
def load(f):
    hd = fits.open(f)[0]; d = hd.data.astype(float); w = WCS(hd.header)
    g = h.g
    ra = g.ra + np.nan_to_num(g.pmra) * (h.EP - 2016) / 3.6e6 / np.cos(np.radians(g.dec)); de = g.dec + np.nan_to_num(g.pmdec) * (h.EP - 2016) / 3.6e6
    x, y = w.all_world2pix(ra, de, 0)
    r0 = w.all_pix2world([[x[0], y[0]], [x[0] + 1, y[0]], [x[0], y[0] + 1]], 0); c = np.cos(np.radians(r0[0, 1]))
    J = np.array([[(r0[1, 0] - r0[0, 0]) * c * 3600, (r0[2, 0] - r0[0, 0]) * c * 3600], [(r0[1, 1] - r0[0, 1]) * 3600, (r0[2, 1] - r0[0, 1]) * 3600]])
    return dict(d=d + 200.0, x=x, y=y, J=J, hdr={'PSF_FWHM': np.nan}, band='')
h.load = load
out = []
for f in sorted(glob.glob(V3 + '/imaging/skymapper/sm_*.fits')):
    try:
        out.append(h.run(f))
    except Exception as e:
        print('fail', f, e)
to_json(out, V3 + '/results/h05_skymapper.json')
