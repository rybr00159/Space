import numpy as np, glob, os
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
import warnings; warnings.filterwarnings('ignore')
V3 = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v3'
g = Table.read('/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v4/gaia/wide_3p6deg_G17.fits').to_pandas()
EP = 2012.64
for f in sorted(glob.glob(V3 + '/imaging/big/*.fz')) + sorted(glob.glob(V3 + '/imaging/vhs_ADP.2014-12-22T14:5[67]*.fz')):
    h = fits.open(f)
    for ext in range(1, len(h)):
        d = h[ext].data.astype(float); hd = h[ext].header
        w = WCS(hd)
        ra = g.ra + np.nan_to_num(g.pmra) * (EP - 2016) / 3.6e6 / np.cos(np.radians(g.dec)); de = g.dec + np.nan_to_num(g.pmdec) * (EP - 2016) / 3.6e6
        x, y = w.all_world2pix(ra, de, 0)
        on = (x > 20) & (x < d.shape[1] - 20) & (y > 20) & (y < d.shape[0] - 20)
        sky = np.median(d)
        rows = []
        for i in np.where(on & (g.phot_g_mean_mag < 14.5))[0]:
            ix, iy = int(round(x[i])), int(round(y[i]))
            pk = d[iy - 3:iy + 4, ix - 3:ix + 4].max() - sky
            dn = np.hypot((g.ra - g.ra[i]) * np.cos(np.radians(g.dec[i])), g.dec - g.dec[i]) * 3600
            dn = np.sort(dn)[1]
            rows.append((g.source_id[i], g.phot_g_mean_mag[i], g.bp_rp[i], pk, dn, ix, iy))
        rows.sort(key=lambda r: -r[3])
        print(os.path.basename(f), ext, d.shape, hd.get('PSF_FWHM'), 'sky', round(sky), 'n', len(rows))
        for r in rows[:14]:
            print('   ', r[0], f'G={r[1]:.2f} C={r[2]:.2f} peak={r[3]:.0f} nn={r[4]:.1f}" pos=({r[5]},{r[6]})', '<-- target' if r[0] == 6571343840265900544 else '')
