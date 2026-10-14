"""v4: Gaia DR3 G < 17 within 3.6 deg of the target (for the FFI lateral-colour calibration), tiled TAP queries."""
import pyvo, time, numpy as np
from astropy.table import vstack
from concurrent.futures import ThreadPoolExecutor
O = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v4/gaia'
svc = pyvo.dal.TAPService('https://gea.esac.esa.int/tap-server/tap')
cols = "source_id,ra,dec,pmra,pmdec,parallax,parallax_error,phot_g_mean_mag,phot_bp_mean_mag,phot_rp_mean_mag,bp_rp,ruwe,ipd_frac_multi_peak,phot_bp_rp_excess_factor"
ra0, de0, R, GL = 330.887948, -42.835511, 3.6, 17.0
tiles = []
for d1 in np.arange(de0 - R, de0 + R, 0.6):
    d2 = d1 + 0.6
    cm = np.cos(np.radians(max(abs(d1), abs(d2))))
    for a1 in np.arange(ra0 - R / cm, ra0 + R / cm, 1.2):
        tiles.append((a1, a1 + 1.2, d1, d2))
def q(t):
    a1, a2, d1, d2 = t
    s = f"SELECT {cols} FROM gaiadr3.gaia_source WHERE ra BETWEEN {a1} AND {a2} AND dec BETWEEN {d1} AND {d2} AND phot_g_mean_mag < {GL}"
    for k in range(4):
        try:
            return svc.run_sync(s, maxrec=500000).to_table()
        except Exception as e:
            time.sleep(5)
    print('FAILED', t); return None
t0 = time.time()
with ThreadPoolExecutor(6) as ex:
    res = list(ex.map(q, tiles))
T = vstack([r for r in res if r is not None])
sep = np.degrees(np.arccos(np.clip(np.sin(np.radians(de0)) * np.sin(np.radians(T['dec'])) + np.cos(np.radians(de0)) * np.cos(np.radians(T['dec'])) * np.cos(np.radians(T['ra'] - ra0)), -1, 1)))
T = T[sep < R]
_, iu = np.unique(np.array(T['source_id']), return_index=True); T = T[iu]
T.write(f'{O}/wide_3p6deg_G17.fits', overwrite=True)
print(len(tiles), 'tiles', len(T), 'rows', f'{time.time() - t0:.0f}s', sum(r is None for r in res), 'failed')
