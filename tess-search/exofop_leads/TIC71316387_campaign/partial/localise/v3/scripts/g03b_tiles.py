import pyvo, time, numpy as np, sys
from astropy.table import vstack, Table
from concurrent.futures import ThreadPoolExecutor
svc = pyvo.dal.TAPService('https://gea.esac.esa.int/tap-server/tap')
cols = "source_id,ra,dec,pmra,pmdec,parallax,phot_g_mean_mag,bp_rp,ruwe"
ra0, de0, R = 330.887948, -42.835511, 2.7
ddec = np.arange(-R, R, 0.45)
tiles = []
for d in ddec:
    d1, d2 = de0 + d, de0 + d + 0.45
    cmin = np.cos(np.radians(min(abs(d1), abs(d2)) if d1*d2>0 else 0)); cmin = np.cos(np.radians(max(abs(d1), abs(d2))))
    wra = R / cmin
    for a in np.arange(-wra, wra, 0.9):
        tiles.append((ra0 + a, ra0 + a + 0.9, d1, d2))
def q(t):
    a1, a2, d1, d2 = t
    s = f"SELECT {cols} FROM gaiadr3.gaia_source WHERE ra BETWEEN {a1} AND {a2} AND dec BETWEEN {d1} AND {d2} AND phot_g_mean_mag < 18.5"
    for k in range(3):
        try:
            return svc.run_sync(s, maxrec=200000).to_table()
        except Exception as e:
            time.sleep(3)
    print('FAILED', t); return None
t0 = time.time()
with ThreadPoolExecutor(6) as ex:
    res = list(ex.map(q, tiles))
T = vstack([r for r in res if r is not None])
c = np.cos(np.radians(de0))
sep = np.hypot((T['ra'] - ra0) * c, T['dec'] - de0)
T = T[sep < R]
T.write('wide_2p7deg_G18p5.fits', overwrite=True)
print(len(tiles), 'tiles', len(T), 'rows', time.time() - t0, 's', sum(r is None for r in res), 'failed')
