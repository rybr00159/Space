import pyvo, time
from astropy.table import vstack
svc = pyvo.dal.TAPService('https://tapvizier.cds.unistra.fr/TAPVizieR/tap')
cols = '"Source" AS source_id, "RA_ICRS" AS ra, "DE_ICRS" AS dec, "pmRA" AS pmra, "pmDE" AS pmdec, "Plx" AS parallax, "Gmag" AS phot_g_mean_mag, "BP-RP" AS bp_rp, "RUWE" AS ruwe'
t0 = time.time()
q = f'SELECT {cols} FROM "I/355/gaiadr3" WHERE 1=CONTAINS(POINT(\'ICRS\',"RA_ICRS","DE_ICRS"),CIRCLE(\'ICRS\',330.887948,-42.835511,2.7)) AND "Gmag" < 18.5'
t = svc.run_sync(q, maxrec=500000).to_table()
t.write('wide_2p7deg_G18p5.fits', overwrite=True)
print(len(t), time.time()-t0)
