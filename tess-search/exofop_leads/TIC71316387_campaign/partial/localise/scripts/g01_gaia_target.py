import pyvo, json, sys
import numpy as np
L='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise'
svc = pyvo.dal.TAPService('https://gea.esac.esa.int/tap-server/tap')
q = """SELECT * FROM gaiadr3.gaia_source WHERE source_id=6571343840265900544"""
r = svc.run_sync(q, maxrec=10).to_table()
row = {c: (r[c][0].item() if hasattr(r[c][0],'item') else str(r[c][0])) for c in r.colnames}
json.dump(row, open(f'{L}/gaia/target_gaia_source.json','w'), indent=1, default=str)
for c in r.colnames:
    print(c, row[c])
q2 = """SELECT source_id, ra, dec, parallax, pmra, pmdec, phot_g_mean_mag, bp_rp, ruwe, ipd_frac_multi_peak, DISTANCE(POINT(330.887948,-42.835511), POINT(ra,dec))*3600 AS sep
FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT(ra,dec), CIRCLE(330.887948,-42.835511, 0.0333)) ORDER BY sep"""
r2 = svc.run_sync(q2, maxrec=1000).to_table()
r2.write(f'{L}/gaia/neighbours_120arcsec.csv', overwrite=True, format='csv')
print(r2[:20])
