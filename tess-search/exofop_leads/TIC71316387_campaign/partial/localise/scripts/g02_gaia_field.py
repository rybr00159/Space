import pyvo
L='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise'
svc = pyvo.dal.TAPService('https://gea.esac.esa.int/tap-server/tap')
q = """SELECT source_id, ra, dec, ref_epoch, parallax, pmra, pmdec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, bp_rp, ruwe,
ipd_frac_multi_peak, DISTANCE(POINT(330.887948,-42.835511), POINT(ra,dec))*3600 AS sep
FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT(ra,dec), CIRCLE(330.887948,-42.835511, 0.105)) AND phot_g_mean_mag < 19.5 ORDER BY sep"""
r = svc.run_sync(q, maxrec=100000).to_table()
r.write(f'{L}/gaia/field_378arcsec_G19.5.csv', overwrite=True, format='csv')
print(len(r)); print(r[:12])
