import pyvo, time, sys
R = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v2'
svc = pyvo.dal.TAPService('https://gea.esac.esa.int/tap-server/tap')
which = sys.argv[1]
if which == 'col':
    q = """SELECT source_id, ra, dec, pmra, pmdec, parallax, phot_g_mean_mag, bp_rp, ruwe, ipd_frac_multi_peak, teff_gspphot
    FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT(ra,dec), CIRCLE(330.887948,-42.835511, 9.0)) AND phot_g_mean_mag < 15.0"""
    out = f'{R}/gaia/cone9deg_G15.fits'
else:
    q = """SELECT source_id, ra, dec, phot_g_mean_mag FROM gaiadr3.gaia_source
    WHERE 1=CONTAINS(POINT(ra,dec), CIRCLE(330.887948,-42.835511, 9.0)) AND phot_g_mean_mag < 16.5"""
    out = f'{R}/gaia/cone9deg_G16.5.fits'
t0 = time.time()
job = svc.submit_job(q, maxrec=3000000); job.run()
job.wait(phases=['COMPLETED', 'ERROR', 'ABORTED'], timeout=1500)
print(which, job.phase, flush=True)
r = job.fetch_result().to_table()
print(len(r), time.time() - t0)
r.write(out, overwrite=True)
