import pyvo, time
R = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v2'
svc = pyvo.dal.TAPService('https://gea.esac.esa.int/tap-server/tap')
q = """SELECT source_id, ra, dec, phot_g_mean_mag, astrometric_params_solved FROM gaiadr3.gaia_source
WHERE 1=CONTAINS(POINT(ra,dec), CIRCLE(330.887948,-42.835511, 3.01))"""
t0 = time.time()
job = svc.submit_job(q, maxrec=3000000); job.run()
job.wait(phases=['COMPLETED', 'ERROR', 'ABORTED'], timeout=900)
print(job.phase)
r = job.fetch_result().to_table()
print(len(r), time.time() - t0)
r.write(f'{R}/gaia/cone3deg_all.fits', overwrite=True)
