import pyvo, time
R = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v2'
svc = pyvo.dal.TAPService('https://gea.esac.esa.int/tap-server/tap')
q = """SELECT a.source_id, b.source_id AS nb_id, DISTANCE(POINT(a.ra,a.dec), POINT(b.ra,b.dec))*3600 AS sep,
b.phot_g_mean_mag AS nb_g, b.astrometric_params_solved AS nb_aps
FROM gaiadr3.gaia_source AS a JOIN gaiadr3.gaia_source AS b
ON 1=CONTAINS(POINT(b.ra,b.dec), CIRCLE(a.ra,a.dec, 8.0/3600.))
WHERE 1=CONTAINS(POINT(a.ra,a.dec), CIRCLE(330.887948,-42.835511, 3.0)) AND a.phot_g_mean_mag BETWEEN 10.75 AND 12.75
AND a.source_id != b.source_id"""
t0 = time.time()
job = svc.submit_job(q, maxrec=200000); job.run()
job.wait(phases=['COMPLETED', 'ERROR', 'ABORTED'], timeout=600)
r = job.fetch_result().to_table()
print(len(r), time.time() - t0)
r.write(f'{R}/gaia/field3deg_neighbours8as.csv', overwrite=True, format='csv')
