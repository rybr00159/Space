import signal, numpy as np
from astroquery.gaia import Gaia
def handler(s,f): raise TimeoutError
signal.signal(signal.SIGALRM, handler); signal.alarm(300)
q="""SELECT source_id, ra, dec, pmra, pmdec, parallax, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, bp_rp, ruwe,
DISTANCE(POINT('ICRS',ra,dec),POINT('ICRS',330.887948,-42.835511))*3600 AS sep
FROM gaiadr3.gaia_source WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',330.887948,-42.835511,0.0667))=1 AND phot_g_mean_mag<20.5"""
j=Gaia.launch_job(q); t=j.get_results(); t.sort('sep')
t.write('/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/data/gaia_240as.csv',overwrite=True)
print(len(t)); print(t['source_id','sep','phot_g_mean_mag','bp_rp','parallax'][:15])
