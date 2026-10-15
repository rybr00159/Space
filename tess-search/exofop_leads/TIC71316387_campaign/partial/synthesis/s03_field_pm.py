# Synthesis s03: proper motions / parallaxes of field M dwarfs that could act as an unrelated EB host
# (BP-RP 1.8-2.7, T in the light-fraction window), to quantify the power of a common-proper-motion test.
import numpy as np, json, sys
from astroquery.gaia import Gaia
Gaia.ROW_LIMIT = -1
q = """SELECT source_id, ra, dec, phot_g_mean_mag, bp_rp, parallax, parallax_error, pmra, pmdec, radial_velocity
FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec), CIRCLE('ICRS',330.887948,-42.835511,1.5))
AND bp_rp BETWEEN 1.8 AND 2.8 AND phot_g_mean_mag BETWEEN 14.0 AND 18.5"""
job = Gaia.launch_job_async(q)
t = job.get_results()
t.write('/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/synthesis/field_Mdwarfs_1p5deg.ecsv', overwrite=True)
print(len(t))
