"""h10a: metallicity distribution of nearby dwarfs as a function of space motion (prior for [Fe/H] given the M3V's
kinematics). Gaia DR3 GSP-Spec dwarfs: parallax > 5 mas, parallax_over_error > 20, 4800 < teff_gspspec < 6300,
4.1 < logg_gspspec < 4.75, (mh upper - lower) < 0.25, first 13 flags_gspspec characters all '0' (best quality),
radial_velocity available. [Fe/H] = mh_gspspec + fem_gspspec (fem if available, else mh).
Output: data/gspspec_dwarfs.csv"""
import numpy as np, pandas as pd, warnings, time
warnings.filterwarnings('ignore')
from astroquery.gaia import Gaia
B = 'tic214311908/analysis/host'
q = """select g.source_id, g.ra, g.dec, g.parallax, g.parallax_error, g.pmra, g.pmdec, g.radial_velocity, g.radial_velocity_error,
 g.phot_g_mean_mag, g.bp_rp, g.ruwe, a.teff_gspspec, a.logg_gspspec, a.mh_gspspec, a.mh_gspspec_lower, a.mh_gspspec_upper,
 a.fem_gspspec, a.alphafe_gspspec, a.flags_gspspec
from gaiadr3.gaia_source as g join gaiadr3.astrophysical_parameters as a on a.source_id = g.source_id
where g.parallax > 5 and g.parallax_over_error > 20 and a.teff_gspspec between 4800 and 6300 and a.logg_gspspec between 4.1 and 4.75
 and a.mh_gspspec_upper - a.mh_gspspec_lower < 0.25 and g.radial_velocity is not null"""
for attempt in range(4):
    try:
        t = Gaia.launch_job_async(q).get_results().to_pandas(); break
    except Exception as e:
        print('retry', str(e)[:100]); time.sleep(20)
t['flag13'] = t.flags_gspspec.astype(str).str[:13]
print(len(t), (t.flag13 == '0000000000000').sum())
t.to_csv(f'{B}/data/gspspec_dwarfs.csv', index=False)
