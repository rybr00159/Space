"""Search for further co-moving companions of the TIC 214311908 system: Gaia DR3 sources within 2 deg (7,800 au
projected at 64 pc... up to ~4.5 pc) with parallax 13-18.5 mas; candidates have tangential velocity differences
< 3 km/s from the M3V (or PM difference < 10 mas/yr for 2-parameter sources). Output: out/comovers.json"""
import json, numpy as np, pandas as pd
from astroquery.gaia import Gaia
Gaia.ROW_LIMIT = -1
B = 'tic214311908/analysis/star'
q = ("select source_id, ra, dec, parallax, parallax_error, pmra, pmdec, pmra_error, pmdec_error, phot_g_mean_mag, bp_rp, ruwe, radial_velocity, radial_velocity_error, "
     "DISTANCE(POINT('ICRS', ra, dec), POINT('ICRS', 260.969674, -51.323821))*3600 as sep "
     "from gaiadr3.gaia_source where 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', 260.969674, -51.323821, 2.0)) and parallax between 13 and 18.5")
t = Gaia.launch_job_async(q).get_results().to_pandas()
pm0 = np.array([-113.671152, 85.54355]); plx0 = 15.5494
dpm = np.hypot(t.pmra - pm0[0], t.pmdec - pm0[1])
t['dv_tan_kms'] = 4.74 * dpm / plx0
t['sep_au'] = t.sep / plx0 * 1000 / 1000 * 1  # arcsec / (plx in arcsec) -> au: sep[as] * 1000/plx[mas]
t['sep_au'] = t.sep * 1000 / plx0
c = t[t.dv_tan_kms < 3].sort_values('sep')
print(len(t), 'sources with plx 13-18.5 mas within 2 deg;', len(c), 'with dv_tan < 3 km/s')
print(c[['source_id', 'sep', 'sep_au', 'parallax', 'parallax_error', 'pmra', 'pmdec', 'dv_tan_kms', 'phot_g_mean_mag', 'bp_rp', 'ruwe', 'radial_velocity']].round(3).to_string())
json.dump(dict(n_plx_window=int(len(t)), candidates=c.round(4).to_dict('records')), open(f'{B}/out/comovers.json', 'w'), indent=1, default=str)
