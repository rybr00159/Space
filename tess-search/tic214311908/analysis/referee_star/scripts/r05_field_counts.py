"""Referee: Gaia DR3 source counts within 6 arcmin of the M3V (excluding < 40" to avoid the Mira/M3V surroundings)
for the density of potential blending stars (T from the Stassun+19 Gaia relation; G-0.43 without colour).
Output: data/gaia_cone6arcmin.csv, out/r05_field_counts.json"""
import json, warnings, time; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from astroquery.gaia import Gaia
Gaia.ROW_LIMIT = -1
R = 'tic214311908/analysis/referee_star'
RA, DEC = 260.969674, -51.323821
q = f"""select source_id, ra, dec, phot_g_mean_mag, bp_rp, parallax, parallax_error, distance(point({RA},{DEC}), point(ra,dec))*3600 as sep
from gaiadr3.gaia_source where 1=contains(point(ra,dec), circle({RA},{DEC},0.1))"""
for k in range(4):
    try: t = Gaia.launch_job_async(q).get_results().to_pandas(); break
    except Exception as e: print('retry', e); time.sleep(20)
t.to_csv(f'{R}/data/gaia_cone6arcmin.csv', index=False)
c = t.bp_rp.values
T = np.where(np.isfinite(c), t.phot_g_mean_mag - 0.00522555 * c ** 3 + 0.0891337 * c ** 2 - 0.633923 * c + 0.0324473, t.phot_g_mean_mag - 0.430)
t['T'] = T
ann = (t.sep > 40) & (t.sep < 360)
area = np.pi * (360 ** 2 - 40 ** 2)
out = {}
for Tl in (16.72, 17.0, 17.47, 18.0, 19.0):
    n = int(((t['T'] <= Tl) & ann).sum()); dens = n / area
    out[f'T<={Tl}'] = dict(n=n, per_arcsec2=dens, per_arcsec2_err=np.sqrt(n) / area,
                          expected_within_1=dens * np.pi, expected_within_2=dens * np.pi * 4, expected_within_3=dens * np.pi * 9)
    print(f'T <= {Tl}: n = {n}, density {dens:.5f} /arcsec2 -> expected within 1/2/3": {dens*np.pi:.4f} / {dens*4*np.pi:.4f} / {dens*9*np.pi:.4f}')
# completeness check: G histogram turnover
h, e = np.histogram(t.phot_g_mean_mag[ann], bins=np.arange(12, 22, 0.5)); print('G hist', list(zip(e[:-1], h)))
json.dump(out, open(f'{R}/out/r05_field_counts.json', 'w'), indent=1)
