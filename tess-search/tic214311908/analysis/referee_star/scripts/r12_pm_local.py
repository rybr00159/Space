"""Referee: local calibration of the DR3 - DR2 proper-motion difference.
Two plain cone queries (0.5 deg around the M3V, G 13.3-15.7): Gaia DR3 (5/6-parameter) and Gaia DR2 (5-parameter);
matched by propagating DR3 positions to 2015.5 (< 0.3"). For matched stars, z = (mu_DR3 - mu_DR2)/sqrt(e3^2+e2^2)
per component; robust offset and scatter (MAD) of z for all stars and for red stars (BP-RP > 2); then the target's
calibrated chi2 and the fraction of comparison stars with a larger calibrated difference.
Output: data/pm_local_dr3.csv, data/pm_local_dr2.csv, out/r12_pm_local.json"""
import json, time, warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from astroquery.gaia import Gaia
from scipy.spatial import cKDTree
from scipy.stats import chi2
Gaia.ROW_LIMIT = -1
R = 'tic214311908/analysis/referee_star'
def run(q):
    for k in range(4):
        try: return Gaia.launch_job_async(q).get_results().to_pandas()
        except Exception as e: print('retry', e, flush=True); time.sleep(20)
c = "1=contains(point(ra,dec), circle(260.969674,-51.323821,0.5)) and phot_g_mean_mag between 13.3 and 15.7"
import os
if not os.path.exists(f'{R}/data/pm_local_dr3.csv'):
    d3 = run(f"select source_id, ra, dec, pmra, pmdec, pmra_error, pmdec_error, phot_g_mean_mag, bp_rp, ruwe, parallax from gaiadr3.gaia_source where {c} and astrometric_params_solved > 3")
    d3.to_csv(f'{R}/data/pm_local_dr3.csv', index=False)
    d2 = run(f"select source_id, ra, dec, pmra, pmdec, pmra_error, pmdec_error, phot_g_mean_mag from gaiadr2.gaia_source where {c} and pmra is not null")
    d2.to_csv(f'{R}/data/pm_local_dr2.csv', index=False)
d3 = pd.read_csv(f'{R}/data/pm_local_dr3.csv'); d2 = pd.read_csv(f'{R}/data/pm_local_dr2.csv')
print('DR3', len(d3), 'DR2', len(d2))
cosd = np.cos(np.radians(-51.32))
ra15 = d3.ra - d3.pmra * 0.5 / 3.6e6 / np.cos(np.radians(d3.dec)); de15 = d3.dec - d3.pmdec * 0.5 / 3.6e6
tr = cKDTree(np.column_stack([d2.ra * cosd * 3600, d2.dec * 3600]))
dd, ii = tr.query(np.column_stack([ra15 * cosd * 3600, de15 * 3600]))
k = (dd < 0.3) & (np.abs(d3.phot_g_mean_mag.values - d2.phot_g_mean_mag.values[ii]) < 0.3)
m = d3[k].copy(); m2 = d2.iloc[ii[k]].reset_index(drop=True); m = m.reset_index(drop=True)
for comp in ('pmra', 'pmdec'):
    m[f'd_{comp}'] = m[comp] - m2[comp]; m[f'e_{comp}'] = np.hypot(m[f'{comp}_error'], m2[f'{comp}_error']); m[f'z_{comp}'] = m[f'd_{comp}'] / m[f'e_{comp}']
print('matched', len(m))
out = {}
for lab, sel in (('all', m.ruwe < 1.4), ('red', (m.ruwe < 1.4) & (m.bp_rp > 2.0))):
    s = m[sel]; o = {}
    for comp in ('pmra', 'pmdec'):
        z = s[f'z_{comp}']; o[comp] = dict(n=int(len(z)), med=float(np.median(z)), mad_sd=float(1.4826 * np.median(np.abs(z - np.median(z)))),
                                       med_d=float(np.median(s[f'd_{comp}'])))
    out[lab] = o
    print(lab, o)
# target
t = m[m.source_id == 5925373272850784768]
print(t[['d_pmra', 'd_pmdec', 'e_pmra', 'e_pmdec', 'z_pmra', 'z_pmdec']])
cal = out['all']
zt = [(float(t[f'z_{c_}'].iloc[0]) - cal[c_]['med']) / cal[c_]['mad_sd'] for c_ in ('pmra', 'pmdec')]
ch = zt[0] ** 2 + zt[1] ** 2
s = m[m.ruwe < 1.4]
zz = np.column_stack([(s[f'z_{c_}'] - cal[c_]['med']) / cal[c_]['mad_sd'] for c_ in ('pmra', 'pmdec')])
chs = (zz ** 2).sum(1)
out['target'] = dict(dpm=[float(t.d_pmra.iloc[0]), float(t.d_pmdec.iloc[0])], z_raw=[float(t.z_pmra.iloc[0]), float(t.z_pmdec.iloc[0])], z_cal=zt, chi2_cal=float(ch),
                     p_chi2_2dof=float(chi2.sf(ch, 2)), frac_comparison_larger=float(np.mean(chs >= ch)))
print('target calibrated z', np.round(zt, 2), 'chi2', round(ch, 2), 'p(chi2,2)', round(chi2.sf(ch, 2), 3), 'fraction of local comparison stars larger', round(float(np.mean(chs >= ch)), 3))
json.dump(out, open(f'{R}/out/r12_pm_local.json', 'w'), indent=1)
