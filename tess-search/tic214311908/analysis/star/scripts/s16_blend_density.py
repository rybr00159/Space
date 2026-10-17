"""Sky density of stars that could mimic the signal by blending with TIC 214311908 (l = 338.83, b = -8.55).
A blended eclipsing source must lose >= 15.4 +- 0.8 e-/s (deblended flux loss, phot/): T <= 17.47 for a 100% eclipse,
T <= 16.72 for 50%, T <= 16.2 for 30%. Densities from Gaia DR3 within 330" of the M3V (data/gaia_330.csv; T from the
Stassun+19 relation) and, as a completeness check at the faint end, DECaPS2 within 5' (data/decaps2_field10.csv;
T from the s08 i,z calibration). Expected numbers of such stars inside circles of 0.08", 0.5", 1", 2", 3" and 5"
of the M3V, and inside the 'hidden' region not covered by Gaia DR3/DECaPS2/IPD diagnostics (see README).
Output: out/blend_density.json"""
import json, numpy as np, pandas as pd
B = 'tic214311908/analysis/star'
g = pd.read_csv(f'{B}/data/gaia_330.csv')
c = g.bp_rp.values; G = g.phot_g_mean_mag.values
T = np.where(np.isfinite(c), G - 0.00522555 * c ** 3 + 0.0891337 * c ** 2 - 0.633923 * c + 0.0324473, G - 0.43)
area = np.pi * 330.0 ** 2
# exclude the three named stars' immediate region? (negligible) ; density per arcsec^2 brighter than T
out = dict(l=338.832, b=-8.548, gaia_n_330=int(len(g)), area_arcsec2=area)
dens = {}
for Tl in (14, 15, 16, 16.2, 16.72, 17.0, 17.47, 18.0, 19.0):
    dens[f'T<={Tl}'] = float(np.sum(T <= Tl) / area)
out['gaia_density_per_arcsec2'] = dens
d = pd.read_csv(f'{B}/data/decaps2_field10.csv')
cc = json.load(open(f'{B}/out/census_summary.json'))['T_decaps_cal']['coef']
iz = np.where(np.isfinite(d.mean_mag_z), d.mean_mag_i - d.mean_mag_z, 0.3)
Td = cc[0] + cc[1] * d.mean_mag_i + cc[2] * iz
aread = np.pi * 300.0 ** 2
out['decaps_density_per_arcsec2'] = {f'T<={Tl}': float(np.sum(Td <= Tl) / aread) for Tl in (16.2, 16.72, 17.47, 18.0, 19.0)}
exp = {}
for lab, Tl in (('100pct', 17.47), ('50pct', 16.72), ('30pct', 16.2)):
    rho = max(dens[f'T<={Tl}'], out['decaps_density_per_arcsec2'][f'T<={Tl}'])
    exp[lab] = {f'r<{r}"': float(rho * np.pi * r ** 2) for r in (0.08, 0.5, 1.0, 2.0, 3.0, 5.0)}
out['expected_number_within_radius'] = exp
print(json.dumps(out, indent=1))
json.dump(out, open(f'{B}/out/blend_density.json', 'w'), indent=1)
