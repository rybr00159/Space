"""s11: validation of specs.py: synthetic DECam/SkyMapper photometry of the M3V (BT-Settl, Teff 3220-3350 K,
R from Fbol) vs observed (DECaPS2 is near saturation for the M3V; SkyMapper DR4) and of Mamajek main-sequence
stars (K4V, G5V) colours; prints the table for the README. Output: out/specs_check.json"""
import sys, json, numpy as np
sys.path.insert(0, 'tic214311908/analysis/shape/scripts')
from specs import *
d = 64.31; DM = 5*np.log10(d/10)
out = {}
for T, R in [(3220, 0.4134), (3350, 0.3813)]:
    m = {b: float(absmag(T, R, b) + DM) for b in ['g', 'r', 'i', 'z', 'smv', 'smg', 'smu']}
    out[f'M3V_{T}'] = m; print(T, R, {k: round(v, 3) for k, v in m.items()})
out['M3V_obs'] = dict(g_decaps=16.177, r_decaps=14.737, i_decaps=13.355, z_decaps=12.761, smv=18.192, smg=15.817, smu_note='SkyMapper u not detected/used')
for lab, T, R in [('K4V', 4600, 0.713), ('G5V', 5660, 0.977), ('A0V', 9700, 2.0), ('WD7000', 7000, 0.0125)]:
    m = {b: float(absmag(T, R, b, compact=(lab.startswith('WD')))) for b in ['g', 'r', 'i', 'z', 'T']}
    out[lab] = m; print(lab, {k: round(v, 3) for k, v in m.items()}, 'g-r', round(m['g']-m['r'], 3), 'r-i', round(m['r']-m['i'], 3), 'i-z', round(m['i']-m['z'], 3))
json.dump(out, open(f'{SD}/out/specs_check.json', 'w'), indent=1)
