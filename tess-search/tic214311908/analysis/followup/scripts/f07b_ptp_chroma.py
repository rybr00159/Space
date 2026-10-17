"""f07b: posterior distribution of the chromatic depth ratios expected if the signal is a planet on an unresolved bound
companion (PTP). Weights: out/ptp_draw_weights.npz (f07; referee2_fpp draws and r05 PTP evidence weights); band ratios
as a function of the companion's TESS contrast dT: out/chromatic_ratio_vs_dT.csv (f04, BT-Settl + filter curves).
Ratios are depth_X / depth_Y for the same event, to compare with 1 (planet on the M3V, after limb darkening).
Output: out/ptp_chroma.json"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flib import *
import numpy as np, pandas as pd
W = np.load(f'{OUT}/ptp_draw_weights.npz')
C = pd.read_csv(f'{OUT}/chromatic_ratio_vs_dT.csv').sort_values('dT')
dT = W['dT'].astype(float)
PAIRS = [('g', 'z_s'), ('r', 'z_s'), ('g', 'i'), ('r', 'i'), ('i', 'z_s'), ('V', 'Ic'), ('ASTEP_blue', 'ASTEP_red'), ('V', 'I+z'), ('r', 'Ic'), ('g', 'TESS'), ('I+z', 'TESS'), ('J', 'TESS')]
out = {}
for case in [k for k in W.files if '__' in k]:
    w = W[case].astype(float)
    ok = (dT >= C.dT.min()) & (dT <= C.dT.max()) & (w > 0)
    res = dict(frac_weight_in_grid=float(w[ok].sum()/w.sum()))
    for a, b in PAIRS:
        r = np.interp(dT[ok], C.dT, C[f'ratio_{a}'])/np.interp(dT[ok], C.dT, C[f'ratio_{b}'])
        o = np.argsort(r); cw = np.cumsum(w[ok][o])/w[ok].sum()
        res[f'{a}/{b}'] = [float(x) for x in np.interp([0.025, 0.16, 0.5, 0.84, 0.975], cw, r[o])]
    out[case] = res
    print(case, {k: [round(x, 3) for x in v] for k, v in res.items() if '/' in k})
jdump(dict(ratios_2p5_16_50_84_97p5=out, note='depth ratio predicted for a planet on the unresolved bound companion; planet on the M3V gives 1 (after LD)'),
      f'{OUT}/ptp_chroma.json')
