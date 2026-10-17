"""h08: published photometric M-dwarf metallicity relations applied to the M3V, with Monte Carlo errors, for three
near-IR choices (2MASS; DENIS raw; DENIS corrected for the frame zero point measured in h06b) and for the
single-star and equal-twin interpretations.
Relations (coefficients checked against the papers in data/lit/):
  B05  Bonfils+05:  [Fe/H] = 0.196 - 1.527 MK + 0.091 MK^2 + 1.886 (V-K) - 0.142 (V-K)^2           (sigma 0.20)
  JA09 Johnson & Apps 09: [Fe/H] = 0.56 dMK - 0.05, dMK = MS(V-K) - MK, MS = sum a_i (V-K)^i,
       a = (-9.58933, 17.3952, -8.88365, 2.22598, -0.258854, 0.0113399) (polynomial from Neves+12 text)  (sigma 0.24, Neves+12)
  SL10 Schlaufman & Laughlin 10: [Fe/H] = 0.79 d(V-K) - 0.17, d(V-K) = (V-K)obs - (V-K)iso(MK),
       (V-K)iso = sum b_i MK^i, b = (51.1413, -39.3756, 12.2862, -1.83916, 0.134266, -0.00382023)          (sigma 0.19)
  N12  Neves+12 recalibration of SL10: [Fe/H] = 0.57 d(V-K) - 0.17                                     (sigma 0.17)
  R21  Rains+21: (BP-Ks)_MS = 0.05385 MK^3 - 1.08356 MK^2 + 7.76175 MK - 14.54705; [Fe/H] = 0.71339 d(BP-Ks) - 0.04301 (sigma 0.19)
V: XP synthetic Johnson V (GaiaXPy JKC_Std, 15.591; h05) with 0.03 mag error; cross-check Riello+21 G-V -> 15.608.
Output: out/h08_photfeh.json"""
import json, numpy as np
B = 'tic214311908/analysis/host'
rng = np.random.default_rng(3)
N = 20000
plx = rng.normal(15.549, 0.036, N)
V = rng.normal(15.591, 0.03, N); BP = rng.normal(15.8307, 0.005, N)
Ks_opts = {'2MASS': (10.482, 0.026), 'DENIS_raw': (10.823, 0.09), 'DENIS_zpcorr': (10.823 - 0.3025, np.hypot(0.09, 0.038))}
a = np.array([-9.58933, 17.3952, -8.88365, 2.22598, -0.258854, 0.0113399])
b = np.array([51.1413, -39.3756, 12.2862, -1.83916, 0.134266, -0.00382023])
def rel(MK, VK, BPK):
    out = {}
    out['B05'] = 0.196 - 1.527 * MK + 0.091 * MK ** 2 + 1.886 * VK - 0.142 * VK ** 2
    MS = sum(a[i] * VK ** i for i in range(6)); out['JA09'] = 0.56 * (MS - MK) - 0.05
    VKiso = sum(b[i] * MK ** i for i in range(6)); out['SL10'] = 0.79 * (VK - VKiso) - 0.17; out['N12'] = 0.57 * (VK - VKiso) - 0.17
    BPKms = 0.05385 * MK ** 3 - 1.08356 * MK ** 2 + 7.76175 * MK - 14.54705; out['R21'] = 0.71339 * (BPK - BPKms) - 0.04301
    return out
SIG = dict(B05=0.20, JA09=0.24, SL10=0.19, N12=0.17, R21=0.19)
res = {}
for lab, (k0, ek) in Ks_opts.items():
    K = rng.normal(k0, ek, N)
    for interp, dm in (('single', 0.0), ('equal_twin', 2.5 * np.log10(2))):
        MK = K + 5 * np.log10(plx / 100) + dm
        r = rel(MK, V - K, BP - K)
        res[f'{lab}_{interp}'] = {k: dict(feh=float(np.median(v)), err_meas=float(np.std(v)), err_tot=float(np.hypot(np.std(v), SIG[k]))) for k, v in r.items()}
        res[f'{lab}_{interp}']['MK'] = float(np.median(MK))
for k, v in res.items():
    print(k, ' '.join(f'{n}={d["feh"]:+.2f}+-{d["err_tot"]:.2f}' for n, d in v.items() if n != 'MK'), f'MK={v["MK"]:.3f}')
json.dump(dict(results=res, V_used=[15.591, 0.03], V_riello=15.608, notes='equal_twin: each component 0.753 mag fainter in all bands, colours unchanged'),
          open(f'{B}/out/h08_photfeh.json', 'w'), indent=1)
