"""r08: empirical validation of the colour-residual twin discriminant with KNOWN binaries.
The discriminant assumes an unresolved near-equal pair has r(G-Ks) (and r(G-W1)) of its components at fixed combined
BP-RP (~normal-star value ~+0.01), whereas metal-rich singles at the same over-luminosity have ~+0.07.
Validation sets (field sample, |dBP-RP| < 0.12 around the M3V colour, r01 local ridges, 2MASS AAA, |C*| < 3):
  (a) Gaia-IPD double-peaked sources (ipd_frac_multi_peak >= 10): pairs at ~0.1-1" whose light is blended in G, BP,
      RP and 2MASS -> photometrically 'unresolved binaries'; brightness ratios mostly near-equal (multi-peak needs a
      comparably bright companion);
  (b) astrometric binaries RUWE 1.4-8 (mostly unequal pairs: a near-equal pair has little photocentre motion);
  (c) Gaia NSS (non_single_star > 0).
For each set: median r(G-Ks), r(G-W1) in bins of dMKs, compared with clean (RUWE<1.1) stars at the same dMKs and with
the normal-star value. The colour residual does not depend on the parallax, so astrometric parallax biases of
binaries do not matter for y (they shift x only).
Output: out/r08_binary_validation.json"""
import json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
H = 'tic214311908/analysis/host'
B = 'tic214311908/analysis/referee2_host'
r01 = json.load(open(f'{B}/out/r01_population.json'))
R = {k: (v['slope'], v['zero_at_target_colour']) for k, v in r01['ridge_local'].items()}
C0 = 2.756149
F = pd.read_csv(f'{H}/data/field_mdwarfs.csv')
x = F.bp_rp; fx = 1.162004 + 0.011464 * x + 0.049255 * x ** 2 - 0.005879 * x ** 3
F['Cs'] = (F.phot_bp_rp_excess_factor - fx) / (0.0059898 + 8.817481e-12 * F.phot_g_mean_mag ** 7.618399)
dm = 5 * np.log10(F.parallax / 100)
F['MK'] = F.ks_m + dm; F['GK'] = F.phot_g_mean_mag - F.ks_m; F['GW1'] = F.phot_g_mean_mag - F.w1mpro; F['MG'] = F.phot_g_mean_mag + dm
for q in ('MK', 'GK', 'GW1', 'MG'): F['d' + q] = F[q] - (R[q][0] * (F.bp_rp - C0) + R[q][1])
wok = F.wise_cc.astype(str).str[:2].eq('00') & F.wise_ph_qual.astype(str).str[:2].eq('AA') & (F.wise_ext == 0)
base = ((F.bp_rp - C0).abs() < 0.12) & (F.ph_qual == 'AAA') & (F.ks_msigcom < 0.05) & (F.Cs.abs() < 3) & (F.parallax / F.parallax_error > 10)
clean = base & (F.ruwe < 1.1) & (F.ipd_frac_multi_peak == 0) & (F.ipd_gof_harmonic_amplitude < 0.03) & (F.parallax / F.parallax_error > 20)
sets = {'IPD multi-peak >= 10': base & (F.ipd_frac_multi_peak >= 10),
        'IPD multi-peak >= 30': base & (F.ipd_frac_multi_peak >= 30),
        'RUWE 1.4-8, multi-peak <= 2': base & F.ruwe.between(1.4, 8) & (F.ipd_frac_multi_peak <= 2),
        'NSS > 0': base & (F.non_single_star > 0)}
bins = [(-1.3, -0.95), (-0.95, -0.55), (-0.55, -0.25), (-0.25, 0.15), (0.15, 0.6)]
out = {}
def med(v): v = v[np.isfinite(v)]; return (float(np.median(v)), float(1.4826 * np.median(np.abs(v - np.median(v)))) if len(v) else None, int(len(v))) if len(v) else (None, None, 0)
cl = {}
for lo, hi in bins:
    k = clean & F.dMK.between(lo, hi); cl[f'{lo}..{hi}'] = dict(rGK=med(F.dGK[k].values), rGW1=med(F.dGW1[k & wok].values))
out['clean'] = cl
print('clean', cl)
for name, s in sets.items():
    d = {}
    for lo, hi in bins:
        k = s & F.dMK.between(lo, hi); d[f'{lo}..{hi}'] = dict(rGK=med(F.dGK[k].values), rGW1=med(F.dGW1[k & wok].values),
                                                                frac_rGK_lt_0p04=float(np.mean(F.dGK[k] < 0.04)) if k.sum() else None)
    # in the M3V's window
    k = s & F.dMK.between(-0.9, -0.6)
    d['window_-0.9..-0.6'] = dict(rGK=med(F.dGK[k].values), frac_rGK_ge_target=float(np.mean(F.dGK[k] >= r01['target']['colour_residuals']['dGK'])) if k.sum() else None)
    out[name] = d
    print(name, {kk: (vv['rGK'], vv.get('rGW1')) for kk, vv in d.items()})
k = clean & F.dMK.between(-0.9, -0.6)
out['clean_window_frac_rGK_ge_target'] = float(np.mean(F.dGK[k] >= r01['target']['colour_residuals']['dGK']))
# slope of r(G-Ks) vs dMK for multi-peak pairs vs clean (binaries should have a shallower slope)
for name, s in (('clean', clean), ('IPD multi-peak >= 10', sets['IPD multi-peak >= 10']), ('RUWE 1.4-8', sets['RUWE 1.4-8, multi-peak <= 2'])):
    k = s & F.dMK.between(-1.0, 0.3) & np.isfinite(F.dGK)
    c = np.polyfit(F.dMK[k], F.dGK[k], 1); out[f'slope_{name}'] = dict(slope=float(c[0]), intercept=float(c[1]), n=int(k.sum()))
    print('slope', name, out[f'slope_{name}'])
json.dump(out, open(f'{B}/out/r08_binary_validation.json', 'w'), indent=1)
