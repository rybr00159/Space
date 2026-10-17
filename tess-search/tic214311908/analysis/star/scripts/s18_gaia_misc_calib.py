"""Calibrations of two further Gaia DR3 quantities of TIC 214311908 against comparison M dwarfs
(data/mdwarf_sample.csv, BP-RP 2.65-2.85, G 13.9-14.7, parallax_over_error > 20):
(a) classprob_dsc_specmod_binarystar = 0.41 (DSC 'binary star' probability from XP spectra);
(b) radial_velocity_error = 2.22 km/s from 21 transits (an SB2 or a fast rotator would broaden the combined CCF and
    inflate the error; Gaia gives no epoch RVs for stars this faint);
(c) RUWE and ipd_gof_harmonic_amplitude percentiles.
Output: out/gaia_misc_calib.json"""
import json, numpy as np, pandas as pd
B = 'tic214311908/analysis/star'
m = pd.read_csv(f'{B}/data/mdwarf_sample.csv')
s = m[(m.bp_rp > 2.65) & (m.bp_rp < 2.85) & (m.phot_g_mean_mag > 13.9) & (m.phot_g_mean_mag < 14.7)]
out = dict(n=int(len(s)))
x = s.classprob_dsc_specmod_binarystar.dropna()
out['dsc_binary'] = dict(n=int(len(x)), median=float(x.median()), frac_ge_0p41=float(np.mean(x >= 0.4147)), pct=[float(np.percentile(x, q)) for q in (16, 50, 84, 95)])
r = s[s.rv_nb_transits.notna() & (s.rv_nb_transits >= 15) & (s.rv_nb_transits <= 30)]
out['rv_error'] = dict(n=int(len(r)), median=float(r.radial_velocity_error.median()), frac_ge_2p22=float(np.mean(r.radial_velocity_error >= 2.217)),
                       pct=[float(np.percentile(r.radial_velocity_error, q)) for q in (16, 50, 84, 95)])
c = s[s.ruwe < 5]
out['ruwe'] = dict(frac_le_1p076=float(np.mean(c.ruwe <= 1.076)), median=float(c.ruwe.median()))
out['gofha'] = dict(frac_le_0p0132=float(np.mean(s.ipd_gof_harmonic_amplitude <= 0.013223)), median=float(s.ipd_gof_harmonic_amplitude.median()))
out['fmp0_frac'] = float(np.mean(s.ipd_frac_multi_peak == 0))
print(json.dumps(out, indent=1))
json.dump(out, open(f'{B}/out/gaia_misc_calib.json', 'w'), indent=1)
