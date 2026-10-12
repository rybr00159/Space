"""Is the Gaia DR3 RV uncertainty of TIC 87055315 (radial_velocity_error 7.34 km/s, faint-star method 2, 7 transits,
3 deblended, rv_expected_sig_to_noise 4.17, grvs_mag 13.23, template Teff 5250 / logg 4.0) unusually large, as a
combined-CCF RV for a star with an orbital RV amplitude of K ~ 20-40 km/s would be?
Comparison sample (Gaia archive, ADQL): stars within 12 deg, rv_method_used = 2, |grvs_mag - 13.23| < 0.2,
rv_template_teff 5000-6000, rv_nb_transits 5-9, rv_expected_sig_to_noise 3.5-5.0; subsamples by rv_nb_deblended_transits.
Output: results/gaia_rv.json
"""
import os, json, time
import numpy as np
from astroquery.gaia import Gaia
ROOT = "/home/user/work/tic87/neighbour"
q = """select source_id, radial_velocity, radial_velocity_error, grvs_mag, rv_nb_transits, rv_nb_deblended_transits,
rv_expected_sig_to_noise, rv_template_teff, rv_template_logg, ruwe, non_single_star, phot_g_mean_mag, bp_rp
from gaiadr3.gaia_source
where 1 = contains(point('ICRS', ra, dec), circle('ICRS', 272.609, -43.994, 12))
and rv_method_used = 2 and grvs_mag between 13.03 and 13.43 and rv_template_teff between 5000 and 6000
and rv_nb_transits between 5 and 9 and rv_expected_sig_to_noise between 3.5 and 5.0"""
df = None
for i in range(4):
    try:
        df = Gaia.launch_job_async(q).get_results().to_pandas(); break
    except Exception as e:
        print("retry", i, repr(e)[:200], flush=True); time.sleep(20)
if df is None:
    json.dump(dict(error="Gaia query failed"), open(os.path.join(ROOT, "results/gaia_rv.json"), "w")); raise SystemExit
df.to_csv(os.path.join(ROOT, "data/cat/gaia_rv_comparison.csv"), index=False)
val = 7.337608
out = dict(n_all=int(len(df)), value=val)
for lab, m in [("all", np.ones(len(df), bool)), ("deblended_ge1", df.rv_nb_deblended_transits.values >= 1),
               ("deblended_2to4", (df.rv_nb_deblended_transits.values >= 2) & (df.rv_nb_deblended_transits.values <= 4)),
               ("nss0_ruwe_lt1p2", (df.non_single_star.values == 0) & (df.ruwe.values < 1.2)),
               ("transits_6to8", (df.rv_nb_transits.values >= 6) & (df.rv_nb_transits.values <= 8))]:
    e = df.radial_velocity_error.values[m]
    e = e[np.isfinite(e)]
    if len(e) < 10: out[lab] = dict(n=int(len(e))); continue
    out[lab] = dict(n=int(len(e)), median=float(np.median(e)), p16=float(np.percentile(e, 16)),
                    p84=float(np.percentile(e, 84)), p95=float(np.percentile(e, 95)),
                    percentile_of_value=float(100 * np.mean(e < val)))
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(ROOT, "results/gaia_rv.json"), "w"), indent=1)
