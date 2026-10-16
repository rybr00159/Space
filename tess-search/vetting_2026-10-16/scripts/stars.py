"""TIC + Gaia DR3 parameters for the five recurring signals."""
import json, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from astroquery.mast import Catalogs
from astroquery.gaia import Gaia
TICS = [90086266, 61332759, 214311906, 55829034, 61141412]
rows = []
for tic in TICS:
    t = Catalogs.query_criteria(catalog="Tic", ID=tic)[0]
    rows.append({k: (t[k].item() if hasattr(t[k], "item") else t[k]) for k in
                 ("ID", "ra", "dec", "Tmag", "Vmag", "Teff", "logg", "rad", "mass", "d", "plx", "GAIA", "contratio",
                  "lumclass", "pmRA", "pmDEC", "Kmag", "Jmag")})
tic = pd.DataFrame(rows)
ids = ",".join(str(int(x)) for x in tic.GAIA)
q = f"""SELECT g.source_id, g.ra, g.dec, g.parallax, g.parallax_error, g.pmra, g.pmdec, g.phot_g_mean_mag, g.bp_rp, g.ruwe,
  g.radial_velocity, g.radial_velocity_error, g.rv_nb_transits, g.rv_amplitude_robust, g.rv_chisq_pvalue,
  g.non_single_star, g.phot_variable_flag, g.ipd_frac_multi_peak, g.ipd_gof_harmonic_amplitude, g.phot_bp_rp_excess_factor,
  g.teff_gspphot, g.logg_gspphot, g.distance_gspphot, g.classprob_dsc_combmod_star,
  a.radius_flame, a.mass_flame, a.lum_flame, a.evolstage_flame, a.teff_gspspec, a.logg_gspspec, a.radius_gspphot
FROM gaiadr3.gaia_source AS g LEFT JOIN gaiadr3.astrophysical_parameters AS a ON a.source_id = g.source_id
WHERE g.source_id IN ({ids})"""
g = Gaia.launch_job(q).get_results().to_pandas()
v = Gaia.launch_job(f"SELECT * FROM gaiadr3.vari_classifier_result WHERE source_id IN ({ids})").get_results().to_pandas()
nss = Gaia.launch_job(f"SELECT source_id, nss_solution_type, period, eccentricity, significance FROM gaiadr3.nss_two_body_orbit WHERE source_id IN ({ids})").get_results().to_pandas()
acc = Gaia.launch_job(f"SELECT source_id, nss_solution_type, significance FROM gaiadr3.nss_acceleration_astro WHERE source_id IN ({ids})").get_results().to_pandas()
ebs = Gaia.launch_job(f"SELECT source_id, frequency, global_ranking FROM gaiadr3.vari_eclipsing_binary WHERE source_id IN ({ids})").get_results().to_pandas()
tic["GAIA"] = tic.GAIA.astype("int64")
m = tic.merge(g, left_on="GAIA", right_on="source_id", how="left")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 60)
for _, r in m.iterrows():
    print(f"\nTIC {int(r.ID)}  Gaia {int(r.GAIA)}  RA {r.ra_x:.5f} Dec {r.dec_x:.5f}  T={r.Tmag:.2f} V={r.Vmag}  K={r.Kmag}")
    print(f"  TIC: Teff {r.Teff} logg {r.logg} R {r.rad} M {r.mass} d {r.d} lumclass {r.lumclass} contratio {r.contratio}")
    print(f"  Gaia: plx {r.parallax:.3f}±{r.parallax_error:.3f} (d≈{1000/r.parallax if r.parallax>0 else np.nan:.0f} pc) G {r.phot_g_mean_mag:.2f} BP-RP {r.bp_rp:.2f} RUWE {r.ruwe:.2f} "
          f"ipd_multi {r.ipd_frac_multi_peak} NSS {r.non_single_star} var {r.phot_variable_flag}")
    print(f"  Gaia RV {r.radial_velocity}±{r.radial_velocity_error} n {r.rv_nb_transits} amp {r.rv_amplitude_robust} p {r.rv_chisq_pvalue}")
    print(f"  GSP-Phot Teff {r.teff_gspphot} logg {r.logg_gspphot} R {r.radius_gspphot} d {r.distance_gspphot}; FLAME R {r.radius_flame} M {r.mass_flame} L {r.lum_flame} evol {r.evolstage_flame}; GSP-Spec Teff {r.teff_gspspec} logg {r.logg_gspspec}")
print("\nvari_classifier:\n", v[["source_id", "best_class_name", "best_class_score"]].to_string() if len(v) else "none")
print("NSS two-body:\n", nss.to_string() if len(nss) else "none")
print("NSS acceleration:\n", acc.to_string() if len(acc) else "none")
print("Gaia EB:\n", ebs.to_string() if len(ebs) else "none")
m.to_csv("stars.csv", index=False)
