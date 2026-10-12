"""Harmonised table of independent stellar-parameter estimates for TIC 87055304 (values read from the rows saved by
catalogs.py in out/catalog_raw and from data/gaia_dr3_astrophysical_parameters_target.csv), plus this work's SED fit.
Output: out/catalog_compilation.csv and out/catalog_compilation.md"""
import json, numpy as np, pandas as pd
R = "/home/user/work/tic87/stellar/out/catalog_raw/"
ap = pd.read_csv("/home/user/work/tic87/data/gaia_dr3_astrophysical_parameters_target.csv", header=None, index_col=0)[1]
sup = pd.read_csv(R + "gaiadr3_ap_supp.csv").iloc[0]
tic = pd.read_csv(R + "tic82_mast.csv").iloc[0]
dr2 = pd.read_csv(R + "gaiadr2_source.csv").iloc[0]
sh19 = pd.read_csv(R + "StarHorse19.csv").iloc[0]; sh21 = pd.read_csv(R + "StarHorse21.csv").iloc[0]
shb = pd.read_csv(R + "SHBoost24.csv").iloc[0]; zh = pd.read_csv(R + "Zhang23.csv").iloc[0]
hu = pd.read_csv(R + "Huson25.csv").iloc[0]; ye = pd.read_csv(R + "Ye25.csv").iloc[0]
hg = pd.read_csv(R + "Huang22.csv").iloc[0]; an = pd.read_csv(R + "external_xgboost_table1.csv").iloc[0]
bj21 = pd.read_csv(R + "BJ21.csv").iloc[0]; bj18 = pd.read_csv(R + "BJ18.csv").iloc[0]
f = lambda x: float(x)
def rng(med, lo, hi, fmt="{:.0f}"): return f"{fmt.format(med)} ({fmt.format(lo)}-{fmt.format(hi)})"
rows = [
 dict(source="TIC v8.2 (Stassun+19, Paegert+21)", method="dereddened colour-Teff spline ('splin', 'dered'); d not from a DR2 parallax (plx field empty)",
      uses_plx="no (plx empty) / unclear", Teff="7276 (no error)", logg="-", feh="-", R="1.79", M="-", d_pc="464", ext="-", note="Teff inconsistent with every other estimate; lumclass DWARF"),
 dict(source="Gaia DR2 Apsis-Priam", method="BP-RP/G colours (ML)", uses_plx="no", Teff=rng(f(dr2.teff_val), f(dr2.teff_percentile_lower), f(dr2.teff_percentile_upper)),
      logg="-", feh="-", R="-", M="-", d_pc=f"1/plx = {1000/dr2.parallax:.0f} (plx {dr2.parallax:.3f}+-{dr2.parallax_error:.3f})", ext=f"A_G {dr2.a_g_val:.3f}", note="DR2 excess noise 2.96 mas"),
 dict(source="Gaia DR3 GSP-Phot MARCS (adopted 'best')", method="XP + G + parallax likelihood, MARCS library", uses_plx="yes (likelihood)",
      Teff=rng(f(ap.teff_gspphot), f(ap.teff_gspphot_lower), f(ap.teff_gspphot_upper)), logg=f"{f(ap.logg_gspphot):.2f}", feh=f"{f(ap.mh_gspphot):.2f}",
      R=f"{f(ap.radius_gspphot):.2f}", M="-", d_pc=rng(f(ap.distance_gspphot), f(ap.distance_gspphot_lower), f(ap.distance_gspphot_upper)), ext=f"A0 {f(ap.azero_gspphot):.3f}",
      note="5.1 sigma from the DR3 parallax; MCMC intervals unrealistically narrow"),
 dict(source="Gaia DR3 GSP-Phot PHOENIX", method="same, PHOENIX library", uses_plx="yes",
      Teff=rng(sup.teff_gspphot_phoenix, sup.teff_gspphot_phoenix_lower, sup.teff_gspphot_phoenix_upper), logg=f"{sup.logg_gspphot_phoenix:.2f}", feh=f"{sup.mh_gspphot_phoenix:.2f}",
      R=f"{sup.radius_gspphot_phoenix:.2f}", M="-", d_pc=rng(sup.distance_gspphot_phoenix, sup.distance_gspphot_phoenix_lower, sup.distance_gspphot_phoenix_upper), ext=f"A0 {sup.azero_gspphot_phoenix:.3f}",
      note="lower log-posterior than MARCS (-1483 vs -1124)"),
 dict(source="Gaia DR3 GSP-Phot A-stars library", method="same, A library", uses_plx="yes",
      Teff=rng(sup.teff_gspphot_a, sup.teff_gspphot_a_lower, sup.teff_gspphot_a_upper), logg=f"{sup.logg_gspphot_a:.2f}", feh=f"{sup.mh_gspphot_a:.2f}",
      R=f"{sup.radius_gspphot_a:.2f}", M="-", d_pc=rng(sup.distance_gspphot_a, sup.distance_gspphot_a_lower, sup.distance_gspphot_a_upper), ext=f"A0 {sup.azero_gspphot_a:.3f}", note="log-posterior -2651"),
 dict(source="Gaia DR3 FLAME", method="L from G, parallax, BC; R from L + GSP-Phot Teff; M, age from tracks", uses_plx="yes (directly)",
      Teff="(GSP-Phot)", logg="-", feh="-", R=rng(f(ap.radius_flame), f(ap.radius_flame_lower), f(ap.radius_flame_upper), "{:.2f}"),
      M=rng(f(ap.mass_flame), f(ap.mass_flame_lower), f(ap.mass_flame_upper), "{:.2f}"), d_pc="(parallax)", ext="-",
      note=f"L {f(ap.lum_flame):.2f} Lsun, age {f(ap.age_flame):.1f} Gyr, evolstage {int(f(ap.evolstage_flame))}"),
 dict(source="Gaia DR3 MSC (two-star model)", method="XP + parallax, unresolved binary of two stars", uses_plx="yes",
      Teff=f"{f(ap.teff_msc1):.0f} + {f(ap.teff_msc2):.0f} (very broad)", logg=f"{f(ap.logg_msc1):.2f} / {f(ap.logg_msc2):.2f}", feh=f"{f(ap.mh_msc):.2f}", R="-", M="-",
      d_pc=rng(f(ap.distance_msc), f(ap.distance_msc_lower), f(ap.distance_msc_upper)), ext=f"A0 {f(ap.azero_msc):.2f}", note="both components dwarfs"),
 dict(source="Gaia DR3 ESP-HS / RVS template", method="XP spectral type; RVS template (rv_atm_param_origin 111)", uses_plx="-",
      Teff="SpT F ; RVS template 5500", logg="RVS template 3.5", feh="RVS template -1.0", R="-", M="-", d_pc="-", ext="-",
      note="RVS template parameters are coarse grid nodes, not an independent determination"),
 dict(source="Bailer-Jones+21 (EDR3)", method="geometric / photogeometric distance prior x parallax", uses_plx="yes", Teff="-", logg="-", feh="-", R="-", M="-",
      d_pc=f"rgeo {rng(bj21.rgeo, bj21.b_rgeo, bj21.B_rgeo)}; rpgeo {rng(bj21.rpgeo, bj21.b_rpgeo, bj21.B_rpgeo)}", ext="-", note="no RUWE-dependent inflation"),
 dict(source="Bailer-Jones+18 (DR2)", method="distance prior x DR2 parallax", uses_plx="yes", Teff="-", logg="-", feh="-", R="-", M="-",
      d_pc=rng(bj18.rest, bj18.b_rest, bj18.B_rest), ext="-", note=""),
 dict(source="StarHorse 2019 (Anders+19, DR2)", method="Bayesian isochrone fit, photometry + DR2 parallax", uses_plx="yes",
      Teff=rng(sh19.teff50, sh19.teff16, sh19.teff84), logg=rng(sh19.logg50, sh19.logg16, sh19.logg84, "{:.2f}"), feh=f"{sh19.met50:.2f}",
      R="-", M=rng(sh19.mass50, sh19.mass16, sh19.mass84, "{:.2f}"), d_pc=rng(1000 * sh19.dist50, 1000 * sh19.dist16, 1000 * sh19.dist84), ext=f"A_V {sh19.AV50:.2f}", note=""),
 dict(source="StarHorse 2021 (Anders+22, EDR3)", method="Bayesian isochrone fit, photometry + EDR3 parallax", uses_plx="yes",
      Teff=rng(sh21.teff50, sh21.teff16, sh21.teff84), logg=rng(sh21.logg50, sh21.logg16, sh21.logg84, "{:.2f}"), feh=f"{sh21.met50:.2f}",
      R="-", M=rng(sh21.mass50, sh21.mass16, sh21.mass84, "{:.2f}"), d_pc=rng(1000 * sh21.dist50, 1000 * sh21.dist16, 1000 * sh21.dist84), ext=f"A_V {sh21.AV50:.2f}", note=f"fidelity {sh21.fidelity:.2f}"),
 dict(source="SHBoost (Khalatyan+24)", method="XGBoost on XP + photometry + 5p astrometry", uses_plx="yes (feature)",
      Teff=f"{10**shb.logTeff:.0f} +- {10**shb.logTeffmean*np.log(10)*shb.s_logTeffmean:.0f}", logg=f"{shb.logg:.2f} +- {shb.s_loggmean:.2f}", feh=f"{shb.Met:.2f} +- {shb.s_Metmean:.2f}",
      R="-", M=f"{shb.Mass:.2f} +- {shb.s_Massmean:.2f}", d_pc="(StarHorse21)", ext=f"A_V {shb.AV:.2f} +- {shb.s_AVmean:.2f}", note=""),
 dict(source="Zhang, Green & Rix 2023", method="data-driven XP forward model + 2MASS/WISE, parallax as likelihood", uses_plx="yes (likelihood)",
      Teff=f"{1000*zh.Teff:.0f} +- {1000*zh.e_Teff:.0f}", logg=f"{zh.logg:.2f} +- {zh.e_logg:.2f}", feh=f"{zh['[Fe/H]']:.2f} +- {zh['e_[Fe/H]']:.2f}", R="-", M="-",
      d_pc=f"plx_post {zh.Plx:.3f} +- {zh.e_Plx:.3f} mas -> {1000/zh.Plx:.0f} pc", ext=f"E {zh.Ext:.3f} +- {zh.e_Ext:.3f}",
      note="XP-only parallax implied ~3.0 +- 0.4 mas if Gaia plx enters as independent Gaussian (see README)"),
 dict(source="Andrae, Rix & Chandra 2023 (XGBoost)", method="XGBoost on XP + CatWISE (+ parallax feature)", uses_plx="yes (feature)",
      Teff=f"{an.teff_xgboost:.0f}", logg=f"{an.logg_xgboost:.2f}", feh=f"{an.mh_xgboost:.2f}", R="-", M="-", d_pc="-", ext="-", note="no per-star errors"),
 dict(source="Huson+25 (Gaia Net)", method="neural network on XP spectra", uses_plx="not stated",
      Teff=f"{10**hu.logTeff:.0f} +- {10**hu.logTeff*np.log(10)*hu.e_logTeff:.0f}", logg=f"{hu['log(g)']:.2f} +- {hu['e_log(g)']:.2f}", feh=f"{hu['[Fe/H]']:.2f} +- {hu['e_[Fe/H]']:.2f}",
      R="-", M="-", d_pc="-", ext="-", note="formal errors likely underestimated"),
 dict(source="Ye+25 (FERRE fit of corrected XP spectra)", method="synthetic-spectrum fit to XP shape", uses_plx="no",
      Teff=f"{ye.Teff:.0f} (sys +-167)", logg=f"{ye.logg:.2f} (sys +-0.40)", feh=f"{ye['[Fe/H]']:.2f} (sys +-0.19)", R="-", M="-", d_pc="-", ext="-", note="parallax-independent; logg precision 0.4 dex"),
 dict(source="Huang+22 (SkyMapper DR2 + Gaia colours)", method="photometric Teff/[Fe/H]; dwarf/giant class", uses_plx="no (RUWE 18.8; distance flag 'TO')",
      Teff=f"{hg.Teff:.0f} +- {hg.e_Teff:.0f}", logg=f"class: {hg.Type}", feh=f"{hg['[Fe/H]']:.2f} +- {hg['e_[Fe/H]']:.2f}", R="-", M="-",
      d_pc=f"{1000*hg.Dist:.0f} +- {1000*hg.e_Dist:.0f} (flag TO; not meaningful)", ext=f"E(B-V)SFD {hg['E(B-V)sfd']:.3f}", note=""),
]
sed = json.load(open("/home/user/work/tic87/stellar/out/sed/noapass.json"))
rows.append(dict(source="This work: SED (CK04 + F99), 14 bands", method="Gaia, SkyMapper DR4, 2MASS, WISE W1/W2; jitter", uses_plx="no",
                 Teff=f"{sed['Teff']['med']:.0f} +{sed['Teff']['hi1']:.0f}/-{sed['Teff']['lo1']:.0f}", logg=f"{sed['logg']['med']:.2f} +{sed['logg']['hi1']:.2f}/-{sed['logg']['lo1']:.2f}",
                 feh=f"{sed['mh']['med']:.2f} +{sed['mh']['hi1']:.2f}/-{sed['mh']['lo1']:.2f}", R=f"{sed['R_per_100pc']['med']:.4f} x d/100pc", M="-",
                 d_pc="-", ext=f"A_V {sed['AV']['med']:.2f} +{sed['AV']['hi1']:.2f}/-{sed['AV']['lo1']:.2f}",
                 note=f"theta = {sed['theta_mas']['med']:.4f} +- {0.5*(sed['theta_mas']['hi1']+sed['theta_mas']['lo1']):.4f} mas; Fbol {sed['fbol_cgs']['med']:.3e} erg/s/cm2"))
df = pd.DataFrame(rows)
df.to_csv("/home/user/work/tic87/stellar/out/catalog_compilation.csv", index=False)
with open("/home/user/work/tic87/stellar/out/catalog_compilation.md", "w") as fh:
    fh.write("| " + " | ".join(df.columns) + " |\n|" + "---|" * len(df.columns) + "\n")
    for _, r in df.iterrows(): fh.write("| " + " | ".join(str(x) for x in r.values) + " |\n")
print(df[["source", "uses_plx", "Teff", "logg", "feh", "R", "M", "d_pc", "ext"]].to_string())
