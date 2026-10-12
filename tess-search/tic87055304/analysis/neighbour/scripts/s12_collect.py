"""Collect all numbers of the neighbour/ investigation into results.json (re-run whenever an input changes).
Inputs: results/*.json, mcmc/mcmc_*.json, tri/tri_*.json, data/novelty/exofop_87055315.json."""
import os, json, glob
import numpy as np
ROOT = "/home/user/work/tic87/neighbour"
J = lambda p: json.load(open(os.path.join(ROOT, p)))
st = J("results/stellar.json"); fr = J("results/fluxratio.json"); pt = J("results/pair_tests.json")
co = J("results/companion.json"); eb = J("results/eb2p.json"); nov = J("results/novelty.json"); rv = J("results/gaia_rv.json")
ep = J("results/ephemeris.json"); d1 = J("results/diag_t14.json"); d2 = J("results/diag_t14b.json")
mL = J("mcmc/mcmc_pair_local.json"); mB = J("mcmc/mcmc_pair_Dfree.json")
m1 = json.load(open("/home/user/work/tic87/transit/mcmc_nbhost.json"))


def pm(p, keys=("med", "lo", "hi")):
    return {k: p[k] for k in keys}


tri = {}
for f in sorted(glob.glob(os.path.join(ROOT, "tri/tri_*.json"))):
    r = json.load(open(f))
    tri[os.path.basename(f)[4:-5]] = r
res = dict(
    object=dict(TIC=87055315, gaia_dr3=6720986208895208960, twomass="J18102616-4359369", ra_2016=272.609023, dec_2016=-43.993633,
                sep_from_TIC87055304_arcsec=9.33, PA_deg=45, G=14.0234, BP_RP=0.9276, Tmag_TIC=13.5466),
    stellar=dict(adopted=st["adopted"], distance=st["distance"], sed_fit={k: v for k, v in st["sed_fit"].items() if k != "residuals_best"},
                 mist=st["mist"], catalogue_consensus=st["catalogue_consensus"], catalogues=st["catalogues"],
                 gaia_dr3_source=st["gaia_dr3_source"], gaia_rv_error_test=rv,
                 script="scripts/s04_stellar.py, scripts/s07_gaia_rv.py"),
    flux_ratio=dict(r_stassun_dr3=fr["A1_stassun_gaiadr3"]["r"], r_tic82=fr["A2_tic82"]["r"], r_xp_sampled=fr["A3_xp_sampled"]["TESS"]["r"],
                    r_adopted=st["adopted"]["r_TESS"], r_err=st["adopted"]["r_TESS_err"], T_synthetic_SED=st["sed_fit"]["T_synthetic"],
                    D_pair=st["adopted"]["D_pair"], D_pair_err=st["adopted"]["D_pair_err"],
                    D_qlp_empirical_S13_S93=[fr["C_empirical"]["D_eff_S13_S93"], fr["C_empirical"]["D_eff_S13_S93_err"]],
                    D_qlp_empirical_S104=[fr["C_empirical"]["D_eff_S104"], fr["C_empirical"]["D_eff_S104_err"]],
                    D_qlp_PRF_predicted={s: v["apertures"][[k for k in v["apertures"] if k.startswith("qlp")][0]]["D_corr"] for s, v in fr["B_qlp_apertures"].items()},
                    script="scripts/s03_fluxratio.py"),
    transit_fit=dict(
        nominal_pair_local=dict(params={k: pm(mL["params"][k]) for k in ["k", "b", "rho", "a", "inc", "T14_h", "T23_h", "depth_ppm", "depth_undiluted_ppm", "D", "P", "T0", "u1", "u2"]},
                                n_transits=mL["n_transits"], n_over_tau=mL["n_over_tau"], script="scripts/nb_mcmc_local.py",
                                note="PSF pair LC (F_t+F_n), per-transit quadratic baselines profiled, D ~ N(5.32, 0.06)"),
        alt_pair_biweight=dict(params={k: pm(mB["params"][k]) for k in ["k", "b", "rho", "a", "T14_h", "T23_h", "depth_ppm", "depth_undiluted_ppm", "D"]},
                               script="scripts/nb_mcmc.py (nbmodel.py)", note="same LC, 1-d biweight detrended (s02_lc.py)"),
        first_round_qlp_sap=dict(params={k: pm(m1["params"][k]) for k in ["k", "b", "rho", "a", "T14_h", "T23_h", "depth_ppm"]},
                                 file="/home/user/work/tic87/transit/mcmc_nbhost.json", note="QLP SAP S13-S93, D = 4.35 fixed"),
        ml_diagnostics_biweight=d1, ml_diagnostics_local_baseline=d2,
        psf_photometry_depths=dict(pair_ppt=[3.44, 0.14], nb_only_cat_norm_ppt=[17.9, 1.6], f_target=[0.16, 0.07, 0.25],
                                   source="/home/user/work/tic87/psf_photometry/README.md")),
    companion=dict(radius=co["radius"], density=co["density"], secondary=co["secondary"], mass_limits=co["mass_limits"],
                   mass_from_radius=co["mass_from_radius"], expectations_vs_mass=co["expectations_vs_mass"],
                   oddeven=co["oddeven"], eb_twin_2P=eb, red_noise_beta_pair=J("results/beta_pair.json"), pair_tests=dict(primary_nb=pt["primary_nb"], secondary_nb=pt["secondary_nb"],
                   secondary_scan=dict(max_abs_snr=pt["secondary_scan"]["max_abs_snr"], phase=pt["secondary_scan"]["phase_max"],
                                       upper2sig_any_ppt=pt["secondary_scan"]["upper2sig_any_ppt"]),
                   phase_curve=pt["phase_curve"]),
                   scripts="scripts/s06_pair_tests.py, s09_eb2p.py, s11_companion.py"),
    triceratops=dict(runs=tri, lc="data/tri_lc_nb_200s.npz (pair LC x (1+D), S66/S93/S104)", inputs=J("data/tri_inputs.json"),
                     runner="/home/user/work/src/run_tri.py (A-C); scripts/run_tri_ontarget.py (D, copy + remove_nearby_stars)"),
    novelty=dict(exofop=dict(tois=nov["exofop"]["tois"], ctois=nov["exofop"]["ctois"], imaging=nov["exofop"]["imaging"],
                             spectroscopy=nov["exofop"]["spectroscopy"], time_series=nov["exofop"]["time_series"],
                             notes="no notes field in the ExoFOP JSON"),
                 vizier={k: {kk: vv for kk, vv in v.items() if kk in ("catalogue_exists", "cone3", "cone21")} if isinstance(v, dict) else v
                         for k, v in nov["vizier"].items()},
                 gaia_dr3_archive_rows=nov["gaia_dr3_archive_rows"], simbad=nov["simbad"],
                 mast_products="only HLSP QLP and GSFC-ELEANOR-LITE for S13; no SPOC 2-min / TESS-SPOC / DV products",
                 web_search=dict(queries=['"TIC 87055315"', '"6720986208895208960"', '"2MASS J18102616-4359369"',
                                          '"Gaia DR3 6720986208895208960" OR "TIC 87055304" OR "TYC 7912-71-1"'],
                                 result="no page mentions any identifier of the star (2026-10-12)"),
                 script="scripts/s05_novelty.py"),
    ephemeris=ep,
)
json.dump(res, open(os.path.join(ROOT, "results.json"), "w"), indent=1, default=float)
print("written; TRICERATOPS runs:", {k: (v.get("FPP"), v.get("NFPP")) for k, v in tri.items()})
