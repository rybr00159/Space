"""Assemble out/summary_ctx.json from the outputs of c04-c17 plus the verified literature values (c14/c14b abstracts,
ext/lit/*.txt full-text extracts) and the known-planet list (c15)."""
import json
import numpy as np, pandas as pd

fp = pd.read_csv("out/fold_products.csv")
def row(prod, star, sector, col):
    r = fp[(fp["product"] == prod) & (fp.star == star) & (fp.sector.astype(str) == str(sector)) & (fp.column == col)].iloc[0]
    f = lambda x: None if not np.isfinite(x) else round(float(x) * 1e6, 0)
    return dict(depth_ppm=f(r.depth), err_formal_ppm=f(r.err), err_empirical_ppm=f(r.err_emp_std), snr_empirical=round(float(r.snr_emp), 2),
                odd_ppm=f(r.odd), odd_err=f(r.odd_err), even_ppm=f(r.even), even_err=f(r.even_err),
                phase05_ppm=f(r.half), phase05_err=f(r.half_err), n_in_transit_points=int(r.n_in),
                p_random_phase_ge=round(float(r.p_rand_ge), 4) if np.isfinite(r.p_rand_ge) else None)
products = {
    "QLP_M3V_DET_ALL": row("QLP", "M3V", "ALL", "DET_FLUX"), "QLP_M3V_SAP_ALL": row("QLP", "M3V", "ALL", "SAP_FLUX"),
    "QLP_Mira_DET_ALL": row("QLP", "Mira", "ALL", "DET_FLUX"), "QLP_pairA_DET_ALL": row("QLP", "pairA", "ALL", "DET_FLUX"),
    "QLP_M3V_DET_S12": row("QLP", "M3V", 12, "DET_FLUX"),
    "TGLC_M3V_psf_S39": row("TGLC", "M3V", 39, "psf_flux"), "TGLC_M3V_aper_S39": row("TGLC", "M3V", 39, "aperture_flux"),
    "TGLC_M3V_psf_S12": row("TGLC", "M3V", 12, "psf_flux"),
    "TGLC_pairA_psf_S39": row("TGLC", "pairA", 39, "psf_flux"), "TGLC_pairB_psf_S39": row("TGLC", "pairB", 39, "psf_flux"),
    "TGLC_Mira_psf_S39": row("TGLC", "Mira", 39, "psf_flux"), "TGLC_ref1_psf_S39": row("TGLC", "ref1", 39, "psf_flux"),
    "TARS_M3V_ALL": row("TARS", "M3V", "ALL", "FLUX"),
    "ELEANORLITE_S12_PCA_blend": row("ELEANOR-LITE", "M3V", 12, "PCA_FLUX"), "ELEANORLITE_S12_CORR_blend": row("ELEANOR-LITE", "M3V", 12, "CORR_FLUX"),
    "TEQUILA_S12_Mira_pixel": row("TEQUILA", "Mira_pixel", 12, "FLUX"),
}
leak = json.load(open("out/tglc_leak.json")); stars = json.load(open("out/tglc_leak_stars.json"))
rob = json.load(open("out/tglc_leak_robust.json")); shape = json.load(open("out/tglc_psf_shape.json"))
tars = json.load(open("out/tars_pixel.json"))
gaia = json.load(open("out/gaia_in_transit.json")); gval = json.load(open("out/gaia_scanlaw_validation.json"))
fut = json.load(open("out/future_tess.json")); futx = json.load(open("out/future_tess_crosscheck.json"))
fu = json.load(open("out/followup.json")); neo = json.load(open("out/neowise_fold.json"))
cat = json.load(open("out/catalog_status.json")); known = json.load(open("out/known_mdwarf_neptunes.json"))
s39 = leak["psf_flux"]["S39"]
summary = {
 "target": "TIC 214311908 = Gaia DR3 5925373272850784768 (M3V, 64.5 pc); ephemeris P = 5.6914968 d, T0 = BTJD 3840.9153",
 "1_independent_light_curves": {
   "products_found": {"QLP": "S12,39,66,93,103,104 for TIC 214311906/08/10 (+ TIC 1515709367 = pair duplicate)",
                      "TGLC": "S12, S39 only (TGLC HLSP covers Sectors 1-55); all six Gaia stars incl. pair B and reference stars",
                      "TARS": "S12,39,66,93, TIC 214311908 only", "GSFC-ELEANOR-LITE": "S12: five TIC IDs, identical RAW/PCA flux (same aperture L_180 centred on col 1588,row 268) = blend",
                      "TEQUILA": "S12, pixel x1515 y0269 (Mira pixel)",
                      "absent": "no SPOC 2-min/20-s, TESS-SPOC, CDIPS, PATHOS, DIAmante, TASOC, T16 products (MAST queries c01/c02)"},
   "fold_results": products,
   "tglc_S39_absolute_flux_loss_e_per_s": dict(observed=s39["observed_dF"], empirical_err=s39["emp_err"]),
   "tglc_S39_leak_model_chi2_dof5": {k: round(v["chi2"], 1) for k, v in s39["hypotheses"].items()},
   "tglc_S39_leak_model_faint_stars_chi2": {"G17.5_3.1asec_NE": rob["V1_baseline"]["G17.5_3.1asec"], "G17.7_5.5asec_N": rob["V1_baseline"]["G17.7_5.5asec"]},
   "tglc_S39_free_position_best_EN_from_Mira_arcsec": stars["psf_flux_all6"]["grid_best_EN"], "tglc_S39_free_chi2_min_dof3": stars["psf_flux_all6"]["grid_chi2_min"],
   "tglc_S39_sim_injected_on_M3V": {k: v for k, v in stars["psf_flux_all6"]["sim"].items() if k != "best_EN_sample"},
   "tglc_S39_robustness_M3V_minus_best_dchi2": {k: round(v["M3V"] - v["free_best_chi2"], 1) for k, v in rob.items()},
   "tglc_S39_psf_shape_variants_M3V_minus_best_dchi2": {k: round(v["M3V"]["chi2"] - v["free_best_chi2"], 1) for k, v in shape.items() if isinstance(v, dict) and "free_best_chi2" in v},
   "tglc_S12": "uninformative (30-min cadence; Mira/background systematics dominate; null chi2 21.9/6)",
   "tars_implied_total_dF_e_per_s": {h: {s: [round(tars[s]["implied_total_dF_if"][h], 1), round(tars[s]["implied_total_dF_err_if"][h], 1)] for s in tars} for h in ("M3V", "pairA", "Mira")},
   "tars_constancy_chi2_S39_S66_S93_dof2": {"M3V": 2.0, "pairA": 10.4, "Mira": 20.0},
 },
 "2_other_archival_time_series": {
   "NEOWISE": dict(n=neo["neowise"]["n"], W1_scatter=neo["neowise"]["W1"]["scatter_frac"], W2_scatter=neo["neowise"]["W2"]["scatter_frac"], n_in_transit=3,
                   W1_depth=[neo["neowise"]["W1"]["depth"], neo["neowise"]["W1"]["depth_err"]], W2_depth=[neo["neowise"]["W2"]["depth"], neo["neowise"]["W2"]["depth_err"]],
                   verdict="not usable: 5-8% scatter, halo-flagged by the Mira (cc_flags hh00), 3 in-transit points"),
   "AllWISE_MEP": "28 epochs (2010), none in transit",
   "Gaia_DR3_epoch_photometry": "exists for the Mira only (has_epoch_photometry=False for the M3V, the pair and the faint neighbours)",
   "Gaia_DR4_Dec2026": dict(n_fov_transits=gaia["DR4"]["n_fov_transits"], in_transit=0, nearest_hours_from_mid=6.7,
                            verdict="DR4 epoch photometry (data to 2020-01-20) contains no in-transit epoch"),
   "Gaia_DR5": dict(n_fov_transits=gaia["DR5"]["n_fov_transits"], expected_in_transit=round(gaia["DR5"]["expected_in_transit"], 2), epochs=gaia["DR5"]["in_or_near_transit"],
                    caveat="commanded scan law; in the DR3 window 19 of 58 predicted crossings of the Mira have no published epoch"),
   "scanlaw_validation": gval,
   "ASAS-SN": "not usable/not accessible: Sky Patrol v2 API timed out (port 9006), service reported down 13 Oct 2026; 8-arcsec pixels blend the Mira (17.8 arcsec, V 12.9-17.0); M3V V~15.6-15.9 gives 0.05-0.1 mag per epoch, ~1% duty cycle -> depth error >~2%",
   "ATLAS": "forced photometry needs an account (not accessed); pixels 1.86 arcsec would resolve the Mira",
   "SuperWASP": "not usable: 13.7-arcsec pixels, 48-arcsec aperture blends Mira+pair+M3V; V 15.6 below useful range",
   "KELT": "not usable: 23-arcsec pixels; target far below KELT's V 8-11 range",
   "MEarth-South": "not a plausible MEarth target (64.5 pc, 0.38 Rsun); target list not checked",
   "ZTF/Pan-STARRS": "Dec -51 outside footprint",
 },
 "3_future_data": {
   "TESS": dict(tess_point_version=fut["tess_point_version"], observed_sectors=[12, 39, 66, 93, 103, 104],
                future_sectors_through_134=[], last_planned_sector_end="2028 Sep 13 (S134)",
                closest_camera_centre_deg=dict(sorted(futx["closest_camera_centre_deg_105_134"].items(), key=lambda kv: kv[1])[:3]),
                heasarc_vs_tess_point_pointing_max_diff_deg=futx["pointing_comparison_max_diff_deg"],
                verdict="no TESS data before at least Sept 2028; any later sector depends on a not-yet-scheduled Cycle 11+"),
   "Gaia_DR4": "Dec 2026: no in-transit epoch (see 2); epoch astrometry/photometry baseline only",
   "Gaia_DR5": "~2030: 1 full + 1 partial (+1 marginal) in-transit epoch expected",
 },
 "4_catalog_status": dict(exofop={k: {kk: vv for kk, vv in v.items() if kk in ("tois", "ctois", "planet_parameters", "files", "imaging", "spectroscopy", "time_series")} for k, v in cat["exofop"].items()},
                          exofop_obsnotes="'Matched Notes for TIC ID ...' empty for all three (checked 2026-10-17)",
                          toi_table_rows=cat["toi_table"]["n_rows"], nearest_toi_arcmin=round(cat["toi_table"]["min_sep_arcmin"], 1),
                          ctoi_table_rows=cat["ctoi_table"]["n_rows"], nearest_ctoi_arcmin=round(cat["ctoi_table"]["min_sep_arcmin"], 1),
                          simbad="Mira = 2MASS J17235464-5119254 (LP*); M3V = Gaia DR3 5925373272850784768 (PM*); refs only Gaia catalogue papers (GCNS 2021, DR3 LPV 2023)",
                          vsx="only the Mira: ASASSN-V J172354.64-511925.3, type M:, V 12.91-17.01, P 294 d",
                          eb_catalogs="none (no TESS EB catalogue, Gaia DR3 vari or VSX entry for the M3V or the pair)",
                          arxiv="no title/abstract/comment hit for any TIC, Gaia or 2MASS identifier; ADS full-text search not available (no token)",
                          wide_binary="El-Badry, Rix & Heintz 2021 catalogue lists M3V + TIC 214311910 (Gaia DR3 ...782208) as a wide pair",
                          xray="eRASS1 1eRASS J172352.9-511928, 3.6 arcsec from the M3V, coronal (p=0.995), M3V most likely counterpart (p_ij=0.73); F_X=4.9e-14 cgs -> L_X=2.5e28 erg/s, log L_X/L_bol=-3.4 (near saturated; pair 11.7 arcsec away could contribute)"),
 "5_priors": {
   "occurrence_4to7Re_Mdwarfs": [
     "Dressing & Charbonneau 2015 (ApJ 807, 45; arXiv:1501.01623), Kepler early M (<4000 K): bins stop at 4 Re; 3.5-4.0 Re: 0.36 (+0.72/-0.18)% at P 1.7-5.5 d and 0.87 (+1.65/-0.43)% at 5.5-18.2 d per star; text: lack of detections where completeness is high 'implies that hot mini-Neptunes and Neptunes are rare around low-mass stars'",
     "Gillis, Cloutier, Pass et al. 2026 (arXiv:2602.23364), TESS mid-to-late M (<0.4 Msun, median 0.22): no Neptunes detected; 95% upper limit 0.0341 planets per star for 3-6.5 Re, P 0.2-30 d; sub-Neptunes (1.8-3 Re) 0.148 +- 0.045",
     "Ment & Charbonneau 2023 (AJ 165, 265; arXiv:2302.04242): 1-sigma upper limit 0.07 planets > 1.5 Re per mid-to-late M dwarf at P 0.5-7 d",
     "Hardegree-Ullman et al. 2019 (AJ 158, 75; arXiv:1905.05900): 0.5-2.5 Re, P 0.5-10 d: 0.86 (+1.32/-0.68) planets per M3V star (small planets only)",
     "Hardegree-Ullman et al. 2025 (Scaling K2 VIII, arXiv:2508.05734): short-period sub-Neptune occurrence peaks at Teff 3750 (+153/-97) K and drops for cooler M dwarfs",
     "Gan, Wang, Wang et al. 2023 (AJ, doi:10.3847/1538-3881/ac9b12; arXiv:2210.08313): hot Jupiters (7 Re-2 RJ, 0.8-10 d) around early M (0.45-0.65 Msun) 0.27 +- 0.09%",
     "Bryant, Bayliss & Van Eylen 2023 (MNRAS, doi:10.1093/mnras/stad626; arXiv:2303.00659): giants 0.6-2 RJ, P<10 d: 0.108 +- 0.083% for 0.26-0.42 Msun; 0.194 +- 0.072% for 0.088-0.71 Msun",
     "Sabotta et al. 2021 (A&A 653, A114; arXiv:2107.03802), CARMENES RV: 10-100 Mearth planets 0.18 (+0.07/-0.05) per star (P<1000 d); 'planets more massive than 10 Mearth become rare' for M < 0.34 Msun",
     "Powers et al. 2023 (AJ 166(2), doi:10.3847/1538-3881/acd8bf; TOI-3785 b): 'rare population of Neptunes (4 < Rp < 7 Re) orbiting ... M-dwarf host stars, of which only ~10 have been confirmed'",
     "Adopted for FPP (our synthesis, not a published number): f(4-7 Re, P 3-10 d | M3V, 0.37 Msun) ~ 0.1-1% per star, best guess ~0.3%; i.e. between the giant-planet rate (0.1-0.3%) and the 3.5-4 Re Kepler rate (~1% over 1.7-18 d)"],
   "known_transiting_4to7Re_P_lt_10d_Teff_lt_4000K": dict(n_total=known["n_4to7_P10_teff_lt_4000"], n_teff_lt_3500=known["n_4to7_P10_teff_lt_3500"],
                                                           n_mass_lt_0p45=known["n_4to7_P10_mass_lt_0p45"], list="out/known_mdwarf_neptunes.csv",
                                                           note="only TOI-3884 b (3180 K, 0.30 Msun, 6.4 Re, 4.54 d) has Teff < 3500 K; TOI-674 b (0.42 Msun) is the only other host < 0.45 Msun. The candidate (3300 K, 0.37 Msun, 5-6 Re, 5.69 d) would be the second such planet."),
   "multiplicity": ["Winters et al. 2019 (AJ, doi:10.3847/1538-3881/ab05dc): M-dwarf multiplicity 26.8 +- 1.4%, companion rate 32.4 +- 1.4%; primaries 0.30-0.60 Msun: 28.2 +- 2.1% (uncorrected, 25 pc); projected separations log-normal, peak 20 au, sigma_log a = 1.16 (fraction of companions within 65 au = 1 arcsec at 64.5 pc: ~0.67 -> ~19% of M dwarfs with a companion inside 1 arcsec; our calculation from their fit); 37 triples among 1120 primaries (3.3%)",
                    "Shan, Johnson & Morton 2015 (ApJ 813, 75; arXiv:1509.04189): 0.11 (+0.02/-0.04) close (P = 1-90 d) M-dwarf companions per apparently single M dwarf from Kepler EBs; power-law period distribution exponent 0.52 +- 0.22",
                    "Moe & Kratter 2021 (MNRAS, doi:10.1093/mnras/stab2328): binaries with a < 1 au fully suppress S-type planets; a = 10 au at 15%; a > 200 au negligible (the M3V's known companions are at ~725 au)"],
   "eb_geometry": "M+M binary at P = 5.69 d (0.3 Msun total): a ~ 0.042 au, eclipse probability (R1+R2)/a ~ 0.04; planet transit probability R*/a = 0.38 Rsun / 0.045 au ~ 0.039",
 },
 "6_followup": dict(csv_check=fu["csv_check"], n_fully_observable_Oct2026_Dec2027=fu["epochs_observable_from_any_site"], per_site=fu["n_observable"],
                    season=fu["season_months_observable"], late_2026="none observable (best: 2026-11-02 23:04 UTC, 16% of the window from CTIO)",
                    list="out/followup_transits.csv",
                    precision_per_event_snr_depth1p5pct={k: dict(nominal_red=v["snr"]["0.015"]["snr_red"], poor_red=v["snr_poor"]["0.015"]["snr_red"], per_exposure_ppt=v["nominal"]["per_exposure_ppt"]) for k, v in fu["precision"].items()},
                    mira_contamination="Moffat beta 2.5, FWHM 2 arcsec, r=3 arcsec: Mira adds 4e-4 (i) to 2.6e-3 (i, Mira at max) of the M3V flux; FWHM 3 arcsec: up to 2-4% at max in i/z (beta 2.5) -> use FWHM < 2.5 arcsec or bluer bands; Mira core saturates in i/z in 60 s on a 1 m (peak ~2e5 e-)"),
}
json.dump(summary, open("out/summary_ctx.json", "w"), indent=1, default=float)
print(json.dumps(summary, indent=1, default=float)[:6000])
