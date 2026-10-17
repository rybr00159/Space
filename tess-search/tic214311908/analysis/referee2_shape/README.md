# Referee 2: independent check of the 'shape' analysis of TIC 214311908

Folder: `tic214311908/analysis/referee2_shape/`
(AI-run analysis, 2026-10-17. Nothing was submitted, uploaded or posted anywhere; no one was contacted. The only
network use was read-only Gaia archive queries: the target's astrometric record and XP sampled spectra.)

This is an adversarial re-derivation of `../shape/` (transit model, stellar-density test, dilution limit,
false-positive fits, secondary eclipse). All code is my own, written in `scripts/`. The only things reused are the
phase-1 light curves, the referee's fake-window noise ratios (`../referee_phot/out/per_sector.csv`), the comparison-M-dwarf
table `../referee_star/data/mdw_cmd.csv` and the TESS response curve file `../shape/data/filt_TESS.Red.dat` (data only).

## Verdict in one paragraph

The planet fit, the derived parameters, the uninformative density test, the dilution profile and the
bound/brown-dwarf exclusions **reproduce**. Two parts of the white-dwarf (WD) exclusion **do not reproduce as stated**:
* the Gaia RUWE likelihood is computed with the wrong geometry and an arbitrary per-CCD error;
* the secondary-eclipse error is underestimated by ~15%.

With both corrected, light curve + secondary + RUWE disfavour a WD eclipsed by the M3V by a minimum dchi2 of
**~9-11, or ~10-12 with my new Gaia XP test**, instead of 13.2-15.2. If one insists on the analyst's own
"conservative" 0.35 mas per-CCD error, the minimum drops to only **6-8**.

The WD case is nevertheless strongly disfavoured by a test the analyst set aside as qualitative, the Gaia RVS
radial-velocity error. A K1 >= 35 km/s orbit would smear the combined CCF and inflate the error by >= 2-4x, giving
-2 ln L >~ 10-15. Its error model is approximate.

On the planet side, the calibrated shape statistics are stronger than the analyst says:
* a V-shaped (grazing) truth reproduces the observed flat-bottom preference in 0/25 REF and 1/25 INV injections;
* the WD-occultation shape is disfavoured by a likelihood ratio of ~4-6 in each photometry (not independent).

## Data and noise (rlib.py, r01)

* INV = `../phot/lc_M3V_all.csv` (quality_ok, no flares, rescaled by 935/979), windows +-0.30 d.
* REF = `../referee_phot/lc_M3V_referee_windows.csv`.
* Per-point sigma = err x per-sector robust scale (window OOT) x beta, with per-window offsets profiled.
* **INV beta re-derived (r01)** from 62 fake ephemerides with the raw flux re-detrended and the fake windows masked
  (single-window depth scatter / formal): 1.15, 1.08, 1.50, 1.48, 1.48, 1.25 (S12...S104), `out/beta_inv.json`.
  The analyst's INV betas (1.00, 1.08, 1.24, 1.21, 1.27, 1.13), from the time-averaging method, are 0-22% lower than
  this transit-timescale calibration. The INV dchi2 values of `../shape` are therefore overstated by ~1.2-1.5x in the
  sectors that dominate (S66-S104; chi2 scales as 1/beta^2). This is why my INV dilution profile never reaches dchi2 = 9.
* REF beta: the referee fake-window ratio (as the analyst).
* **Exposure integration**: batman's `supersample_factor` uses linspace(-0.5, 0.5, n) offsets, i.e. an effective
  exposure n/(n-1) too long (1.25x for the 200-s data, n = 5; demonstrated in r03a).
  * I integrate with the midpoint rule.
  * My own uniform-disk occultation code (`rwd.py`) agrees with batman to 1e-11 instantaneous, `out/validate_rwd.json`.
  * The effect on the analyst's results is small (it smooths the WD box slightly, in the WD's favour).

## Results (my value vs the analyst's)

| Quantity | ../shape | referee 2 (method, file) | verdict |
|---|---|---|---|
| Rp/R* free INV / REF | 0.119 (+0.008/-0.005) / 0.119 +- 0.005 | 0.1203 (+0.0093/-0.0063) / 0.1198 (+0.0061/-0.0048) (emcee 32x8000, r02, `out/mcmc_*_free_fake.json`) | upheld |
| rho_circ INV / REF (g/cm3) | 10.3 (+5.4/-5.9) / 12.3 (+4.0/-5.8) | 9.0 (+5.8/-5.5) / 11.6 (+4.5/-6.0) | upheld |
| T14 INV / REF | 1.69 (+0.17/-0.12) / 1.65 (+0.13/-0.09) h | 1.72 (+0.21/-0.14) / 1.66 (+0.15/-0.09) h | upheld |
| grazing posterior fraction | 4.9% / 2.2% | 7.3% / 1.5% | upheld |
| ln B single/twin | -0.43 / -0.53 | -0.16 / -0.60 (rho-only marginal likelihood, r09) | upheld (no preference) |
| Rp single / twin | 4.94 +- 0.26 / 5.23 +- 0.25 R_E | 4.94 +- 0.27 / 5.25 +- 0.27 R_E (r09, `out/derived.json`) | upheld |
| b, a/R*, i, T14 (single) | 0.53 +- 0.10, 25.3 +- 0.8, 88.81 +- 0.24, 1.70 +- 0.07 | 0.52 +- 0.11, 25.4 +- 0.8, 88.82 +- 0.26, 1.70 +- 0.08 | upheld |
| Teq single / twin | 471 +- 10 / 420 +- 9 K | 470 +- 11 / 420 +- 9 K | upheld |
| depth (model, mid-transit) | 15.7-15.8 ppt | 15.9 (14.7-17.3) INV, 16.1 (15.0-17.1) REF | upheld |
| V-shape - planet dchi2 | +3.0 / +5.4, "only moderately supported" | +2.8 / +5.7. **Calibrated (r07b, 25 shift injections of the V-shaped ML model): V truth gives 0.43 +- 1.26 (INV, 1/25 >= 2.76) and 0.66 +- 0.94 (REF, 0/25 >= 5.68, max 4.2); planet truth gives 4.5 +- 3.3 / 5.4 +- 3.1** | number upheld; interpretation revised: a grazing / V-shaped source is disfavoured at p ~0.04 (INV) and p < 0.04 (~5 sigma Gaussian, REF) |
| WD (box) - planet, LC only | +0.7 / +2.9 (s06), +0.79 / -0.33 (s13b) | +1.1 (INV) / +1.5 (REF), cool branch (r03) | upheld |
| shape calibration of dchi2(WD - planet) | INV only: planet +0.8 +- 1.4, WD -2.9 +- 1.6 (8 epochs); "REF none" | 30 shift injections per truth (r07): INV planet +1.8 +- 2.2, WD -2.3 +- 1.6, real +1.1; REF planet +3.9 +- 2.1, WD -2.6 +- 2.0, real +1.5; 0/30 WD truths reach the real value in either. LR(planet/WD) = 6.5 / 4.3 (1-D), 4.4 / 5.9 (2-D with dchi2(V - planet)) | upheld, extended to REF; the shape gives a modest LR ~4-6, not ~15 |
| secondary at phase 0.5 | -0.16 +- 0.79 ppt (INV, delivered flux) | **+0.20 +- 0.90 ppt** (INV, re-detrended with the window masked, empirical std of 62 phases; beta-weighted: -0.12 +- 0.88), 2-ppt injection recovered 100% (r04; also 100% without the mask, r04b) | revised: error +14%; the value difference is within noise (template / weighting) |
| 3-sigma secondary limit / T_b | 2.2 ppt / < 2460 K | **3.0 ppt (white weights) or 2.5 ppt (beta weights) / < 2540 K** (r14) | revised; M8V (2570 K) only marginally excluded |
| k on the M3V with density prior | k=0.14 +12.2, 0.17 +38.9, 0.20 +45.3 (INV) | 0.14: +14.9/+19.8, 0.17: +40/+55, 0.20: +44/+61, 0.25: +47/+65 (INV/REF, r14) | upheld (BD on M3V excluded) |
| shape dilution limit, transit geometry | dchi2 = 9 at f_src 0.031-0.04 (dT 3.5-3.8) | REF: +7.1 at dT 3.37, +10.6 at dT 3.64 -> dT_max ~3.5 (3 sigma). **INV: never reaches 9 (+8.3 at dT 4.1) with the re-calibrated INV noise** (r05) | partly revised: INV gives no 3-sigma limit |
| occultation geometry (k > 1) | +1.8 to +7.8 at any f_src | INV +2.6 to +3.7, REF +5.7 to +7.0 | upheld |
| per-sector depth chi2 | 8.4/5 INV, 3.4/5 REF | 7.7/5, 3.4/5 (r12) | upheld |
| odd - even | INV +5.0 +- 2.1, REF +1.3 +- 1.9 ppt | INV +4.7 +- 2.2, REF +1.3 +- 1.9 | upheld |
| late (S93-104) - early (S12-66) depth | "1.5-2 sigma" | INV +5.5 +- 2.2 ppt (2.5 sigma), REF +2.8 +- 1.9 (1.5 sigma); post-hoc split | flag |

## The white-dwarf scenario in detail (r03, r04, r06, r06c, r08, r10, r11, r13)

**Light curve + secondary + M-dwarf priors** (r10, WD radius on the cool branch, joint fit; own secondary
+0.20 +- 0.90 ppt).
* The minimum is **+2.2 (INV) and +2.6 (REF), at M_WD ~0.45**, where the predicted secondary is 0.9 ppt.
  The analyst found +3.1 at 0.45.
* Predicted secondary, including LD and self-lensing (rwd.secondary_depth): 3.3, 2.4, 1.9, 1.6, 0.9, 0.2 ppt at
  M_WD = 0.17, 0.25, 0.30, 0.35, 0.45, 0.60. This agrees with the analyst.

**Gaia RUWE** (observed 1.076; astrometric_excess_noise 0.036 mas at 0.5 sigma; 37 transits, 328 CCD obs).
Problems in `../shape/scripts/s16_wd_checks.py`:
1. Random orbital inclinations. The system eclipses, so the orbit is edge-on, and a 1-D wobble can hide along the
   cross-scan direction.
2. Uniform scan angles and uniform random epochs. This source's scans are strongly anisotropic:
   scan_direction_strength_k2 = 0.41. Gaia epochs also come in clumps.
3. 0/400 passing simulations are converted into a fixed -2 ln L = 12.3 through a 1e-3 floor.
4. sigma_CCD = 0.22 / 0.35 mas is assumed.

What I did instead:
* **Scanning law.** I built a nominal-scanning-law model (r06). With precession phase 60 deg it reproduces the
  catalogue scan statistics: k1..k4 strengths 0.51/0.44/0.32/0.23 vs 0.50/0.41/0.32/0.26, mean angles -164/8/-42/37
  vs -163/10/-34/37 deg.
* **Per-CCD error.** The five catalogue astrometric errors through this geometry imply a per-CCD sigma of
  **0.19-0.25 mas** (r06b, `out/sigma_ccd_check.json`).
* **Validation.** r06c reproduces the analyst's numbers under the analyst's assumptions: 6.8 vs 6.4 at 0.3 Msun,
  sigma 0.35. With the edge-on orbit and the NSL scans the same case gives 2.9.

-2 ln L(RUWE) (r06, density at 1.076, edge-on, NSL, WD G-band share 2.1%):

| M_WD | 0.2 | 0.3 | 0.45 | 0.6 |
|---|---|---|---|---|
| sigma 0.18 (from parallax_error) | 4.6 | 8.2 | 14.4 | 18.5 |
| sigma 0.22 | 3.0 | 5.9 | 9.7 | 15.4 |
| sigma 0.25 | 2.1 | 4.5 | 7.8 | 11.4 |
| sigma 0.35 | 0.8 | 2.1 | 4.0 | 5.9 |
| ../shape (sigma 0.35 / 0.22) | 2.0 / 8.0 | 6.4 / 12.3 | 12.3 (floor) | 12.3 (floor) |

**New test, Gaia XP blue excess** (r08; 228 comparison M dwarfs with BP-RP 2.66-2.86).
* A WD supplying 1.57% of the TESS light at the Teff it needs (4800-6100 K for 0.25-0.6 Msun) supplies 8-13% of the
  420-540 nm flux.
* The blue index is regressed on red-only indices, which the WD barely changes. The target sits at +0.017 (71st
  percentile, robust sd 0.032).
* dchi2(WD - planet) = +0.8, +2.8 and +6.7 at 0.3, 0.45 and 0.6 Msun.
* Regressing on BP-RP instead absorbs the WD signal: the WD makes BP-RP bluer too, giving dchi2 -1 to -1.6. Use the
  red-only version.

**New semi-quantitative test, Gaia RVS error** (r13).
* For grvs_mag 12.55 the DR3 RV and its error come from the summed CCF of 21 transits.
* With K1 = 35-73 km/s (0.2-0.6 Msun) the Zucker (2003) error grows by 2.1-3.0x (5th percentile of random phases)
  to 2.9-4.4x (median).
* The intrinsic error would then have to be <~1.0 km/s, against a minimum of ~1.3 km/s for 1642 M dwarfs of the
  same grvs_mag. That is -2 ln L ~ 10-16 at the most favourable smearing.
* Caveats: the Gaia error model, the CCF width (12-16 km/s assumed), and KDE tail extrapolation.

**Combined minimum over M_WD** (`out/wd_combined.json`, own secondary; with the REF r18 secondary +0.4 to +0.8 higher):

| per-CCD sigma | LC+sec+RUWE INV / REF | + XP |
|---|---|---|
| 0.18 | 13.2 / 13.7 | 14.0 / 14.5 |
| 0.22 | 10.9 / 11.3 | 11.7 / 12.1 |
| 0.25 | 9.3 / 9.8 | 10.3 / 10.8 |
| 0.35 | 5.9 / 6.3 | 7.6 / 8.1 |

Adding the RV-smearing term (>= ~12 at the minimum, 0.3-0.45 Msun) gives >= ~20-25.

## Missed tests and recommendations

1. The Gaia RVS error (above) and the Gaia XP blue excess (above) should enter the WD likelihood.
2. Calibrate flat-vs-V with a V-shaped truth (done here: p ~0.04 INV, p < 0.04 REF). The raw dchi2 of 3-5 understates
   it, because a V-shaped truth gives dchi2(V - planet) of only 0.4-0.7 +- 1.
3. Decisive follow-up for the remaining WD / blend space:
   * multi-colour photometry: a WD occultation is ~6-8x deeper in g than in TESS, while a planet is achromatic;
   * one high-cadence transit: 2-3 min ingress (WD) vs ~14 min (planet);
   * two low-resolution RVs: K1 35-106 km/s vs < 0.01 km/s.
4. The bound-EB grid (../shape s08) has no occulter larger than 0.105 Rsun and fixed LD. It never tests a giant
   planet on a bound 0.2-0.3 Msun companion, which still makes a planet (Rp ~7 R_E), not a false positive.
5. The late-vs-early depth difference (2.5 sigma in INV) should be checked against the Mira's phase and the
   deblending templates.

## Scripts (`scripts/`; /home/user/work/venv/bin/python, OMP_NUM_THREADS=1)

| Script | Purpose | Output |
|---|---|---|
| rlib.py | data loading, noise model, batman planet model with midpoint exposure integration, profiled offsets, Kepler/density helpers | (library) |
| rwd.py | uniform-source occultation by a large disk (own code), exposure integration, WD secondary with LD + microlensing | (library) |
| r00_gaia.py | Gaia DR3 astrometric record of the target | out/gaia_astrom.json |
| r01_beta_inv.py | INV red-noise factors from 62 re-detrended fake ephemerides | out/beta_inv.json, out/beta_inv_rows.csv |
| r02_planet_mcmc.py | emcee fits free / single / twin, both photometries | out/mcmc_{ds}_{model}_fake.{json,npz}, out/r02_*.log |
| r03a_validate.py | own occultation model vs batman; batman supersampling quirk | out/validate_rwd.json |
| r03_scenarios.py | ML fits: planet free/single/twin, V-shape, WD cool/hot/twin, WD mass profile (LC only) | out/scen_{ds}_fake.json |
| r04_secondary.py | phase-0.5 depth with re-detrending, empirical errors, injection check | out/secondary.json, out/secondary_rows.csv |
| r04b_unmasked_bias.py | recovery of a 2-ppt secondary when phase 0.5 is NOT masked in the detrending (as in the delivered flux): 100%, so ../shape s10 is not biased by this | out/secondary_unmasked_bias.json |
| r05_dilution.py | shape-only dilution profile, free LD, both geometries | out/dilution_{ds}.json |
| r06_ruwe.py | NSL scan model fitted to the catalogue scan statistics; per-CCD sigma; RUWE likelihood vs M_WD (edge-on) | out/ruwe_nsl_{plx,0.22,0.25,0.35}.json |
| r06b_sigma_check.py | per-CCD sigma implied by each of the five catalogue errors | out/sigma_ccd_check.json |
| r06c_reproduce_analyst.py | analyst's RUWE assumptions reproduced, then changed one at a time | out/ruwe_reproduce.json |
| r07_inject.py | shift (INV, REF) injection calibration of dchi2(WD - planet), dchi2(V - planet) | out/inject_{ds}_shift.{csv,json} |
| r07b_inject_vshape.py | V-shaped truth calibration of dchi2(V - planet) | out/inject_{ds}_vshape.{csv,json} |
| r08_xp.py | Gaia XP blue-excess test (downloads XP sampled spectra) | out/xp_test.json, data/xp_sampled.csv |
| r09_derived.py | Vanderburg T12/T13 dilution bound, rho-only Bayes factor, adopted parameters | out/derived.json |
| r10_wd_total.py | WD: joint LC + secondary profile in M_WD, combined with RUWE | out/wd_total_{ds}.json |
| r11_combine.py | WD: LC+sec + RUWE (4 sigma choices) + XP + RV table | out/wd_combined.json |
| r12_per_sector.py | per-sector / odd-even / early-late depth, per-sector duration | out/per_sector.json |
| r13_rv_smear.py | Gaia RVS combined-CCF error under orbital smearing | out/rv_smear.json |
| r14_kprofile_tb.py | k-profile on the M3V with the density prior; occulter brightness-temperature limit | out/kprofile_tb.json |
| r15_summary.py | collects the key numbers | out/summary.json |

Logs: `out/*.log`. A sign error in my first r04 version (secondary reported with the opposite sign) was found and
fixed before r10/r11/r14 were re-run. The superseded r10 outputs were deleted.
