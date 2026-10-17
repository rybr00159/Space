# Photometric systematics next to the Mira, cleaned light curves and final ephemeris for TIC 214311908 (`clean`)

Folder: `tic214311908/analysis/clean/`
(AI-run analysis, 2026-10-17. Nothing was posted, uploaded or submitted anywhere, and no one was contacted. No
network access was used.)

**Inputs (read only):** the TESScut cubes in `/home/user/work/tesscut5/214311906/`; the code and calibration products
of `../phot/` (photlib, p03_scene, WCS shifts, barycentric caches) and `../referee_phot/` (r05_local, WCS offsets,
Mira b images, barycentric corrections); `../search/out/gaia_variables.json` (Gaia DR3 variables in the cutout);
`../referee2_search/out/r01_rrab_leak.json` (predicted per-transit leak bias); `../referee_phot/out/invscale.json`
(per-sector PRF/TIC flux scale). Nothing was written outside this folder.

Notation: **INV** = the phase-1 phot "model B" deblending method (global pixel regression + Mira empirical image +
M3V PRF, per cadence); **REF** = the referee's local-window method (per-pixel baselines in +-0.35 d windows, PRF fit of
the residual frames). Depths are in ppt of the M3V's flux in PRF units, 979 e-/s (referee_phot r09/r17).

## Bottom line

1. **The RRab leak is removed at its source.** Excluding the catalogued Gaia variables from the field-wide estimators
   (pointing, flux scale, faint-pixel background plane and level) of the INV pipeline removes the inverted RRab signal.
   * First-harmonic amplitude at the Gaia RRab period, before -> after: 1.93 -> 0.45, 2.51 -> 0.35, 3.24 -> 0.78,
     3.79 -> 0.23 and 2.08 -> 0.69 ppt in S39, S66, S93, S103 and S104. S12 is 0.76 -> 1.05 ppt, not significant
     either way.
   * After cleaning, p_emp is 0.10-0.88 in five sectors, and the 95% upper limits are 0.75-1.8 ppt.
   * The route is confirmed as the pointing regressor. dx(t) carried the RRab period at p < 0.01 in S39, S93, S103
     and S104; after masking, p = 0.12-0.70.
   * **One residual:** S93 keeps a 0.78-ppt signal at the RRab period (p < 1/120). Its phase differs from the leak's.
     Its pixel pattern fits an M3V-PRF or Mira-pointing-dipole shape. Larger masks do not remove it.
   * **Using the variables' own light curves as explicit pixel regressors makes things worse.** It re-injects the
     RRab (1.6-2.8 ppt) and the EBs (Delta-chi2 up to 488) and raises the noise by 3-37%. It was rejected.
   * The REF photometry never carried the leak. On a contiguous REF-method curve, the equivalent amplitude is
     0.23-0.63 ppt (p >= 0.18) before and after; the transfer of 0.55 is calibrated by injection.
2. **What cleaning does to the transit** (identical per-transit fits on every version, `c11`):
   * INV per-transit depths change by 2.7 ppt rms (maximum 7.4 ppt).
   * The change tracks the leak bias predicted by referee2_search, with slope -0.82 +- 0.11 (-1 means full removal).
     The phase-1 depths regressed on that bias with slope +1.13 +- 0.33; the cleaned ones give +0.31 +- 0.35.
   * The global depth barely moves: INV 15.34 +- 1.03 -> 15.10 +- 0.89 ppt; REF 15.86 -> 15.82 ppt.
3. **Final transit parameters (cleaned, batman, single-star density prior; fake-epoch-calibrated noise):**
   * INV: depth 14.50 +- 0.98 ppt, k 0.1131 +- 0.0047.
   * REF: depth 16.06 +- 1.00 ppt, k 0.1199 +- 0.0043.
   * **Adopted depth 15.2 +- 0.9 (stat) +- 1.1 (sys) ppt**, i.e. a flux loss of 14.9 e-/s.
   * **Rp/R* = 0.117 +- 0.004 (stat) +- 0.004 (sys) and Rp = 4.90 +- 0.22 (stat, incl. R*) +- 0.17 (sys) R_earth.**
   * b 0.49-0.54, a/R* 25.3 +- 0.8, T14 1.69-1.73 h, flat fraction 0.71-0.74, Teq 471 +- 11 K. The grazing posterior
     fraction is 0 with the density prior and 0.5-1.4% with free density.
   * With free density, rho* = 8.4 (+5.0/-4.4) g/cm3 (INV) and 12.3 (+3.9/-5.7) (REF). Both agree with the single star
     (9.3 +- 0.9).
4. **Final ephemeris** (both pipelines combined, injection-calibrated timing errors; `out/ephemeris_final.json`):
   * **T0 = 2460840.91488 +- 0.00093 BJD_TDB (3840.91488 BTJD_TDB, epoch 0).**
   * **P = 5.6915020 +- 0.0000078 d.** At the reference epoch E = -26, T0 = 3692.93583 +- 0.00091 and corr(T0, P) = -0.002.
   * TTV chi2 9.9/16 (INV) and 15.4/19 (REF), Monte Carlo p = 0.56 / 0.63 against the empirical pull distributions.
     The combined fit gives 15.9/19. No TTVs. O-C rms 8-9 min.
   * The 1-sigma prediction error in 2027 (E = 100-160) is 1.9-2.5 min.
   * The changes from phase 1 (3840.9141 +- 0.0013, 5.6915087 +- 0.0000090 d) are +1.1 min and -6.7e-6 d (0.6 sigma).
5. **Depth versus epoch: statistical, not instrumental bias and not astrophysical** (`c07`, `c05`, `c09`, `c11`,
   `c16`).
   * After cleaning, S39+S66 versus S103+S104 is +3.3 +- 2.1 ppt (INV, 1.6 sigma), +3.5 +- 2.1 (REF, 1.7 sigma) and
     +2.9 +- 2.1 (REF without the scattered-light event E58, 1.35 sigma).
   * The look-elsewhere-corrected p of this post-hoc split is 0.31-0.66 across three families of possible splits.
     The omnibus sector chi2 is 7.4/5 (p 0.19, INV) and 4.8/5 (p 0.44, REF). The linear trend is
     +0.5 to +0.7 +- 0.4 ppt/yr.
   * **Pixel-level injection** of a constant 15.7-ppt transit into the raw cubes, 10 fake ephemerides per sector,
     recovers 0.998-1.000 (matched templates) in every sector in both pipelines.
   * Plausible template errors (0.07-px shifts, 0.15-px blur, the other pipeline's PRF) change the recovery by
     -7% to +6%. Even the most adverse combination per sector gives a late/early ratio of at most 1.11 (INV) and
     1.07 (REF), against an observed ~1.25.
   * The pipelines do not manufacture a sector-dependent depth of the observed size.
   * No correlation with the Mira's flux, its brightness slope, its pulsation phase, the background or the jitter
     (all |z| <= 1.7, permutation p >= 0.03 before multiple-test correction).
   * The flux loss is uncorrelated with the Mira's flux: slope -0.0002 +- 0.0004 ppt per e-/s. If the eclipsed body
     orbited the Mira, the slope would be +0.0010; that is excluded at 2.7 sigma (INV) and 3.4 sigma (REF).
   * The contrast moves by -3.7 +- 1.9 ppt (2.0 sigma) when the nuisance-regressor set is changed (`c16`). This
     photometry's systematics are as large as the contrast itself.
6. **Phase 0.490-0.495 bump and phase-0.57 excursions** (cleaned scans, `c08`):
   * The bump persists in both pipelines at 1.2-1.8 ppt (S/N 1.3-2.2 after calibration; 1.5-2.0 ppt, S/N 1.7-2.4
     without S12). The global p is 0.6-0.9 for ~56 independent phases.
   * The phase-0.57 excursions no longer have opposite signs. Both pipelines show a brightening of -1.9 to -2.8 sigma
     at phase 0.575-0.58, and the earlier REF dip at 0.5725 is gone (-0.8 sigma).
   * The largest REF scan value (S/N 4.0 at phase 0.5375) is one 30-min S12 window at a gap edge, a single scan step
     wide. It disappears without S12.
   * Secondary at phase 0.5: +0.3 +- 0.8 ppt (INV) and +0.3 +- 0.9 ppt (REF); 3-sigma limits 2.8 / 3.1 ppt.
7. **Weightings:**
   * Within INV, apertures of 1, 2, 4 px and r <= 1.5 px agree with the PRF fit sector by sector (chi2 1.6-4.4/5).
     They are deeper on average by +2.2 to +3.0 +- 1.8-2.2 ppt.
   * Within REF, the aperture-PRF offsets average +1.3 to +1.8 +- 0.8-0.9 ppt. They scatter more (chi2 9-19/5),
     driven by S12 and two S103 events (E55, E58).
   * The late-early contrast appears in every weighting, so it is not a PRF-weighting artefact. Shared noise is
     not ruled out.
8. **Odd/even** (cleaned, fixed shape, calibrated errors):
   * INV +3.7 +- 1.8 ppt (2.1 sigma). It was +5.3 +- 2.1 before cleaning with the same code; cleaning removes only
     part of it.
   * REF +1.3 +- 1.8 ppt (0.7 sigma). Combined, allowing for the correlation of the two pipelines' errors
     (rho = 0.56): +2.5 +- 1.6 ppt (1.6 sigma).
   * The INV excess remains a property of the INV photometry of a few events (S104 E60/E63, the missing E58), as
     the referee found.

**For the planet hypothesis (from this folder):**
* The transit survives the removal of the field-variable coupling unchanged in depth, shape and timing, in both
  independent photometries.
* The flux loss is constant in time within the errors and uncorrelated with the Mira's brightness.
* The light curve is flat-bottomed. The density is consistent with the single M3V. There are no TTVs and no
  secondary eclipse.

**Against or still open:**
* The cleaned INV odd-even is +3.7 +- 1.8 ppt (2.1 sigma). REF does not confirm it.
* There is a ~25% (1.4-1.7 sigma) late/early depth contrast in all weightings. The injections exclude a pipeline
  bias of that size, so if it were real it would point to an unmodelled PSF effect in the S103/S104 geometry, or to
  an intrinsically variable eclipsed source.
  * S103/S104 are on different cameras/CCDs, with the M3V on the other side of the Mira. The static Mira PRF
    residual at the M3V is ~0.9 x F_M3V there, against -0.4 to +0.3 in S12-S93.
  * The REF per-sector depth correlates with that residual (z = 2.1, permutation p = 0.008 before correction for
    20 tested variables; ~0.15 after).
* The adopted depth carries a 1.1-ppt systematic. The two methods differ by 1.6 ppt on the same pixels.
* A 0.8-ppt residual signal at the RRab period remains in the INV S93 curve.

## 1. Removing the field-variable coupling (c00, c01, c02, c12, c15)

**Variables masked** (`cleanlib.var_mask`, `out/varpos.json`):
* Every Gaia DR3 phot_variable_flag = VARIABLE source with G < 18.6 inside the cutout, except the Mira, which both
  pipelines model separately.
* Mask radius: 2.5 px (G >= 14), 3.0 px (12 <= G < 14), 3.5 px (G < 12).
* This includes the RRab 5925372615693467904, the EBs 5925384852082638720 (P 0.7503 d) and 5925378392451781120
  (P 1.5953 d), the bright LPV 5925373341570270208 (G 11.6, 113"), an RS CVn and further ECL/LPV/SOLAR_LIKE
  sources.
* In S12-S93 the RRab sits at the cutout edge. In S103/S104 it is inside, at (18.9, 16.6) and (17.4, 17.2).

**INV** (`cleanlib.jitter_clean`, `cleanlib.photometry`, a marked copy of p04/p05/p06):
* The masked pixels are removed from the pointing / flux-scale fit and from the faint-pixel background (level and
  plane).
* With no masks the copy reproduces the phase-1 model-B fluxes exactly (max |diff| = 0.0 e-/s, `scripts/t_repro.py`).
* Versions: none, rr3 (RRab and the 2 EBs only), all (**adopted**), and allreg (all, plus the aperture light curves
  of the RRab and both EBs as explicit pixel regressors).

**REF** (`refclean.aux_clean`, `refclean.SectorC`): the same masks are removed from the referee's faint-pixel
B(t) and field-star pointing pixel sets. With no masks it reproduces the referee's aux file exactly.

**Leak statistic:**
* First-harmonic semi-amplitude of a 3-harmonic fit at the Gaia period, phased to the Gaia epoch of maximum, on
  the out-of-transit, quality-filtered, detrended light curve.
* p_emp from 120 fake periods of 0.45-0.95 d (the RRab and EB periods and their halves are avoided).
* 95% upper limit = sqrt(max(A^2 - 2 sigma^2, 0)) + 1.645 x 1.2 sigma, with sigma the per-quadrature noise from the
  fake periods.
* The EBs are tested with 4-harmonic Delta-chi2.

| Sector | INV phase-1 (none) A1, ppt (p) | INV cleaned (all) A1 (p), UL95 | INV allreg A1 (p) | dx RR p none -> all | REF orig -> all, equivalent A1 (p) |
|---|---|---|---|---|---|
| 12 | 0.76 (0.34) | 1.05 (0.15), 1.80 | 1.19 (0.18) | 0.52 -> 0.70 | 0.33 (0.95) -> 0.62 (0.73) |
| 39 | 1.93 (<0.01) | 0.45 (0.23), 0.78 | 1.58 (<0.01) | <0.01 -> 0.44 | 0.54 (0.53) -> 0.72 (0.25) |
| 66 | 2.51 (<0.01) | 0.35 (0.73), 0.75 | 2.76 (<0.01) | 0.12 -> 0.48 | 0.63 (0.31) -> 0.48 (0.43) |
| 93 | 3.24 (<0.01) | **0.78 (<0.01)**, 1.21 | 0.92 (<0.01) | <0.01 -> 0.44 | 0.52 (0.34) -> 0.71 (0.18) |
| 103 | 3.79 (<0.01) | 0.23 (0.88), 0.75 | 0.57 (0.33) | <0.01 -> 0.64 | 0.31 (0.60) -> 0.25 (0.71) |
| 104 | 2.08 (<0.01) | 0.69 (0.10), 1.16 | 1.16 (0.02) | <0.01 -> 0.12 | 0.23 (0.78) -> 0.26 (0.69) |

(`out/leak_summary.json`, `out/leak_inv.json`, `out/leak_ref.json`, `figs/rrab_leak.png`.)

* **EBs.** In the M3V curves the EB signals were not significant before cleaning (p >= 0.07) and remain so after,
  except S103 EB 0.75 d, where Delta-chi2 = 20 gives p = 0.01; this is expected once in 12 tests.
* **REF tiled curve.** REF windows tiled every 0.3 d, keeping the central +-0.15 d, give a contiguous curve.
  * A 3 e-/s sinusoid at P_RR injected at the M3V is recovered with a transfer of 0.54-0.56 (S39, S93, S103;
    `out/tile_transfer.json`). The REF amplitudes above are divided by 0.55.
  * Equivalent REF 95% upper limits are 0.65-1.39 ppt.
* **S93 residual** (`scripts/e01_residual_rr.py`, `scripts/e02_rr_pixmap.py`, exploratory).
  * With 1.6x mask radii, the faint-pixel plane gradients are free of the RRab, but the M3V amplitude stays at
    0.79 ppt.
  * The per-pixel RRab-period amplitude images within 3 px of the M3V give chi2 117 for 58 values. Fitted by the
    M3V PRF or by the Mira gradient images, this drops to 75-87, so the two cannot be distinguished.
  * It is a small pointing-like or M3V-like term from a route not captured by the masked estimators.
  * Its effect on transit depths is <= 0.8 ppt per event and <= 0.2 ppt on the global depth, which is in the
    systematic budget.
  * An optional per-sector RRab model is provided in the INV CSV (`rrab_resid_model`).

## 2. Cleaned-light-curve transit analysis (c03, c04, c05, c06, c14)

**Data sets** (`c03`, `data/win_*.npz`):
* INV: undetrended flux 1 + dF/979, +-0.30 d windows, per-window offset + slope. The detrended curve gives the same
  fake-epoch scatter.
* REF: strict windows.
* Acceptance (both): >= 2 in-transit points and >= 25% of the out-of-transit cadences on each side.
  * INV keeps 18 windows. It loses S12 E-388/E-387 (too few in-transit points), S104 E62 (no data before
    ingress) and S103 E58 (background cut).
  * REF keeps 21 windows, including S103 E58 (a scattered-light event at 35 ppt, with a per-point scatter 1.5x the
    sector norm) and S12 E-388/E-387.

**Noise calibration** (`out/noise_calib.json`, 150 fake epochs per sector):
* beta = std of fake single-window depths / median formal error.
* INV: 1.22, 1.02, 1.25, 1.13, 1.03, 1.38 (S12-S104).
* REF: 1.37, 1.54, 1.62, 1.58, 1.02, 1.33.
* Cleaning lowered the REF S93 fake-epoch scatter from 6.45 to 4.26 ppt per transit.

**Global fits** (`c04`, emcee 40 x 6000; `out/global_*.json`, `out/chains_*.npz`).
* Model: batman with Claret LD priors and per-window baselines profiled.
* The autocorrelation time is 62-82 steps with the density prior and 90-176 free. The free-density chains are
  shorter than 50 tau: medians are reliable, tails are approximate.

| Fit | k | depth (ppt) | b | rho* (g/cm3) | T14 (h) | T23/T14 | Rp (R_E) | T0 (BTJD_TDB) | P (d) |
|---|---|---|---|---|---|---|---|---|---|
| INV, rho prior | 0.1131 +- 0.0047 | 14.50 +- 0.98 | 0.49 (+0.09/-0.13) | 9.4 +- 0.9 | 1.73 +- 0.07 | 0.74 +- 0.04 | 4.73 +- 0.24 | 3840.9146 +- 0.0017 | 5.6915035 +- 0.0000090 |
| REF, rho prior | 0.1199 +- 0.0043 | 16.06 +- 1.00 | 0.54 (+0.06/-0.08) | 9.4 +- 0.9 | 1.69 +- 0.05 | 0.71 +- 0.03 | 5.01 +- 0.24 | 3840.9156 +- 0.0012 | 5.6915051 +- 0.0000080 |
| INV, free rho | 0.1151 (+0.0076/-0.0053) | 14.74 +- 1.03 | 0.55 +- 0.29 | 8.4 (+5.0/-4.4) | 1.75 +- 0.14 | 0.72 +- 0.13 | 4.83 +- 0.31 | 3840.9149 +- 0.0017 | 5.6915039 +- 0.0000093 |
| REF, free rho | 0.1194 (+0.0053/-0.0046) | 16.00 +- 0.98 | 0.40 +- 0.27 | 12.3 (+3.9/-5.7) | 1.65 +- 0.10 | 0.75 +- 0.08 | 5.00 +- 0.26 | 3840.9158 +- 0.0012 | 5.6915050 +- 0.0000080 |

* With the density prior: a/R* = 25.3 +- 0.8, i = 88.8-88.9 deg, Teq = 471 +- 11 K (Teff 3350 +- 60 K, zero albedo,
  full redistribution). Rp uses R* = 0.383 +- 0.012 Rsun.
* Grazing (b > 1 - k) posterior fraction: 0 (prior), 1.4% (INV free), 0.5% (REF free).

**Adopted depth** (`c14`, `out/final_depth.json`):
* GLS combination of INV and REF with their noise correlation rho = 0.56 (pooled no-injection depths at the same
  fake transits, c05).
* Result: **15.23 +- 0.87 (stat) +- 1.08 (sys) ppt** (total 1.39). Flux loss 14.9 e-/s.

| Systematic term | ppt | basis |
|---|---|---|
| INV vs REF method | 0.78 | half the difference of the two global depths |
| PRF template | 0.52 | rms of (recovery - 1) for 0.07-px shifts and 0.15-px blur, both pipelines, all sectors (c05) |
| F_M3V normalisation (979 e-/s) | 0.46 | +-3% (XP statistics 2%; route spread -1.3/+5.7%) |
| window selection (E58) | 0.17 | REF with / without E58 |
| residual RRab | 0.20 | upper bound from the S93 residual |
| baseline / detrending version | 0.05 | INV raw vs detrended |
| **total** | **1.08** | quadrature |

* k = 0.1170 +- 0.0039 (stat) +- 0.0041 (sys); Rp = 4.90 +- 0.22 (stat, incl. R*) +- 0.17 (sys) R_E.
* Agrees with phase 1 (15.7 +- 0.8 +- 1.0 ppt; shape: Rp 4.94 +- 0.27 R_E).

**Per-transit depths and times** (`c06`, `out/per_transit_clean.csv`):
* Fixed shape (k 0.1195, b 0.52, a/R* 25.4), free depth and mid-time on a grid with parabolic refinement, linear
  baseline.
* Errors are calibrated by the pixel-level injections (c05: 23-42 fake transits per sector and pipeline).
  * Timing: sigma = formal x std(pulls). INV pull std 1.85, 3.20, 2.70, 1.42, 2.17, 1.44 (robust 0.8-1.9).
    REF 1.97, 1.85, 1.53, 1.27, 2.49, 1.37. 3-15% of pulls are > 3 sigma.
  * Depth: the no-injection scatter at the same fake epochs (2.7-8.6 ppt per transit).
* INV-REF error correlation at identical fake transits: timing 0.24-0.74, depth 0.33-0.70.
* Single-transit comparison on 18 common events: depth difference rms 3.9 ppt (correlation 0.57), time difference
  rms 7.2 min.

**TTV chi2** (`out/timing_clean.json`):

| | INV | REF |
|---|---|---|
| Calibrated chi2 | 9.9/16 (p 0.87; MC p 0.56) | 15.4/19 (p 0.70; MC p 0.63) |
| Formal errors | 40.8 | 48.1 |
| Robust pull scale | 27.1 (0 events > 3 sigma) | 32.2 (0 events > 3 sigma) |
| O-C rms | 8.3 min | 9.0 min |

The phase-1 S66 E-128 "outlier" is -20 min in INV and +1.8 min in REF; the calibrated INV error is 17 min.

**Ephemeris** (`out/ephemeris_final.json`):
* Per event, the inverse-variance average of the INV and REF times (non-negative weights; the variance includes the
  per-sector correlation), followed by a weighted linear fit.
* A full GLS of the 39 measurements gives P = 5.6915000 +- 0.0000075 and T0 = 3840.91504 +- 0.00091. It is reported
  as a variant only, because with rho_dt > sigma_REF/sigma_INV it puts negative weight on INV and its T0 falls
  outside both pipelines' values.

| Solution | T0 at E = 0 (BTJD_TDB) | P (d) | chi2/dof |
|---|---|---|---|
| INV | 3840.91425 +- 0.00140 | 5.6915106 +- 0.0000151 | 9.9/16 |
| REF | 3840.91476 +- 0.00096 | 5.6914997 +- 0.0000079 | 15.4/19 |
| **Combined (adopted)** | **3840.91488 +- 0.00093** (BJD_TDB 2460840.91488); T0(E = -26) = 3692.93583 +- 0.00091, corr -0.002 | **5.6915020 +- 0.0000078** | 15.9/19 |

* Covariance matrices at both reference epochs are in the JSON.
* All times are target-barycentred BJD_TDB (phot p00 / referee r02 corrections, which agree within 0.6 s).
* The ephemerides of `../ctx/out/followup_transits.csv` (phase 1) are shifted by +1.1 min at E = 0 and by -0.1 min
  at E = 130 (2027); their uncertainty was ~2 min.

## 3. Depth versus epoch (c07, c05, c09, c11, c16)

**Per-sector depths** (fixed shape of each pipeline's posterior median, beta-scaled errors; `out/depth_epoch.json`,
`figs/depth_epoch.png`):

| | S12 | S39 | S66 | S93 | S103 | S104 | all | chi2/5 | S39+66 vs S103+104 | early vs late |
|---|---|---|---|---|---|---|---|---|---|---|
| INV cleaned | 19.7 +- 6.3 | 13.1 +- 1.6 | 11.9 +- 2.5 | 15.8 +- 1.9 | 13.7 +- 2.1 | 19.9 +- 2.6 | 14.6 +- 0.9 | 7.4 (p 0.19) | +3.3 +- 2.1 (1.6 sigma) | +2.9 +- 1.8 (1.6 sigma) |
| REF cleaned | 12.2 +- 3.8 | 14.6 +- 2.0 | 14.5 +- 2.7 | 14.7 +- 2.3 | 17.5 +- 1.9 | 18.7 +- 1.9 | 16.1 +- 0.9 | 4.8 (p 0.44) | +3.5 +- 2.1 (1.7 sigma) | +3.0 +- 1.9 (1.6 sigma) |
| REF without E58 | 12.2 +- 3.8 | 14.6 +- 2.0 | 14.5 +- 2.7 | 14.7 +- 2.3 | 16.0 +- 2.0 | 18.7 +- 1.9 | 15.8 +- 0.9 | 4.1 (p 0.53) | +2.9 +- 2.1 (1.35 sigma) | +2.5 +- 1.9 (1.3 sigma) |

* Before cleaning, the identical c11 measurement gave INV S39+66 vs S103+104 = +5.1 +- 2.4. Cleaning lowered it to
  +4.0 +- 2.1 (c11 weighting). The RRab leak explains about a quarter of the phase-2 contrast.
* **Look-elsewhere.** 200,000 null draws of six sector depths with the measured errors.
  * For S39+S66 vs S103+S104: global p = 0.31-0.39 against the 5 contiguous early/late splits, 0.55-0.64 against
    all 31 bipartitions, and 0.57-0.66 against all 130 pairs of disjoint sector groups (each of >= 2 sectors).
  * The largest |z| of any split is 2.1-2.6, with global p 0.13-0.39.
  * Trend in time: +0.52 +- 0.41 (INV) and +0.72 +- 0.40 (REF) ppt/yr.
* **Injection test** (c05; 23-42 fake transits per sector and pipeline; `out/inject_pix_s*.csv`,
  `out/inject_summary.csv`):
  * Matched templates recover 0.998-1.000 in all six sectors.
  * Cross-pipeline templates: the INV PRF injected into REF gives 0.98-1.00; the REF PRF into INV gives 0.99-1.02.
  * Shifted by 0.07 px: 0.93-1.06 (INV), 0.95-1.03 (REF). Blurred by 0.15 px: 0.96-0.99.
  * Most adverse per-sector choice of template error: late/early recovery ratio <= 1.11 (INV) and <= 1.07 (REF).
  * The noisy recovered depths per sector (matched): INV 15.2-17.5, REF 13.7-16.1 ppt, all consistent with 15.7.
  * **The pipelines cannot produce the observed contrast from a constant transit unless the true PSF error is
    several times larger than the field-star WCS/PRF calibration allows.**
* **Correlates.**
  * (a) Mira aperture flux (12,500-18,800 e-/s at the transits): z = -0.4 (INV), -1.2 (REF).
    * The Mira-host scaling of +0.00095 ppt per e-/s is excluded at 2.7 / 3.4 sigma.
    * Mira flux slope: z = -0.8 / -1.6 (p_perm 0.38 / 0.04).
    * Pulsation phase (VSX: P 294 d, epoch JD 2457654.0) cos/sin: |z| <= 1.7 (p_perm >= 0.03).
  * (b) Static Mira PRF residual near the M3V (free-position PRF fit to the mean image):
    * Projection on the M3V PRF: -0.40, +0.19, +0.31, -0.06, +0.88, +0.89 x F_M3V (S12-S104); rms 60-310 e-/s.
    * Correlation of the M3V PRF with the Mira b image over the fit pixels: 0.50-0.56 (S12-S93) against 0.68-0.71
      (S103/S104).
  * (c) Geometry:
    * S12-S93: camera 1, CCD 2/1, M3V at PA ~180 deg from the Mira.
    * S103: camera 2 CCD 3, PA -53 deg. S104: camera 1 CCD 4, PA -39 deg.
    * Pixel phases vary randomly. Separation 0.90-0.94 px.
  * Sector-level regressions (6 points, 10 variables x 2 pipelines): the only |z| >= 2 are REF on (b), z = 2.08
    (p_perm 0.008), and on the template correlation, z = 2.04 (p 0.04). After 20 tests this is p ~0.15.
  * (d) Background level and maximum in the window: |z| <= 0.5. Jitter: z <= 1.1.
  * (e) Per-sector PRF/TIC flux scale of field stars (1.064-1.141): renormalising the depths per sector leaves the
    contrast at 1.7-1.8 sigma. The out-of-transit M3V flux itself cannot be measured next to the Mira (static fits
    are degenerate, as in phase 1).
* **Regressor sensitivity** (c16).
  * The allreg variant (the variables' light curves as pixel regressors; rejected for its leaks and noise) has no
    contrast (-0.2 +- 2.5).
  * Paired against the adopted version at 60 fake epochs per sector, the contrast change is -3.7 +- 1.85 ppt
    (2.0 sigma).
  * The contrast is therefore not robust to nuisance-regression choices at its own size.
* **Verdict.**
  * Statistical: 1.35-1.7 sigma before look-elsewhere, p ~0.3-0.7 after.
  * Possibly with a minor instrumental component specific to the S103/S104 scene.
  * There is no evidence for an astrophysical change. A real 25% change in flux loss would need the eclipsed source
    to vary by 25%; the M3V cannot (spots < 2%), and a Mira host is excluded by the flux test above.

## 4. Phase 0.490-0.495 and phase 0.57 (c08)

* **Scan:** 0.15-0.85, step 0.0025, fixed planet-occultation shape (T14 1.70 h, T23/T14 0.72), per-event +-0.30 d
  windows.
  * INV on the cleaned undetrended curve.
  * REF with dedicated referee-method windows at every scan phase and epoch.
  * S/N is calibrated by the scan's robust std: 1.15 (INV) and 1.56 (REF).
  * Results: `out/phase_scan_*.csv`, `out/phase_scan.json`, `figs/phase_scan.png`.

| | INV | REF | INV without S12 | REF without S12 |
|---|---|---|---|---|
| depth at phase 0.490 / 0.4925 / 0.495 (ppt) | 1.77 / 1.63 / 1.28 | 1.20 / 1.18 / 1.51 | 2.04 / 1.96 / 1.40 | 1.52 / 1.55 / 1.79 |
| calibrated S/N | 2.15 / 1.91 / 1.56 | 1.28 / 1.26 / 1.59 | 2.38 / 2.23 / 1.68 | 1.75 / 1.78 / 2.03 |
| minimum near 0.57 (S/N, phase) | -1.98 (0.5775) | -2.84 (0.58) | -2.07 (0.575) | -1.83 (0.575) |
| at 0.5725 (phase-2 REF dip) | -0.70 | -0.80 | -0.83 | -0.34 |
| secondary at 0.5 (ppt), 3-sigma limit | +0.32 +- 0.83, 2.8 | +0.30 +- 0.93, 3.1 | +0.38 +- 0.86, 3.0 | +0.80 +- 0.90, 3.5 |
| largest |S/N| anywhere (phase) | 2.97 (0.190) | 4.01 (0.5375; one S12 window) | 2.39 (0.490) | 2.69 (0.345) |

* The 0.49 bump is common to both pipelines, which share pixels; the INV-REF scan correlation is 0.45.
* At ~2 sigma with ~56 independent phases its global p is ~0.6-0.9. A test restricted a priori to |phase - 0.5| <
  0.01 would still give only p ~0.07.
* If real, it would need e cos(omega) ~ -0.012 and a 1.3-2 ppt occultation, i.e. T_b ~2400-2500 K. Not significant.
* The "phase-0.57 excursions" are now a shallow brightening in both curves. Their opposite signs in phase 2 came
  from the uncleaned photometry. Nothing at phase 0.57 calls for an explanation.
* The REF spike at 0.5375 is S12 E-385 at BTJD 1652.75: 21 points, 14 out of transit, formal S/N 16.6 at one scan
  step only (`scripts/e03_scan_spike.py`).

## 5. Photometric weightings (c09)

* Same cleaned pixel data in both pipelines.
  * INV: model-A residual cube (Mira b, pair PRF and plane subtracted) summed over the aperture and divided by its
    PRF fraction.
  * REF: residual frames of the per-pixel baselines.
* Apertures:
  * px1, px2, px4: the 1, 2 and 4 M3V-PRF-brightest pixels, holding 19-27%, 31-39% and 50-60% of the PRF.
  * r <= 1.5 px: 72-81% of the PRF.
* Errors are paired fake-epoch differences (60 per sector). Results: `out/apertures.json`,
  `out/apertures_per_transit.csv`.

| Aperture - PRF (ppt) | S12 | S39 | S66 | S93 | S103 | S104 | mean | chi2/5 |
|---|---|---|---|---|---|---|---|---|
| INV px1 - B | -1.3 +- 7.0 | +2.6 +- 4.6 | -2.5 +- 6.2 | +5.0 +- 3.4 | +3.8 +- 6.5 | -2.6 +- 15 | +2.6 +- 2.2 | 1.6 |
| INV px4 - B | -0.6 +- 6.6 | +1.2 +- 4.2 | -6.2 +- 6.4 | +7.2 +- 3.3 | +4.2 +- 4.6 | -1.2 +- 12 | +3.0 +- 2.0 | 4.3 |
| INV r1.5 - B | -0.5 +- 6.9 | +2.0 +- 3.8 | -5.5 +- 5.6 | +6.8 +- 3.2 | +3.3 +- 3.8 | +0.3 +- 10 | +2.8 +- 1.8 | 4.1 |
| INV A - B | -1.3 +- 6.2 | +1.9 +- 4.1 | -4.7 +- 5.6 | +6.0 +- 3.2 | +3.4 +- 3.7 | +0.4 +- 10 | +2.6 +- 1.9 | 3.3 |
| REF px1 - M | +12.5 +- 3.7 | +0.2 +- 1.5 | +4.2 +- 3.8 | -1.3 +- 1.6 | +8.3 +- 3.1 | -5.7 +- 7.1 | +1.4 +- 1.0 | 19.1 |
| REF px2 - M | +4.5 +- 2.0 | -0.3 +- 1.7 | +4.1 +- 6.5 | +0.5 +- 1.2 | +9.5 +- 2.6 | -3.1 +- 4.3 | +1.8 +- 0.8 | 14.9 |
| REF px4 - M | +7.5 +- 4.2 | +1.3 +- 2.3 | +2.5 +- 3.4 | +0.3 +- 1.2 | +4.3 +- 1.4 | -2.1 +- 2.2 | +1.6 +- 0.8 | 10.0 |
| REF r1.5 - M | +11.7 +- 5.5 | +3.4 +- 3.1 | +2.8 +- 3.0 | -0.3 +- 1.6 | +2.5 +- 1.2 | -1.6 +- 1.7 | +1.3 +- 0.8 | 9.0 |
| REF 3-star T - M | -13.8 +- 6.7 | -1.9 +- 2.5 | -1.0 +- 2.4 | +2.6 +- 3.3 | +1.7 +- 2.2 | +3.1 +- 5.6 | -0.2 +- 1.2 | 6.6 |

* INV weightings agree sector by sector.
* REF apertures deviate in S12 (30-min data) and in S103, where E55 and E58 are 10-25 ppt deeper in small apertures
  than in the PRF fit.
* Apertures are on average 1.3-3.0 ppt (8-20%) deeper than the PRF fits, at 1.2-2.2 sigma. This may reflect a PSF
  core more compact than the SPOC PRF near this red star.
  * It is not added to the budget: the 979 e-/s normalisation is itself defined in PRF units from PRF fits to field
    stars.
  * It would raise the depth if apertures were the reference.
* The late-early contrast appears in every weighting: INV apertures +5 to +8 +- 4-6, REF apertures +4 to +11 +- 3-4
  ppt. Those c09 contrasts are unweighted sector means and are inflated by E58.
* So the contrast is not produced by PRF weighting. Since all weightings share the same pixels and noise, this
  does not make it more significant.

## 6. Final light curves

### `lc_M3V_INV_clean.csv`

All six sectors, every finite cadence, 50,560 rows. INV method with all Gaia variables excluded from the field-wide
estimators.

| Column | Meaning |
|---|---|
| time_btjd_tdb | BJD_TDB - 2457000, barycentred for the M3V |
| time_bjd_tdb | BJD_TDB |
| sector, cadence_s | sector and exposure time (1800 / 600 / 200 s) |
| quality_ok | 1 = use: QUALITY in {0, 2048}, background excess < 60 e-/s/px, pointing within 6 robust sigma, model-B chi2 < 3x median, not flare-flagged |
| flux_norm | detrended normalised flux: wotan biweight, 0.75 d, transits +-0.1 d masked; 1 = out of transit; normalised to 979 e-/s |
| flux_err | per-cadence error: sector robust sigma x sqrt(local chi2 ratio >= 1) |
| beta_transit | fake-epoch red-noise factor of the sector for transit-length fits; multiply flux_err by it |
| flux_raw_norm | undetrended 1 + dF/979 (used for the transit fits here, with per-window linear baselines) |
| trend | biweight trend that was divided out |
| dflux_e_s | M3V flux change from model B, e-/s |
| flare_flag | > 4 robust sigma and its neighbours |
| rrab_resid_model | per-sector 3-harmonic model of any residual RRab-period signal, out of transit (fractional) |
| flux_norm_modelA, flux_err_modelA | all-free model A (noisier) |
| epoch | from the final ephemeris |
| in_transit_window | within +-0.1 d of a final-ephemeris transit |
| tess_quality | TESS QUALITY flag |
| bkg_e_s_px | faint-pixel background, variables excluded |
| jitter_dx_px, jitter_dy_px | cleaned pointing offsets |
| time_btjd_tesscut | original TESScut TIME |

The optional `rrab_resid_model` can be subtracted from flux_norm - 1. It is tiny except in S93.

### `lc_M3V_REF_clean.csv`

One strict window per transit, 21 windows and 4,663 rows; REF method with cleaned auxiliary regressors.

| Column | Meaning |
|---|---|
| time_btjd_tdb, time_bjd_tdb | as above |
| sector, epoch, window_id, cadence_s | |
| dflux_e_s | M3V flux change, e-/s, after the per-pixel baselines |
| flux_norm | 1 + dflux/979 |
| flux_err | window out-of-transit robust std |
| beta_transit | as above (REF values) |
| oot | 1 = outside +-0.07 d |
| T_M3V_e_s, T_pair_e_s, T_Mira_e_s | three-star fit |
| bkg_e_s_px | cleaned B(t) |

S103 E58 (BTJD 4171.02) lies in scattered light. Treat it with care.

### `lc_M3V_REF_tiled.csv`

The REF method run as contiguous tiles every 0.3 d, keeping +-0.15 d of each, for out-of-transit work.

* Columns: time, sector, tile_id, dflux_e_s, flux_norm, flux_err, in_transit_window.
* The per-tile quadratic baselines suppress slow signals. A signal at 0.65 d keeps 55% of its amplitude.
* Not for transit fitting.

## Key numbers

| Quantity | Value | Method | Uncertainty | File |
|---|---|---|---|---|
| RRab leak in INV, before -> after | 1.9-3.8 -> 0.2-0.8 ppt (S39-S104) | 3-harmonic fit at Gaia P, 120 fake periods | UL95 after 0.75-1.8 ppt | out/leak_summary.json |
| RRab in REF | <= 0.6 ppt equivalent, p >= 0.18 | tiled REF curve, transfer 0.55 | UL95 0.65-1.4 ppt | out/leak_summary.json, out/tile_transfer.json |
| Per-transit depth change from cleaning (INV) | rms 2.7, max 7.4 ppt; slope vs predicted leak -0.82 +- 0.11 | same fit on all versions | regression | out/modes_compare.json |
| Adopted depth | 15.23 ppt (flux loss 14.9 e-/s) | GLS of INV and REF batman posteriors (rho prior), rho 0.56 | +- 0.87 stat +- 1.08 sys | out/final_depth.json |
| Rp/R*, Rp | 0.1170, 4.90 R_E | as above | +- 0.0039 +- 0.0041; +- 0.22 +- 0.17 R_E | out/final_depth.json |
| b, a/R*, T14, T23/T14 | 0.49-0.54, 25.3, 1.69-1.73 h, 0.71-0.74 | emcee, rho prior | +-0.07-0.13, +-0.8, +-0.05-0.07, +-0.03-0.04 | out/global_*_single.json |
| rho* free | 8.4 (+5.0/-4.4) INV, 12.3 (+3.9/-5.7) REF | emcee, log-uniform rho | posterior | out/global_*_free.json |
| Ephemeris (adopted) | T0 = 2460840.91488 BJD_TDB, P = 5.6915020 d | per-event combined times, WLS | +- 0.00093 d, +- 0.0000078 d | out/ephemeris_final.json |
| TTV chi2 | INV 9.9/16, REF 15.4/19 | injection-calibrated errors | MC p 0.56 / 0.63 | out/timing_clean.json |
| Odd - even | INV +3.7 +- 1.8, REF +1.3 +- 1.8, combined +2.5 +- 1.6 ppt | fixed-shape group fit | beta-calibrated | out/final_depth.json |
| S39+66 vs S103+104 | +3.3 +- 2.1 (INV), +3.5 +- 2.1 (REF) ppt | fixed-shape sector fits | LEE p 0.31-0.66 | out/depth_epoch.json |
| Injection recovery (constant 15.7 ppt) | 0.998-1.000 matched; 0.93-1.06 mismatched | 10 fake ephemerides x 6 sectors, both pipelines | adverse late/early <= 1.11 | out/inject_summary.csv, out/depth_epoch.json |
| Mira-host flux-scaling test | slope -0.0002 +- 0.0004 vs +0.00095 expected | per-transit regression on Mira flux | excluded 2.7-3.4 sigma | out/depth_epoch.json |
| Secondary at phase 0.5 | +0.3 +- 0.8 (INV), +0.3 +- 0.9 (REF) ppt | fixed occultation shape | 3-sigma 2.8 / 3.1 ppt | out/phase_scan.json |
| Phase 0.49 bump | 1.2-2.0 ppt, S/N 1.3-2.4 | phase scan | global p ~0.6-0.9 | out/phase_scan.json |
| Aperture - PRF depth | INV +2.2 to +3.0 +- 1.8-2.2; REF +1.3 to +1.8 +- 0.8-1.0 ppt | paired fake-epoch differences | chi2 INV 1.6-4.4/5, REF 9-19/5 | out/apertures.json |

## Scripts (`scripts/`; /home/user/work/venv/bin/python, OMP_NUM_THREADS = 1)

| Script | What it does | Outputs |
|---|---|---|
| `cleanlib.py` | INV pipeline copy with variable masks (jitter_clean, photometry, make_lc), helpers | - |
| `refclean.py` | REF auxiliary regressors with variable masks; SectorC (r05 windows on cleaned aux); tiled light curve | - |
| `tfit.py` | batman model with midpoint exposure integration, window container, linear depth fits, single-transit fit | - |
| `c00_varpos.py` | pixel positions of the Gaia variables per sector | out/varpos.json |
| `t_repro.py` | check: unmasked copy = phase-1 phot fluxes (0.0 e-/s) | stdout |
| `c01_build_inv.py` | INV photometry for none/rr3/all/allreg; leak statistics in M3V/Mira/pair curves and regressors | data/lcinv_*_s*.npz, data/jitter_*_s*.npz, out/leak_inv.json |
| `c02_build_ref.py` | REF photometry for orig/none/rr3/all: strict + 150 fake windows, tiled curve, leak statistics | data/refwin_*_s*.npz, data/reftile_*_s*.npz, data/refaux_*_s*.npz, out/leak_ref.json |
| `c03_windows.py` | transit-window data sets and fake-epoch noise calibration (beta) | data/win_*.npz, out/noise_calib.json |
| `c04_global.py` | emcee batman fits (rho prior / free) for INV and REF | out/global_*.json, out/chains_*.npz |
| `c05_inject.py` | pixel-level injection of a constant 15.7-ppt transit; both pipelines; 7 template variants | out/inject_pix_s*.csv |
| `c06_timing.py` | per-transit depths and times, calibration, TTV, odd/even, final ephemeris | out/per_transit_clean.csv, out/timing_clean.json, out/ephemeris_final.json |
| `c07_epoch.py` | per-sector depths, splits, look-elsewhere, injection summary, geometry, correlates | out/depth_epoch.json, out/per_transit_correlates.csv, out/inject_summary.csv |
| `c08_phasescan.py` | phase scans (with and without S12), secondary | out/phase_scan*.csv, out/phase_scan.json |
| `c09_apertures.py` | aperture vs PRF weightings (`--from-csv` re-summarises) | out/apertures.json, out/apertures_per_transit.csv, out/aperture_prf_fractions.json |
| `c10_export.py` | final CSVs | lc_M3V_INV_clean.csv, lc_M3V_REF_clean.csv, lc_M3V_REF_tiled.csv |
| `c11_modes.py` | same per-transit measurement on every light-curve version; leak regression | out/modes_compare.json, out/modes_per_transit.csv |
| `c12_tile_transfer.py` | transfer of the tiled REF curve at the RRab period | out/tile_transfer.json |
| `c13_figs.py` | figures | figs/*.png |
| `c14_final.py` | adopted depth, sys budget, k/Rp, odd/even | out/final_depth.json |
| `c15_leak_summary.py` | leak table | out/leak_summary.json |
| `c16_allreg_check.py` | contrast change between regressor sets vs fake-epoch null | out/allreg_check.json |
| `e01_residual_rr.py`, `e02_rr_pixmap.py`, `e03_scan_spike.py` | exploratory: S93 residual, its pixel pattern, the REF scan spike | stdout |

Run order: c00, c01, c02, c03, c05, c06, c04, c07, c08, c11, c12, c09, c14, c15, c16, c10, c13. c06 writes the
ephemeris used by c08-c11 and c10; after it changed, c08, c09, c10 and c11 were rerun with the final ephemeris. Logs
of every run are in `logs/`.

## Caveats

* **Masks:** 'all' masks 168-193 of 441 pixels per sector. The pointing fit then rests on 222-245 pixels,
  against 411-414 in phase 1.
  The pointing noise did not grow measurably; per-cadence scatter changed by -1.5 to +2.4%.
* **Injections:** they use PRF models. They cannot test a true-PSF error larger than the 0.07-px shift and
  0.15-px blur explored, or a scene-specific error of the Mira model near the M3V. The static Mira PRF residual is
  ~0.9 x F_M3V near the M3V in S103/S104. That is the one instrumental difference between the early and late
  sectors.
* **Window acceptance:** for the 30-min S12 data it depends on the exact ephemeris at the edge (18 vs 19 INV
  windows between c03/c04 and c11). S12 carries little weight.
* **REF E58:** the event is retained in REF (strict by the referee's background estimator) but fails the INV
  background cut. It is deep (35 ppt) and noisy; its influence is shown separately (REF without E58).
* **Free-density chains:** these are shorter than 50 autocorrelation times (90-176 steps). Medians are reliable,
  the tails approximate.
* **Not tested:** chromatic (colour-dependent) PSF differences between the M3V and field stars beyond the shifts
  and blur above. The VSX Mira period (294 d) and epoch were taken as given.
