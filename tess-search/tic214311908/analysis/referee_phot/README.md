# Referee check of the 'phot' investigation (deblended photometry of TIC 214311908)

Folder: `tic214311908/analysis/referee_phot/`
(AI-run analysis, 2026-10-17. Nothing was posted, uploaded or submitted anywhere.)

The investigator's work is in `../phot/`. Everything here was re-derived with separate code. Only the
investigator's deliverables (`lc_M3V_all.csv`, `out/per_transit_M3V.csv`, `out/timing_table.csv`, `out/wcs_s*.json`,
`out/timecorr.json`) were read, to audit them. The referee's photometry is a different method: local windows around
each transit, with per-pixel baselines and PRF fits to the residual frames, and no global Mira template.

The scripts r00-r16 were written in an earlier run of this referee session, which was interrupted by a machine
restart. In this run they were all re-read and re-run from r01 to r15 (`scripts/run_chain.sh`, logs in `logs/`).
The outputs reproduced to better than 1e-6. r16 was re-run after its variant list was fixed, and r17-r23 are new.

## Verdict on the investigation (short)

**Upheld:**
* The dip is on the M3V or within about 3" of it, and not on the pair or the Mira.
* It is present in all six sectors at a constant flux loss.
* No secondary eclipse.
* Flat-bottomed.
* No significant TTVs.
* The QLP S103/S104 time-stamp offsets.
* The Mira's image is not shifted by 16-17".

**Revised or refuted:**
1. **Depth normalisation is inconsistent (revised).** The flux loss (15.4 e-/s) is measured in units of the PRF fit,
   but it was divided by the *catalogue* flux of the M3V (935 e-/s).
   * The investigator's own PRF-fitted field-star fluxes are 1.107 +- 0.011 x the catalogue flux (r17).
   * The M3V's flux in PRF units is 979 e-/s from Gaia XP synthetic photometry calibrated on field stars (r09), and
     968-1035 e-/s from the other routes.
   * Depth = **15.7 +- 0.8 (stat) +- 1.0 (sys) ppt**, not 16.4 ppt. The quoted alternative of 18.4 ppt (837 e-/s) is
     wrong, because 837 e-/s is also on the catalogue scale.
   * The XP and SED-synthetic Tmag of the M3V is 13.08-13.10, against 13.0135 in the TIC.
2. **The odd/even difference is not reproduced (refuted as a signal).**
   * Independent photometry at fixed ephemeris: odd - even = -0.6 +- 1.5 e-/s (-0.4 sigma). The 12 method variants
     all give |odd - even| <= 1.1 sigma (r16).
   * On the investigator's own event list: +1.1 +- 1.6 e-/s (0.7 sigma).
   * The investigator's +4.3 ppt (2.3-3.2 sigma) is a property of their model-B light curve. Measured with the
     referee's code on their CSV it gives +5.5 +- 2.4 ppt (2.2 sigma) with empirical errors, and 3.8 sigma with white
     errors (r20).
   * Per-event depths from the two photometries of the same pixels differ by 4.5 e-/s rms (correlation 0.59). The
     differences are largest for S104 E60/E63 and S103 E59, the events that make up the investigator's odd/even
     difference.
   * Their quoted per-transit depth errors are about 1.3x smaller than the scatter of their own light curve at fake
     epochs: 4.3-6.0 ppt per transit.
3. **The S66 E-128 "timing outlier" is an artefact (refuted).**
   * The binned data in the two light curves are nearly identical. The referee photometry gives O-C +1.4 min, and
     the investigator's own model A gives +7.2 min.
   * Only their model B, fitted with the same code, gives -21.6 min. Single-event timing is not robust below
     ~10-20 min (r21).
4. **"Cross-talk <= 2%" holds only by construction.**
   * The investigator injected exactly the PRF that the fit uses, so with matched templates linear least squares
     returns the injected flux exactly. The referee's matched injections show the same.
   * Their own template-mismatch test moves model-A pair cross-talk to -54...+40%.
   * The exclusion of the pair and the Mira should rest on the difference-image centroid. That exclusion is robust:
     the pair is 11" (> 6 sigma) and the Mira 18.6" from the dimming position.
5. **Minor:** QLP S103/S104 times were corrected only by the FFI TIMECORR median (+5.17 / +7.47 min).
   * QLP S103/S104 TIME equals TESScut (TIME - TIMECORR) - 4.1 s.
   * The full correction to target-barycentric BJD_TDB is therefore +5.74 / +7.27 min, leaving residuals of
     +0.57 / -0.20 min. Negligible for timing.
6. **Mira variable-image offsets** are model-dependent. The PRF fits to the Mira have chi2/dof 40-1300, and the
   referee finds offsets of 0.6-7.6" (S66: 10.6 +- 3.4"). The claim "within 0-3.4"" is overstated, but there is no
   16-17" shift.

**Open (agreed with the investigator, quantified further):**
* The background star Gaia DR3 5925373272827260544 (bg3, 3.1" NE, F = 30 e-/s in PRF units) would need **~50%** deep
  eclipses. It is disfavoured only at about 2-3 sigma:
  * Referee centroid against bg3: chi2 25 (stat) or 5.3 (with 1.5" systematic), for 2 dof.
  * Injection-calibrated: Mahalanobis^2 14.9 for bg3, against 2.4 for the M3V.
  * Photometric host test: delta-chi2 = 9.9.
  * A 0.1 px PRF/WCS error towards bg3 would erase the photometric discrimination (r19). **bg3 is not excluded by
    TESS.**
* A second faint neighbour, Gaia DR3 5925373272827262208 (G 17.68, 5.4" NNE), would need 68% eclipses. It is
  disfavoured at p = 0.001 with 1.5" systematic (r22). The investigator did not list it.
* Sectors S103/S104 are deeper (17.0 +- 1.3 e-/s) than S39/S66 (13.7 +- 1.4 e-/s), 1.7 sigma. This is not
  significant (all-sector chi2 3.3/5 to 7.1/5) but is worth watching.
* Transit shape depends on the photometry. The same fitting code gives T14 = 1.83 +- 0.12 h, flat fraction 0.47 on
  the investigator's light curve, and 1.57 +- 0.07 h, flat fraction 0.73 (flat 0.6-0.9 within delta-chi2 2.4) on the
  referee's. Both reject a V shape: delta-chi2 = 18 and 38.

## Key referee numbers

| Quantity | Referee value | Method | Script / file |
|---|---|---|---|
| Flux loss on M3V (M3V-only fit) | 14.9 +- 0.6 (global fit); 15.3 +- 0.9 e-/s (sector mean, empirical errors), sector chi2 3.3/5 | trapezoid on local-window PRF photometry; errors from 60-104 fake windows per sector | r06 `out/transit_fit.json`, `out/per_sector.csv` |
| Same, 12 method variants | 14.1-15.5 e-/s (std 0.37) | fit radius, baseline order, window, mask, PRF shift | r16 `out/variants.json` |
| PRF-template mismatch | recovery 0.91-1.06 (0.1 px shift + 0.15 px blur), 0.98 (blur only) | pixel injections | r07 `out/inject_timing.json`, r19 |
| M3V flux in PRF units | 979 +- 20 (stat) e-/s (966 without > 1020 nm tail); TIC x 1.107 = 1035; SED Tmag x 1.107 = 968 | Gaia XP synthetic photon rate calibrated on 27 field stars' PRF fluxes | r09 `out/xpflux.json`, r17 `out/invscale.json` |
| PRF flux / TIC flux, field stars | investigator's own fits 1.107 +- 0.011 (1.06-1.14 per sector); referee 1.10 | isolated T < 14 stars | r17, r01b |
| Depth on M3V | **15.7 +- 0.8 (stat) +- 1.0 (sys) ppt**; investigator's 15.4 e-/s / 979 = 15.7 ppt | flux loss / F(PRF units) | r17 |
| Three-star fit at transit | M3V 13.8 +- 1.6, pair -0.6 +- 4.0, Mira 2.1 +- 2.7 e-/s | M3V PRF + pair PRF + Mira b-image | r06 |
| Host test (6 sectors x 4 observables) | M3V chi2 27.8/23 (p 0.22); bg3 37.7/23 (p 0.028); pair 224/23; Mira 387-394/23 | injection responses + empirical covariance from fake windows | r19 `out/crosstalk.json` |
| Difference-image fits | chi2 M3V 232.5, bg3 247.2, pair 357.8, Mira 477.7 (PRF) / 506.6 (b), free 208.8 for 189 dof | stacked flat-part residual frames | r10 `out/diffimg.json` |
| Dimming centroid from M3V | E +0.52 +- 0.89", N -1.36 +- 0.76" (sector scatter chi2 8.7, 10.0 / 5) | weighted mean of per-sector free PRF fits | r10 |
| Centroid, injection-calibrated | Mahalanobis^2 M3V 2.4, mismatched-PRF M3V 3.3, bg3 14.9, pair 49.2 (2 dof) | 40 trials x 6 sectors of injected stacks | r12 `out/centroid_inject.json` |
| Odd - even | -0.62 +- 1.47 e-/s (fixed ephemeris, empirical errors); +1.45 +- 1.69 (per-transit, calibrated errors) | | r08 `out/oddeven_timing.json`, r16, r20 |
| Secondary at phase 0.5 | +0.28 +- 0.75 e-/s = 0.29 +- 0.77 ppt | windows at phase 0.5, empirical errors | r18 `out/phase_scan_clean.json` |
| Phase scan (clean windows) | max S/N 3.27 at phase 0.5725 (2.4 e-/s); scan S/N robust std 1.16 | windows overlapping a transit removed; the uncleaned r08 maximum (3.9 at phase 0.94) is an overlap artefact | r18 |
| Ephemeris (20 transit times, calibrated errors) | T0 = 3840.91564 +- 0.00118 BTJD_TDB, P = 5.6915051 +- 0.0000092 d | linear fit | r08 |
| Timing chi2 | 47.8/18 formal; 14.0/18 calibrated; Monte Carlo p = 0.54 against the injection pull distribution | | r08, r15 `out/timing_dist.json` |
| Timing pulls (injections) | std 1.6-2.6, robust 0.7-1.15, 3-19% beyond 3 sigma | 47-70 injections per sector | r07 |
| Referee - investigator event times | rms 8.1 min (robust 5.6) on identical pixels; depths rms 4.5 e-/s | | r11 `out/compare_events.json` |
| Barycentric correction | target - FFI: -18.4, -12.6, -6.7, -7.4, +30.9, -13.1 s (S12...S104); agrees with the investigator within 0.6 s | astropy Earth ephemeris | r02 `out/timecheck.json` |
| QLP time stamps | S12-S93 already target-barycentred (+1.3 to +2.5 s); S103/S104 = TIME - TIMECORR - 4.1 s (no barycentric correction) | cadence matching | r02b `out/qlp_stamps.json` |

## Scripts (`scripts/`, run with /home/user/work/venv/bin/python, OMP_NUM_THREADS=2)

* `rlib.py`: cutout loader, Gaia propagation, SPOC PRF via `PRF.TESS_PRF.locate`, helpers. It is independent of the
  investigator's photlib.
* `r00_gaia.py`: Gaia DR3 and TIC within 330" (network). Output: `data/gaia330.csv`, `data/tic330.csv`.
* `r01_wcs.py`: per-sector WCS offset from isolated normal-colour field stars; red-star offsets. Output:
  `out/wcs_s*.json`.
* `r01b_fluxscale.py`: global field fit; ratio of PRF-fit flux to TIC flux. Output: `out/fluxscale.json`,
  `out/fluxscale_stars.csv`.
* `r02_time.py`: target barycentric correction (astropy) and comparison with the investigator. Output:
  `out/timecheck.json`, `data/tcorr_s*.npy`.
* `r02b_qlpstamps.py`: QLP time stamps against TESScut FFI, spacecraft and target times. Output:
  `out/qlp_stamps.json`.
* `r03_aux.py`: background, pointing and scale per cadence. Output: `data/aux_s*.npz`.
* `r04_mirab.py`: Mira variable-light image and its PRF position. Output: `out/mirab_s*.json`, `data/mirab_s*.npz`.
* `r04b_mirastatic.py`: static Mira position, with the M dwarfs fixed or free. Output: `out/mirastatic.json`.
* `r05_local.py`: core local-window photometry (functions are reused by later scripts; option linbase added in this
  run). Output: `data/win_s*.npz`.
* `r06_transit.py`: global trapezoid, flat-fraction profile, per-sector depths, three-star amplitudes. Output:
  `out/transit_fit.json`, `out/per_sector.csv`.
* `r07_pertransit.py`: per-transit times and depths, with injection-calibrated errors and a mismatched-template
  test. Output: `out/per_transit.csv`, `out/inject_timing.json`.
* `r08_oddeven_timing.py`: odd/even, timing chi2, phase 0.5 and phase scan. Output: `out/oddeven_timing.json`,
  `out/phase_scan*.csv`. The scan near phases 0.06/0.94 is biased; see r18.
* `r09a_getxp.py` (network) and `r09_xpflux.py`: Gaia XP synthetic TESS flux of the M3V in PRF units. Output:
  `data/ext/xp_sampled.csv`, `out/xpflux.json`. The TESS response is in `data/ext/tess_red*.dat`.
* `r10_diffimg.py`: in-transit difference images, host fits, centroid. Output: `out/diffimg.json`.
* `r11_compare.py`: event-by-event comparison of referee, investigator and QLP. Output: `out/compare_events.csv`,
  `out/compare_events.json`.
* `r12_centroid_inject.py`: injection calibration of the centroid (M3V, bg3, pair, mismatched M3V). Output:
  `out/centroid_inject.json`.
* `r13_csvcheck.py`: audit of `../phot/lc_M3V_all.csv`. Output: `out/csvcheck.json`.
* `r14_shape_crosscheck.py`: the same shape fit on the investigator's and the referee's light curves. Output:
  `out/shape_crosscheck.json`.
* `r15_timing_dist.py`: Monte Carlo of timing chi2 against the injection pulls. Output: `out/timing_dist.json`.
* `r16_variants.py`: 12 method variants of the flux loss and odd-even. Output: `out/variants.csv`,
  `out/variants.json`.
* `r17_invscale.py`: the investigator's own PRF/TIC flux scale; depth under each normalisation. Output:
  `out/invscale.json`.
* `r18_phasescan_clean.py`: phase scan with windows that overlap a transit removed. Output:
  `out/phase_scan_clean.csv`, `out/phase_scan_clean.json`.
* `r19_crosstalk_host.py`: injection cross-talk matrix (M3V, pair, Mira b and PRF images, bg3, three mismatched M3V
  images) and the forward host test. Output: `out/crosstalk.json`.
* `r20_inv_oddeven.py`: odd/even on the investigator's light curve and table, against the referee light curve on the
  same events. Output: `out/inv_oddeven.json`.
* `r21_E128.py`: S66 E-128 and the odd/even-driving events, compared in both light curves. Output: `out/E128.json`.
* `r22_neighbours.py`: census of possible hosts within 40" (required depth, centroid chi2). Output:
  `out/neighbours.csv`.
* `r23_export.py`: export of the referee local-window M3V light curve. Output: `lc_M3V_referee_windows.csv`.
* `run_chain.sh`: re-runs r01-r15 in order (no network).

## Output for downstream use

`lc_M3V_referee_windows.csv` has one window per transit (strict quality; the S12 E-388 window is a non-detection and
is kept).

Columns:
* time_btjd_tdb (target-barycentred)
* sector, epoch, window_id, cadence_s
* dflux_e_s
* flux_norm = 1 + dflux/979
* flux_err, oot
* T_M3V_e_s, T_pair_e_s, T_Mira_e_s (three-star fit)
* bkg_e_s_px

Its per-transit noise is about 1.25-1.7x lower than that of the investigator's model-B curve: the empirical
single-transit error is 3.1-3.3 e-/s, against 4.0-5.6 e-/s (`out/inv_oddeven.json` b_fake_single_err x 935). It covers only +-0.35 d around each transit (S39 +-0.4 d, S12 +-0.6 d).
