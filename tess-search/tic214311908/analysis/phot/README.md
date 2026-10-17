# Deblended photometry of TIC 214311908 (M3V) next to the Mira TIC 214311906

Folder: `tic214311908/analysis/phot/`
(AI-run analysis, 2026-10-16. Nothing was submitted, uploaded or posted anywhere.)

## Bottom line

* A difference-image PSF photometry of the 21x21 px TESScut cubes (S12, S39, S66, S93, S103, S104) separates the
  M3V (TIC 214311908), the M-dwarf pair (TIC 214311910 A+B) and the Mira. The 5.6915 d dip is found in the M3V
  light curve in **every sector** with a constant flux loss of **15.4 +- 0.8 e-/s** (model B; sector chi2 6.6/5);
  depth **16.4 +- 0.9 ppt** if the M3V gives 935 e-/s (TIC Tmag 13.0135), **18.4 ppt** for the Gaia-colour estimate
  (837 e-/s). The model-choice spread is +-0.8 e-/s (8 variants).
* The deblended pair light curve shows **-0.3 +- 3.9 e-/s** at the transit times, the Mira **-0.8 +- 1.5 e-/s**.
  Injections into the pixel data show that this pipeline returns 100 +- 1% of a dip on the M3V, and 0-2% of it in the
  pair (model A) and Mira light curves. A forward model of six light curves per sector, using injection-calibrated
  responses and an empirical noise covariance, gives chi2 = 38.9/35 for the M3V as host, 230.8/35 for the pair and
  407.6/35 for the Mira.
* The in-transit difference image is fitted by the M3V PRF (total chi2 172.5 for 191 dof; pair 290.1, Mira 382.9).
  The free-position centroid is (-1.7 +- 1.4, +0.0 +- 1.1)" from the M3V, from noise alone. The pair lies 11.6" E
  and the Mira 18.6" E.
* **Open:** the faint background star Gaia DR3 5925373272827260544 (G 17.48, T ~16.8, ~28 e-/s), ~3" NE of the M3V,
  cannot be separated by TESS. It would need ~55% deep eclipses. The centroid disfavours it only weakly: chi2 = 12.4
  (2 dof) with statistical errors alone, 4.5 after adding 1.5" of systematics.
* **Marginal odd/even difference:** odd transits are deeper than even ones by 4.3 +- 1.9 ppt (2.3 sigma, per-transit,
  model B) to 5.0 +- 1.6 ppt (3.2 sigma, phase-folded). Model A gives 1.5 sigma. Per sector the difference comes almost
  entirely from S103, where the only even event (E56) is shallow in model B (4.9 +- 5.8 ppt) but normal in model A
  (14-18 ppt). It also comes partly from S12. S39, S66, S93 and S104 give 0.6, -6.0, 3.5 and 5.3 e-/s
  (+-3-4 e-/s). Not robust, but it should be re-tested.
* **No secondary eclipse:** -0.17 +- 0.93 ppt at phase 0.5. The phase scan's maximum is S/N 2.6, which noise
  explains (scan S/N rms 1.15). Flat-bottomed: a V shape is worse by delta-chi2 = 18.6.
* **Timing:** 19 deblended transit times fit a linear ephemeris with chi2 = 11.5/17 when the errors are calibrated by
  injection, or 28.6/17 with formal or prayer-bead errors. Either way the O-C rms is 8.0 min. The QLP excess
  (53.4/11) comes from underestimated QLP per-event errors (median 4.5 min against a true scatter of ~8 min), not from
  TTVs. The QLP time stamps agree with TESScut to within 0.25 min after the S103/S104 correction. Without it they are
  -5.17 and -7.41 min off, an independent confirmation of that correction. One event, S66 E-128, is 21 min early
  (2.3 sigma with calibrated errors, 3.8 sigma without). QLP puts the same event 19 min later than we do.
* **Corrections to earlier statements:**
  1. The Mira's TESS image is **not** shifted by 16-17" by its colour. Its variable-light image
     (d pixel / d F_Mira) sits within 0-3.4" of its Gaia position in five sectors and 7" in S66, where it hardly
     varies (0.3%). The earlier wrong WCS came from a degenerate direct-image fit (Mira PSF mismatch plus the pair 6"
     away; see p03_scene: static fits give negative fluxes), not from a real colour shift. The conclusion that
     calibrating on the Mira was wrong still stands; its explanation should change.
  2. QLP's strongest event (S66 E-129, 2.7% deep, S/N 6.6) falls in strong scattered light, with background
     +300 to +700 e-/s/px. With relaxed cuts we measure 28.8 +- 5.7 ppt there, about 1.8x the mean. It should not
     anchor any depth estimate.
  3. "Odd/even within 1.9 sigma" (QLP) becomes 1.5-3.2 sigma here, depending on the model (see above).

## Method (pipeline)

1. `p01_gaia.py`: Gaia DR3 within 360" (8916 sources); TESS magnitudes from Stassun+19 colour relation.
   `p00_timecorr.py`: target-specific barycentric correction of TESScut times from JPL Horizons
   (-19 s ... +31 s per sector, applied).
2. `p02_wcs.py`: per-sector WCS shift from a global PRF fit to the field (all Gaia stars; Mira +-3 px masked).
   Shifts are -0.064...+0.048 px (x) and -0.013...+0.103 px (y), +-0.003-0.006 px. Individual normal-colour stars
   agree to 0.01-0.04 px. PRF extra blur is 0-0.03 px. Field residual rms is 1.4-1.8% of the flux.
3. `p03_scene.py`: static scene fits of the median image. They are **degenerate**: Mira, pair and M3V fluxes swing
   to unphysical, even negative, values, because the SPOC PRF does not describe the Mira (chi2/dof 6-800 against its
   empirical image). Not used for photometry; it documents why.
4. `p03b_mirapsf.py`: empirical Mira image b = d pixel / d F_Mira from the Mira's 0.3-4.6% in-sector variability
   (pixel regression on the Mira flux proxy, background and pointing). Its position offset is in
   `out/mirapsf_s*.json`.
5. `p04_jitter.py`: per-cadence pointing offsets from the field (rms 0.006-0.05 px) and the background level.
6. `p05_phot.py` (core): the faint-pixel background plane is subtracted with unit response. A pixel regression in
   1-h bins over quality==0 or 2048 cadences, with background excess < 150 e-/s/px and the real (and injected)
   transit windows masked, fits a + b * dF_Mira + c * dB + d * dx + e * dy. The difference image
   D = pixel - a - c dB - d dx - e dy keeps the Mira's variation and all short-term changes. A per-cadence weighted
   linear fit of D, over pixels within 3 px of the Mira or the M3V, uses three template sets:
   * **B** (primary for the M3V and Mira): Mira empirical image b + M3V PRF + plane.
   * **Bp** (pair light curve, M3V not fitted): b + pair PRF + plane.
   * **A** (all free): b + M3V PRF + pair PRF + plane. The templates are strongly correlated (pair-Mira -0.95,
     pair-M3V -0.86), so A is 2-3x noisier.

   The pair A+B are merged with B/A = 0.367 from Gaia G, at their flux-weighted position (0.83" apart).
7. `p06_lc.py`: times -> BJD_TDB for the target. quality_ok = finite & QUALITY in {0, 2048} & background excess
   < 60 e-/s/px & pointing within 6 robust sigma & per-cadence chi2 < 3x median; the background cut is set from
   `e09_bkgcut.py`. The flux is normalised to the reference fluxes (M3V 935 e-/s; pair 895 e-/s; Mira = mean
   aperture flux minus the M dwarfs). It is detrended with a wotan biweight (0.75 d) with the transit windows
   +-0.1 d masked. NB: the wotan biweight is biased for series centred on 0, so detrending is always done on the
   normalised flux (~1).

Per-cadence scatter of the M3V (model B): 8.6 ppt (S12, 30 min), 7.8 ppt (S39, 10 min) and 13.9 / 14.7 / 17.7 /
16.2 ppt (S66 / S93 / S103 / S104, 200 s), about 3.5-4 ppt per hour. Model A: 14-35 ppt per cadence.

## Key numbers

| Quantity | Value | Method | Uncertainty | File / script |
|---|---|---|---|---|
| Flux lost in transit, M3V (model B), all sectors | 15.4 +- 0.8 e-/s | linear trapezoid fit at the global ephemeris and shape | max(white, phase-scramble 150 random phases); model spread +-0.8 e-/s | out/per_sector_depths.csv (p08), out/variants.json (p13) |
| Same, per sector S12/39/66/93/103/104 | 14.8+-3.9, 13.0+-1.6, 12.9+-2.6, 15.9+-1.7, 17.5+-2.3, 18.8+-2.2 e-/s | as above | chi2 about the mean 6.6/5 | out/per_sector_depths.csv |
| Depth on the M3V | 16.4 +- 0.9 ppt (F = 935 e-/s) | above / reference flux | +5% systematic per 0.05 mag of Tmag; 18.4 ppt if F = 837 e-/s (Gaia colour) | p06/p08 |
| Model A M3V | 16.1 +- 1.8 e-/s | as above | | out/per_sector_depths.csv |
| Pair (model A) at transit | -0.3 +- 3.9 e-/s | as above | | out/per_sector_depths.csv |
| Mira (model B / A) at transit | -0.8 +- 1.5 / -1.0 +- 3.1 e-/s | as above | | out/per_sector_depths.csv |
| Injection recovery, dip on M3V -> B:M3V | 1.000 (sector means 0.996-1.005) | 15 e-/s trapezoids in the pixel cube, 7 fake ephemerides x 6 sectors, full pipeline rerun | scatter 0.003-0.019 | out/inject_summary.json (p07, p07b) |
| Cross-talk M3V -> A:pair, A:Mira, B:Mira | -0.017, +0.011, -0.001 | same | +-0.07, +-0.05, +-0.01 | same |
| Cross-talk pair -> B:M3V, B:Mira, A:M3V | 0.40, 0.52, 0.02 | same | +-0.07, +-0.14, +-0.04 | same |
| Cross-talk Mira -> B:M3V, A:M3V | 0.005, 0.012 | same | +-0.013, +-0.027 | same |
| Template-mismatch injection (M3V PSF blurred 0.15 px, shifted 0.14 px) | B:M3V 0.88-1.14, A:M3V 0.74-1.31, A:pair -0.54...+0.40 | same | systematic | same |
| Host test chi2 (6 LCs x 6 sectors) | M3V 38.9/35 (p 0.30); pair 230.8/35; Mira 407.6/35; flux loss 14.8 +- 0.7 e-/s | forward model with injection responses, empirical 6x6 covariance from 300 fake epochs | | out/host_test.json (p09) |
| In-transit difference image chi2 (6 sectors) | M3V 172.5, pair 290.1, Mira 382.9, free position 167.5 / 191 dof | PRF fits to the stacked in-out image; errors from 200 fake epochs | | out/validate_diffimg.json (p14) |
| Dimming centroid relative to the M3V | E -1.68 +- 1.36", N +0.01 +- 1.06" | weighted mean of six free-position fits | statistical only; add ~1-2" systematic | out/centroid_combined.json (p14b) |
| Global ephemeris (deblended) | T0 = 3840.91414 +- 0.00134 BTJD_TDB, P = 5.6915087 +- 0.0000090 d | trapezoid least squares on all model-B data | formal | out/transits.json (p08) |
| Ephemeris from the 19 transit times | T0 = 3840.91272 +- 0.00175, P = 5.6915092 +- 0.0000122 d | weighted linear fit, injection-calibrated errors | | out/timing.json (p10) |
| T14 / flat fraction | 1.83 +- 0.12 h / 0.47 +- 0.09 (another start gives 1.65 h / 0.63: shallow chi2 surface) | trapezoid | formal | out/transits.json |
| Flat vs V | V worse by delta-chi2 = 18.6 | refit with flat = 0 | | p08 |
| Secondary at phase 0.5 | -0.17 +- 0.93 ppt | linear fit | scan-calibrated (robust rms of the S/N scan = 1.15) | out/transits.json |
| Phase scan max | S/N 2.62 at phase 0.300 (1.9 ppt) | all phases, step 0.002 | consistent with noise | out/phase_scan.npy |
| Odd - even (model B) | per-transit +4.3 +- 1.9 ppt (2.3 sigma; jackknife 1.8-2.7 sigma); folded +5.0 +- 1.6 ppt (3.2 sigma) | | per-transit depth chi2 28.1/18 | out/transits.json |
| Odd - even (model A) | +5.0 +- 3.4 ppt (1.5 sigma) | | | out/transits.json |
| Per-transit timing error scale (injection pulls) | std 1.16-2.32 (robust 0.82-1.19), with 3-7% outliers > 3 sigma | 60 synthetic transits per sector in the final LC | | out/timing.json |
| Timing chi2 (ours) | 11.5/17 calibrated; 28.6/17 formal or prayer-bead; O-C rms 8.0 min | | | out/timing.json, figs/oc.png |
| QLP vs ours, 12 common events | mean -6.8 min, rms 10.1 min, chi2 19.4/12 | | | out/timing_table.csv |
| QLP - TESScut time-stamp lag (corrected QLP) | -0.25...+0.05 min (all sectors) | cross-correlation of high-passed systematics | ~0.1 min | out/public_timing_lag.json (p11) |
| QLP - TESScut lag, uncorrected S93 / S103 / S104 | -0.04 / -5.17 / -7.41 min | same | | out/public_timing_lag_uncorrected.json (p11b) |
| Mira image offset from Gaia (E, N), S12/39/66/93/103/104 | (1.2,1.2), (2.7,-2.8), (-1.5,-7.1), (0.0,-1.4), (-0.7,1.3), (-3.4,-0.2)" | PRF fit to the variable-light image | +-0.5-1.3" (S66 weakest) | out/mirapsf_s*.json (p03b) |

## Per-transit results (model B, strict quality; `out/per_transit_M3V.csv`)

Errors listed are max(formal, prayer-bead); multiply the timing errors by the sector's injection factor
(`out/timing.json`, injection_calibration.pull_std) for honest values (column err_cal in `out/timing_table.csv`).
O-C is relative to the global ephemeris above.

| Sector | Epoch | O-C (min) | sigma_t (min) | depth (ppt) | S/N | note |
|---|---|---|---|---|---|---|
| 12 | -387 | +1.9 | 18 | 23 +- 12 | 1.9 | 30-min data, 2 in-transit points |
| 12 | -385 | +7.7 | 13 | 19 +- 6 | 3.0 | |
| 39 | -259 | +4.1 | 7 | 14.4 +- 3.0 | 4.9 | |
| 39 | -258 | +2.5 | 8 | 12.9 +- 3.0 | 4.3 | |
| 39 | -257 | -5.9 | 7 | 14.6 +- 4.1 | 4.8 | |
| 66 | -130 | +1.1 | 8 | 11.5 +- 3.1 | 3.8 | |
| 66 | -128 | -23.2 | 5 | 18.4 +- 3.4 | 5.8 | timing outlier |
| 66 | -127 | +2.8 | 20 | 11.7 +- 5.2 | 2.2 | 60% coverage |
| 93 | -1 | -5.9 | 8 | 16.8 +- 4.3 | 4.7 | |
| 93 | 0 | +3.3 | 9 | 15.8 +- 4.6 | 4.3 | |
| 93 | 1 | -2.0 | 5 | 20.2 +- 3.8 | 6.0 | |
| 93 | 2 | -3.3 | 10 | 14.3 +- 3.4 | 4.2 | |
| 103 | 55 | +8.5 | 5 | 22.3 +- 3.9 | 5.7 | |
| 103 | 56 | +11.9 | 29 | 4.9 +- 5.8 | 1.2 | shallow in B; 14-18 ppt in model A |
| 103 | 59 | -0.5 | 6 | 28.0 +- 5.8 | 6.7 | |
| 104 | 60 | +7.1 | 7 | 15.9 +- 3.7 | 4.3 | |
| 104 | 61 | +9.7 | 7 | 23.4 +- 4.9 | 6.3 | |
| 104 | 62 | -8.7 | 5 | 25.9 +- 5.1 | 5.1 | |
| 104 | 63 | -5.6 | 6 | 25.7 +- 4.2 | 6.9 | |
| relaxed | 39 / -256 | +3.7 | 8 | 14.6 +- 3.9 | 4.2 | partly in elevated background |
| relaxed | 66 / -129 | +14.7 | 4 | 28.8 +- 5.7 | 7.0 | strong scattered light; QLP's "best" event |
| relaxed | 12 / -386 | +0.6 | 15 | 17 +- 7 | 2.4 | |

Relaxed fits (all cadences with background excess < 800 e-/s/px, raw flux with a local linear baseline) are flagged
`relaxed=True` in the CSV. Further relaxed entries (S12 E-388, S39 E-255, S103 E58) are dominated by systematics and
should not be used.

## Public light curves (p11; `out/public_depths.csv`)

* No TGLC, TESS-SPOC, SPOC 2-min, CDIPS or PATHOS light curves exist for these stars; MAST was queried 2026-10-16
  (`data/mast/mast_timeseries.csv`).
* **QLP** (S103/S104 times corrected), fractional depth in QLP's own normalisation, M3V DET:
  12.9+-3.7, 14.9+-2.0, 19.1+-1.8, 13.0+-2.2, 10.7+-1.6, 14.0+-1.5 ppt (weighted 14.1 +- 0.8, chi2 13.3/5).
  QLP's pair light curve (17.9 ppt) and Mira light curve (2.4 ppt) show the dip too, because QLP does not deblend.
* **TARS** (S12-S93, M3V): 2.75-3.35 ppt, constant (chi2 1.1/3). Its normalisation includes the Mira's light.
  TARS lists periods of 8.1-8.7 d.
* **GSFC-ELEANOR-LITE** S12: the three stars' light curves are essentially identical (the same pixels), with no
  detectable dip (S/N < 0.2). Not usable.

## Variability (p12; `out/variability.json`, `figs/variability.png`)

* No credible rotation period. The largest Lomb-Scargle peaks of the systematics-regressed M3V flux (4.5-10 d) change
  from sector to sector and match the peaks in the pair's light curve (S12 4.97/4.97 d, S39 6.6/6.8, S66 10.1/9.7,
  S93 9.4/9.8, S103 9.43/9.43). That points to residual Mira light or systematics common to both. A control star
  shows 2-6 ppt systematics at 5.7-9.7 d. Slow (> 1 d) variability below ~1-2% cannot be measured with this
  photometry. The TARS 8.1-8.7 d periods come from a Mira-dominated aperture and are not confirmed.
* Flares: 5 candidates in 110.8 d of good data. Four are in S12 (single 30-min points of 5-6 sigma, 44-57 ppt, plus
  one 4-point event at BTJD 1634.41), and one single 200-s point in S104 (4203.016, 82 ppt, 5.1 sigma). None is
  convincing. The 5-sigma per-cadence threshold is ~40-90 ppt.

## Evidence summary for the planet hypothesis

For:
* The dip appears in the deblended M3V light curve in all six sectors (2019-2026) with a constant flux loss
  (chi2 6.6/5).
* It is absent from the deblended pair and Mira light curves at the level the injections predict for an M3V host.
* The forward-model host test and the in-transit difference images favour the M3V over the pair and the Mira by
  delta-chi2 > 100.
* The dip is flat-bottomed, with no secondary eclipse (< 1.7 ppt at 2 sigma at phase 0.5) and a linear ephemeris
  without significant TTVs.

Against or open:
* The marginal odd/even difference (1.5-3.2 sigma, model-dependent, driven by S103 E56 and S12).
* The ~3" background star Gaia DR3 5925373272827260544 is not excluded by TESS.
* The depth normalisation rests on the catalogue Tmag of the M3V; the static flux of the M3V cannot be measured next
  to the Mira.
* Model B assumes the pair is constant: a pair flare or variability leaks 40% into the M3V. Model A avoids this but
  is noisier.

## Output files

Top level (for downstream use):
* `lc_M3V_all.csv`: deblended M3V light curve, all six sectors, all finite cadences.
  * Primary columns: time_btjd_tdb (BJD_TDB - 2457000, barycentred for the M3V); flux_norm (detrended, model B,
    1 = out of transit); flux_err (sector robust sigma x sqrt(local chi2 ratio >= 1)); sector; cadence_s;
    quality_ok (1 = use).
  * Further columns: flux_raw_norm (= 1 + dF/935 e-/s, not detrended); trend (biweight trend divided out);
    dflux_e_s (dF of the M3V in e-/s, model B); flare_flag; flux_norm_modelA, flux_err_modelA,
    flux_raw_norm_modelA (all-free model); in_transit_window (+-0.1 d of the reference ephemeris); tess_quality;
    bkg_e_s_px (faint-pixel background); jitter_dx_px, jitter_dy_px; time_btjd_tesscut (original TESScut TIME).
* `lc_pair_all.csv` and `lc_Mira_all.csv`: the same layout. The pair's flux_norm is model Bp, which absorbs the M3V
  transit by construction (1.25-2.3x). **Use flux_norm_modelA for the pair.** The Mira's flux_norm is model B and
  flux_norm_modelA is model A. Both have dflux_e_s and dflux_e_s_modelA.
* `lc_allstars_s{012,039,066,093,103,104}.npz` (copies of data/lc_s*.npz): per sector, all stars and models.
  * Common keys: t_tdb, t_tess, q, cadn, quality_ok, fin, bkg, dx, dy, intransit_window, exptime_s, sector.
  * Per star and model, with suffix M3V_B, M3V_A, pair_Bp, pair_A, Mira_B or Mira_A: dF_* (e-/s), raw_*
    (normalised), trend_*, det_* (detrended), err_*, flare_*.

`out/`:
* `wcs_s*.json`: WCS shift, blur, field star fluxes, per-star offsets.
* `scene_s*.json`: static scene fits (degenerate; documentation only).
* `mirapsf_s*.json`: Mira variable-image fits.
* `phot_s*.json`: per-sector photometry metadata (model template correlations, formal errors, static fits).
* `lc_stats.json`: per-sector scatter and reference fluxes.
* `transits.json`: global fit, odd/even, secondary, variability first pass, flare first pass.
* `per_sector_depths.csv`: every light curve and sector, white and phase-scramble errors.
* `per_transit_M3V.csv`: per-transit fits, models B and A, strict and relaxed.
* `inject_s*.json`, `inject_all.csv`, `inject_summary.json`: pixel injections and cross-talk matrices.
* `host_test.json`: forward-model host test.
* `timing.json` and `timing_table.csv`: timing chi2, injection calibration, QLP comparison.
* `public_depths.csv`, `public_timing_lag*.json`: public light curves and the time-stamp check.
* `variability.json`: Lomb-Scargle with control stars; flare list.
* `variants.csv` and `variants.json`: systematics of the depth from model choices.
* `validate_diffimg.json` and `centroid_combined.json`: in-transit difference images and centroid.
* `timecorr.json`: barycentric corrections.
* `phase_scan.npy`: columns phase, depth, error.
* `*.log`: console output of each script.

`data/`:
* `gaia/gaia_360arcsec.csv`, `gaia/tic_three.csv`: catalogues.
* `horizons/`: JPL Horizons TESS vectors.
* `jitter_s*.npz`: pointing and background per cadence.
* `mirapsf_s*.npz`: Mira variable image and its error.
* `phot_s*.npz`: raw per-cadence fluxes of every model, chi2, quality masks, a and b images.
* `lc_s*.npz`: light curves.
* `mast/`: MAST product list and downloaded TARS / eleanor-lite files.

`figs/`:
* `lc_sectors.png`: 1-h binned, detrended light curves of the three stars per sector.
* `fold_by_sector.png`: per-sector folds in e-/s for every light curve.
* `fold_all.png`: all-sector folds plus the phase scan.
* `transits_M3V.png`: every transit with its fit.
* `oc.png`: O-C, ours against QLP.
* `diffimg_s*.png`: in-transit difference images and residuals per host.
* `mirapsf_s*.png`: Mira variable image and PRF residuals.
* `wcs_s*.png`: field residuals and per-star offsets against colour.
* `scene_s*.png`: degenerate static fits.
* `variability.png`: periodograms and control stars.

`explore/`: exploratory plots from early, superseded pipeline versions. Not for use: e03/e06 show the wotan bias
around zero and the earlier per-pixel background regression.

## Scripts (`scripts/`; run with /home/user/work/venv/bin/python, OMP_NUM_THREADS=2)

Libraries:
* `photlib.py`: Cutout loader, Gaia propagation, a fast SPOC-PRF evaluator (PRF package model files, optional
  Gaussian blur), Horizons barycentric correction, helpers.
* `tlib.py`: exposure-integrated trapezoid, linear depth fits, and per-transit fit with formal, prayer-bead and
  chi2-profile errors. Note that it fits offsets dtc rather than absolute times, because scipy's diff_step is
  relative.

Pipeline, in run order:
* `p01_gaia.py`: Gaia query.
* `p00_timecorr.py`: barycentric corrections.
* `p00_mastquery.py`: public light curves list.
* `p02_wcs.py`: WCS from field stars.
* `p03_scene.py`: static scene fits (shows the degeneracy).
* `p04_jitter.py`: pointing and background.
* `p03b_mirapsf.py`: Mira variable image.
* `p05_phot.py`: core photometry. Options: inject, extra_mask_times, prf_perturb, rfit, mira_template,
  free_others, controls, return_D, jitter_order. A second-order jitter test (`e11_jitter_order.py`) changed nothing.
* `p06_lc.py`: light curves and CSVs.
* `p07_inject.py`, `p07b_summary.py`: pixel injections and cross-talk.
* `p08_transits.py`: global fit, per-sector, per-transit, odd/even, secondary, scan, variability first pass.
* `p09_host.py`: host test.
* `p10_timing.py`: timing, injection calibration, QLP comparison.
* `p11_public.py`, `p11b_lagcheck.py`: public products and the time-stamp lag.
* `p12_variability.py`: rotation and flares with control stars.
* `p13_variants.py`: model-choice systematics.
* `p14_validate.py`, `p14b_centroid.py`: in-transit difference images and centroid.
* `p15_figs.py`: figures.

Exploratory (kept for the record):
* `e01_look.py`: images.
* `e02_ab.py`: mean against variable image.
* `e03_quick.py`, `e06_quick2.py`, `e04_dbg.py`: early light-curve looks, affected by the wotan bias.
* `e05_cond.py`: template correlations.
* `e07_rfit.py`: fit radius.
* `e08_qual.py`, `e09_bkgcut.py`: quality cuts.
* `e10_depths.py`: first depth table.
* `e11_jitter_order.py`: quadratic pointing terms.

## Caveats

* The SPOC PRF fits field stars to 1.4-1.8% rms but the Mira very badly. The Mira is modelled by its empirical
  variable image, which also absorbs any colour-dependent change of its PSF over the sector. In S66 the Mira varies
  by only 0.3%, so its template there is noisy (that sector's model-A errors are larger).
* Model B (primary) does not fit the pair. Any pair variability leaks into the M3V curve at 33-52%. Model A is free
  of this but has 2-3x the noise.
* The normalisation (935 e-/s) is catalogue-based. The static M3V flux could not be measured (p03 / p05 static fits
  are degenerate).
* Errors: phase-scramble and injection calibration are used where formal errors are unreliable. Timing errors are
  non-Gaussian, with 3-7% outliers.
* Quality cuts drop 45% of S66 and 50% of S12 (scattered light). S66 E-129 and S39 E-255 are only in relaxed data.
