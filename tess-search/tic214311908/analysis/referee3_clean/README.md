# Referee 3 of the `clean` analysis (photometric systematics, cleaned light curves, final ephemeris) — TIC 214311908

Folder: `tic214311908/analysis/referee3_clean/`
(AI-run adversarial referee, 2026-10-17. Nothing was posted, uploaded or submitted, and nobody was contacted. The
only network use was read-only JPL Horizons queries for the TESS spacecraft ephemeris (q01).)

Inputs (read only): the TESScut cubes; the analyst's delivered light curves `../clean/lc_M3V_INV_clean.csv`,
`../clean/lc_M3V_REF_clean.csv`; `../clean/out/*` (for comparison only); the analyst's REF fake windows
`../clean/data/win_ref_fake.npz` (noise calibration of REF); `../phot/data/gaia/gaia_360arcsec.csv`, the photlib SPOC
PRF (only for encircled energies), `../search/out/gaia_variables.json`. q11 imports the analyst's pipeline code
(read-only) to run a new injection test. All outputs are in this folder (~7 MB).

## Bottom line

1. **The late-versus-early depth contrast is a property of the PRF/scene-model photometry, not of the pixels.**
   PSF-independent raw-pixel apertures (sum of raw pixels around the M3V, local time baselines, no Mira model, no
   field-wide estimators; encircled energy EE 0.76-0.96 changes by <= 2 % for 0.1-px shifts or 0.15-px blur) give
   S103+S104 minus S39+S66 = **-1.0 +- 2.4 ppt (3x3 box) and -4.9 +- 3.2 ppt (r <= 2.5 px)**, i.e. the opposite sign
   of the PRF fits (INV +3.3 to +4.1, REF +3.4 without E58). Paired at identical real and fake epochs, the INV-minus-
   aperture contrast is +6.6 to +7.4 +- 3.1-3.3 (box3) and +10.6 to +12.2 +- 3.9-4.0 ppt (r <= 2.5) in all 6 baseline/
   background/window/quality variants (2.0-3.2 sigma in the post-hoc split; split-independent omnibus chi2 4.7-10.7/5,
   p 0.06-0.45). The apertures are consistent with a constant depth (sector chi2 1.5/5 and 3.9/5). The analyst's
   statements "present in every photometric weighting" and "pipeline injections cannot produce it" are not supported:
   their weightings were all applied to Mira-model residual cubes, and their matched-template injections recover
   1.000 by construction (both pipelines are linear least squares with the injected image in the design matrix).
2. **Realistic PSF-mismatch injections give sector-dependent responses.** Injecting the transit with the empirical
   Mira PSF image (shifted to the M3V) instead of the SPOC PRF, both cleaned pipelines recover 0.92-1.53 (INV) and
   0.96-1.55 (REF) depending on the sector (late/early ratio 1.16 INV); a PRF taken 400 px away on the detector gives
   0.94-1.03 (`q11`). The response to a PSF-shape error is not bounded by the analyst's 0.07-px/0.15-px tests.
3. **Depth is probably a little deeper than adopted.** Refitting the delivered light curves with my own windows and
   noise calibration: INV 15.47 +- 1.05 ppt (analyst 14.50), REF 16.08 +- 0.97 (analyst 16.06), REF without E58
   15.69. GLS combination (rho 0.56): **15.8 +- 0.9 ppt** (15.6 without E58). Raw apertures: 16.2-16.9 ppt on the INV
   windows, 16.7-17.5 on all windows. INV depends on the fake-epoch beta set (S39 1.02 vs 1.43, S104 1.38 vs 1.12) by
   ~1 ppt and on the LD prior (wide q1/q2: INV 16.06 ppt, k 0.1136) by ~0.6 ppt; neither is in the analyst's budget,
   nor is the PRF-vs-aperture term. Suggested: **15.8 +- 0.9 (stat) +1.5/-1.1 (sys) ppt; k ~0.1185; Rp ~4.95 +- 0.24
   (stat) +0.25/-0.18 (sys) R_E.** Planet-sized either way.
4. **Ephemeris: per-transit times reproduced (0.1 min rms), but P depends on the timing-error model at ~1 sigma.**
   With sector-uniform LC-level injection-calibrated errors: combined T0 = 3840.91500-3840.91510 +- 0.0012-0.0016,
   P = 5.6915097-5.6915111 +- 0.0000080-0.0000108 d (analyst 3840.91488 +- 0.00093, 5.6915020 +- 0.0000078); global
   MCMC fits give P 5.6915045-5.6915067. The 2027 predictions differ by +1.5 to +2.2 min between error models; my
   1-sigma for 2027 is 2.6-3.5 min (analyst 1.9-2.5). Suggested for planning: T0 = 3840.9150 +- 0.0013 BTJD_TDB
   (2460840.9150 BJD_TDB), P = 5.691506 +- 0.000011 d; 1 sigma ~3 min in 2027; use +-10 min (3 sigma) margins.
5. **Upheld:** time system (target-barycentred BJD_TDB to 0.01 s, own Horizons check); RRab leak removed (0.26-0.78
   ppt residual, S93 0.78 ppt p 0.02); no TTVs (chi2 8.9/17, 12.6/19 with my errors); no odd/even difference in any
   photometry (INV +2.9 +- 2.1, REF -0.6 +- 1.8 / +2.7 +- 1.9 without E58, apertures -1.1 +- 2.2 and -0.1 +- 2.8);
   no secondary (apertures -1.7 +- 1.1 and -2.3 +- 1.4 ppt at phase 0.5); LEE p of the post-hoc split 0.31-0.64
   (recomputed with the analyst's sector numbers); shape b ~0.51-0.53.
6. **Phase 0.490-0.495 bump not in raw pixels:** -0.04, -0.36, -1.23 +- 1.08 ppt (box3) and -1.1, -1.5, -2.8 +- 1.4
   (r <= 2.5 px), against +1.2-2.0 ppt in the PRF photometries.

Net effect on planet vs false positive: favourable. The only time-variable-depth hint (late > early) does not exist
in PSF-independent photometry, odd/even and secondary are null in every photometry, and the depth is, if anything,
slightly larger (Rp ~5 R_E). The adopted ephemeris precision is somewhat optimistic.

## Scripts (`scripts/`, /home/user/work/venv/bin/python, OMP_NUM_THREADS=1)

| Script | What it does | Outputs |
|---|---|---|
| `q01_timesys.py` | TESS position from JPL Horizons (id -95), Roemer delay to the M3V vs TIMECORR; compares with the analyst's time column | out/q01_timesys.json, out/q01_addtab.npz |
| `q02_bigap.py` | raw-pixel large-aperture flux losses per real and fake transit (apertures box3, r 2.5/3/3.5 px on the M3V, r 4/5 px on the group light centre), own background and pointing, EE from the SPOC PRF; variants via env (Q2_BKG, Q2_ORDER, Q2_HALF, Q2_QUAL, Q2_TAG, Q2_PHRANGE) | out/q02_bigap{tag}.json, out/q02_bigap_per_transit{tag}.csv, out/q02_bigap_fake{tag}.csv |
| `run_q02_variants.sh` | runs q02/q03/q04 for 6 variants (no/wide background, linear baseline, +-0.40 d, +-0.25 d, QUALITY 0 or 2048) | logs/q03_*.log, logs/q04_*.log |
| `q03_bigap_summary.py` | per-sector EE-corrected depths, global depth, contrast, comparison with the analyst's per-transit depths | out/q03_bigap_summary{tag}.json |
| `q04_paired.py` | INV delivered curve vs apertures at identical real and fake epochs: per-sector differences with paired errors, contrast difference | out/q04_paired{tag}.json, out/q04_paired_per_transit{tag}.csv |
| `q05_timing.py` | own per-transit mid-times (INV, REF delivered curves), LC-level injection timing calibration, WLS ephemerides (INV, REF, combined independent/correlated) | out/q05_timing.json, out/q05_per_transit.csv |
| `q05b_eph_variants.py` | ephemeris sensitivity to the timing-error model; sector jackknife | out/q05b_eph_variants.json |
| `q06_rrab.py` | residual RRab amplitude in the delivered INV curve, 300 fake periods | out/q06_rrab.json |
| `q07_lee.py` | look-elsewhere p of the post-hoc split from the analyst's per-sector depths | out/q07_lee.json |
| `q08_global.py` | own emcee batman fits of the delivered INV/REF curves (density prior; Gaussian or wide LD; REF with/without E58) | out/q08_global.json |
| `q09_followup.py` | 2026-2027 predictions from the adopted ephemeris vs phase-1 CSV, BJD_TDB->UTC check, S/N requirements | out/q09_followup.json, out/q09_predictions.csv |
| `q10_oddeven.py` | odd/even for INV, REF, apertures, with fake-epoch errors and a permutation null | out/q10_oddeven.json |
| `q11_inject_altpsf.py` | analyst's cleaned pipelines re-run with a NON-template injected PSF (empirical Mira image shifted to the M3V; PRF 400 px away) | out/q11_inject_s*.csv, out/q11_diag_s*.json, out/q11_summary.json |
| `q12_secondary_ap.py`, `q12b_phasescan_ap.py` | secondary eclipse and 0.40-0.60 phase scan in the raw apertures | out/q12_secondary_ap.json, out/q12b_phasescan_ap.json |
| `q13_wcscheck.py` | sanity check of WCS pixel positions with bright-star centroids | out/q13_wcscheck.json |
| `q14_contrast_all.py` | contrast and sector chi2 for every photometry with this referee's errors | out/q14_contrast_all.json |
| `q15_omnibus.py` | split-independent chi2 of per-sector INV-minus-aperture differences, all variants | out/q15_omnibus.json |

## Key numbers

| Quantity | Value | Method | Uncertainty | File |
|---|---|---|---|---|
| Barycentric correction check | analyst time column = TIME - TIMECORR + Roemer(M3V): max diff 0.012 s | Horizons TESS vectors | - | out/q01_timesys.json |
| RRab residual in delivered INV | 0.70, 0.40, 0.36, 0.78, 0.26, 0.74 ppt (S12-S104); p 0.32, 0.30, 0.77, 0.02, 0.85, 0.16 | 3-harmonic fit, 300 fake periods | UL95 0.6-1.1 ppt | out/q06_rrab.json |
| Late-early contrast, raw apertures | box3 -1.0 +- 2.4; r<=2.5 -4.9 +- 3.2 ppt (variants: -0.2..-3.1, -4.0..-7.5) | fixed shape, quadratic baseline, fake-epoch errors | as quoted | out/q14_contrast_all.json, out/q03_bigap_summary*.json |
| Contrast INV minus aperture (paired) | +6.6..+7.4 +- 3.1-3.3 (box3); +10.6..+12.2 +- 3.9-4.0 (r<=2.5) | same epochs, paired fake differences | omnibus chi2 4.7-10.7/5, p 0.06-0.45 | out/q04_paired*.json, out/q15_omnibus.json |
| Non-template injection response | mirab: INV 0.92-1.53, REF 0.96-1.55 by sector; altpos 0.94-1.03 | analyst pipelines, 10 fake ephemerides/sector, paired | SEM 0.01-0.05 | out/q11_summary.json |
| Global depth, delivered curves | INV 15.47 +- 1.05, REF 16.08 +- 0.97 (no E58 15.69 +- 1.01); wide LD 16.06 / 16.45 | emcee, density prior | posterior | out/q08_global.json |
| Global depth, raw apertures | 16.7 +- 1.1 (box3), 17.5 +- 1.4 (r<=2.5) all windows; INV-AP paired -0.6..-2.0 +- 1.2-1.7 | as above | | out/q14_contrast_all.json, out/q04_paired.json |
| Ephemeris (combined, my errors) | T0 3840.91510 +- 0.00116, P 5.6915097 +- 0.0000080 (indep.); 3840.91500 +- 0.00156, 5.6915111 +- 0.0000108 (corr.) | WLS of per-event times | 2027 1-sigma 2.6-3.5 min | out/q05_timing.json |
| Timing calibration (per transit) | INV 8.9-14.1 min (200-600 s sectors), 25 min S12; REF 6.1-10.5, 17.7 min S12 | LC-level injections at 57-102 fake epochs/sector | std (robust 4.8-7.3 min) | out/q05_timing.json |
| Odd - even | INV +2.9 +- 2.1 (perm p 0.10); REF -0.6 +- 1.8, +2.7 +- 1.9 without E58; box3 -1.1 +- 2.2; r<=2.5 -0.1 +- 2.8 ppt | fixed-shape per-transit, fake errors | as quoted | out/q10_oddeven.json |
| Secondary, phase 0.5 | box3 -1.7 +- 1.1, r<=2.5 -2.3 +- 1.4 ppt | raw apertures | 3-sigma UL ~3.2-4.1 ppt | out/q12_secondary_ap.json |
| Phase 0.490/0.4925/0.495 | box3 -0.04/-0.36/-1.23 +- 1.08; r<=2.5 -1.1/-1.5/-2.8 +- 1.4 ppt | raw-aperture scan | | out/q12b_phasescan_ap.json |
| LEE p (post-hoc split) | INV 0.37/0.63/0.64, REF 0.31/0.59/0.62, REF no E58 0.52/0.79/0.81 (5/31/130-split families) | 200,000 null draws | | out/q07_lee.json |
| Follow-up | V 15.6, I_C 12.9, SkyMapper i 13.39, z 12.77; 5-sigma single transit needs <= 14 ppt per 2-min bin (white); 10-sigma <= 7 ppt | depth 15.2 ppt, T14 1.70 h | ctx CSV utc_mid truncated to the minute (0-1 min early) | out/q09_followup.json, out/q09_predictions.csv |

## Caveats
* The raw-aperture method has no Mira model; its per-transit noise is 1-2x the PRF fits' (box3 similar, r <= 2.5 px
  ~2x), so it cannot measure sector depths better than ~2.5-4 ppt; it is used to test, not to replace, the PRF depths.
  Its only PSF dependence is the EE (<= 2 %); a PSF with more light outside the aperture than the SPOC PRF would make
  the aperture depths deeper still.
* Background for the apertures uses 5-12 faint pixels per sector (variant 'wide' uses more); without background
  subtraction the apertures are too noisy to be useful.
* The 'mirab' injected image is the Mira's empirical image (a much redder star, and noisy in S66/S104: L2 difference
  from the PRF 102 % and 41 %); it is a stress test of the pipelines' response to a PSF-shape error, not an estimate
  of the M3V's true PSF.
* The post-hoc split was chosen from the PRF data; the split-independent omnibus test (q15) is the fairer statistic.
* No background jobs are left running.
