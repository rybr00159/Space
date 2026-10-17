# Referee check of the `ctx` investigation (TIC 214311908)

AI-run adversarial referee report, 2026-10-17. Nothing was posted, uploaded or submitted. All code here was written
independently. The `ctx` scripts were read only to locate errors; no `ctx` code was imported. Inputs that were reused
are listed: ctx's downloaded light curves, ctx's `tglc_leak.json` observed vector, ctx's `gaia_fov_transits.csv`
crossing list, and loc's `wcs_cal.json` WCS offset.

## Verdict in one paragraph

Most of ctx's archival, catalogue, future-data and follow-up numbers reproduce. The one new result that counted
against the M3V does not survive. ctx's claim that the TGLC Sector 39 "leak model" prefers a source about 7″ north
of the M3V (Δχ² ≈ 13, p ≈ 0.003) is a modelling artefact, not evidence:

- **Where the Δχ² comes from.** It sits almost entirely in near-degenerate directions of a near-singular covariance.
  The dominant one is the 0.3–0.4 e⁻/s difference between the TGLC light curves of pair A and pair B, two stars 0.83″
  apart (σ ≈ 0.09 e⁻/s).
- **TGLC does not match the model in exactly those directions.** Real TGLC values differ from a linear TGLC emulation
  applied to the same pixels by χ² ≈ 10 in ctx's metric.
- **ctx's geometry is wrong for two stars.** The TGLC HLSP did not propagate proper motions. A noise-free 1″ (0.05 px)
  geometry error alone produces a Δχ² of 4–48.
- **The pixels themselves disagree with ctx.** A direct PRF fit to the same S39 pixels puts the source 3.6″ *south*
  of the M3V, as expected for a red star's lateral-colour shift. It excludes the G = 17.7 star 5.5″ north (0/400
  injections) and disfavours the G = 17.5 star 3.1″ away (6/400).

Smaller issues:

- The vetting file `predicted_transits_2027.csv` has a `utc_mid` column with no barycentric correction (−3.4 to
  +7.3 min). ctx checked only `bjd_tdb` and called the file correct.
- ctx's prior comparison with eclipsing binaries is numerically wrong: it confuses 0.11 companions per star with a
  rate of about 0.1%.
- The QLP depth depends on the estimator: 1.34% (ctx) against 1.51% (referee), both unbiased in injections.
- The odd/even "systematics" argument is overstated. The difference is simply not significant: 1.7σ with
  empirical errors.

## Scripts (`scripts/`) and outputs (`out/`, `figs/`)
| Script | What it does | Output |
|---|---|---|
| `rlib.py` | Ephemeris; trapezoid template (T14 1.63 h, flat 70 %, exposure-integrated); joint depth fit with an independent quadratic baseline per ±0.3-d window; fake-window helpers | – |
| `leak.py` | Referee re-implementation of TGLC's linear operators (psf: fit of [ePSF/Σ11×11, x, y, 1] on the 5×5 cut, from tglc 0.8.1 `fit_lc`; aperture: 3×3 sum / portion). TESScut WCS + loc's S39 offset; cut grids on the detector grid from TGLC `STAR_X/Y`; `tglc_positions()` gives the positions TGLC actually used | – |
| `r01_tglc_dF.py` | TGLC S12/S39 flux lost per star (e⁻/s) with the local-baseline estimator; fake-triplet noise | `out/r01_tglc_dF.json`, `out/r01_fake_*.npy` |
| `r02_leak_gls.py` | Leak-model GLS: ctx's and the referee's dF; covariances A (independent single windows), B (ctx-style random phases, re-implemented), diagonal, Ledoit–Wolf; without pair B / Mira; leave-one-transit-out; noise-free model-mismatch tests | `out/r02_leak_gls.json` |
| `r02b_tglc_geometry.py` | TGLC-faithful geometry: no proper-motion propagation for the M3V and pair A; refit | `out/r02b_tglc_geometry.json` |
| `r02c_eigen.py` | Eigen-decomposition of covariance B; χ² of each hypothesis broken down by eigenvector | `out/r02c_eigen.json` |
| `r03_pixels_s39.py` | **Same S39 pixels, analysed directly.** Background plane removed; per-pixel amplitude image over the 4 complete S39 transits; PRF fits for each hypothesis plus a free grid; TGLC emulation; injections on the M3V, the colour-shifted M3V, G17.7 and G17.5 (400 each) through both the direct fit and the leak model | `out/r03_pixels_s39.json`, `out/r03_log.txt`, `figs/r03_amp_s39.png` |
| `r03b_variants.py` | Direct-fit robustness: WCS calibration on/off, fit radius 2.5/3.5/5 px, background plane on/off, window ±0.2/0.3/0.45 d, linear baseline | `out/r03b_variants.json` |
| `r03c_tglc_vs_emulation.py` | Real TGLC dF minus emulation from the same pixels, in ctx's covariance metric | `out/r03c_tglc_vs_emulation.json` |
| `r04_qlp_fold.py` | QLP M3V / pair A / Mira depth, phase 0.5 and odd/even with the local-baseline estimator; empirical errors from fake window sets; S103/S104 time shift applied | `out/r04_qlp_fold.json` |
| `r04b_depth_injection.py` | 1.5 % injections into the M3V QLP light curves: bias of the referee estimator and of a re-implementation of ctx's (biweight, masked) | `out/r04b_depth_injection.json` |
| `r05_followup.py` | Transit times 2026-10-17 to 2027-12-31 with σ from the full covariance; check of both columns of `predicted_transits_2027.csv`; observability at CTIO, La Silla, SAAO and Siding Spring (alt > 30°, Sun < −12° over mid ± 1.815 h, and ± 2σ) | `out/r05_followup.json`, `out/r05_observable.csv` |
| `r06_gaia_dr3_epochs.py` | Mira DR3 G epochs (the same field-of-view crossings as the M3V) phased on the ephemeris; recomputes Δt and σ for ctx's crossing list | `out/r06_gaia.json` |
| `r07_rotation.py` | Test ctx did not run: rotation and synchronisation, using a Lomb–Scargle periodogram of TARS (1-px aperture) and TGLC S39 | `out/r07_rotation.json` |
| (curl) | ExoFOP target JSON for 214311906/08/10 and the full TOI table; NASA Exoplanet Archive `pscomppars` query | `ext/exofop/`, `ext/nea_mdwarf_neptunes.csv` |
| (astroquery) | eRASS1 (J/A+A/682/A34) and its counterparts (J/A+A/684/A121); MAST product inventory within 40″ | numbers quoted below |

## Key re-derived numbers

### TGLC S39 leak model (ctx's main new claim)

- **TGLC dF (e⁻/s), referee estimator** (M3V, pair A, pair B, Mira, ref1, ref2): 13.4, 12.4, 12.1, 9.7, −4.5, −3.8.
  ctx's values are 15.0, 13.9, 13.5, 12.2, −4.5, −2.2. They are consistent within the 2–3 e⁻/s errors.
- **Reproduction.** With ctx's dF and a re-implemented ctx-style covariance: χ²(M3V) = 19.7, G17.7 = 10.9,
  free fit = 6.1, Δχ²(M3V − free) = 13.6. The best position is now (−14, +13)″ from the Mira, against ctx's
  (−15, +7)″.
- **Which covariance directions carry it.** Eigenvalue σ's are 0.089, 0.42, 0.87, 1.17, 2.30 and 8.5 e⁻/s.
  - Of the M3V's χ² of 19.7, 8.9 comes from the pair A − pair B direction (σ 0.089) and 6.5 from the Mira-vs-rest
    direction.
  - With a diagonal covariance, Δχ² = 1.2. With Ledoit–Wolf shrinkage it is 6.4; without pair B, 4.95.
- **TGLC against a linear emulation of the same pixels.**
  - Difference vector: (0.0, −0.8, −1.1, −0.6, −2.2, −1.8) e⁻/s, which is χ² = 10.1 in ctx's metric.
  - Of that χ², 5.9 lies in the pair A − pair B direction (TGLC +0.30 against emulation −0.03 e⁻/s) and 3.1 in
    M3V − pair. The leak model therefore does not describe TGLC to the precision its χ² requires.
- **TGLC's own positions** (`STAR_X/Y`) for the M3V and pair A are offset by (+0.61, −0.49)″ from the
  proper-motion-propagated positions. All other stars are within 0.03″, so the HLSP did not propagate proper motions.
  - Refit with this geometry: χ²(M3V) = 30.4 and the free best fit jumps to (−12, +18)″, the edge of the grid.
- **Noise-free mismatch.** A true M3V source with random per-star geometry errors of 0.05 px (≈ 1″) gives
  Δχ²(M3V − free) = 3.8, 5.1, 5.6, 10.3, 36.9 and 47.8 (median 8). At 0.10 px the median is 27.
- **Direct S39 pixel fit** (39 px, 37 dof):
  - χ² for each hypothesis: M3V 30.1; M3V with a 3.5″ colour shift 24.9; G17.5 36.5; G17.7 46.3; pair A 48.8;
    Mira 75.8.
  - The free fit lands 0.9″ E and 3.6″ S of the M3V, with χ² 24.8.
  - Injections:
    - Truth G17.7: P(Δχ²(G17.7 − M3V) ≥ +16.1) = 0/400.
    - Truth G17.5: P(Δχ²(G17.5 − M3V) ≥ +6.3) = 6/400.
    - Truth M3V: p = 0.21 and 0.34 for those two statistics.
  - Robust across the 7 variants: G17.7 +13 to +31, G17.5 +5 to +13, pair A +17 to +33. The free fit always lies
    1.6–5.1″ south.
- **Leak model applied to the emulated (self-consistent) dF:** M3V 5.3, G17.7 4.9, pair A 6.6. It cannot even
  separate the pair from the M3V. ctx's "TGLC excludes the pair (Δχ² +34) and the Mira (+67)" is not a valid test
  either. The exclusion of the pair and the Mira rests on the direct pixel fits (loc, and r03 here).

### QLP and other products

- **QLP M3V DET depth.** The referee gets 15139 ± 967 ppm (S/N 15.7); ctx got 13351 ± 754.
  - Both estimators are unbiased in injections: −0.01 % for the referee and −0.01 % for ctx's masked biweight.
    ctx's estimator would be −5.5 % if the transit were not masked.
  - The 13 % difference is estimator-dependent baseline noise, mainly from S93 and S103. Quote 1.3–1.5 % ± 0.1 %
    (diluted QLP depth).
- **Odd − even (M3V DET):** +2889 ± 1690 ppm, 1.7σ with empirical errors; SAP gives 0.9σ. Not significant.
- **Phase 0.5:** −133 ppm (DET) and +69 ppm (SAP). No secondary eclipse.
- **TARS:** the header keyword `APERTURE = 1` confirms the 1-pixel aperture that ctx called an unverified
  assumption.

### Timing, observability, Gaia

- **`predicted_transits_2027.csv`:**
  - `bjd_tdb` matches to 0.07 min and σ to 0.05 min.
  - **`utc_mid` = BJD_TDB → UTC with no barycentric correction.** It is off by −3.4 to +7.3 min, up to 1.6σ of the
    timing uncertainty. ctx's own `followup_transits.csv` is correct.
- **Observability:** 17 epochs fully observable (15 with ±2σ margins), 5–6 per site, April–August 2027. The first is
  2027-04-11 07:36 UTC from CTIO/La Silla. Late 2026 at best covers 14 % of the window. ctx is confirmed.
- **Gaia:**
  - The Mira's 39 DR3 epochs are all ≥ 404 min from a predicted transit; the nearest has σ 12.7 min. DR4 has no
    in-transit epoch, as ctx says.
  - DR5 crossing 2024-10-29: Δt = +13.9 ± 2.5 min. Δt and σ were recomputed from ctx's crossing list; the crossing
    times themselves were not re-derived.
- **Future TESS:** tess-point 0.9.5 (pointing table to Sector 134) lists only Sectors 12, 39, 66, 93, 103 and 104 for
  the Mira and the M3V. ctx is confirmed.

### Catalogues and literature

- **ExoFOP:** 0 TOIs, CTOIs, imaging, spectroscopy, time-series entries or files for any of the three TICs. The
  nearest TOI is 4398.01, 19.1′ away. ctx is confirmed.
- **NASA Exoplanet Archive:** 18 transiting planets of 4–7 R⊕ with P < 10 d around hosts below 4000 K. Only
  TOI-3884 b has Teff < 3500 K. ctx is confirmed.
- **eRASS1 1eRASS J172352.9−511928:**
  - Position error 5.4″; F_X = 4.93 × 10⁻¹⁴ erg/s/cm².
  - Counterparts: M3V p_ij = 0.73 at 4.2″; pair A 0.27 at 9.0″. ctx is confirmed.
- **Priors:**
  - These quotes are confirmed in the source text: Winters et al. 2019 (log-normal peak 20 au, σ = 1.16 in the text;
    the abstract says a "broad peak at 4–20 AU"), Shan et al. 2015 (0.11), Gillis et al. 2026 (≤ 0.0341), Bryant et
    al. 2023, Gan et al. 2023 and Dressing & Charbonneau 2015.
  - The fetched "Moe & Kratter 2021" abstract is actually arXiv:1706.09894 (Moe & Kratter 2018), so that citation was
    not checked against its own text.
  - **Wrong statement in ctx:** "prior 0.1–1 % … not much above the per-star rate of eclipsing close companions
    (0.11 per M dwarf)".
    - 0.11 is the number of P = 1–90 d companions per M dwarf, not an eclipsing rate.
    - Eclipsing companions with P = 3–10 d are about 0.11 × (0.16–0.27) × 0.04–0.05 ≈ 0.07–0.15 % per M dwarf.
    - On the target star these cannot produce a 1.5 %, flat-bottomed, 1.6-h dip, which needs a body of about 5 R⊕.
    - The relevant false-positive priors are for blended eclipsing binaries, background or bound, and ctx did not
      compute them.

### Rotation (new; ctx did not test)

- **TARS (1-pixel aperture) M3V periodogram:** in all four sectors the strongest peak is at 8.5–9.2 d (S12 8.84,
  S39 8.87, S66 8.51, S93 9.16 d). The LS power at 5.69 d is 0.001–0.18.
- **TGLC S39:** every star, including ref1, peaks at about 7.2 d, so this is a TGLC systematic.
- **Reading:** if the ~8.8-d signal is the M3V's rotation, it is not synchronised with the 5.69-d period. It would
  also fit the near-saturated X-ray activity (Ro ≈ 0.1). This is not decisive, because the 1-px aperture also takes
  in the pair and the Mira's wings.

## Files sound for downstream use (from ctx)
- **Sound:**
  - `lc/` (all downloaded light curves).
  - `out/fold_products.csv` and `out/fold_events.csv` (diluted product depths; estimator spread ±1–2 ppt).
  - `out/followup_transits.csv` (correct barycentric UTC) and `out/followup.json`.
  - `out/future_tess*.json`, `out/gaia_fov_transits.csv`, `out/gaia_in_transit.json`.
  - `out/catalog_status.json`, `out/vizier_extra.json`, `out/known_mdwarf_neptunes.csv`, `out/mags.json`.
  - `out/neowise_fold.json` (non-constraining).
- **Do not use as evidence:**
  - `out/tglc_leak.json`, `tglc_leak_stars.json`, `tglc_leak_robust.json` and `tglc_psf_shape.json`.
  - `fig/tglc_leak_map_s39.png`.
  - The "TGLC prefers 7″ N" and "TGLC excludes the pair and the Mira" statements.
  - The prior sentence comparing 0.1–1 % with "0.11 eclipsing companions".
- **Fix before use:** the `utc_mid` column of the vetting file
  `/home/user/Space/tess-search/vetting_2026-10-16/results/214311906/predicted_transits_2027.csv`.
