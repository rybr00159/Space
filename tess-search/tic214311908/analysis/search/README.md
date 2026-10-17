# TIC 214311908: additional planets, aliases, rotation and flares (`search`)

Folder: `tic214311908/analysis/search/`
(AI-run analysis, 2026-10-17. Nothing was posted, uploaded or submitted anywhere, and no one was contacted.)

**Inputs (read-only):**
* From `../phot/`: `lc_M3V_all.csv`, `lc_pair_all.csv`, `lc_Mira_all.csv`, `lc_allstars_s*.npz`, `data/jitter_s*.npz`,
  `out/wcs_s*.json`, `out/timing_table.csv` and `out/lc_stats.json`. The phot modules `photlib.py`, `p05_phot.py` and
  `p03_scene.py` were imported, and nothing was written into `../phot`.
* The TESScut cubes in `/home/user/work/tesscut5/214311906/`.
* From other folders: `../ctx/lc/` (QLP, TARS), `../referee_phot/out/per_transit.csv`, `../referee_phot/data/ext/`
  (TESS response and Gaia XP spectrum) and `../star/star_params.json`.
* Fractional fluxes were rescaled to the PRF-unit M3V flux of 979 e-/s: f = 1 + 0.955 (flux_norm - 1).

## Bottom line

### 1. A new systematic: an RR Lyrae star 200" away leaks into the phot M3V light curve

**The signal.** The strongest periodicity in the deblended M3V light curve (phot model B) is at P = 0.648528 d. It has
a strongly non-sinusoidal shape. Its first-harmonic semi-amplitude is 1.8-3.9 ppt (1.8-3.8 e-/s) in S39-S104; S12
shows less than 1.8 ppt.

**Its source.** A search of the whole cutout (all 441 pixels, 6 sectors; `s05b`) puts the source 198-204" S of the
M3V, within 0.09-0.13 px of **Gaia DR3 5925372615693467904**.
* That star has G 15.67 and lies at the cutout edge.
* Gaia DR3 classifies it as an **RRab star with P = 0.6485275 d** and a 1.0-mag peak-to-peak G amplitude.
* Placing the source on the M3V is worse by Delta-chi2 = 2,800-7,100 per sector.
* A coherent fit over 2021-2026 picks P = 0.648528 d, the Gaia period. The neighbouring alias teeth are worse by
  Delta-chi2 = 106 and 133.

**How it gets in.** The leak enters with inverted sign: the folded M3V curve shows a fast drop and a slow recovery.
* It does not come through the faint-pixel background. The RRab's 3x3 neighbourhood holds 0-1 of the faint pixels,
  and B(t) shows no signal.
* The most likely route is the pointing regressor. The field-star pointing offset dx(t) carries the RRab period
  (p_emp <= 0.09 in 4 of 6 sectors; amplitude 3-4e-4 px). The Mira's steep PSF converts that into a few e-/s around
  the M3V.

**Consequences:**
1. It is **not** rotation of the M3V.
2. It was the top peak of every BLS run on the light curve as delivered, together with its harmonics.
3. It changes **individual transit depths by up to 5.7 ppt (rms 2.8 ppt)**.
   * Removing it lowers the chi2 of the per-transit depths from 21.7 to 16.8 (19 dof).
   * S103 E55 / E56 / E59 move from 20.4 / 3.4 / 24.8 ppt to 14.8 / 8.1 / 19.1 ppt.
   * The model-B odd-even difference falls from 3.1 +- 1.8 to 2.3 +- 1.8 ppt.
   * The global depth does not change (14.96 -> 14.87 ppt, fixed shape).
   * **Transit fits on the phot model-B curve should remove it.** Use `data/m3v_clean.npz` (column `det_clean`), or
     refit 3 harmonics at 0.648528 d per sector.

**The same coupling elsewhere.**
* A control star's light curve shows Delta-chi2 = 378 at 1.595 d, which matches the Gaia DR3 EB
  5925378392451781120 (G 17.1, 168" away; matched by period only).
* The model-A pair and Mira curves peak at 0.7499 d. That is close to the Gaia DR3 EB 5925384852082638720
  (P = 0.7503 d, 167" away).

### 2. The 5.69-d transit is not such a leak

The same whole-cutout search for the transit signal (`s05c`) finds:
* The dimming lies within 1.3-4.8" of the M3V in all six sectors, with a flux loss of 12.3-17.9 e-/s.
* The M3V position is worse than the best position by only Delta-chi2 = 0.2-7.3.
* The best position more than 2 px away is worse by Delta-chi2 = 24-121.

No other source in the 7.3' x 7.3' cutout carries the signal.

### 3. No additional transiting planet

**The empirical threshold.** On the light curve with the RRab removed and the 5.69-d transits masked, the coherent BLS
(0.3-30 d, 947,276 frequencies, all 6 sectors) was calibrated on 60 null versions of the same light curve. These were
30 versions with each sector rotated in time, 16 that were also inverted, and 14 inverted maxima from an earlier run.
* The null maxima have a median of Delta-chi2 = 107 and range 87-138.
* A Gumbel fit gives a **global 1 % false-alarm level of Delta-chi2 = 147** (0.1 %: 171).

**The real data.**
* The highest peak is Delta-chi2 = 94.3 at 26.36 d (3 events; one supplies 59 % of the Delta-chi2), with
  p_emp = 0.80. **Nothing is near the threshold.**
* The per-sector searches (null p-values 0.02-0.94) and TLS per sector (SDE 5.1-6.8) find nothing either.
* The only per-sector excess is S104 (Delta-chi2 = 92, p_emp = 0.02 before allowing for the 7 sector groups). It is
  one event at BTJD 4197.97 that fails vetting (see 7).
* The S103+S104 TLS peak at 16.79 d (SDE 12.6; TLS's white-noise FAP 8e-5) is not significant against 12 time-rolled
  copies of the same data (p = 0.42).
  * It rests on two events, E0 and E2. E0 sits in scattered light, with 55 % good cadences and the background 25 %
    above the sector median. E2 is the S104 event, where the model-A pair takes -221 and the M3V +86 e-/s, a
    degenerate split.
  * The middle event (BTJD 4181.18) is absent: 2.7 +- 4.3 e-/s.

**Completeness** (800 injections through the full chain, recovered at Delta-chi2 >= 147):

| P range (d) | R50 (R_earth) | R90 (R_earth) |
|---|---|---|
| 0.5-1 | 2.8 | 3.1 |
| 1-2 | 3.0 | 3.3 |
| 2-4 | 3.6 | 4.1 |
| 4-8 | 4.2 | 4.8 |
| 8-15 | 4.5 | 5.2 |
| 15-30 | 5.2 | 6.0 |

* Planets with R >= 4.5 R_earth are recovered 88 % of the time (95 % at P < 15 d).
* Nothing below 2.5 R_earth is recovered.
* A planet like the known one (5 R_earth at 5.7 d) would be found at about 100 %.

**What this means.** A second transiting planet would have strongly supported the planet interpretation, but the
search is only sensitive to Neptune-size companions. Most companions in M-dwarf systems are 1-2.5 R_earth, so the
null result carries almost no information either way.

### 4. The 5.69-d period is unique

**BLS without the mask.**
* The highest peak is 5.691505 d (Delta-chi2 = 402).
* Its simple aliases are lower: P/2 184, 2P 250, P/3 116, 3P/2 126, 3P 168.
* The next unrelated peak is Delta-chi2 = 124 (28.0 d).
* Alias comb (teeth every P^2/baseline = 0.0126 d): the next tooth is worse by Delta-(S/N)^2 = 366.

**Phases 0.5, 1/3 and 2/3.** The depths there are -0.92, -0.29 and -0.90 +- 0.76 ppt, against 14.7 +- 0.76 ppt at
phase 0. **P/2 and P/3 are excluded at more than 19 sigma.**

**2P is not excluded by photometry.** An equal-depth EB at 11.38 d would look the same. Odd - even =
+3.3 +- 1.8 ppt (1.8 sigma, empirical errors) on the cleaned model-B curve. The referee's independent photometry gives
-0.6 +- 1.5 e-/s.

**Cycle count.**
* From S103+S104 alone, the number of cycles back to S93 is fixed to +-0.013 cycles (1 sigma).
* S66, S39 and S12 then follow to +-0.002-0.005 cycles.
* The best alternative cycle assignment of the 19-20 transit times is worse by Delta-chi2 >= 3.2e4, with calibrated
  errors, for both the phot and the referee times.
* Period: 5.6915092 +- 0.0000122 d (phot times), 5.6915051 +- 0.0000092 d (referee times).

### 5. Rotation: no period can be tied to the M3V

* **0.6485 d** is the RRab (item 1).
* **The ~4.5-10.7-d peaks** of the deblended M3V curve change from sector to sector. They also appear at ten empty
  "phantom" PSF positions next to the Mira (`s22`), at the M3V's own period:
  * In every sector, 10-100 % of the phantoms are at least as strong as the M3V.
  * In S66 and S103 most phantoms peak at 9.3-10.3 d themselves.
  * These are systematics of the Mira's surroundings.
* **The TARS ~8.0-8.8-d modulation** (`s18`, `s11`) is not on the M3V.
  * If it came from the M3V, it would imply 16 e-/s (1.6 %) in S12 and 7.9 e-/s in S93.
  * The deblended M3V shows 6.0 e-/s (95 % upper limit 9.8) and 3.0 e-/s (upper limit 4.9) at the TARS periods.
  * Pixel injections recover 70-100 % of 8.7-d signals on the M3V.
  * The raw 1-px light curves show no PRF-shaped source at that period.
* **Synchronisation test** (`s17`). This tests whether the star rotates in step with the 5.69-d orbit, as a tidally
  locked binary would.
  * No modulation at 5.6915 d: per-sector 95 % upper limits 0.4-2.3 %, best 0.37 % in S93.
  * None at 11.383 d either: upper limits 0.4-3 %. But pixel injections recover only 20-50 % of a 12-d signal.
  * So there is no hint of a synchronised binary, but the 11.4-d test is weak.
* **Sensitivity** (`s14`, pixel-level injections through the full phot photometry):
  * A 1 % semi-amplitude sinusoid at P = 1-8.7 d is detected in 50-67 % of sector-trials; at 2 %, in 83-100 %.
  * P >= 12 d cannot be measured (recovered amplitude 0.2-0.5 of the injected one).
* **X-ray and age.**
  * L_X/L_bol = 4.2e-4 implies P_rot ~11 d if the star is unsaturated (star agent: Ro = 0.17).
  * Within the ~0.3 dex scatter and the eRASS counting error (7 +- 3 counts), it is also consistent with saturation,
    i.e. P_rot below ~8.5 d at Ro below ~0.13.
  * At a typical amplitude of 0.5-1 % (Newton+16), the photometry cannot exclude either range.
  * Gyrochronology priors: Medina+22 (0.1-0.3 M_sun) give 2.0 +- 1.2 Gyr for P_rot < 10 d and 5.6 +- 2.7 Gyr for
    10-90 d; Newton+16 give < 2 Gyr on average for P < 10 d.
  * The old kinematics (W = +48 km/s, star agent) favour P_rot of tens of days, which TESS cannot measure here.
  * The tension between the X-ray level and the kinematics remains open. Rotation does not resolve it.

### 6. Flares: none securely on the M3V, and the test is not constraining

**The detector** (`s10`). Its threshold was set to at most 1 false event on the inverted light curve. Completeness is
50 % at an equivalent duration of 75 s and 90 % at 237 s, i.e. E_TESS = 1.0e33 and 3.3e33 erg
(L_TESS = 1.4e31 erg/s).

**4 candidates in 110.8 d:**
* two in S12 (30-min cadences);
* one 2-point event in S104, which the model-A pair light curve takes (+267 e-/s);
* one single 200-s point in S93 exactly at the threshold (z = 4.1). Its pixel fit lies 9.4" from the M3V; the M3V is
  worse by Delta-chi2 = 11 and the pair by 31.

None passes the pixel attribution.

**Rate.** No flares above E_TESS = 1.05e33 erg gives a rate below 0.036 per day (95 %).
* Medina+20 active mid-M dwarfs (log R31.5 = -1.30, alpha = 1.98) predict 0.18 such flares in these data.
* Inactive ones predict fewer than 0.0005.
* The data cannot tell active from inactive.

### 7. Single transits: none credible

**Above the inverted-curve threshold (S/N 7.1)** two events exceed it, and both fail the diagnostics:

| Event | S/N | Depth | Duration | Model-A pair | Pixel image |
|---|---|---|---|---|---|
| S104, BTJD 4197.97 | 9.0 | 24 ppt | 4 h | +379 e-/s | not PRF-like at the M3V; best position 223" away |
| S66, BTJD 3117.31 | 7.5 | — | — | +64 e-/s | not PRF-like at the M3V; best position 166" away |

**Above the phantom/control-based threshold (S/N 11.3)**: none.

**Sensitivity.** Single transits are recovered at more than 64 % only for depths of 15-25 ppt or more with T14 > 2.5 h,
at the S/N-7.1 threshold. At the S/N-11.3 threshold the limit rises to depths of 25 ppt or more, i.e.
>= 6.6 R_earth.

## Evidence summary for the planet hypothesis

**For** (from this folder):
* The period is unique: aliases, alias comb and cycle count (item 4).
* The dip is not a leak from anywhere else in the cutout (item 2). This rules out the coupling mechanism exposed by
  the RRab as an explanation of the transit.
* There is no modulation at 5.69 d (upper limit 0.4-2.3 %), so no sign of a tidally synchronised host.

**Against or open:**
* **No second transiting planet.** This is uninformative: completeness is 0 % below 2.5 R_earth (item 3).
* **The deblended light curve has a cutout-wide coupling to distant variables** (RRab, Gaia EBs). The transit is
  unaffected, but the phot per-transit depths, odd/even and S103 anomalies are partly this systematic.
* **2P (an equal-eclipse EB at 11.38 d) cannot be excluded by photometry.** Odd - even is +3.3 +- 1.8 ppt on the
  cleaned model-B curve; the referee's photometry gives -0.6 +- 1.5 e-/s.
* **No rotation period and no flares are measurable.** The activity and age tension (X-ray active, kinematically old)
  is not resolved here.

## Key numbers

| Quantity | Value | Method | Uncertainty / calibration | File |
|---|---|---|---|---|
| Leak period in the M3V light curve | 0.648528 d | coherent 3-harmonic fit, S39-S104 | next alias teeth worse by Delta-chi2 = 106 and 133 | out/periodic_model.json (s07) |
| Leak amplitude (1st harmonic) per sector | S12 0.73+-0.55; S39 1.84+-0.19; S66 2.46+-0.42; S93 3.23+-0.23; S103 3.90+-0.41; S104 2.03+-0.28 ppt | per-sector harmonic fit | 0.5-d block bootstrap; fake-period p99 0.8-1.9 ppt | out/periodic_model.json |
| Leak source | 198-204" S of the M3V, 0.09-0.13 px from Gaia DR3 5925372615693467904 (RRab, P = 0.6485275 d, G 15.67) | whole-cutout single-PRF fit of 6 harmonic amplitude images | M3V position worse by Delta-chi2 = 2,806-7,132 per sector | out/pixfull_periodic_p0p6485.json (s05b), out/gaia_variables.json (s08) |
| Leak route | dx(t) carries P_RR (p <= 0.09 in 4/6 sectors); not B(t) | harmonic fit to the phot regressors | fake periods | out/rrlyr_mechanism.json (s23) |
| Leak effect on transit depths | rms 2.8 ppt, max 5.7 ppt per transit; depth chi2 21.7 -> 16.8 (19 dof); global 14.96 -> 14.87 ppt | fixed trapezoid, local linear baselines | empirical per-sector errors (100 fake epochs) | out/rrlyr_transit_effect.json (s16) |
| 5.69-d dimming position | 1.3-4.8" from the M3V; M3V position worse by Delta-chi2 = 0.2-7.3; next minimum > 2 px away worse by 24-121 | whole-cutout PRF fit of the per-pixel transit depth | 60 fake-epoch error images | out/pixfull_transit.json (s05c) |
| Global BLS false-alarm level | Delta-chi2 = 147 (1 %), 171 (0.1 %) | Gumbel fit to 60 null maxima | null median 107, range 87-138 | out/summary.json (s02, s20) |
| Highest real BLS peak | Delta-chi2 = 94.3, P = 26.36 d, 3 events | coherent BLS on the RRab-cleaned light curve | p_emp = 0.80 | out/bls_real_clean.json (s01), out/summary.json |
| Per-sector BLS | maxima 33-102; p_emp 0.02 (S104, one systematic event) to 0.94 | as above, per sector | 46 null realisations each | out/summary.json |
| TLS per sector | SDE 5.1-6.8 (inverted 4.8-11.6); S103+S104 12.6 at 16.79 d | TLS 2.0.1, R 0.381, M 0.370 | rolled null SDE 7.5-18.6, p = 0.42 | out/tls.json, out/tls_null_S103_104.json (s19, s19c) |
| TLS sanity check (transits kept) | SDE 12.0 (S39), 18.9 (S93), 25.1 (S103+S104) at 5.686-5.693 d | | | out/tls_sanity.json (s19b) |
| Completeness R50 / R90 | 2.8/3.1 (0.5-1 d), 3.0/3.3 (1-2), 3.6/4.1 (2-4), 4.2/4.8 (4-8), 4.5/5.2 (8-15), 5.2/6.0 (15-30 d) R_earth | 800 batman injections through the full chain; logistic fit | binomial (n = 115-146 per P bin) | out/inject_*.csv (s09), out/summary.json, figs/completeness.png |
| Unmasked BLS peak | 5.691505 d, Delta-chi2 = 402; 2P 250, P/2 184 | coherent BLS | next unrelated 124 | out/alias.json (s12) |
| Alias comb | next tooth worse by Delta-(S/N)^2 = 366 | fixed trapezoid, T0 refit | | out/alias.json |
| Depth at phases 0.5 / 1/3 / 2/3 | -0.92 / -0.29 / -0.90 ppt (primary 14.72 ppt) | fixed trapezoid | +-0.76 ppt (300 random phases) | out/alias.json |
| Cycle-count uncertainty | S103+104 -> S93 +-0.013; -> S66 0.005; -> S39 0.002; -> S12 0.004 cycles | progressive linear ephemeris | calibrated per-transit errors | out/alias.json |
| Next-best cycle assignment | Delta-chi2 >= 3.2e4 (referee), 1.8e5 (phot) | timing chi2 over all assignments within +-5 teeth | | out/alias.json |
| Odd - even (cleaned model B) | +3.3 +- 1.8 ppt (1.8 sigma) | fixed trapezoid | 300 random phases | out/alias.json |
| Modulation at 5.6915 d | <= 0.37-2.3 % (95 %, per sector, model B) | sinusoid + quadratic trend, 1-h bins | 200 fake periods | out/sync_limits.json (s17) |
| Rotation sensitivity | 1 % at P = 1-8.7 d detected in 50-67 % of sector-trials; 2 % in 83-100 %; P = 12 d: recovery 0.2-0.5 | pixel injections through phot p05 | 6 sectors x 4 amplitudes x 6 periods | out/rotation_inject.csv, out/rotation_inject_summary.json (s14) |
| ~5-10-d peaks: phantom test | fraction of phantoms >= the M3V at its period: 0.2, 0.3, 0.6, 1.0, 0.4, 0.1 (S12-S104) | same fit at 10 empty positions next to the Mira | | out/phantom_rotation.json (s22) |
| TARS 8.7-d amplitude | 2.49 (S12), 1.37 (S39), 0.47 (S66), 1.44 (S93) ppt; implied on the M3V 16, 7, 2.3, 7.9 e-/s | sinusoid fit; scaled by the TARS transit depth | p_emp 0.00, 0.00, 0.15, 0.18 | out/tars_amp.json (s18) |
| Flare completeness | ED50 = 75 s (E_TESS 1.0e33 erg), ED90 = 237 s | 1,350 Davenport-template injections | | out/flares.json, out/flare_injections.csv (s10) |
| Flare rate on the M3V above 1.05e33 erg | 0 observed; < 0.036 / d (95 %) | 110.8 d of good data | Medina+20 active: 0.18 expected | out/flares.json |
| L_TESS of the M3V | 1.4e31 erg/s | Gaia XP spectrum x TESS response, d = 64.31 pc | +-20 % (zero-point 1.16e31, 3350-K blackbody 1.73e31) | s10 docstring |
| Single-transit threshold | S/N 7.1 (inverted curve), 11.3 (median of phantom/control maxima) | box matched filter, 0.75-6 h | | out/single_transits.json (s15) |

## Scripts (`scripts/`; run with /home/user/work/venv/bin/python, OMP_NUM_THREADS <= 2)

| Script | What it does | Outputs |
|---|---|---|
| `slib.py` | Shared helpers: loading and rescaling, transit mask, binning, wotan detrending (as phot p06), duration-aware frequency grid, numba coherent BLS (box, five durations 0.3-1.25 x the central T14 for rho = 9.4 g/cm3), batman injection model | - |
| `t_timing.py`, `t_ls.py`, `t_fold.py`, `t_ctrl.py` | First looks: timing, Lomb-Scargle and folds where the 0.6485-d signal was found; control-star census | stdout only |
| `s01_bls_real.py` | Coherent and per-sector BLS (real and inverted) of the as-delivered (model B, model A) and, with argument `clean`, the RRab-cleaned light curve | out/bls_real.json, out/bls_real_clean.json, data/bls_real_*.npz |
| `s02_bls_null.py` | Null realisations (rolled, inverted+rolled, pair/Mira, controls/phantoms), same search; incremental saving | out/null_*_clean_*.json |
| `run_nulls.sh`, `run_inject.sh` | Background runners for s02 and s09 | logs/ |
| `s03_controls.py` | Control-star light curves via phot p05 `controls`, and 10 phantom-position light curves (Mira b + M3V + phantom PRF + plane) | data/ctrl_lcs.npz, out/controls.json |
| `s04_pixcube.py`, `s04b_pixfull.py` | Cache the phot difference-image cubes D (11x11 crop; full 21x21) | data/pix_s*.npz (pixfull deleted after use; regenerate with s04b in about 1 min) |
| `s05_pix_periodic.py` | Crop-level localisation of the 0.6485-d signal with injections (superseded by s05b: the source is outside the crop) | out/pix_periodic_p0p6485.json, figs/pix_periodic_p0p6485.png |
| `s05b_pixfull_periodic.py` | Whole-cutout localisation of a periodic signal (harmonic amplitude images, PRF scan) | out/pixfull_periodic_p0p6485.json, figs/pixfull_periodic_p0p6485.png |
| `s05c_pixfull_transit.py` | Whole-cutout localisation of the 5.69-d transit | out/pixfull_transit.json |
| `s06_public_periodograms.py` | Lomb-Scargle of the TARS and QLP products | out/public_periodograms.json |
| `s07_periodic_model.py` | Per-sector harmonic model of the 0.6485-d signal, coherence test, residual periodograms, cleaned light curve | out/periodic_model.json, data/m3v_clean.npz, data/coherence_scan.npy, figs/periodic_fold.png |
| `s08_gaia_variables.py` | Gaia DR3 variability tables within 4' (network) | data/gaia/gaia_vari_4arcmin.csv, out/gaia_variables.json |
| `s09_inject.py` | Injection-recovery of the coherent BLS (full chain re-run for each injection) | out/inject_[1-8].csv |
| `s10_flares.py` | Flare detector, injection calibration, attribution, rate comparison | out/flares.json, out/flare_candidates.csv, out/flare_injections.csv, figs/flares.png |
| `s11_tars_pixels.py` | 1-px light curves from the cube; amplitude maps at the TARS periods | out/tars_pixels.json, figs/tars_pixels.png |
| `s12_alias.py` | Aliases, alias comb, cycle count, alias phases, odd/even | out/alias.json, figs/alias.png |
| `s13_lit.py` | Crossref abstracts (Medina+20/22, Ilin+19/21, Newton+16, Pass+22, Popinchalk+21 and others) | data/lit/crossref_abstracts.json |
| `s14_rotation.py` | Real periodograms of the undetrended deblended curves and pixel-level sinusoid injections | out/rotation_real.json, out/rotation_inject.csv, out/rotation_inject_summary.json |
| `s15_single.py` | Single-transit matched filter, thresholds, vetting, injections | out/single_transits.json, out/single_events.csv, out/single_inject.csv |
| `s16_rrlyr_effect_on_transits.py` | Per-transit depths with and without the RRab leak | out/rrlyr_transit_effect.json |
| `s17_sync_limits.py` | Modulation limits at P_orb, 2 P_orb and the TARS period | out/sync_limits.json |
| `s18_tars_amp.py` | TARS 8.7-d amplitude against the TARS transit depth | out/tars_amp.json |
| `s19_tls.py`, `s19b_tls_sanity.py`, `s19c_tls_null.py` | TLS per sector (real and inverted), sanity check with transits, S103+S104 rolled null | out/tls.json, out/tls_sanity.json, out/tls_null_S103_104.json |
| `s20_summary.py` | Thresholds, p-values, completeness map, figures | out/summary.json, figs/bls_clean.png, figs/null_dist.png, figs/completeness.png |
| `s21_vet_candidate.py` | Per-event vetting of a periodic candidate (run for P = 16.79 d) | out/vet_P16p79.json |
| `s22_phantom_rotation.py` | M3V long-period power against 10 phantom positions | out/phantom_rotation.json |
| `s23_rrlyr_mechanism.py` | How the RRab couples in (faint-pixel set, B, dx, dy) | out/rrlyr_mechanism.json |

Logs of every run are in `logs/`. `logs/s02_inv_roll_1_killed.log` holds the 14 global maxima of an interrupted run,
which are used in s20.

## Caveats

* **The pipeline is phot's, not re-derived.** The search uses the phot model-B photometry, which is deblended against
  the Mira but carries the cutout-wide coupling found here. Model A is noisier and has long-period single-event
  systematics (its top BLS peaks have 76-94 % of the Delta-chi2 in one event); it was not used for the planet search.
* **The nulls are imperfect.**
  * The rolled and inverted M3V nulls keep the same systematics. Inverted copies turn positive excursions into dips;
    one inverted event, Delta-chi2 = 119.5, recurs in half of them.
  * The pair, Mira, control and phantom curves contain leaked real variables and larger systematics, so they were
    reported but not used for the threshold.
  * The threshold relies on 60 maxima and a Gumbel extrapolation to 1 %.
* **The combined TLS was not run.** The 2,574-d baseline makes the TLS period grid too large (more than 1e6 periods).
  TLS was run per sector and for S103+S104, and the combined search uses the numba BLS.
* **Rotation tests are limited by Mira systematics**, which reach 10-50 e-/s at 5-10 d at phantom positions, and are
  blind at P >= 12 d. The TARS argument assumes TARS's processing preserves an 8.7-d signal as it preserves the
  transit.
* **Flare energies are TESS-band energies** (+-20 % from L_TESS). The Medina+20 comparison sample is 0.1-0.3 M_sun,
  slightly less massive than this 0.37 M_sun star.
* **The RRab coupling route is likely but not fully proven.** dx(t) is significant in 4 of 6 sectors, but not in S93
  and S103, where the leak is largest.
