# Independent pixel-level localisation of the P = 4.136 d transit near TIC 87055304

Folder: `/home/user/work/tic87/loc_independent/`. All code here was written from scratch. I did not open any file in
`../transit/` until my own numbers were final (see "Comparison with ../transit/" at the end).

## Verdict

**The transit comes from the 9.4" neighbour TIC 87055315 (Gaia DR3 6720986208895208960, T = 13.54), not from the target.**

- **Verdict scale:** *supports false positive*. The signal is not on TIC 87055304, so it is not a planet of the target.
- **Joint 5-sector sky position of the transit source** (TESS PRF):
  - offset from the target: (ΔE, ΔN) = (+6.37 ± 0.42, +6.08 ± 0.38)", errors from a bootstrap over transits;
  - add about 0.5–1" per axis of registration systematics (see the registration section);
  - the source lies **8.80" from the target** and **0.66" from TIC 87055315**, which sits at (+6.64, +6.68)" in 2023.
  - Gaussian PSF: (+5.93 ± 0.45, +5.96 ± 0.41)", 1.0" from the neighbour.
- **Injection-recovery calibration.** Synthetic transits were injected at fake epochs into the real cutouts and recovered with
  the identical pipeline. The real combined statistic s_par is the projection on the target→neighbour axis: 0 means the
  target, 1 means the neighbour.
  - Real value: s_par = 0.91 (PRF) or 0.86 (Gaussian).
  - Target host predicts −0.02 ± 0.10, so **z_T = 8.6–10.1σ** over the four combinations of injection PSF and analysis PSF.
  - Neighbour host predicts 0.99 ± 0.09, so z_N = −0.8 to −1.6.
  - **ln LR(target / neighbour) = −36 to −50**, using Gaussian approximations of the calibration distributions.
  - None of the 20,000 combined target-host simulations reaches the real value.
  - A conservative Fisher combination of per-sector empirical p-values gives p_T ≤ 1.2e-7 (≥ 5.3σ). This value is limited
    only by the 600 simulations per sector.
- **Every one of the 26 individual transits** has s_par > 0.5, i.e. lies on the neighbour side.
  - Under the target hypothesis this has probability 10^-31.6 (nominal), or < 10^-17.8 (conservative: per-sector 95%
    upper limits and half the effective number of fake windows).
- **Robustness.** The result holds across:
  - both PSF models (TESS PRF and a fitted Gaussian), including cross-model injections;
  - all five sectors, including S104 on CCD 3, where the pixel grid is rotated by about 135° relative to CCDs 1–2;
  - 16 analysis variants (fit radius, registration method, baseline model, windows, pointing decorrelation, planar
    background, dropping scattered-light windows). Combined s_par stays within 0.89–1.00 in every variant.
- **Third candidate.** Gaia DR3 6720986208880830336 (G = 17.96, T ≈ 17.35) lies 3.4" from the neighbour. It is the only
  other star any single sector prefers (S13 alone). Jointly it is disfavoured by Δχ² = 97 and by ln LR(neighbour/it) = 24.8
  in a three-host injection calibration. It would also need a ~65–75 % deep eclipse.

## Data and conventions

- **Cutouts.** TESScut 21×21 calibrated FFI cutouts for S13 (cam1/ccd1, 1800 s), S39 (1/2, 600 s), S66 (1/2, 200 s),
  S93 (1/2, 200 s) and S104 (1/3, 200 s). I used QUALITY == 0 cadences only, and excluded the S104 scattered-light episode at
  BTJD 4178.8–4179.2.
  - The FLUX in these cutouts is not background-subtracted (FLUX_BKG = 0). There is a floor of about 180–240 e/s/px
    (sky plus unresolved stars), which is fitted as a constant.
- **Ephemeris.** P = 4.136195 d, T0 = 1658.41325 BTJD, T14 = 0.1875 d.
- **WCS.** Taken from the cutout APERTURE extension.
- **Gaia DR3 positions:**
  - `data/gaia_dr3_90arcsec.csv` within 90";
  - a VizieR copy of Gaia DR3 (I/355/gaiadr3), G < 17.5 within 330", for modelling the whole stamp
    (`gaia_wide/`, `fetch_gaia_wide_vizier.py`). It matches the CSV exactly for the 47 common stars;
  - positions propagated with proper motions to each sector's mid-epoch.
- **T magnitudes** from G and BP−RP (Stassun et al. 2019 relation).
- **Geometry check.** The target–neighbour separation is **9.33" at epoch 2016.0** and 9.40–9.44" in 2019–2026, at position
  angle about 45° (NE).
  - The "8.94"" quoted in CONTEXT.md is the distance from the TIC J2000 position, not the star-to-star separation.
- **PSF models:**
  - (1) the TESS PRF from the `PRF` package (`PRF.TESS_PRF`, Bell & Higgins), with SPOC PRF files (start_s0004) for cam1
    ccd1–3 downloaded to `prf_files/` and interpolated to the stamp's CCD position from `1CRV4P`/`2CRV4P`. I verified that the
    model centroid tracks the input position with slope 1.
  - (2) a pixel-integrated elliptical Gaussian, with widths fitted to each sector's mean image: σ ≈ 0.47–0.64 px.
  - The PRF fits the mean image far better: residual rms within 2 px of the target is 0.3–1.0 % for the PRF versus 1.3–3.4 %
    for the Gaussian.

## Method

### Registration of the direct image (2a) — `s01_registration.py` → `registration.json`, `fig_registration_s*.png`

The mean out-of-transit image (good cadences more than T14 from any transit) is modelled as

  Σ F_i·PSF(x_i+dx, y_i+dy) + b

- x_i, y_i are the Gaia positions projected with the WCS.
- **Free fluxes** (F ≥ 0) for all stars with T < 17 within 60", plus T < 15 stars near the fit region.
- **Catalogue flux × one free scale** for all other T < 17.5 stars in the stamp.
- A free constant b and a **free global offset (dx, dy)**.
- dx, dy are solved by variable projection over a circle of radius R around the target.

Variants:

- R = 2.5, 3, 4 and 5 px;
- neighbour/target flux ratio fixed to the Gaia value;
- "external" registration: target region masked (annulus 2.5–8 px), so only the other stars register the frame;
- both PSF models.

Results:

- PRF offsets are 0.02–0.11 px (≤ 2.5"). The six PRF variants, including the external registration, agree to ≤ 0.064 px
  (≤ 1.3") in every sector.
- **WCS check.** Bright stars (T < 13) were fitted individually for (x, y, F) with all others fixed. Residuals are 0.2–0.9"
  rms per axis, with means of ≤ 0.3".
- **Blend photocentre.** The mean image minus all other stars, fitted with a single PSF, sits at s_par = 0.02–0.17 relative to
  the registered target. The Gaia flux ratio predicts 0.16. So the target position in the pixel frame is pinned to about ≤ 1.3".
- With free fluxes, the Gaussian model drives the neighbour flux to zero. For the Gaussian analyses I therefore use the
  registration with the flux ratio fixed (`local_r4_ratiofix`).

### Difference images (1, 2b) — `pipeline.py`, `s02_diffimg.py` → `diff_<tag>.json/.npz`

**Windows (per transit):**

- In-transit cadences: |t − tc| < 0.35·T14, the central 70 %.
- Out-of-transit cadences on **both** sides: T14 ≤ |t − tc| ≤ 2.5·T14, a gap of 0.5·T14 after the contact points.
- Required coverage: ≥ 50 % of the nominal in-transit cadences and ≥ 25 % of each side.
- Result: 26 transits (6, 5, 5, 5, 5). One-sided and largely flagged events are dropped. The S39 and S66 stray-light flag
  2048 removes several events.

**Difference image:**

- Per pixel, a linear-in-time baseline is fitted to both OOT sides and evaluated at the in-transit times.
- D = baseline − in-transit mean, so a dimming star gives a positive PSF.
- Per-pixel σ comes from the OOT residual rms.
- Transits are stacked with one scalar weight per transit, which preserves the PSF shape.

**Fits** inside r ≤ 3 px of the target:

- (i) **One PSF at a free position** plus a constant. This gives the offset from the target and the neighbour in arcsec, via
  the local WCS Jacobian, and s_par / perpendicular offset.
- (ii) **Amplitudes at fixed positions** plus a constant:
  - target + neighbour, giving f_N = a_N/(a_T + a_N);
  - target + neighbour + ...830336;
  - all T < 17.4 stars within 30" (the stars that can physically produce the depth), linear fit and NNLS;
  - all G < 20 stars within 30" (NNLS);
  - single-host χ² for every Gaia star within 60" that could produce the depth.
- The multi-star fits within 30" are strongly degenerate and are reported only for completeness, in the JSON. The two- and
  three-source fits and the single-host χ² scan are the meaningful ones.

**Uncertainties:**

- Bootstrap over transits: 2000 resamples per sector.
- Per-transit fits: their scatter, and a sign test.
- Bootstrap with only 5–6 transits underestimates errors, so all significances below use the **injection calibration**.

### Injection–recovery (3) — `s03_inject.py`, `s03b_inject_hosts.py` → `inj_*.json/.npz`, `inj3_prf.json`

**What is injected:**

- Trapezoid transits (T14 = 0.1875 d, T12 = 0.045 d, integrated over the exposure) with amplitude equal to each sector's real
  difference-image amplitude, measured with the injection PSF model.
- Each is injected as a scaled PSF dip on the target or on the neighbour.
- Epochs: tc + φ·P, φ = 0.25–0.75 in steps of 0.05. The fake windows never overlap real transit windows.
- 51–59 usable fake windows per sector.

**How it is analysed:**

- The identical pipeline is applied.
- D is linear in the cube, so D(cube + injection) = D(cube) + D(injection) exactly. I verified this on a literally injected
  cube: maximum difference about 1e-12 e/s.
- 600 "resampled" data sets per sector and host: one random fake window per real transit, with no repeats.
- 11 "literal" data sets: the whole ephemeris shifted by φ.
- Null data sets with no injection.

**PSF combinations:** PRF→PRF, Gauss→Gauss, PRF→Gauss and Gauss→PRF. The cross-model runs test sensitivity to PRF-model
mismatch.

**Three-host run:** target, neighbour and ...830336, PRF only. It uses 2-D Gaussian likelihoods of the fitted (E, N) offset.

### Combination (4) — `s04_combine.py`, `s05_joint.py`, `s11_final.py`

- Per-sector likelihoods use Gaussian fits and KDEs of the calibration distributions for s_par, Δχ²(T−N) and (E, N).
- The combined s_par is an inverse-variance-weighted mean; its null distributions are built by drawing one simulation per
  sector.
- Fisher combination of empirical p-values.
- Joint sky fit: one common (E, N) for all sectors, with per-sector amplitude and constant, χ² rescaled per sector, and a
  bootstrap over transits. It also produces a χ² map.

## Results

### Per sector, TESS PRF

Offsets are measured from the target, with bootstrap errors. s_par uses the bootstrap 68 % interval. z_T and z_N are the
injection-calibrated z-scores of s_par against each hypothesis.

| Sector | n_tr | ΔE (") | ΔN (") | dist. to target / neighbour (") | s_par | f_N (2-src) | Δχ²(T−N) | s_par if target host (inj.) | s_par if neighbour host (inj.) | z_T | z_N |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 13 | 6 | +11.7 ± 3.2 | +6.6 ± 0.4 | 13.4 / 5.0 | 1.37 [1.22, 1.65] | 1.62 ± 0.23 | 97 | 0.27 ± 0.49 | 1.20 ± 0.34 | 2.2 | 0.5 |
| 39 | 5 | +7.4 ± 0.4 | +5.8 ± 0.6 | 9.4 / 1.1 | 0.99 [0.95, 1.04] | 1.02 ± 0.03 | 85 | −0.04 ± 0.16 | 1.01 ± 0.19 | 6.3 | −0.1 |
| 66 | 5 | +6.9 ± 1.2 | +4.4 ± 0.7 | 8.1 / 2.3 | 0.84 [0.69, 0.99] | 0.88 ± 0.13 | 60 | −0.07 ± 0.25 | 0.93 ± 0.29 | 3.7 | −0.3 |
| 93 | 5 | +6.7 ± 1.0 | +7.1 ± 0.7 | 9.7 / 0.4 | 1.03 [0.92, 1.15] | 0.99 ± 0.08 | 132 | −0.11 ± 0.25 | 0.91 ± 0.41 | 4.6 | 0.3 |
| 104 | 5 | +4.4 ± 0.6 | +6.6 ± 0.9 | 8.0 / 2.3 | 0.83 [0.75, 0.91] | 0.78 ± 0.09 | 63 | −0.04 ± 0.15 | 0.98 ± 0.12 | 5.8 | −1.3 |

- The neighbour is at (+6.6, +6.7)".
- The Gaussian PSF gives s_par = 1.39, 0.90, 0.83, 0.98, 0.78 and Δχ²(T−N) = 100, 66, 63, 115, 59.
- Every sector puts the source NE of the target, at the neighbour. This includes S104, whose pixel axes are rotated by about
  135°, so a pixel-frame systematic (PRF asymmetry, CCD orientation) cannot produce a sky-fixed offset like this.

**Fraction of the signal on the target** (2-source fits, f_T = 1 − f_N; a target host would give 1):

| PSF | f_T, bootstrap errors | f_T, injection-calibrated errors |
|---|---|---|
| TESS PRF | 0.00 ± 0.03 | 0.12 ± 0.10 |
| Gaussian | 0.12 ± 0.03 | 0.19 ± 0.09 |

The 3-source fit (target, neighbour, ...830336) puts a_T = −5.5 ± 2.2, 0.4 ± 1.7, 2.7 ± 1.8, 0.1 ± 1.2 and 1.7 ± 1.5 e/s on
the target out of about 13 e/s in total.

**Model-light check** (`s07_centroid.py`): the flux-weighted centroid of D, after subtracting the annulus median, minus the
centroid of the direct image in an r = 1.5 px aperture. No PSF fit is involved. Projected on the target→neighbour axis:

| Sector | Observed | Predicted, target host | Predicted, neighbour host |
|---|---|---|---|
| 13 | 0.92 ± 0.16 | −0.17 to +0.07 | 0.55–0.99 |
| 39 | 0.66 ± 0.04 | −0.17 to +0.07 | 0.55–0.99 |
| 66 | 0.49 ± 0.06 | −0.17 to +0.07 | 0.55–0.99 |
| 93 | 0.58 ± 0.09 | −0.17 to +0.07 | 0.55–0.99 |
| 104 | 0.58 ± 0.07 | −0.17 to +0.07 | 0.55–0.99 |

The predicted ranges span both PSF models.

### Joint sky fit (`s05_joint.py`, `fig_joint_chi2map_*.png`)

| PSF | Source (E, N) (") | Offset from target | Offset from neighbour | Δχ² at target | Δχ² at ...830336 | Δχ² at neighbour |
|---|---|---|---|---|---|---|
| PRF | (+6.37 ± 0.42, +6.08 ± 0.38) | 8.80" | 0.66" (1.6σ, bootstrap) | 559 | 97 | 3.1 |
| Gauss | (+5.93 ± 0.45, +5.96 ± 0.41) | 8.41" | 1.01" (2.1σ) | 451 | 100 | 6.1 |

- s_par = 0.935 ± 0.048 (PRF); the perpendicular offset is −0.24 ± 0.35".
- The 0.7–1.0" residual offset from the neighbour is within the registration systematics, which are about 1" (see
  registration).

### Injection-calibrated significance (`combine_*.json`, `fig_injcal_*.png`)

| Inject → analyse | Real s_par (combined) | Target host | Neighbour host | z_T | z_N | ln LR(T/N) combined s_par | Σ sectors ln LR (E, N) | Fisher p_T (floor-limited) |
|---|---|---|---|---|---|---|---|---|
| PRF → PRF | 0.911 | −0.02 ± 0.10 | 0.99 ± 0.09 | 9.1 | −0.9 | −41.0 | −58.7 | 8e-8 (5.4σ) |
| Gauss → Gauss | 0.865 | −0.02 ± 0.09 | 0.99 ± 0.09 | 9.7 | −1.4 | −45.6 | −58.3 | 1.3e-7 (5.3σ) |
| PRF → Gauss | 0.876 | −0.08 ± 0.10 | 0.95 ± 0.10 | 10.1 | −0.8 | −50.3 | −66.5 | 9e-8 (5.3σ) |
| Gauss → PRF | 0.908 | +0.05 ± 0.10 | 1.04 ± 0.08 | 8.6 | −1.6 | −35.8 | −49.4 | 1.2e-7 (5.3σ) |

- The summed Δχ²(T−N) over sectors is +438 (PRF). The target-host simulations give −435 ± 103 (8.4σ away); the
  neighbour-host simulations give +441 ± 114.
- Three-host calibration (PRF): Σ ln L = −76.9 (target), −19.2 (neighbour), −44.0 (...830336). So ln LR(neighbour/target) =
  57.7 and ln LR(neighbour/...830336) = 24.8.
- Caveats on these numbers:
  - They rely on Gaussian approximations of the calibration distributions, which is optimistic in the far tails.
  - The fake windows overlap, so the simulations are correlated.
  - The Fisher numbers are lower bounds set by the number of simulations.
  - The headline I would defend is "**> 8σ** (injection-calibrated), conservatively ≥ 5σ".

### Adversarial checks

**Arguments for the target, and why each fails:**

- **Registration error.** About 9" (0.45 px) of registration error would be needed. The registration variants agree to
  ≤ 1.3", and registration from other stars only (target masked) gives the same answer: combined s_par 0.91.
- **PRF-model mismatch.**
  - Gaussian PSF: combined s_par 0.86–0.90.
  - Cross-model injections recover the correct host in both directions.
  - Model-free centroid shifts point the same way.
- **Pointing jitter** (`s08_jitter_pertransit.py`, `fig_jitter_pertransit.png`). This was the most dangerous systematic.
  - The leverage is F_target/A ≈ 260–295, so an in-transit pointing shift of only about 1.6 mpx would mimic the full
    target→neighbour offset.
  - I measured pointing per cadence from stars 2.5–9 px from the target. The in-transit shifts are 0.1–0.9 mpx, with one
    S93 transit at 3 mpx that has negligible weight.
  - The predicted spurious s_par is −0.31 to +0.61 per transit, with mixed sign; the zero-weight S93 epoch 531 reaches +0.95.
    The mean over transits is +0.10. Per-transit s_par is 0.61–1.95. After subtracting the prediction, 23 of 26 stay above
    0.5 and all stay above 0.36.
  - Decorrelating every pixel against the measured pointing changes nothing: s_par = 1.04, 0.97, 0.78, 1.03, 0.89 (PRF).
  - No momentum dumps fall inside the windows, except one in S13 epoch 5.
- **Scattered light.**
  - Eight windows have uniform background changes larger than 2 e/s/px (up to 351 e/s/px in S93 epoch 531). Their
    stacking weights are small (≤ 0.07).
  - Dropping them, and/or adding a planar background, gives combined s_par 0.92–0.95.
- **Crowding and other stars.**
  - Single-host χ² over all 28 Gaia stars within 60" that could produce the depth: the neighbour is best in S39, S66, S93
    and S104. ...830336 is best only in S13 (by Δχ² = 14). Jointly the neighbour wins by Δχ² = 94 over ...830336 and by 556
    over the target.
  - No other Gaia star can host the signal: Δχ² ≥ 388 relative to the neighbour.
- **Wrong Gaia positions or epochs.**
  - Proper-motion propagation changes the target–neighbour vector by < 0.1".
  - The target's RUWE of 18.8 implies mas-level astrometric noise, which is irrelevant at 9".
  - The VizieR copy matches the shared CSV exactly.
- **CCD orientation.** S104 (CCD 3, rotated about 135°) gives the same sky position.

**Arguments for the neighbour, and the caveats:**

- S13 alone, with only six 30-min transits and one transit carrying 70 % of the weight, scatters to s_par 1.37 and slightly
  prefers ...830336.
- The combined s_par is about 1σ short of 1, and the joint fit sits 0.66–1.0" from the neighbour.
- S104 retains a target amplitude of f_T = 0.22 ± 0.09. The earlier investigator independently finds 0.25 ± 0.08.

None of these changes the host verdict.

### Light-curve amplitude sanity check (`s09_depth_consistency.py`)

- The total deficit is A = 12.2–14.1 e/s and is consistent across sectors.
- If the target were the host, its intrinsic depth would be A/F_target = 3.4–3.8 ppt.
- If the neighbour is the host, its intrinsic depth is about **1.8–2.0 %**, using F_N = F_target × the Gaia flux ratio. The
  free-flux fits give a less stable 2.0–3.4 %.
- A background-subtracted r = 2.5 px aperture would show about 1.9–2.4 ppt. QLP reports 2.9–4.2 ppt. I attribute the
  difference to QLP's background estimate in this crowded field, but I have not verified that (see open issues).

### Robustness table (`s10_variants.py` → `variants_table.json`, `fig_variants.png`)

Combined s_par with bootstrap errors (0 = target, 1 = neighbour). Values in brackets are S13, S39, S66, S93, S104.

| Variant | Combined s_par | Per sector |
|---|---|---|
| prf_main | 0.97 | [1.37, 0.99, 0.84, 1.03, 0.83] |
| gauss_main | 0.90 | [1.39, 0.90, 0.83, 0.98, 0.78] |
| prf_pointing | 0.94 | [1.04, 0.97, 0.78, 1.03, 0.89] |
| gauss_pointing | 0.89 | |
| fit radius 2 px | 0.93 | |
| fit radius 2.5 px | 0.97 | |
| fit radius 4 px | 0.96 | |
| registration: external | 0.91 | |
| registration: flux ratio fixed | 0.99 | |
| registration: R = 2.5 px | 0.99 | |
| constant baseline | 1.00 | |
| OOT 1.5–3 T14 | 0.94 | |
| in-transit central 50 % | 0.94 | |
| planar background | 0.94 | |
| clean windows only | 0.95 | |
| clean windows + plane | 0.92 | |

## Comparison with ../transit/ (read only after the numbers above were final)

- **`localize_twosource.json`** (2-source amplitudes): f_target = 0.00, −0.10, 0.01, −0.12, +0.25 (±0.08–0.09).
  - Mine (1 − f_N, PRF): −0.62, −0.02, 0.12, 0.01, 0.22. These agree in S39–S104, including the ~0.2 target amplitude in
    S104.
  - My S13 is dominated by one 30-min transit with a negative target lobe. With pointing decorrelation it becomes −0.12.
  - The quoted combined f_target = 0.02 ± 0.04 agrees with my bootstrap value (0.00 ± 0.03). However, the
    injection-calibrated error is about 3× larger (± 0.10).
- **`localize_checks.json`, `localization_ratio_maps.json`, `pixel_depth.json`.** Centroid and ratio-map χ² strongly favour
  the neighbour, and dF(5×5) = 13.4 ± 0.5 e/s. This is consistent with my A = 12–14 e/s and with my Δχ²(T−N) = 60–132 per
  sector (CONTEXT quotes 75–145).
- **`prf_localize.json`** (free single-PRF position) gives combined (E, N) = (4.18 ± 0.25, 3.30 ± 0.30)": 5.3" from the
  target and "7σ from the neighbour", i.e. between the two stars. **This disagrees with my (6.37, 6.08)".**
  - The raw pixel positions of the fitted difference-image source agree with mine to ≤ 0.1 px in S39–S104 (0.26 px in S13).
  - The disagreement is therefore entirely in the registration. Their direct-image fit allowed **negative fluxes**: the
    neighbour came out at −75 to −1059 e/s, with residual rms 5–9 %.
  - This shifted their registered target 0.05–0.24 px toward the neighbour, by 3.5–5" in S104 for example. That shrinks the
    apparent offset.
  - My registration has non-negative fluxes, residuals of 0.3–1 %, a fitted N/T flux ratio of 0.11–0.21 (Gaia predicts 0.19),
    and agrees with the target-masked "external" registration. I consider the `prf_localize.json` position a registration
    artefact.
- **The earlier `src/centroid.py` claim** ("on target, 9" neighbour excluded at > 6σ") is **wrong**. Its errors were so large
  (S66 1.1" at 0.3σ means σ ≈ 4") that it could not exclude anything.

## Open issues

1. **Unresolved companions.** A star unresolved from TIC 87055315 within about 1–2" (e.g. a bound companion) cannot be
   separated by TESS pixels.
2. **The neighbour's own light curve.** As an independent check, the neighbour's S13 QLP and eleanor light curves (`data/mast_dl`)
   should show a dip of about 2 % at the same ephemeris. I did not do this here.
3. **QLP depth.** QLP depths of 2.9–4.2 ppt exceed what my amplitude predicts for a background-subtracted aperture
   (1.9–2.4 ppt). This needs to be understood before deriving the neighbour-host k from QLP depths.
4. **Minor tensions.** The 1–2σ pull of the combined position toward the target (S104 f_T ≈ 0.2) and the S13 preference for
   ...830336 are modest and unexplained.
5. **PRF version.** The SPOC start_s0004 PRF was used for all sectors. The real PSF may have changed, which motivates the
   Gaussian and cross-model tests; both agree.

## Files

**Code:**

- `common.py`: I/O, WCS, Gaia, windows, difference image, PSF models, fitting.
- `pipeline.py`: sector context, injections, analysis.
- `s01_registration.py`
- `s02_diffimg.py` (with environment options `LOC_IN_HALF`, `LOC_OOT_IN`, `LOC_OOT_OUT`, `LOC_BASELINE`, `LOC_PLANE`,
  `LOC_CLEAN`)
- `s03_inject.py`
- `s03b_inject_hosts.py`
- `s04_combine.py`
- `s05_joint.py`
- `s06_figures.py`
- `s07_centroid.py`
- `s08_jitter_pertransit.py`
- `s09_depth_consistency.py`
- `s10_variants.py`
- `s11_final.py`
- `run_variants*.sh`
- `fetch_gaia_wide_vizier.py`

**Results:**

- `results.json`: master summary.
- `registration.json`
- `diff_<variant>.json/.npz`
- `inj_<inj>_<ana>.json/.npz`
- `inj3_prf.json`
- `combine_*.json`
- `joint_*.json`
- `centroid_shift.json`
- `jitter_pertransit.json`
- `depth_consistency.json`
- `variants_table.json`

**Figures:**

- `fig_diffimg_{prf,gauss}_main.png`: per-sector mean image, difference image, S/N, and residuals after target-only and
  neighbour-only fits, with Gaia stars and the fitted position.
- `fig_per_transit_offsets_*.png`
- `fig_joint_chi2map_*.png`
- `fig_injcal_*.png`: injection calibration.
- `fig_jitter_pertransit.png`
- `fig_variants.png`
- `fig_registration_s*.png`

**Reproduce** (with `OMP_NUM_THREADS=1`, about 30 min in total):

1. `s01`
2. `s02 prf 3 local_r4 prf_main`
3. `s02 gauss 3 local_r4_ratiofix gauss_main`
4. `s03 prf prf prf_main` (and the other three PSF combinations)
5. `s04 <real> <inj>`
6. `s05`, `s06`, `s07`, `s08`, `s03b`, `s09`
7. `run_variants*.sh`, then `s10`, then `s11`
