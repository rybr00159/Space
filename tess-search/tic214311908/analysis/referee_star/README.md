# Referee check of the 'star' investigation of TIC 214311908

Folder: `tic214311908/analysis/referee_star/`

*AI-run adversarial referee check, 2026-10-17. Nothing was submitted, uploaded or posted anywhere.*

The investigation under review is in `../star/`. Every number below was recomputed with the referee's own code (`scripts/`). It uses
the referee's own archive queries (Gaia DR3/DR2, VizieR eRASS1, Mamajek table). The investigator's downloaded images (DECaPS2 cutout)
and its stored DSS measurement table were used **as data only**.

## Verdict in one paragraph

The parameters that matter for the transit-shape test are confirmed:
* R = 0.381 ± 0.012 R_sun, M = 0.370 ± 0.012 M_sun, rho = 9.4 (8.6–10.3) g/cm3;
* twin case: rho = 15.5 (14.1–17.0) g/cm3;
* d = 64.31 ± 0.15 pc.

The over-luminosity is confirmed: ΔMKs = −0.75 mag, ΔMG = −0.68 mag; 5.2–5.9% of comparison stars are as bright. Among stars with RUWE < 1.1, like
the target, **only 1.4%** are as bright, so the target is a rarer outlier than the investigator states. The twin or metal-rich
alternatives remain open. The planet numbers are confirmed: Rp = 5.12 R_earth at 1.64% and b = 0.6; a = 0.0448 au; Teq = 471 K; K = 6.95 m/s
per 10 M_earth; TSM = 129 at 23 M_earth.

Three exclusion claims are overstated:
1. **DECaPS2 contrast at ≤ 1.25".** The 'injection test' in the investigator's code is vacuous, and a single-star fit absorbs 57–66% of a companion's
   flux at 1". The corrected limit is ΔT ≈ 4.0 (z) to 4.2 (i) at 1", not 5.2. The M3V core is also saturated in DECaPS2.
2. **DSS look-back limits.** A sign/logic error in the limit rule assigns T_lim = 24 to grid points that have positive residuals, including
   residuals above 3σ. The corrected I-plate limits are T > 18.27 (within 3") and > 18.28 (within 2"), not 18.92 / 19.30. The R-plate limit is
   18.67 for red sources. The conclusion still holds, because a blend needs T ≤ 17.47.
3. **Gaia IPD "excluded".** Likelihood ratios of 0.04–0.09, with 84% upper bounds of up to 0.18, are called "excluded". They are
   factor-10–25 disfavourings, not exclusions.

Two numbers were mis-specified without changing any conclusion:
* **eRASS1 positional σ.** Per-axis σ is POS_ERR/√2. The M3V is then favoured over the pair by 8.6:1, not ~3:1.
* **Wide-companion velocity test.** The 0.57 km/s vs 1.43 km/s comparison is not meaningful, because A's PM is perturbed by B (RUWE 13.7). The
  binding rests on the El-Badry+21 R = 2e-5.

## Scripts and outputs

| Script | What it does | Output |
|---|---|---|
| `scripts/r01_gaia.py` | Gaia DR3 cone 45" around the M3V (all columns); DR2 counterparts via dr2_neighbourhood; DR3 astrophysical_parameters | `data/gaia_cone45.csv`, `data/gaia_dr2_xm.csv`, `data/gaia_ap_target.csv` |
| `scripts/r02_samples.py a` | Comparison sample: plx > 12 mas, plx/err > 20, BP−RP 2.55–2.95, 2MASS best neighbour (9,239 stars) | `data/mdw_cmd.csv` |
| `scripts/r03_dss_limits_recompute.py` | Recomputes the DSS look-back T limits from the investigator's stored per-grid-point measurements. It uses a correct 3σ upper-limit rule and lists > 3σ positive residuals | `out/r03_dss_limits.json` |
| `scripts/r04_decaps_inject.py` | Image-level injection–recovery in the DECaPS2 i and z cutout. Uses an empirical PSF and a single-PSF refit of the M3V after injection. Gives the recovered (non-absorbed) flux fraction vs separation, corrected limits, and a saturation check | `out/r04_decaps_inject.json` |
| `scripts/r05_field_counts.py` | Gaia DR3 counts in a 40–360" annulus. Gives the density of stars with T ≤ 16.72 / 17.47 etc. and the expected number within 1/2/3" | `data/gaia_cone6arcmin.csv`, `out/r05_field_counts.json` |
| `scripts/r06_star_params.py` | Monte Carlo R, M, rho and twin case from the Mann+15/19 relations, with the zero_point package run independently. Includes a [Fe/H] = +0.2 variant, a DENIS-Ks sensitivity case, and L from the Mamajek BC_Ks | `out/r06_star_params.json` |
| `scripts/r07_xray.py` | eRASS1 Main (J/A+A/682/A34) cone query. Separations to candidate counterparts at the eRASS1 epoch, association likelihoods for both σ conventions, L_X, L_X/L_bol, Ro, P_rot | `data/erass1_cone.csv`, `out/r07_xray.json` |
| `scripts/r08_kin_wide.py` | UVW (astropy, Monte Carlo), Bensby+14 TD/D; wide-system separations, ΔPM, Δv_tan with and without A's RUWE-inflated errors, escape speeds | `out/r08_kin_wide.json` |
| `scripts/r09_cmd.py` | Over-luminosity relative to a KDE-mode ridge (MKs, MG vs BP−RP); colour–colour consistency (G−Ks, G−RP, J−Ks at fixed BP−RP); over-luminous fraction split by RUWE | `out/r09_cmd.json`, `figs/r09_cmd.png` |
| `scripts/r10_ruwe_check.py` | Independent simplified RUWE simulation at spot points (M2 0.08–0.33 M_sun, a 0.03–6 au) for σ_CCD = 0.22 and 0.35 mas | `out/r10_ruwe_check.json` |
| `scripts/r11_planet.py` | batman: k, Rp and T14 for depths 1.4–1.9% and b = 0–0.8 (single and twin; two LD laws); a, S, Teq, K, TSM | `out/r11_planet.json` |
| `scripts/r12_pm_local.py` | Local (0.5° cone) DR3−DR2 proper-motion calibration: plain DR3 and DR2 cone queries matched by position. Gives the calibrated z and percentile of the target | `data/pm_local_dr3.csv`, `data/pm_local_dr2.csv`, `out/r12_pm_local.json` |

Other files:
* `data/EEM_dwarf.txt`: Mamajek EEM table, downloaded 2026-10-17.
* `out/r02a.log`, `out/r12.log`: query logs.

## Key recomputed numbers (method, uncertainty, file)

### Stellar parameters (`out/r06_star_params.json`)

* **Inputs:** parallax ZP −0.0521 mas (zero_point package) gives d = 64.31 ± 0.15 pc; MKs = 6.441 ± 0.026.
* **Single star:** R = 0.3813 ± 0.0116, M = 0.3702 ± 0.0118, rho = 9.41 (8.59–10.35) g/cm3. Agrees with the investigator.
  With [Fe/H] = +0.2 ± 0.15, R = 0.385 and rho = 9.16.
* **Equal twin:** R = 0.289, M = 0.264, rho = 15.45 (14.09–16.99) g/cm3.
* **Luminosity cross-check:** the Mamajek BC_Ks = 2.716 gives L = 0.0171 L_sun, 4% above the SED value of 0.0165. This corresponds to Teff ≈ 3380 K with R_Mann.
* **2MASS photometry check:**
  * The 2MASS PSC flags the M3V cc_flg = 'ddd' (possible diffraction-spike contamination).
  * DENIS J and K are 0.41 / 0.34 mag fainter than 2MASS, and DENIS I is 0.25 mag fainter than SkyMapper/DECaPS i. The investigator did not discuss any of this.
  * 2MASS Ks agrees with unWISE W1 (Ks − W1 = 0.15 against 0.15–0.17 for M2.5–M3.5V) and with the optical SED. The 2MASS value is therefore adopted; DENIS is
    discrepant, probably because of the Mira.
  * With DENIS Ks the density would be 11.8 g/cm3. This sensitivity case is in the json.

### CMD (`out/r09_cmd.json`)

* **Offsets:** ΔMKs = −0.753 and ΔMG = −0.677 mag. 5.2% and 5.9% of 1,059 stars with |ΔBP−RP| < 0.03 are as bright.
* **Split by RUWE:** among RUWE < 1.1 stars, 1.4% are as bright (n = 840); among RUWE 1.1–1.4, 7.3% (n = 946).
* **Mamajek table:** offsets −0.56 (G) and −0.61 (Ks). Twin case: MKs 7.19 and MG 10.99 per star, against M3.5V 7.10 / 10.87.
* **Colour consistency at fixed BP−RP:**
  * G−RP residual −0.001 (z = −0.25): no sign of RP contamination by the Mira.
  * G−Ks +0.088 (z = +1.9; 96th percentile).
  * J−Ks +0.045 (z = +1.3).
  * Over-luminous comparison stars show the same G−Ks excess (+0.056 for ΔMKs −0.55 to −0.9). The colours look like those of other over-luminous stars,
    whether binary or metal-rich.

### Wide system and kinematics (`out/r08_kin_wide.json`)

* **Kinematics:** UVW = (22.0, −5.6, 40.4) ± (2.1, 0.8, 0.35) km/s; v_LSR = 58.4 km/s; P_thick = 0.30. Confirmed.
* **Pair A:**
  * Separation 11.728" = 754 au; A–B separation 0.834" = 53.7 au. Δv_tan = 0.57 km/s (formal).
  * With A's PM errors scaled by its RUWE (13.7), Δv_tan = 1.0–3.2 km/s (68%), against v_esc = 1.43–1.56 km/s.
  * A's photocentre also moves at ~1 km/s about the A–B barycentre. The velocity test is therefore not diagnostic.
  * The El-Badry+21 R = 2e-5, together with the common parallax and PM, is what supports a bound pair.

### X-ray (`out/r07_xray.json`)

* **Source:** 1eRASS J172352.9-511928; DET_LIKE 6.24; 7.1 ± 3.3 counts; F = 5.4e-14; log L_X = 28.43; L_X/L_bol = 4.2e-4.
* **Rotation implied:** Ro = 0.16 (Wright+11, β = −2.7) or 0.19 (Wright+18), so P_rot ≈ 10.5–12.7 d.
* **Separations** at epoch 2020.23: 4.17" (M3V), 8.9" (A), 8.8" (B), 5.5" (background star 3.1" NE), 15.5" (Mira).
* **Association:** with per-axis σ = POS_ERR/√2 = 3.81", the likelihoods are 0.55 (M3V), 0.064 (A), 0.35 (background star), 2e-4 (Mira).
  The investigator's per-axis σ = 5.38" understates the preference for the M3V over the pair.

### Background density (`out/r05_field_counts.json`)

* T ≤ 17.47: 0.00367 /arcsec², so 0.0115 / 0.046 / 0.104 within 1 / 2 / 3". Confirmed.
* T ≤ 16.72 (50% eclipse): 0.0021 /arcsec², so 0.0066 / 0.026 / 0.059 within 1 / 2 / 3".

### Planet (`out/r11_planet.json`)

* **Single star:** a = 0.0448 au, a/R = 25.26, S = 8.22 S_earth, Teq = 471 K.
  * Rp(1.64%, b = 0.6) = 5.12 R_earth (5.11 with LD (0.15, 0.45)); range 4.6–5.8 for 1.4–1.9% and b = 0–0.8.
  * T14 = 1.93 / 1.63 / 1.37 h for b = 0 / 0.6 / 0.8.
  * K = 6.95 m/s per 10 M_earth; TSM 129 (23 M_earth), 297 (10 M_earth).
* **Twin:** a/R = 29.8; T14 ≤ 1.71–1.73 h at b = 0; Rp = 5.3–5.5 at 1.64%.

### DECaPS2 (`out/r04_decaps_inject.json`)

* **Recovered fraction** of an injected companion after a single-PSF refit, i / z:
  * 0.5": −0.02 / −0.07
  * 0.75": 0.10 / 0.13
  * 1.0": 0.34 / 0.43
  * 1.25": 0.57 / 0.66
  * 1.5": 0.73 / 0.80
  * 2.0": 0.86 / 0.89
* **Corrected 5σ limits:**
  * i: 4.24 (1"), 5.07 (1.25"), 5.80 (1.5"), 5.58 (2").
  * z: 3.96 (1"), 4.53 (1.25"), 4.96 (1.5"), 5.24 (2").
  * No limit below 0.75".
* **Saturation:** the M3V peak is 24–30% below the PSF scaled to the wings. The PSF fit recovers 53–59% of the catalogue flux. DECaPS2 flags most M3V
  detections: ndet_ok 13 of 37, with 2 good i and 2 good z magnitudes.

### DSS (`out/r03_dss_limits.json`)

* **Corrected minimum T limits** (BP−RP 1.0 / 2.5), within 3":
  * SERC-I 1979: 18.27 / 18.51.
  * AAO-SES R 1992: 19.31 / 18.67.
  * Within 2": I 18.28 / 18.52; R 19.88 / 19.25.
* **Positive-residual points:**
  * I plate: 8 grid points at 3.0–3.8σ.
  * R plate: 3 grid points at 4.5–9.1σ, at 2.2–2.9" from the star's position at the plate epoch (PSF-wing residuals). The investigator gave these T_lim = 24.

### Astrometric spot check (`out/r10_ruwe_check.json`)

* **σ_CCD = 0.22 mas** (the investigator-equivalent):
  * M2 = 0.08–0.2 M_sun is excluded (L ≤ 0.03) at 0.1–4 au; at 6 au L = 0.02–0.06.
  * M2 = 0.3 is excluded only at ~0.5–3 au (L = 0.23 at 0.3 au, 0.13 at 4 au).
* **σ_CCD = 0.35 mas:** the inner edge moves to ~0.3 au for M2 ≤ 0.2.
* **Near-equal companions** (q ≈ 0.9) are unconstrained.
* The investigator's "0.08–0.3 M_sun at ~0.1–4 au" is approximately right for 0.08–0.2 M_sun. For 0.3 M_sun it should read ~0.5–3 au.

### DR3−DR2 PM (`out/r12_pm_local.json`)

* **Comparison sample:** 4,786 stars within 0.5° matched between DR3 and DR2 (RUWE < 1.4).
* **Calibration on all comparison stars:** z offsets (−0.12, −0.16) and MAD scatter (0.97, 1.04).
* **Calibration on the 64 red stars (BP−RP > 2):** offsets (+0.48, −0.42) and scatter (1.34, 1.37). The offsets lie in the same direction as the
  target's.
* **Target:** ΔPM = (+0.396, −0.236) mas/yr, raw z = (3.36, −1.79).
  * With the all-star calibration, χ² = 15.2.
  * With the red-star calibration, χ² ≈ 5.6 (p ≈ 0.06).
  * Empirically, 4.2% of local comparison stars have a larger calibrated difference. The investigator reported 6.6% from an all-sky M-dwarf sample.
* **Interpretation:** this is a marginal (~2σ) anomaly. It mildly favours, and does not disfavour, a companion with an orbit of a few au to ~10 au. In the
  investigator's own grid, L_PM reaches 2.4 at 0.24–0.56 au and 5–10 au. The investigator reported the PM test only as an exclusion.
