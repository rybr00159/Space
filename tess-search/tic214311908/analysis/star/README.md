# TIC 214311908: the host star, its companions and its neighbourhood

Folder: `tic214311908/analysis/star/`

*AI-run analysis, 2026-10-16/17. Nothing was submitted, uploaded or posted anywhere. All numbers below were
re-derived here from archive data. Where they rely on another investigation's outputs, this is stated.*

## Bottom line

* **Star.** Taken as a single main-sequence star, TIC 214311908 has the following parameters (Mann+15/19 M_Ks
  relations; Fbol from the SED):
  * R = 0.381 +- 0.012 R_sun, M = 0.370 +- 0.011 M_sun, rho = 9.4 +- 0.9 g/cm3 (6.7 +- 0.6 rho_sun);
  * Teff = 3350 +- 55 K, L = 0.0165 +- 0.0003 L_sun;
  * d = 64.31 +- 0.15 pc, after the parallax zero point is applied.

  These agree with the earlier values (0.38 R_sun, 0.37 M_sun, ~3300 K, 64.5 pc).
* **New and important: the star is over-luminous.** It lies 0.72 mag above the empirical M_Ks–(BP−RP) main-sequence
  ridge, and 0.63 mag above it in M_G. Only 5.6% of 3,061 comparison M dwarfs at the same colour are as bright.
  Two explanations remain open:
  * **(a) a metal-rich single star**, with [Fe/H] about +0.2 to +0.4. A calibration on 186 stars with spectroscopic
    metallicities gives +0.19 +- 0.16 for the offset alone. About +0.37 would be needed to explain the whole offset.
  * **(b) an unresolved near-equal-mass companion within about 0.08" (5 au).** Such a pair, two M3.5V stars of
    0.29 R_sun each, reproduces both the colour and the luminosity exactly.

  Youth is not the explanation:
  * BANYAN Σ gives field membership (P = 1.00).
  * The space motion is hot: W_LSR = +48 km/s and |v_LSR| = 58 km/s. This is typical of an old thin-disk star
    (P_thick = 0.30).
* **The wide system is real.** TIC 214311910 A is a bound companion at 754 au projected separation:
  * El-Badry+21 list the pair with a chance-alignment probability R = 2e-5.
  * Their tangential velocities differ by 0.57 km/s; the escape velocity at that separation is 1.4 km/s.
  * A has its own companion B 0.83" (53 au) away, with M_G = 11.7, about 0.20 M_sun.

  So the system is a hierarchical triple, or a quadruple if (b) holds:
  * The outer orbit is ~2e4 yr and the inner A–B orbit ~550 yr. The ratio a_out/a_in is ~14, so the system is
    dynamically stable.
  * No other co-moving star lies within 2 deg.
* **Unresolved companions to the M3V that could produce the signal.** A companion can supply the 15.4 e-/s flux loss
  only if ΔT <= 4.46 (for a 100% eclipse) or ΔT <= 3.71 (for 50%). The archive data exclude most of the space where
  such a companion could sit; `figs/companion_exclusion.png` shows the full map.
  * **Gaia IPD/RUWE diagnostics** (calibrated on ~3,200 high-resolution-imaged red stars) exclude companions with
    ΔI < 2 at 0.08–2.5" (5–160 au). ΔI = 2–3.5 is disfavoured at 0.08–0.8" and excluded at 0.8–2.5".
  * **Gaia astrometry.** RUWE 1.076 and the DR3−DR2 proper-motion difference were simulated with orbit fits. They
    exclude 0.08–0.3 M_sun companions at ~0.1–4 au, and brown dwarfs of ≥ 40 M_Jup at 0.24–2.4 au.
  * **DECaPS2 2018 residual images** exclude unmodelled sources with ΔT < 5.2–8.6 at 1–6".
  * **Not excluded:**
    1. any companion with a < ~0.1 au, that is, a spectroscopic binary;
    2. near-equal-mass companions (q ≳ 0.9) at any a < 5 au, because their photocentre does not move;
    3. companions with ΔT 3.5–4.5 at 0.08–0.8" (5–50 au), where the calibration has few or no pairs.
* **Background (unbound) blends are excluded close to the star.** On the SERC-I 1979 and AAO-SES R 1992 plates the
  M3V lay 3.8–6.7" from its 2019–2026 positions. Injection tests on those plates show no star brighter than T = 18.9
  (I) and T = 18.35–19.0 (R) within 3" of the TESS-epoch positions, except the catalogued Gaia stars. A background
  eclipsing source needs T ≤ 17.47.
  * The two known stars that could still supply the flux are both catalogued:
    1. the background star 3.1" NE (G 17.48, T 16.8, 1.5 kpc), which needs a 55% eclipse;
    2. a star 5.7" away (T 17.2), which needs a 76% eclipse.
  * Excluding them is left to the TESS localisation.
  * The expected number of uncatalogued stars with T ≤ 17.47 within 1" is 0.012.
* **Activity.**
  * **X-ray.** eRASS1 has a faint X-ray source (DET_LIKE 6.2, 7 ± 3 counts), 4.2" (0.8σ) from the M3V. If all of it
    comes from the M3V, log L_X = 28.4 and L_X/L_bol = 4e-4, which corresponds to Ro ~0.17 (P_rot ~11 d).
    * The association is not unique. The M-dwarf pair is 1.7σ away and could contribute.
    * An active, young-looking X-ray level sits in tension with the old kinematics. It is reconciled if the
      emission comes from the pair, if the detection is spurious (it is near the catalogue threshold), or if the
      M3V is a short-period binary.
  * **Rotation and flares.** TESS gives no usable rotation period: only sinusoids above 1.5–5% are detectable. No
    convincing flares are seen.
  * **No other coverage.** GALEX has no coverage of this field, and ROSAT and XMM have no source here.
* **The planet the signal would imply.** For a single host, a depth of 1.4–1.9% gives Rp = 4.6–5.8 R_earth (5.13 +-
  0.15 R_earth at 1.64%, b = 0.6). Other properties:
  * orbit a = 0.0448 au; insolation S = 8.2 S_earth; Teq = 471 K (A = 0) or 431 K (A = 0.3);
  * the QLP duration T14 = 1.63 h fits b = 0.60 (0.54–0.65), with the density implied by the duration equal to the
    star's;
  * the deblended T14 = 1.83 +- 0.12 h needs b = 0.36 (0–0.53);
  * K = 16 (9–29) m/s for the Chen & Kipping mass of 23 (13–42) M_earth;
  * TSM = 130 (71–234), ESM = 10.2. TSM is above the Kempton+18 threshold of 90 for 1.5–10 R_earth planets.
  * In the twin-binary case, Rp = 5.5 R_earth, still planetary. That case is disfavoured if T14 is really 1.8 h
    (84% of draws exceed the central-transit duration).

## Key numbers

| Quantity | Value | Method | Uncertainty | File / script |
|---|---|---|---|---|
| Parallax (zero point corrected) | 15.549 +- 0.036 mas | Gaia DR3 15.4973 + 0.0521 (Lindegren+21, `zero_point`) | error ×1.3 (El-Badry+21) ⊕ 0.015 mas ZP | `star_params.json` (s11) |
| Distance | 64.31 +- 0.15 pc | 1/parallax, MC | MC | s11 |
| M_Ks / M_G | 6.441 +- 0.027 / 10.233 | 2MASS Ks 10.482 +- 0.026 | | s11 |
| Radius | 0.381 +- 0.012 R_sun (0.385 with [Fe/H] = +0.19) | Mann+15 M_Ks relation | 2.9% relation scatter | s11 |
| Mass | 0.370 +- 0.011 M_sun | Mann+19 (n = 5) | 2% scatter + 2% systematic | s11 |
| Density | 9.42 +0.92/−0.83 g/cm3 (6.68 +0.65/−0.59 rho_sun) | M/R^3 | MC | s11 |
| log g | 4.844 +- 0.028 | | MC | s11 |
| Fbol | 1.275e-10 erg/s/cm2 | BT-Settl CIFIST SED fit (11 bands, χ² 4.4/9) | 0.9% statistical, 2% adopted | `out/sed.json` (s06) |
| L | 0.0165 +- 0.0003 L_sun | 4πd²Fbol | MC | s11 |
| Teff (adopted) | 3350 +- 55 K | Stefan–Boltzmann (L, R_Mann) | MC | s11 |
| Teff checks | SED 3220 +- 100; Mann+15 V−J 3280; BP−RP 3382; TIC 3301 +- 157 | | systematic ±100 K | s06, s11 |
| T (TESS) of the M3V | 13.014 (TIC) / 13.086 (SED synthetic) / 13.134 (Gaia colour) | | M3V depth 1.65 / 1.75 / 1.84% for 15.4 e-/s | s06, s08 |
| CMD offset | ΔM_Ks = −0.72, ΔM_G = −0.63 mag; 5.6% / 6.5% of comparison stars as bright | 8,632 Gaia+2MASS M dwarfs (< 100 pc, RUWE < 1.4) | ridge sd 0.44 mag | `out/cmd.json` (s07) |
| [Fe/H] if single | +0.19 +- 0.16 | ΔM_Ks calibrated on Mann+15, Newton+14, Rojas-Ayala+12 (n = 186) | scatter 0.16 dex | `out/metallicity.json` (s10) |
| Equal-mass twin case | 2 × (0.289 R_sun, 0.264 M_sun, 3237 K), rho 15.5 g/cm3 | M_Ks + 0.753 | | `star_params.json` binary_twin_scenario |
| UVW (heliocentric) | (22.0 +- 2.1, −5.6 +- 0.8, +40.4 +- 0.3) km/s; v_LSR 58 km/s | Gaia DR3 astrometry + RV 16.28 +- 2.22 | MC | `out/kinematics.json` (s09) |
| Population / young groups | P_thin 0.70, P_thick 0.30; BANYAN Σ field 1.00 | Bensby+14; BANYAN Σ (Gagné+18) | | s09 |
| Wide companion A | 11.73" = 754 au; Δμ = 1.88 mas/yr → 0.57 km/s (v_esc 1.43 km/s); El-Badry+21 R = 2e-5 | Gaia DR3 | A has RUWE 13.7 | `star_params.json` wide_system (s20) |
| A–B | 0.83" = 53 au; M_G(A) 10.61 (~0.30 M_sun), M_G(B) 11.70 (~0.20 M_sun); chance alignment ~1e-3 | Gaia DR3 | | s20 |
| Gaia IPD exclusion | ΔI < 2 at 0.08–2.5": LR 0.02–0.09; ΔI 2–3.5 at 0.08–0.8": LR ~0.10 | 425 imaged red pairs, 3,166 red single stars | Jeffreys intervals | `out/ipd_limits.json` (s12b) |
| Astrometric exclusion | 0.08–0.3 M_sun at ~0.1–4 au; ≥ 40 M_Jup at 0.24–2.4 au; ≥ 0.34 M_sun at no separation | RUWE + DR3−DR2 PM simulation | L < 0.05 contour | `out/astrometric_limits.json` (s12c) |
| DR3−DR2 PM difference | (+0.396 +- 0.118, −0.236 +- 0.132) mas/yr: 2.4σ after calibration; 6.6% of comparison stars larger | 1,236 comparison M dwarfs (z scatter ×1.50/×1.60) | | s12c |
| DECaPS2 contrast (5σ) | ΔT ≈ 5.2 (1"), 5.7 (1.5–2"), 7.1 (3"), ≥ 8 (4–6") | residual image − crowdsource model, matched filter | 16 PAs per radius | `out/decaps_contrast.json` (s13c) |
| DSS look-back (unbound) | no star with T < 18.9 (SERC-I 1979) / 18.35–19.0 (AAO-SES R 1992) within 3" of the 2019–2026 positions | scene fits + injection, Gaia-calibrated plates (0.06–0.07 dex) | 3σ local noise | `out/dss_limits.json` (s13b) |
| Blending stars expected | 0.0037 /arcsec² with T ≤ 17.47 → 0.012 (1"), 0.046 (2"), 0.10 (3") | Gaia DR3 within 330" (DECaPS2 consistent) | Poisson | `out/blend_density.json` (s16) |
| X-ray | 1eRASS J172352.9-511928; F(0.2–2.3 keV) = 5.4e-14 +- 2.5e-14 cgs; log L_X = 28.43; L_X/L_bol = 4.2e-4 (if all from M3V) | eRASS1 Main | ±46% (counts) | `out/activity.json` (s14) |
| X-ray position | 4.2" (0.78σ) from M3V; 8.9" (1.66σ) from pair; 15.5" (2.9σ) from Mira | per-axis σ = POS_ERR 5.38" | | s14 |
| Rotation limit | 90% detection for 1.5–5% semi-amplitude sinusoids (P 1–12 d) | injection into the deblended LC | per sector | s14 |
| Rp (single host) | 5.13 +- 0.15 R_earth (1.64%, b 0.6); 4.6–5.8 for 1.4–1.9%, b 0–0.8 | batman, quadratic LD (0.20, 0.38) | | `out/planet_params.json` (s15) |
| a, S, Teq | 0.0448 au, 8.2 S_earth, 471 K (A = 0) | | | s15 |
| b from T14 | 0.60 (0.54–0.65) for 1.63 h; 0.36 (0–0.53) for 1.83 h | circular | | s15 |
| K | 6.9 m/s per 10 M_earth; Chen & Kipping 23 M_earth → 16 (9–29) m/s | | | s15 |
| TSM / ESM | 130 (71–234) / 10.2 | Kempton+18 | mass-dependent | s15 |

## 1. Stellar parameters

* **Gaia DR3** (`data/gaia_120.csv`, `data/gaia_ap.csv`, `data/gaia_nss_vari.json`):
  * Astrometry and image diagnostics: parallax 15.497 +- 0.025, RUWE 1.076, 5-parameter solution,
    ipd_frac_multi_peak 0, ipd_gof_harmonic_amplitude 0.013, BP blended in 2/38 transits.
  * Photometry: C* = 0.008 (0.7σ), so the BP/RP photometry is clean.
  * Spectroscopy: RV 16.28 +- 2.22 km/s from 21 transits; the RV error is at the 40th percentile of comparison stars.
  * Catalogue flags: no NSS solution, no variability processing, not flagged as an emission-line star.
  * Astrophysical parameters: GSP-Phot gives none; MSC gives Teff 3298 K. The DSC 'binary' probability is 0.41, but
    48% of comparison M dwarfs score at least that, so it carries no information (`out/gaia_misc_calib.json`).
  * Zhang+23 XP parameters are unusable (χ² 20944).
* **Other photometry** (`data/vizier_*.csv`, positions propagated to each catalogue epoch):
  * 2MASS JHKs, quality AAA.
  * SkyMapper DR4 griz plus v. SkyMapper r is excluded from the fit. SkyMapper r and DECaPS2 r are both ~0.3 mag
    brighter than the best model and agree with each other once the filter difference is allowed for, so the model's
    TiO-dominated r band is at fault, not the survey. Including SkyMapper r gives Teff 3250 K with χ² 80/10.
  * DECaPS2 grizY. The M3V is near saturation; it is used only as a check.
  * unWISE W1/W2. AllWISE W1–W4 are flagged as contaminated by the Mira's halo (ccf = hhh0): W1/W2 are 0.25/0.30 mag
    too bright and W3/W4 come from the Mira, so they are not used.
  * DENIS I; VHS (saturated); APASS (blended with the Mira); no CatWISE source.
* **SED** (s06, `figs/sed.png`):
  * Fit: BT-Settl CIFIST, logg 5, with synthetic photometry. Teff = 3220 K, χ² = 4.4 for 9 dof.
  * Fbol = 1.275e-10 erg/s/cm2. It changes by < 0.5% for logg 4.5/5.5 or when band groups are dropped.
  * The θ-based radius is 0.413 R_sun at that Teff. It equals the radius from Fbol and Teff, so it is not
    independent of Teff.
  * It is 8% above R(M_Ks). The difference reflects the ~100 K BT-Settl Teff scale for M dwarfs, and is also exactly
    what an equal-mass twin would give (sqrt(2) × 0.289 = 0.409). This is a weak hint only.
* **Adopted values** (s11): the single-star values in the table above. The density, 9.4 +- 0.9 g/cm3, is what the
  transit-shape test should use. The twin case gives 15.5 +- 1.5 g/cm3.
* **Metallicity** (s10): there is no spectroscopic value. The photometric value from the CMD offset is +0.19 +- 0.16
  if the star is single.

## 2. The wide system

* **TIC 214311910 A**: plx 16.41 +- 0.39 (RUWE 13.7), pm (−112.80, +87.21) against the M3V's (−113.67, +85.54).
  * The PM difference is 1.88 mas/yr, or 0.57 km/s, below the 1.43 km/s escape speed at 754 au for ~0.87 M_sun in
    total.
  * The 2.3σ parallax difference is not meaningful given A's RUWE.
  * El-Badry+21 catalogue the pair (sepAU 757, R = 2e-5, MSMS).
* **Companion B** (Gaia DR3 5925373272827253632, 2-parameter solution, G 15.74) is 0.83" from A.
  * It explains A's RUWE.
  * A chance alignment within 0.83" has a probability of ~1e-3.
  * At the system's distance M_G = 11.70, an M4 dwarf of ~0.20 M_sun.
* **A's own luminosity cannot be checked cleanly.** A's BP/RP light is blended with B. After removing B's assumed
  colour, A's BP−RP comes out at 2.71–2.87, so its CMD offset is anywhere from −0.1 to −0.6 mag. This cannot
  separate metallicity from binarity for the M3V.
* **What Gaia does not give:** no RV for A or B, no NSS solution, no variability flag. No other co-mover lies within
  2 deg (s17).

## 3. Unresolved companions to the M3V

See `figs/ipd_limits.png`, `figs/astrometric_limits.png` and `figs/companion_exclusion.png`.

* **IPD calibration (s12a, s12b).**
  * Sample: the TIC 71316387 speckle/AO compilation (Tokovinin, Ziegler, Furlan) plus the Lamman+20 Robo-AO M-dwarf
    survey (3,084 stars with Gaia DR3 matches).
  * 'Clean like the target' means fmp ≤ 1, gofha ≤ 0.025, RUWE ≤ 1.2 and no Gaia neighbour within 3".
  * 54% of all single stars are clean, and 37.6% of red ones.
  * Red pairs with ΔI < 2 at 0.08–2.5" are clean in only 1/99, 1/44, 0/65, 0/52 and 0/78 cases, so they are
    excluded.
  * Translation to the G band: ΔG is larger than ΔI by ~0.3–1 mag for redder companions.
* **Astrometry (s12c).**
  * Setup: simulated Gaia-like epochs (36 transits, 18 visibility periods). The per-transit AL error, 0.072 mas, is
    set so that the simulated parallax error equals DR3's.
  * Companion light: β uses the companion's G flux (Mamajek table).
  * Likelihoods: RUWE_0 is drawn from 4,840 comparison M dwarfs (median 1.099; the target's 1.076 is at the 37th
    percentile). The DR3−DR2 PM difference is in calibrated units. L < 0.05 counts as excluded.
* **Why near-equal-mass companions are not excluded:** a twin moves the photocentre by β ≈ 0. The over-luminosity
  allows such a companion at < 5 au.
* **RUWE among similar stars** (s18): only 10% of comparison M dwarfs of similar G and colour have RUWE ≤ 1.076. The
  star is astrometrically quieter than most.

## 4. Neighbour census and blend probability

The census covers every Gaia DR3 source within 30" of the M3V or of the Mira (92), plus DECaPS2 sources within 30"
of the M3V that have no Gaia match (252, nearly all T > 19.5).
* **Files:** `neighbour_census.csv` (copy of `out/neighbour_census.csv`; s08).
* **Columns:** offsets; G, BP−RP, parallax, PM; T and its source; TESS flux (15000 e-/s at T = 10); the eclipse depth
  needed to remove 15.4 e-/s (with a lower bound for −2σ); the maximum depth with 1 and 2 px apertures centred on the
  M3V for the S12 and S104 PRFs (all catalogued light included); a distance class; notes.
* **T calibration for DECaPS2-only sources:** T = −0.215 + 0.998 i + 0.185 (i−z), with 0.024 mag scatter from 289
  Gaia matches.

Sources that could supply the flux, within 15" of the M3V:

| Source | Separation | T | Eclipse depth needed |
|---|---|---|---|
| Gaia ...260544 (background, 1.5 kpc) | 3.13" | 16.8 | 55% |
| Gaia ...262208 | 5.72" | 17.2 | 76% |
| Gaia ...256448 | 10.6" | 17.5 | 99% |
| pair B | 11.35" | 14.5 | 7% |
| pair A | 11.73" | 13.5 | 2.5% |
| Gaia ...256576 | 12.0" | 16.4 | 38% |
| Gaia ...931072 | 13.5" | 17.3 | 85% |

The Mira, at 17.8", would need 0.33%.

In a 1 px aperture on the M3V, the M3V gives only 12.7% of the light; the Mira dominates.

Deeper archives:
* **DECaPS2 2018 image residuals** (s13c): unmodelled sources at 1–6" are limited to ΔT > 5.2–8.6. Below 1" the
  crowdsource model may absorb a companion.
* **DSS look-back** (s13a, s13b; `figs/dss_resid.png`): stars that do not share the M3V's motion are excluded at
  T < 18.35–19.1 within 3" of the 2019–2026 positions (T < 18.35–19.5 within 2"), in the I (1979) and R (1992) plates. The J (1975) and QV
  (1987) plates have poor calibration (0.49 and 0.23 dex) and are not relied on.
* **Not covered:** VVV/VVVX and Legacy Survey DR10. VHS is saturated on the M3V and its catalogue shows only the Gaia
  stars.

**Star counts** (s16; l = 338.83, b = −8.55): 0.0037 stars/arcsec² with T ≤ 17.47 (Gaia, confirmed by DECaPS2).
The expected number within 1" / 2" / 3" is 0.012 / 0.046 / 0.10. The look-back test removes most of this for unbound
stars. TRILEGAL was not used, because Gaia is complete at these magnitudes.

## 5. Activity and age

* **X-ray** (s14). The emission comes from the M-dwarf system, not from the Mira or an unrelated source:
  * Chance coincidence with an unrelated eRASS1 source within 10" has a probability of 1e-3.
  * The background star 3.1" NE has the same positional likelihood as the M3V. At 1.5 kpc it would need
    L_X ~ 1.5e31 erg/s, a very low prior.
  * The detection is near the DET_LIKE = 6 threshold, so it may be spurious.
  * If real and from the M3V, Wright+11 gives Ro ≈ 0.17 (0.23 if half of it), so P_rot ≈ 11–15 d.
* **Kinematics** (s09): W_LSR = 48 km/s. This would be a > 3σ outlier for a ≲ 2 Gyr population (σ_W ≲ 15 km/s), so
  the star is probably several Gyr old.
* **The tension between them is unresolved.** If the X-ray source is real and from the M3V, a tidally spun-up close
  binary (case b, at a ≲ 0.1 au) would reconcile both. Spectroscopy could test this.
* **Implications for the planet:**
  * Spots or faculae cannot be ruled out below the 1.5–5% variability limit.
  * Unocculted spots deepen transits by a factor of ~1/(1 − f_spot × (1 − contrast)). That is ≤ 1–3% for low spot
    coverage, which is negligible next to the depth spread from the Tmag/deblending choice (1.65–1.84%).

## 6. Planet parameters

`out/planet_params.json` (s15) gives, for depths 1.4/1.64/1.75/1.9% and b = 0–0.9: k, Rp, T14 and T23/T14 for the
single and twin cases.
* **Durations.** At 1.64%, T14 = 1.93 h (b = 0), 1.64 h (b = 0.6) and 1.37 h (b = 0.8).
* **Implied density.** Circular orbit, from T14: 9.5 g/cm3 for T14 = 1.63 h and b = 0.6, equal to the star's
  9.4 +- 0.9.
* **RV semi-amplitudes:** 7 / 14 / 21 / 35 / 70 m/s for 10 / 20 / 30 / 50 / 100 M_earth.

## Caveats and corrections to earlier statements

* **Distance.** 64.5 pc was the uncorrected parallax. With the Lindegren+21 zero point it is 64.31 +- 0.15 pc. The
  difference is negligible.
* **Spectral type.** "M3V" by luminosity; by BP−RP (2.76) it is M3.5V. The two disagree because of the
  over-luminosity.
* **Depth normalisation.** It is uncertain by ~10% from the M3V's own T: 13.01 (TIC), 13.09 (SED), 13.13 (Gaia colour).
  This propagates as ±5% into Rp.
* **Calibration samples are mostly brighter stars.** The IPD calibration has few pairs at ΔI 3.5–5. The astrometric
  simulation uses an approximate scanning law, random scan angles and a Mamajek mass–luminosity relation that ignores
  the star's own over-luminosity.
* **Photographic plates.** The DSS limits come from photographic plates with non-linear response. The calibrations
  for the I/R plates are tight (0.06–0.07 dex), and limits are from injections at the actual positions. Positions
  within 1.5" of catalogued Gaia stars (including the 3.1" background star) were not tested.
* **Twin versus metal-rich single.** The two cases cannot be separated with archival data. High-resolution
  spectroscopy, giving [Fe/H], SB2 lines or an RV trend, and AO/speckle imaging below 0.08" are needed.

## Scripts (run in this order; outputs in `data/`, `out/`, `figs/`)

| Script | What it does | Outputs |
|---|---|---|
| `s01_gaia.py` | Gaia DR3 cone 120" around the M3V; astrophysical_parameters of key stars. Its later steps are superseded by s01b | `data/gaia_120.csv`, `data/gaia_ap.csv` |
| `s01b_gaia_extra.py` | NSS/variability tables, DR2 counterparts, 2MASS/AllWISE best neighbours | `data/gaia_nss_vari.json`, `data/gaia_dr2.csv`, `data/gaia_xm.csv` |
| `s01c_gaia_wide.py` | Gaia DR3 within 330" (plate calibration, star counts) | `data/gaia_330.csv` |
| `s02_phot.py` | VizieR photometry/X-ray/UV/variability catalogues within 20", epoch-propagated offsets | `data/vizier_*.csv`, `data/phot_raw.json` |
| `s02b_more_cats.py` | El-Badry+21, Zhang+23 XP parameters, full eRASS1 records | `data/vizier2_*.csv` |
| `s03a_models.py` | BT-Settl CIFIST spectra (SVO) and filter curves | `data/models/` |
| `s04_mdwarf_sample.py` | Gaia+2MASS+DR2+DSC comparison sample (parallax > 10 mas, BP−RP 2.45–3.05) | `data/mdwarf_sample.csv` |
| `s05_decaps.py` | DECaPS2 catalogue (Data Lab TAP) 90" and 5' | `data/decaps2_90.csv`, `data/decaps2_field10.csv` |
| `s06_sed.py` | SED fit, Fbol, synthetic T | `out/sed.json`, `figs/sed.png` |
| `s07_cmd.py` | Main-sequence ridges and the target's offsets | `out/cmd.json`, `figs/cmd.png` |
| `s08_census.py` | Neighbour census with PRF aperture dilution | `out/neighbour_census.csv`, `out/census_summary.json` |
| `s09_kinematics.py` | UVW, Bensby+14 populations, BANYAN Σ (local copy of the public code in `data/ext/banyan`, numpy-2 shim only) | `out/kinematics.json` |
| `s10_metallicity.py` | [Fe/H]-offset calibration (Mann+15, Newton+14, Rojas-Ayala+12); pair-A deblending | `out/metallicity.json`, `data/feh_calibrators.csv`, `figs/metallicity.png` |
| `s11_star_params.py` | Adopted parameters (MC) and the twin case | `star_params.json`, `out/star_params_mc.npz` |
| `s12a_mdwarf_imaging_sample.py` | Lamman+20 Robo-AO (and Vrijmoet+22) samples with Gaia DR3 diagnostics | `data/mdw_imaging_sample.csv` |
| `s12b_ipd_limits.py` | IPD/RUWE likelihood map vs separation/contrast | `out/ipd_limits.json`, `figs/ipd_limits.png` |
| `s12c_astrometric_sim.py` | RUWE + DR3−DR2 PM companion simulation | `out/astrometric_limits.json`, `figs/astrometric_limits.png` |
| `s13a_dss_look.py` | DSS quick-look | `figs/dss_look.png` |
| `s13b_dss_limits.py` | DSS look-back limits (scene fits, injections) | `out/dss_limits.json`, `out/dss_calib_*.npy`, `figs/dss_resid.png` |
| `s13c_decaps_contrast.py` | DECaPS2 residual-image contrast curve | `out/decaps_contrast.json`, `figs/decaps_resid.png` |
| `s14_activity.py` | X-ray, GALEX coverage, rotation limits, flares (uses `../phot/lc_*_all.csv`) | `out/activity.json`, `figs/rotation.png` |
| `s15_planet.py` | Planet parameters, durations, K, TSM/ESM | `out/planet_params.json` |
| `s16_blend_density.py` | Star densities and expected blends | `out/blend_density.json` |
| `s17_comovers.py` | Co-moving companion search within 2 deg | `out/comovers.json` |
| `s18_gaia_misc_calib.py` | DSC binary probability, RV error, RUWE, gofha percentiles | `out/gaia_misc_calib.json` |
| `s19_exclusion_figure.py` | Combined exclusion map | `figs/companion_exclusion.png`, `out/exclusion_summary.json` |
| `s20_finalize.py` | Adds the system, kinematics, activity, limits, blend and planet summaries to `star_params.json` | `star_params.json` |

Other files:
* `data/ext/EEM_dwarf.txt`: Mamajek dwarf table, downloaded 2026-10-16.
* `data/ext/mann2015_table.ecsv`.
* `data/img/dss/*.fits`: DSS2 3' and 8' scans.
* `data/img/decaps/*.fits`: DECaPS2 coadd and model cutouts, 63".
* `out/s0*.log`: console logs.

## For the next phase

* **Stellar parameters:** `star_params.json`.
  * Top-level single-star values: each is a dict with median, err_lo and err_hi.
  * `binary_twin_scenario` covers the equal-mass case.
  * `teff_checks` lists the alternative Teff estimates.
  * Use rho_star = 9.42 +0.92/−0.83 g/cm3 (6.68 rho_sun), R = 0.381 +- 0.012 R_sun, M = 0.370 +- 0.011 M_sun,
    Teff = 3350 +- 55 K, quadratic TESS LD ≈ (0.20, 0.38).
* **MC samples:** `out/star_params_mc.npz` (R, M, L, Teff, rho, d; 20,000 draws).
* **Neighbour census:** `neighbour_census.csv`. Read source_id as a string.
* **Exclusion maps:** `out/ipd_limits.json` (cells), `out/astrometric_limits.json` (M2 × a grid of likelihoods),
  `out/decaps_contrast.json`, `out/dss_limits.json` and `out/exclusion_summary.json`.
