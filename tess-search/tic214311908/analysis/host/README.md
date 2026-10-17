# TIC 214311908: single metal-rich star or unresolved twin?

Folder: `tic214311908/analysis/host/`

*AI-run analysis, 2026-10-17. Nothing was posted, uploaded or submitted anywhere. All numbers were derived here from
public archive data (Gaia DR3 archive and VizieR, CDS XMatch, arXiv PDFs for published coefficients). Phase-1 files
(`../star`, `../referee_star`) were read as inputs only.*

## Bottom line

* **The over-luminosity is real.**
  * ΔMG = −0.68 mag comes from Gaia G alone, and ΔMKs = −0.76.
  * 2MASS Ks is right. DENIS is fainter only because its frame has a zero-point error (+0.30 in K and +0.385 in J for
    every star in that image). VHS saturated-star photometry and unWISE agree with 2MASS.
  * The XP spectrum shows no contamination from the Mira: at most 0.7% of the 800-nm flux.
* **Every luminosity-independent metallicity indicator says the M3V is metal-rich, at [Fe/H] ≈ +0.2 to +0.3.** That
  is exactly what a single star needs (+0.25 ± 0.15). A near-equal twin would need each star at about −0.13 ± 0.15
  (essentially solar).
  * **Colour excess at fixed BP−RP:** r(G−Ks) = +0.089 gives +0.26; r(G−W1), independent of 2MASS, = +0.092 gives
    +0.33. A twin's colours are those of a single star of its type, so it would show about −0.015.
  * **Gaia XP spectral shape:**
    * a regression on 738 calibrators gives +0.13 raw, i.e. +0.16 to +0.21 once calibrated;
    * it lies above 97% of normal stars of the same colour;
    * its residual pattern correlates with that of metal-rich calibrators at r = 0.71 (98.8th percentile).
  * **The bound wide pair** TIC 214311910 A+B, from its blended colour: +0.13 ± 0.11.
  * The published XP value from **Gaia Net (−0.03) is uninformative for M dwarfs** (r ≈ 0 against calibrators) and is
    not evidence for a twin.
* **The Gaia RV error** (61st percentile of comparable stars) excludes SB2 twins at a ≲ 0.1 au and disfavours them out
  to ~0.4 au. **Gaia astrometry** excludes q ≤ 0.93 companions at ~0.4–2.4 au. **Only q ≥ 0.95 twins at ~0.3–5 au
  remain hidden.** These are about 2% of M dwarfs from literature multiplicity statistics, and 0.4% from the field
  population mixture.
* **Posterior odds single : twin ≈ 100 : 1.** P(twin) ≈ 1%, with a range of 0.2–2.5% over priors, calibrations and
  photometric systematics; the CMD alone gives 6–10%.
  * Adopt a **single host with [Fe/H] ≈ +0.2 ± 0.1**: R = 0.385 ± 0.012 R_sun, M = 0.370 ± 0.012 M_sun,
    ρ ≈ 9.2 g/cm³.
  * **Planet radius:** Rp = 5.16 R_earth if single (5.11 for the solar-metallicity radius); 5.4–5.6 R_earth in the
    surviving twin case. The twin case is not 5.5–6 R_earth: only +6–8%.
* **For the false-positive analysis:**
  * The near-equal twin case drops to ~1% prior weight.
  * Even if real, it would not be a false-positive channel: an equal bound companion would host a ~3% transit
    (Rp/R* ≈ 0.18), which is still planetary.
  * Nothing here constrains fainter bound companions (ΔT 3.5–4.5), the 4–5 au companions, hierarchical EBs, or the
    background stars at 3.1" and 5.5".

## Key numbers

| Quantity | Value | Method | Uncertainty / calibration | File |
|---|---|---|---|---|
| CMD offsets | ΔMKs −0.761, ΔMG −0.680 | Field-ridge KDE mode (50,583 M dwarfs) | ±0.03 (meas.) | `out/h13_colour_feh.json` |
| DENIS−2MASS frame offset (image 1724791) | J +0.385, Ks +0.302 | Same-frame field stars (n = 17 / 30) | robust sd 0.040 / 0.038; M3V +0.411 / +0.341 | `out/h06b_nir_local.json` |
| DENIS Ks, zero-point corrected | 10.52 ± 0.10 (2MASS 10.482 ± 0.026) | Frame offset applied | – | h06b |
| VHS Ksap3 (saturated) | 11.167 vs 11.159 for stars with the same 2MASS Ks | Forward test, n = 7 | sd 0.021; Ks good to ~±0.1 | h06b |
| unWISE W1 / W2 (AllWISE scale) | 10.318 / 10.199 | Local unWISE−AllWISE offsets +0.014 / +0.030 | ±0.013 | h06b |
| Photometric [Fe/H] (single; twin) | B05 +0.07 (−0.15); JA09 +0.47 (+0.04); SL10 +0.29 (−0.08); N12 +0.16 (−0.10); R21 +0.29 (+0.02) | Published relations, V(XP) = 15.591 | ±0.17–0.24 | `out/h08_photfeh.json` |
| [Fe/H] a single star needs | +0.25 ± 0.15 (+0.33 high-precision set) | ΔMKs = −0.282 − 1.880 [Fe/H] (n = 1521) | sd 0.29 mag | h13 |
| [Fe/H] each twin component needs | −0.15 ± 0.15 | Same relation, offset −0.753 | – | h13, h15 |
| r(G−Ks), r(G−W1) at fixed BP−RP | +0.089, +0.092 | Field medians | ±0.027, ±0.03 | h13 |
| Colour [Fe/H] | +0.26 (G−Ks), +0.33 (G−W1), +0.24 (G−H) | 3σ-clipped calibrator fits | sd 0.03 mag → ±0.12–0.15 dex | h13 |
| XP-shape [Fe/H] | +0.13 raw (+0.21 de-biased; calibrators with the same prediction: +0.16, 16–84% +0.04/+0.27) | PCA + ridge, 738 calibrators | CV 0.11 dex; noise ±0.08 | `out/h12b_xp.json` |
| XP residual vs metal-rich pattern | r = 0.71 (98.8th percentile) | Correlation with calibrator means | Null: 84 normal field stars | `out/h12e_residual_corr.json` |
| Mira light in the XP | ε = 0.0014 ± 0.0012 (≤ 0.7% at 800 nm) | Template + Mira fit | Null from field stars | `out/h12c_xp_population.json` |
| XP twin LR (twin/single) | 0.40–0.47 | Over-luminous field mixture (twin fraction 0.16–0.18) | Model-dependent | h12c |
| Gaia Net / SHBoost [Fe/H] | −0.03 (uninformative, r ≈ 0) / −0.16 (r = 0.4, uses parallax) | Calibrated on Mann-scale stars | – | `out/h16_lit_xp_feh.json` |
| Wide pair A+B [Fe/H] | +0.13 ± 0.11 | Blend colour given separate G | Single needs +0.29, twin −0.11 (±0.09) | `out/h14_wide_pair_feh.json` |
| Gaia rv_error percentile | 61% (n = 928 comparison stars) | Empirical | – | `out/h07_gaia_rv.json` |
| Twins flagged by the RV error | 96% (0.04 au), 82% (0.1), 54% (0.25), 34% (0.4), 12% (1 au) | CCF simulation | σ_CCF 11–20 km/s | h07 |
| Twin astrometric exclusion | q 0.80–0.85 at 0.3–2.4 au; q 0.90–0.93 at 1.0–1.3 au; q ≥ 0.95 unconstrained | RUWE + DR3−DR2 PM simulation | L < 0.05 | `out/h09_twin_astrometry.json` |
| P([Fe/H] > +0.2 / > +0.3), kinematic analogues | 7.3% / 3.1% (n = 763; tight analogues 3.7% / 2.4%) | GSP-Spec dwarfs < 200 pc | Poisson | `out/h10b_feh_prior.json` |
| Twin prior (q ≥ 0.95, a < 5 au) | 2.4% before the Gaia tests; 1.9% after them (q ≥ 0.8: 3.9% / 2.35%) | Winters+19, El-Badry+19, Susemiehl & Meyer 22 | ×0.5–1.5 over variants | h15 |
| Twin fraction among RUWE < 1.1 field stars | 0.43% (0.23–0.69%) in the (ΔMKs, r(G−Ks)) plane; 0.73% (0.34–1.04%) in (ΔMG, r(G−W1)) | Structured population mixture | Bootstrap (12) | `out/h17b_structured_mixture.json` |
| **P(twin)** | **≈ 1% (0.2–2.5%)**; population-based 0.5–1.9% before the XP and wide-pair factors | Prior-based (h15) and population-based (h17b) | Variants | h15, h17b |
| Rp | single 5.16 R_earth (5.11); twin 5.40–5.56 R_earth | Depth 15.7 ppt, LD factor 1.04 | ±0.2 R_earth from the depth error | `out/h15_posterior.json` |

## 1. Near-IR photometry: which is right, 2MASS or DENIS? (h06a, h06b, h08)

**The 2MASS photometry is right. The DENIS values are fainter because the DENIS frame has a zero-point error. The
over-luminosity does not come from the near-IR.**

* **2MASS (1999.5):**
  * J 11.363 ± 0.028, H 10.740 ± 0.031, Ks 10.482 ± 0.026.
  * Qflg AAA; Cflg 'ddd'. The 'd' (diffraction-spike) flag is plausible: the M3V lies almost due west of the Mira
    (17.8", PA ~84°).
* **DENIS (1998.8):** J 11.774 ± 0.08, Ks 10.823 ± 0.09, I 13.094 ± 0.04.
* **The DENIS frame is offset** (h06b; local 8' catalogues; 2MASS-matched stars with Qflg A, Cflg 000).
  * In the DENIS image that measured the M3V (1724791), every star is fainter than in 2MASS:
    * J by **+0.385 ± 0.040** (robust sd; n = 17, J 10–11.5);
    * Ks by **+0.302 ± 0.038** (n = 30).
  * The M3V's own offsets, +0.411 (J) and +0.341 (Ks), are 0.7σ and 1.0σ from the frame values.
  * Corrected for the frame offset, DENIS gives J = 11.39 ± 0.09 and Ks = 10.52 ± 0.10. Both agree with 2MASS to
    0.03–0.04 mag.
  * DENIS I shows the same offset. For 25 local stars, DENIS I is 0.21 ± 0.08 mag fainter than I_C predicted from
    Gaia G and BP−RP (Riello+21). For the M3V the offset is +0.15.
* **VHS DR5 (2014.6)**
  * VHS saturates at the M3V's brightness. Saturated aperture magnitudes still increase monotonically with true
    brightness, so they can be compared with those of other stars.
  * Local stars whose 2MASS magnitudes are within ±0.12 of the M3V's have:
    * Ksap3 = 11.159 (sd 0.021, n = 7); the M3V has 11.167 (71st percentile);
    * Jap3 = 11.894 (sd 0.055); the M3V has 11.940.
  * The saturated response is shallow (slope 0.2–0.35), so VHS only constrains Ks to about ±0.1 mag at 1σ. It shows
    no sign of a contamination of 2MASS Ks larger than ~0.1 mag at a third epoch and with different optics.
* **unWISE (Mira modelled; flags 'halo')**
  * W1 = 10.332 and W2 = 10.229.
  * The local unWISE−AllWISE offsets are +0.014 (W1) and +0.030 (W2) ± 0.013 (n = 25/22). On the AllWISE scale this
    gives W1 = 10.318 and W2 = 10.199.
  * At fixed BP−RP, the Ks−W1 residual is +0.002, the H−Ks residual +0.002 and the J−Ks residual +0.045. 2MASS and
    unWISE are mutually consistent.
  * AllWISE itself is contaminated by the Mira's halo (W1 10.112, ccf 'hhh0') and is not used.
* **Size of the over-luminosity.** It is
  * ΔMKs = −0.761 with 2MASS;
  * ΔMKs = −0.72 with zero-point-corrected DENIS;
  * **ΔMG = −0.680 from Gaia G alone.**

  The raw DENIS Ks would give −0.42, but those values are wrong. A Ks error therefore cannot remove the
  over-luminosity.
* **G and BP−RP are not contaminated by the Mira.** C* is 0.7σ, the G−RP residual at fixed BP−RP is −0.001, and the XP
  fit allows less than 0.7% Mira light at 800 nm (section 2).
* **Published photometric [Fe/H] relations** (h08)
  * Inputs: V = 15.591 ± 0.03 from XP synthetic photometry (Riello+21: 15.608). The coefficients were checked against
    the papers in `data/lit/`.

  | Relation (scatter) | single star, 2MASS | equal twin, 2MASS | single, DENIS corrected | single, DENIS raw |
  |---|---|---|---|---|
  | Bonfils+05 (0.20) | +0.07 | −0.15 | +0.03 | −0.21 |
  | Johnson & Apps 09 (0.24) | +0.47 | +0.04 | +0.41 | −0.02 |
  | Schlaufman & Laughlin 10 (0.19) | +0.29 | −0.08 | +0.24 | −0.14 |
  | Neves+12 (0.17) | +0.16 | −0.10 | +0.13 | −0.15 |
  | Rains+21 BP−Ks (0.19) | +0.29 | +0.02 | +0.25 | −0.07 |

  * These relations only restate the CMD offset. They cannot separate a metal-rich single star from a solar twin.
    Neves+12 found that JA09 over-estimates [Fe/H] by 0.14 dex and B05 under-estimates it.

## 2. Gaia XP spectrum (h01, h05, h11b, h12b-e, h16)

**The XP spectrum is clean, and its shape is that of a mildly metal-rich M dwarf. It shows no composite signature.**

An equal-mass twin has exactly the XP shape of a single star of the same type, so XP cannot see a q ≈ 1 twin
directly. What XP can do is measure [Fe/H] without using the parallax. The CMD then implies:
* about +0.25 if the star is single;
* about −0.13 for each component if it is a twin.

* **Data quality.**
  * Spectra come from the Gaia DR3 archive DataLink (h01) and the official sampled spectra in VizieR I/355/xpsample
    (h11b). The two agree to 0.04% rms.
  * BP: 36 transits, 1 blended. RP: 39 transits, 1 contaminated.
  * Standard deviations 1.07 (BP) and 1.17 (RP); S/N 50–100 per 2-nm sample in RP.
  * The A(+B) entry is blended (31/33 blended transits). It has a spurious bump at 1000 nm and is unusable for shape
    work.
* **Contamination by the Mira** (h12b, h12c).
  * The M3V's shape deviates from the median of 84 colour-matched normal field stars by an rms of 1.07σ. That is the
    61st percentile of the field stars themselves.
  * A fit of template + ε × (Mira XP flux) gives ε = 0.0014 ± 0.0012. That is 1.2σ from the null distribution measured
    on field stars, and means at most 0.7% (1σ) of the M3V's 800-nm flux comes from the Mira. **No contamination.**
* **XP synthetic photometry** (GaiaXPy JKC_Std):
  * V = 15.591 and I_C = 12.859;
  * DECam-system g/r/i/z/Y match DECaPS2 within the saturation effects (`out/xp_synphot_target.csv`).
* **XP-shape [Fe/H]** (h12b). This regression uses no parallax and no magnitudes.
  * Model: log-flux shape over 420–1000 nm, 30 PCA components plus BP−RP, ridge regression.
  * Training set: 738 Mann-scale calibrators with BP−RP 2.45–3.10 and |ΔMKs| < 1.3.
  * Performance: 10-fold CV robust scatter 0.11 dex, r = 0.78, CV slope 0.64 (the estimator regresses toward the
    mean).
  * **M3V: raw prediction +0.13 ± 0.08 (flux noise) ± 0.11 (CV), de-biased +0.21.** Calibrators with the same raw
    prediction have true [Fe/H] of +0.16, with a 16–84% range of +0.04 to +0.27.
  * Comparison with field stars of the same colour:
    * normal-CMD field stars: median −0.04 (16–84%: −0.15 to +0.05); the M3V lies above 97% of them;
    * over-luminous field stars (ΔMKs −1.0 to −0.55): median +0.13; the M3V is at their median.
* **Residual pattern** (h12e). The M3V's residual from the normal-star median correlates with:
  * the mean residual of metal-rich calibrators ([Fe/H] > +0.25): r = 0.71, the **98.8th percentile** of normal field
    stars;
  * the mean residual of metal-poor calibrators: r = 0.13 (62nd percentile);
  * the mean pattern of over-luminous field stars: r = 0.71 (97.6th percentile).
* **XP twin null** (h12c). The XP-[Fe/H] distribution of over-luminous field stars was fitted as twins (shaped like
  normal stars) plus metal-rich singles, with the width of the metal-rich component fixed.
  * Twin fraction: **0.16 (0.11–0.20)** for RUWE < 1.2 and 0.18 (0.06–0.27) for RUWE < 1.1.
  * M3V likelihood ratio twin/single: **0.40–0.47**.
  * P(twin | XP, over-luminous): 0.07–0.09.
  * With a free width the mixture is degenerate (twin fraction → 0). The fixed-width model is the generous one for
    twins.
* **Single vs composite templates** (h12d). The templates are colour-binned medians of normal stars; companions are
  weighted by their M_G.
  * A two-star composite improves χ² over a single template by 4.2 for the M3V.
  * Normal field stars improve by a median of 26 (the M3V is below 83% of them). Metal-rich calibrators improve by a
    median of 2.7.
  * There is no composite signature, as expected for a single star or a q ≈ 1 twin.
* **Published XP-based values** (h16; calibrated on Mann-scale calibrators within 0.1 mag of the M3V's colour).
  * **Gaia Net** (Huson+25) uses XP coefficients only and gives [Fe/H] = −0.03. **It carries no information for M
    dwarfs:** its correlation with spectroscopic [Fe/H] is r ≈ 0, and it returns −0.04 to −0.03 for almost every M
    dwarf, normal or over-luminous. Its "solar" value must not be read as evidence for a twin.
  * **SHBoost** gives −0.16. It uses parallax and photometry, so it is not luminosity-free. On M dwarfs it has r = 0.4,
    slope 0.28 and offset −0.21.
  * **StarHorse** (+0.29) uses photometry and parallax, so it is just the CMD offset.
  * **Zhang+23** failed (χ² 20944). The star is not in **Andrae+23**, which requires CatWISE, absent here because of
    the Mira.

## 3. Colour-only (luminosity-independent) metallicity and the CMD (h13, h17, h17b)

**The M3V's colours at fixed BP−RP are those of a metal-rich star ([Fe/H] +0.26 to +0.33), not those of a solar-metallicity
twin (−0.15).**

This is the most direct discriminant. At a given BP−RP:
* metallicity makes an M dwarf brighter **and** redder in G−IR;
* an unresolved near-equal companion makes it brighter only. With the Mamajek main sequence at fixed combined BP−RP,
  Δr(G−Ks) = 0.000 (q = 1) to −0.015 (q = 0.8).

So the colour residual measures the composition, whichever hypothesis holds.

* **Ridges** (h13): KDE mode of M_Ks and M_G, and running medians of colours, against BP−RP. They come from 50,583
  field M dwarfs (RUWE < 1.2, 2MASS AAA, |C*| < 3).
* **M3V offsets:**
  * ΔMKs = −0.761 and ΔMG = −0.680;
  * colour residuals: r(G−Ks) = **+0.089**, r(G−H) = +0.089, r(G−J) = +0.044, r(J−Ks) = +0.045, r(H−Ks) = +0.002,
    r(G−W1) = **+0.092** (unWISE, AllWISE scale), r(Ks−W1) = +0.002, r(G−RP) = −0.001, r(BP−G) = +0.001.
* **Calibration** (1521 Mann-scale calibrators; RUWE < 1.1, |ΔMKs| < 1.3, 3σ clipping; the sd includes [Fe/H]
  label noise).

  | Quantity | Relation | sd | [Fe/H] implied for the M3V |
  |---|---|---|---|
  | ΔMKs | −0.282 − 1.880 [Fe/H] | 0.29 mag | +0.25 (the [Fe/H] a single star needs) |
  | ΔMG | −0.255 − 1.629 [Fe/H] | 0.29 | +0.26 |
  | r(G−Ks) | +0.024 + 0.250 [Fe/H] | 0.030 | **+0.26** |
  | r(G−H) | +0.029 + 0.242 [Fe/H] | 0.032 | +0.24 |
  | r(J−Ks) | +0.017 + 0.182 [Fe/H] | 0.033 | +0.15 |
  | r(G−W1) (n = 1317) | +0.016 + 0.230 [Fe/H] | 0.028 | **+0.33** |

  * With the 62 high-precision (non-Birky) calibrators:
    * ΔMKs = −0.138 − 1.877 [Fe/H] (sd 0.20), so a single star needs +0.33;
    * r(G−Ks) = +0.022 + 0.193 [Fe/H], giving +0.35.
* **Colour-only ridge regression** (10 colour residuals and BP−RP, no absolute magnitude).
  * Performance: CV scatter 0.08 dex, r = 0.79, forward slope 0.63.
  * M3V: raw +0.15 ± 0.13, de-biased ≈ +0.25.
  * The largest weight is on r(W1−W2), which is the least trustworthy here: unWISE W2 fracflux is 0.82 in the Mira's
    halo.
* **What each hypothesis predicts at the observed ΔMKs:**
  * single star ([Fe/H] = +0.25): r(G−Ks) = +0.087 and r(G−W1) = +0.074;
  * twin (components at −0.15): −0.013 and −0.018.
  * Observed: +0.089 ± 0.027 (meas.) and +0.092 ± 0.03.

  The twin prediction is 2.5–3σ off once the 0.03 mag intrinsic scatter is added; the single-star prediction fits
  within 0.5σ.
* **ΔMG/ΔMKs.** The M3V has 0.894; metallicity predicts 0.866 (0.90 with the high-precision set); a twin predicts 1.00.
  This is the same G−Ks information, not an independent test.
* **Population test** (h17b; no [Fe/H] labels, no literature priors).
  * Sample: 3807 field stars like the M3V (|ΔBP−RP| < 0.08, RUWE < 1.1, ipd_frac_multi_peak = 0, gof_harm < 0.03).
  * Structured mixture: singles lie along a metallicity line, r(G−Ks) = c + r ΔMKs; twins are the same population
    shifted by −0.753 mag; plus unequal binaries and outliers.
  * The singles' slope from the data, r = −0.085, is close to the high-precision calibrator ratio (0.193/−1.877 =
    −0.103). The Birky-label ratio is −0.133. The twin weight is **0.43% (0.23–0.69%)**.
  * At the M3V's ΔMKs the singles line predicts r(G−Ks) = 0.064 and the twins line 0.000. The M3V has 0.089, so
    **P(twin) = 0.005 (0.003–0.008)**.
  * Variants: 
    * Ks 0.04 mag fainter (in case the 2MASS diffraction spike contaminated it): P(twin) = 0.012 (0.007–0.021).
    * **2MASS-free plane (ΔMG, r(G−W1)):** singles slope −0.064, twin weight 0.73% (0.34–1.04%). The M3V has
      r(G−W1) = 0.091 against 0.040 on the singles line and −0.009 on the twins line, so **P(twin) = 0.011
      (0.005–0.016)**. With W1 0.03 mag fainter, P(twin) = 0.019.
    * v_tan 30–70 km/s subsamples (n = 1707 and 1325): the twin component collapses to zero and P(twin) → 0. At these
      sample sizes the model is not identifiable, so these runs are not used.
  * A free (unstructured) 2D Gaussian mixture (h17) is not identifiable: it absorbs every star as 'single'. It does
    reproduce the observed 2.7% of stars brighter than the M3V.

## 4. The wide companions as a metallicity probe (h14)

**The wide pair favours a mildly metal-rich system ([Fe/H] ≈ +0.13 ± 0.11): 1.1σ from the single-star requirement
(+0.29) and 1.7σ from the twin requirement (−0.11). The likelihood ratio single:twin is ≈ 2.**

* **Inputs.** TIC 214311910 A (G 14.654) and B (G 15.741) lie 0.83" apart, so they share one Gaia BP/RP measurement.
  * The blend has BP−RP = 2.878 at M_G(A+B) = 10.231.
  * The excess factor computed with G_A + G_B is 1.446, which is consistent with a complete blend
    (C* = −1.7σ, i.e. at most ~1% of B's flux lost from the window).
* **Model.** Mann-scale calibrators with RUWE < 1.1 and |ΔMKs| < 1.5 (n = 2486) give
  M_G = poly3(BP−RP) + (c0 + c1(BP−RP − 2.75))[Fe/H], with dM_G/d[Fe/H] = −1.88 at BP−RP 2.75 and 0.17 mag intrinsic
  scatter.
  * For each [Fe/H], each component's colour is set by its M_G, with that scatter included.
  * BP−G and G−RP come from robust field relations.
  * The likelihood compares the predicted blend colour with the observed one (σ = 0.022 mag).
* **Result** (flat prior): **[Fe/H](A+B) = +0.13 ± 0.11**.
* **The same model applied to the M3V** requires +0.29 if single or −0.11 if a twin (±0.09).
* **Caveats.**
  * A or B could themselves be unresolved binaries. A's RUWE of 13.7 is attributed to B.
  * B's G magnitude comes from a 2-parameter solution.
  * DECaPS2 resolves the pair only in g and r (Δr = 1.45 against ΔG = 1.09, as expected for a redder B). It gives no
    separate colours.
  * This is a weak test. In h15 it is used as a likelihood in [Fe/H].

## 5. Gaia RV and astrometry against an SB2/astrometric twin (h07, h09)

**The Gaia RV error rules out twins inside ~0.1 au and disfavours them inside ~0.4 au. Astrometry rules out q ≤ 0.93
at ~0.4–2.4 au. Near-equal twins (q ≥ 0.95) at ~0.3–5 au remain allowed.**

* **What Gaia RV provides.**
  * DR3 measured this faint star with rv_method_used = 2 (grvs_mag 12.55): the RV comes from the combined CCF of
    21 transits.
  * rv_chisq_pvalue, rv_renormalised_gof and rv_amplitude_robust are not published for such stars.
  * Template: 3500 K / logg 3.5 / [Fe/H] −0.5. 12 of the 21 transits were deblended (A and the Mira are nearby).
* **Empirical comparison** (h07; field stars with method 2, G_RVS 12.30–12.80, expected S/N within ±25%, 15–30
  transits, 2MASS AAA).
  * 928 normal stars (RUWE < 1.2, |ΔMKs| < 0.25) have rv_error median 2.11, p84 2.69 and p95 3.58 km/s.
  * **The M3V's 2.22 km/s is at the 61st percentile.** Among normal stars with more than 30% deblended transits it is
    at the 64th.
  * Binaries do inflate this error empirically. 17.8% of over-luminous low-RUWE stars have an error above 1.5× the
    median, against 7.9% of normal stars (and 13% of stars with RUWE > 2).
* **Simulation** (h07).
  * Setup: equal twin of 2 × 0.264 M_sun; per-transit CCF = sum of two Gaussians (σ_CCF = 15 km/s, with 11 and 20
    as checks); 21 transits in 8 visibility periods over 974 d; error inflation from Zucker's curvature formula at the
    highest peak; 0.5 km/s floor.
  * Fraction of twins whose error would exceed the singles' 95th percentile, i.e. be flagged:

    | a (au) | 0.01 | 0.04 | 0.1 | 0.25 | 0.4 | 0.63 | 1.0 | ≥ 1.6 |
    |---|---|---|---|---|---|---|---|---|
    | P (d) | 0.5 | 4 | 16 | 63 | 126 | 252 | 503 | ≥ 1000 |
    | Flagged | 99% | 96% | 82% | 54% | 34% | 20% | 12% | ~5% (false-alarm level) |
    | L(twin)/L(single) | 0.009 | 0.033 | 0.10 | 0.28 | 0.50 | 0.65 | 0.89 | ≈ 1 |

* **Astrometry** (h09). The same physics as phase 1, re-coded for twins with q = 0.80–1.00, flux ratios from the
  Mamajek M_G(mass) sequence, RUWE 1.076 and the DR3−DR2 proper-motion anomaly.

  | q | Excluded (L < 0.05) | Disfavoured (L < 0.3) |
  |---|---|---|
  | 0.80 | 0.32–2.4 au | — |
  | 0.85 | 0.42–2.4 au | — |
  | 0.90 and 0.93 | ~1.0–1.3 au | 0.56–2.4 au |
  | 0.95 | — | 1.3–1.8 au only |
  | ≥ 0.97 | — | — (unconstrained at all a) |

  * Twins at a > 5 au (ρ > 0.08") are disfavoured by the Gaia IPD statistics (phase 1, LR ~0.05). Beyond ~160 au
    they would be resolved.
* **Fraction of twins that survive all these tests**, with the Winters+19 separation prior: 0.22–0.24 for q ≥ 0.97,
  0.15 for q = 0.95 and 0.05–0.13 for q = 0.80–0.93 (`out/h15_posterior.json` detail).
* **What survives:** a q ≥ 0.95 twin at ~0.3–5 au (orbital period ~0.2–15 yr). Gaia can neither see nor exclude it.

## 6. Priors: metallicity given kinematics, and twin frequency (h10a, h10b, literature)

* **[Fe/H] given the kinematics** (h10a/h10b).
  * Sample: best-quality Gaia GSP-Spec dwarfs within 200 pc (all 13 first flag characters 0; 4800–6300 K; n = 8629).
    [Fe/H] = mh + fem, uncalibrated; the Recio-Blanco+23 logg term is ~0.02 dex for dwarfs.
  * The M3V's UVW_LSR is (33.1, 6.6, 47.6) km/s and v_tot 58.4 km/s; Bensby TD/D = 0.25.

  | Sample | n | median [Fe/H] | P(> +0.2) | P(> +0.3) |
  |---|---|---|---|---|
  | All dwarfs | 8629 | −0.05 | 9.2% | 3.7% |
  | \|W_LSR\| 30–65 km/s (adopted) | 763 | −0.12 | 7.3% | 3.1% |
  | W within 10 km/s and v_tot within 15 km/s of the M3V | 82 | −0.20 | 3.7% | 2.4% |

  * **The [Fe/H] a single M3V would need:**
    * **+0.25 ± 0.15** from the CMD calibration (h13: ΔMKs = −0.282 − 1.880[Fe/H], sd 0.29 mag, n = 1521; the sd
      includes label noise);
    * +0.33 from the 62 high-precision (non-Birky) calibrators (ΔMKs = −0.138 − 1.877[Fe/H], sd 0.20);
    * +0.29 in the h14 model.
  * About **4–9%** of kinematically similar dwarfs are that metal-rich (> +0.2), and **2.5–3.7%** exceed +0.3.
* **Twin frequency.**
  * Inputs:
    * M-dwarf multiplicity rate 26.8 ± 1.4% (Winters+19, 25 pc);
    * separations log-normal, peak 20 au projected (~25 au physical), σ = 1.16 dex, so P(a < 5 au) = 0.27;
      Susemiehl & Meyer 22 give a peak of 48 au and σ 0.97, so 0.16;
    * mass ratios p(q) ∝ q^0.5 (El-Badry+19, primaries 0.1–0.4 M_sun);
    * twin excess F_twin = 0.19 ± 0.05 at 50–350 au, growing toward smaller separations. Adopted 0.25, range 0.15–0.35.
  * Result: P(companion with q ≥ 0.95 and a < 5 au) = **2.4%** (3.9% for q ≥ 0.8; about 1.3–3% across the
    alternative priors).
  * After the Gaia survival factors above (RUWE, PM, RV, IPD; they remove mainly a < 0.3 au and q < 0.95 at 0.3–3 au):
    **1.9%** (q ≥ 0.95) and 2.35% (q ≥ 0.8) of M dwarfs (h15 detail).
  * Empirically (h17b), only **0.43% (0.23–0.69%)** of RUWE < 1.1, IPD-clean field M dwarfs at this colour belong to
    a 'twin' population in the CMD/colour plane. That is ~4× below the literature-based value. Possible reasons:
    * the singles' x-distribution absorbs part of the twins;
    * the twin excess may be smaller among old field M dwarfs.

    The prior-based odds in section 7 use the larger, literature value. That is conservative for the twin case.
  * Winters+19 also note that fast-moving (older) M dwarfs are somewhat less often multiple.
  * Membership of a wide triple may raise the inner-binary probability. This is tested in h15 as multiplicity × 1.5.
* **Base rate of over-luminosity.** 2.7% of RUWE < 1.1, IPD-clean field stars at this colour are as over-luminous as
  the M3V (h17; the referee found 1.4% with their ridge). Metal-rich singles and twins share that tail.

## 7. Posterior odds (h15, h17b) and consequences

### Prior-based model (h15)

* **Hypotheses.**
  * S: single star with [Fe/H] = f.
  * T: q = 0.80–1.00 companion at a, with both stars on the same sequence; binary offsets −0.738 (q = 0.8) to −0.753
    (q = 1).
* **Likelihoods:**
  * ΔMKs given f, from the h13 calibration;
  * the luminosity-free indicators as functions of f: XP shape (forward-calibrated), r(G−W1) or r(G−Ks), and the
    wide-pair likelihood (h14).
* **Priors:**
  * f from the kinematic analogues (h10b);
  * T from multiplicity × q-distribution × separation distribution;
  * T's survival from the Gaia RUWE/PM (h09), RV error (h07) and IPD.

| Run | P(twin) | odds T:S |
|---|---|---|
| CMD only (no metallicity indicator) | 6.3% (7.9% with intrinsic sd 0.15; 10.5% high-precision calibration) | 0.068 |
| CMD + XP shape | 2.5% | 0.026 |
| CMD + r(G−W1) | 1.6% | 0.016 |
| CMD + r(G−Ks) | 1.7% | 0.017 |
| CMD + wide pair | 2.0% | 0.020 |
| CMD + Gaia Net (uninformative check) | 6.3% | 0.068 |
| **Baseline: CMD + XP + r(G−W1) + wide pair** | **0.9%** | **0.0088** |
| Baseline, intrinsic ΔMKs sd 0.15 mag (label noise removed) | 0.2% | 0.0021 |
| Baseline with r(G−Ks) instead of r(G−W1) | 0.9% | 0.0093 |
| CMD + XP + wide (no photometric colour) | 1.4% | 0.014 |
| Baseline, [Fe/H] prior of all local dwarfs / of tight analogues | 0.8% / 1.0% | |
| Baseline, Susemiehl & Meyer separations; F_twin 0.15 / 0.35 | 0.6%; 0.7% / 1.1% | |
| Baseline, multiplicity × 1.5 (triple system) | 1.4% | |
| Twin-friendly: all-dwarf prior, multiplicity × 1.5, F_twin 0.35, S&M separations | 1.2% | |
| High-precision calibration of ΔMKs and r(G−Ks): CMD + XP + r(G−Ks) + wide | 2.3% | 0.024 |

* Posterior mean [Fe/H] under S is +0.16 to +0.21. Under T the components would need ~−0.13 before the indicators
  are applied.

### Population-based model (h17b)

P(twin | ΔMKs, colour), using only the field population of stars selected like the M3V: 
* (ΔMKs, r(G−Ks)): **0.5% (0.3–0.8%)**; 1.2% if Ks is 0.04 mag too bright.
* (ΔMG, r(G−W1)): **1.1% (0.5–1.6%)**; 1.9% if W1 is 0.03 mag too bright.
* Kinematically matched subsamples: ~0 (not identifiable).

Multiplying in the independent XP likelihood ratio (twin/single 0.40–0.47; h12c) and the wide-pair ratio (~0.5) lowers
these values by a further factor of ~4–5.

### Verdict

* **Posterior odds single (metal-rich) : twin ≈ 100 : 1** (range ≈ 40 : 1 to 500 : 1). **P(twin) ≈ 1% (0.2–2.5%).**
* The CMD alone gives P(twin) ≈ 6–10%. Every luminosity-independent metallicity indicator points the same way:
  * G−Ks and G−W1 colour excesses;
  * the XP spectral shape and its residual pattern;
  * the wide companions.

  Each lowers the twin probability by a factor of 2–4, and together by ~7.
* **Adopted host composition: [Fe/H] ≈ +0.2 ± 0.1.** The indicators alone give +0.15 to +0.33; the posterior including
  the kinematic prior gives +0.17 to +0.21.

### Consequences

* **Host parameters.** Use the single-star values with the [Fe/H] term (phase 1 `star_params.json`, radius_rsun_with_feh):
  * R = 0.385 ± 0.012 R_sun;
  * M = 0.370 ± 0.012 M_sun (the Mann+19 mass–M_Ks relation depends only weakly on [Fe/H]);
  * ρ = 9.2 (8.4–10.1) g/cm³ for the transit-shape/density test;
  * Teff ≈ 3330 K (Stefan–Boltzmann with the larger R).

  The twin values (2 × 0.289 R_sun, ρ 15.5 g/cm³) carry only ~1% weight.
* **Planet radius** (h15; depth 15.7 ppt on the combined light, limb-darkening factor 1.04; the ±0.8 (stat) ±1.0 (sys)
  ppt depth error adds ±0.2 R_earth):
  * single, metal-rich: **Rp = 5.16 R_earth** (5.11 with the solar-metallicity radius);
  * twin, q = 1: 5.47 R_earth;
  * twin, q = 0.9–0.95: 5.40–5.56 R_earth, depending on which star is the host.

  The twin case changes Rp by only +6–8%, so it is not the 5.5–6 R_earth range assumed earlier. The planet stays
  Neptune-sized either way.
* **Transit duration.** It cannot separate the two cases:
  * single host (ρ 9.2): T14 = 1.93 h at b = 0, so the observed 1.5–1.6 h needs b ≈ 0.6–0.7;
  * twin host (ρ 15.5): T14 ≤ 1.71–1.73 h, so it needs b ≈ 0.3–0.5.
* **False-positive analysis.**
  * The "near-equal twin" entry in the open false-positive space can be down-weighted to ~1% prior. Even if real it is
    not a false-positive channel. It is bound and equal-mass, so an eclipse of it would need the same ~3% depth on an
    M3.5V. That is still planetary (Rp/R* ≈ 0.18), and it only changes Rp by +7%.
  * A hierarchical EB in which the companion is itself an eclipsing pair is not addressed by this test. Such a
    companion must be fainter (ΔT ≥ 3.5) to avoid adding over-luminosity, and the over-luminosity is now explained by
    metallicity.
  * Fainter bound companions (ΔT 3.5–4.5 at 5–65 au, any companion at 4–5 au) are untouched by this analysis. A
    ΔK = 3 companion adds only ~0.07 mag of over-luminosity, inside the scatter. They remain open for the blend
    analysis.
  * The two catalogued background stars (3.1" and 5.5") are unaffected.
  * The metal-rich host is consistent with the planet–metallicity correlation (mild for Neptunes around M dwarfs). It
    is neither evidence for nor against the planet.
* **Activity.** If the eRASS1 source is the M3V, L_X/L_bol = 4e-4 implies Ro ≈ 0.17, i.e. P_rot ≈ 11–13 d and a
  moderately young age. That sits in tension with the hot kinematics.
  * A tidally synchronised twin (a < 0.1 au) would explain the activity, but the Gaia RV error excludes such twins
    (L ≤ 0.1).
  * The activity therefore does not favour the twin. It is consistent with a single star of intermediate age (or with
    emission from the A+B pair).

## Caveats

* **[Fe/H] labels.** Most calibrators (~2100 of 2275) have Birky+20 APOGEE-Cannon labels.
  * Those labels look compressed: the G−Ks slope is 0.25 mag/dex with Birky and 0.19 with the 62 high-precision
    stars.
  * The robust 3σ-clipped fits use the Birky-dominated set. Variants with the high-precision set raise P(twin) by about
    a factor of 2 (2.3% instead of 0.9%).
  * LAMOST (Ding+22) [Fe/H] disagrees with the Mann scale (offset −0.33, r = 0.39) and is not used for training.
  * Without the |ΔMKs| < 1.3 sanity cut, giants and cross-match failures contaminate the Birky set. The first-pass
    fits included them; all numbers above use the cut.
* **The colour tests depend on photometry at the ~0.03 mag level.**
  * r(G−Ks) relies on 2MASS Ks, which carries the 'd' (diffraction-spike) flag.
  * r(G−W1) relies on unWISE W1 inside the Mira's halo (flag 32; fracflux 0.88).
  * The two agree (+0.089 and +0.092), and VHS confirms Ks to ±0.1. The 'Ks 0.04 fainter' and 'W1 0.03 fainter'
    variants are reported in h17b.
* **XP regression** uses only stars with G < 15, because sampled XP spectra are published only for those. It
  regresses toward the mean (slope 0.64), which is handled through its forward calibration. Gaia Net was checked and
  found uninformative.
* **The twin model assumes both components lie on the same single-star sequence.** A pair of two metal-rich stars
  cannot reproduce ΔMKs = −0.76: it would be over-luminous by ~1.2 mag.
* **The prior-based odds (h15) depend on literature multiplicity statistics.** Their spread (×0.5–2) is covered by the
  variants. The population-based odds (h17b) depend on the mixture model being identifiable. In the v_tan 30–70 km/s
  subsample the twin component collapsed to 0, so treat those values as "≲ 1%", not as precise numbers.
* **The Gaia RV simulation** assumes that DR3 faint-star errors scale with the curvature of the combined CCF (Zucker
  2003). The empirical excess of large errors among over-luminous stars supports that binaries inflate the error, but
  the exact inflation is model-dependent (σ_CCF 11–20 km/s were tested).
* **Not tested here:** high-resolution spectroscopy, which would settle the question at once.
  * An SB2 at a < 0.5 au would show double lines or RV variation of 5–25 km/s.
  * [Fe/H] can be measured directly from NIR spectra (Mann+13 / Newton+14 indices).
  * At 64 pc, 0.3–5 au corresponds to 5–80 mas. 8-m speckle imaging (Zorro / 'Alopeke, ~20 mas) would resolve a
    q ≈ 1 pair at ≳ 1.3 au. Closer pairs need interferometry, aperture masking or spectroscopy.

## Scripts and outputs

All scripts are in `scripts/` and run with `/home/user/work/venv/bin/python -I`. GaiaXPy 2.1.5 was installed into the
venv with uv.

| Script | What it does | Outputs |
|---|---|---|
| `h00_xputil.py` | Helper: Gaia DataLink XP_CONTINUOUS download (DATAMODEL_GAIA, ≤50 ids per anonymous request) into GaiaXPy-ready arrays | – |
| `h01_xp_target.py` | XP continuous spectra of the M3V and of the A(+B) entry; GaiaXPy calibration (336–1020 nm, 2 nm) | `data/xp/xp_cont_target.pkl`, `out/xp_target_sampled.csv`, `out/xp_target_meta.csv`, `figs/xp_target_raw.png` |
| `h02_vizier_allcats.py` | Every VizieR catalogue with the M3V within 3" (finds Gaia Net, SHBoost, StarHorse, Verberne+24, …) | `data/vizier_allcats.txt`, `data/vizier_allcats_rows.csv` |
| `h03_field_sample.py` | Gaia archive: all M dwarfs with plx > 10 mas, plx/err > 20, 2.30 < BP−RP < 3.25; 2MASS and AllWISE best neighbours; RV statistics (n = 50,583) | `data/field_mdwarfs.csv`, `out/h03.log` |
| `h04_calibrators.py` | Literature spectroscopic [Fe/H] (Mann+15, Mann+13 wide binaries, Newton+14 (+cpm), Rojas-Ayala+12, Terrien+15, Gaidos+14, Birky+20; Ding+22 LAMOST kept separate) and CDS XMatch to Gaia DR3. Its Gaia-archive stage failed (HTTP 500/503) | `data/calibrators_raw.csv`, `data/calibrators_xmatch.csv` (reduced), `out/h04.log` |
| `h04b_calibrators_vizier.py` | Rebuilds the calibrator table from the cached Gaia XMatch plus CDS XMatch to 2MASS and AllWISE | `data/calibrators_gaia.csv` |
| `h05_xp_synphot_target.py` | XP synthetic photometry (JKC, SDSS, PS1, DECam) for the M3V and A(+B) | `out/xp_synphot_target.csv` |
| `h06a_nir_field.py` | Local (8') DENIS, 2MASS, unWISE, AllWISE, VHS DR5 and Gaia catalogues | `data/local_*.csv` |
| `h06b_nir_local.py` | DENIS−2MASS frame zero points; unWISE−AllWISE offsets; VHS saturated-star forward and inverse calibration; DENIS I vs Gaia I_C | `out/h06b_nir_local.json`, `figs/h06b_nir_local.png` |
| `h07_gaia_rv.py` | Empirical rv_error percentiles for normal, over-luminous and RUWE > 2 stars; SB2-twin CCF simulation (flagged fraction and likelihood ratio vs a) | `out/h07_gaia_rv.json`, `figs/h07_gaia_rv.png`, `out/h07.log` |
| `h08_photfeh.py` | B05, JA09, SL10, Neves+12 and Rains+21 photometric [Fe/H] for 2MASS, DENIS raw and DENIS corrected; single and twin | `out/h08_photfeh.json` |
| `h09_twin_astrometry.py` | RUWE + DR3−DR2 PM likelihood for twins (q 0.80–1.00, a 0.01–30 au) | `out/h09_twin_astrometry.json`, `figs/h09_twin_astrometry.png` |
| `h10a_gspspec_kin.py` | Gaia GSP-Spec dwarfs within 200 pc (best flags) | `data/gspspec_dwarfs.csv` |
| `h10b_feh_prior.py` | UVW and Bensby TD/D; [Fe/H] distributions of kinematic analogues; Gaussian-mixture priors | `out/h10b_feh_prior.json`, `figs/h10b_feh_prior.png` |
| `h11_xp_batch.py` | DataLink XP batch download. Abandoned: too slow (~50 spectra per 4 min); kept for the record | – |
| `h11b_xp_vizier.py` | Official sampled XP spectra (VizieR I/355/xpsample, multi-position cone queries) for calibrators, field stars, the M3V, A(+B) and the Mira | `data/xp/xps_chunk_*.npz` (1646 spectra) |
| `h12a_xp_ids.py` | Selection of XP sources (calibrators; field normal and over-luminous at the M3V's colour) | `data/xp/xp_ids.csv` |
| `h12b_xp_analysis.py` | XP peculiarity and contamination check, XP-shape [Fe/H] regression with CV, LAMOST check, field normal vs over-luminous predictions, composite fit | `out/h12b_xp.json`, `out/h12b_xp_pred_field.csv`, `figs/h12b_xp_main.png`, `figs/h12b_xp_twin_null.png` |
| `h12c_xp_population.py` | Twin fraction among over-luminous stars from XP [Fe/H]; the M3V's LR; Mira-flux fit calibrated on field stars | `out/h12c_xp_population.json`, `figs/h12c_xp_population.png` |
| `h12d_xp_composite.py` | Single vs two-star composite χ² for the M3V, normal, over-luminous and metal-rich/poor stars | `out/h12d_xp_composite.json`, `out/h12d.log` |
| `h12e_residual_corr.py` | Correlation of the M3V's XP residual with metal-rich, metal-poor and over-luminous patterns (null from normal stars) | `out/h12e_residual_corr.json` |
| `h13_colour_feh.py` | Field ridges (M_Ks, M_G, colours vs BP−RP); calibrator fits of ΔM and colour residuals vs [Fe/H] (3σ-clipped; Birky-dominated and high-precision sets); colour-only ridge regression | `out/h13_colour_feh.json`, `data/calib_sel.csv`, `figs/h13_colour_feh.png`, `out/h13.log` |
| `h14_wide_pair_feh.py` | System [Fe/H] from the A+B blend colour given separate G magnitudes | `out/h14_wide_pair_feh.json`, `figs/h14_wide_pair_feh.png`, `out/h14.log` |
| `h15_posterior.py` | Prior-based posterior odds S vs T with all variants; planet radii | `out/h15_posterior.json` |
| `h16_lit_xp_feh.py` | Calibration of Gaia Net and SHBoost [Fe/H] on M dwarfs; field normal vs over-luminous | `out/h16_lit_xp_feh.json`, `figs/h16_lit_xp_feh.png`, `out/h16.log` |
| `h17_population_mixture.py` | Free 2D Gaussian-mixture population model. Not identifiable (all stars absorbed as 'singles'); kept as a negative result | `out/h17_population_mixture.json`, `figs/h17_population_mixture.png`, `out/h17.log` |
| `h17b_structured_mixture.py` | Structured, identifiable population mixture (singles on a metallicity line, twins shifted by −0.753 mag) in the (ΔMKs, r(G−Ks)) and (ΔMG, r(G−W1)) planes, with variants | `out/h17b_structured_mixture.json`, `figs/h17b_structured_mixture.png`, `out/h17b.log` |

Other files:
* `data/lit/`: arXiv texts used for coefficients and priors:
  * Johnson & Apps 09, Schlaufman & Laughlin 10, Neves+12, Rains+21;
  * Winters+19, Susemiehl & Meyer 22, El-Badry+19;
  * Huson+25 (Gaia Net), Khalatyan+24 (SHBoost).

  PDFs were deleted after text extraction.
* `out/h03.log` … `out/h17b.log`: run logs.
