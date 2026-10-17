# Figure captions: TIC 214311908 transit candidate

These figures come from an AI-run analysis of public TESS data. No professional astronomer has reviewed them. Every number is taken from the analysis folders named in each caption; scripts are in `../figscripts/` (`run_all.sh` rebuilds all seven figures in about 20 s). INV and REF are two separate deblending photometries of the same TESS pixels (INV: global pixel regression with a Mira image and an M3V PRF; REF: local per-pixel baselines plus PRF fits). Depths are in ppt of the M3V's TESS flux (979 e-/s). Times are target-barycentred BJD_TDB. The adopted ephemeris is T0 = 2460840.9150 ± 0.0013 BJD_TDB, P = 5.691506 ± 0.000011 d (referee3_clean).

## fold.png
**(a)** The transit folded on the adopted ephemeris, in 15-min bins, for three photometries (offset vertically):
* the cleaned INV curve (18 transits);
* the cleaned REF curve (21 transits);
* a raw-pixel 3x3-box aperture with no scene model (21 transits).

The black curve is the adopted model: depth 15.8 ppt, Rp/R* 0.119, b 0.52, a/R* 25.3, T14 1.70 h. All three photometries show the same flat-bottomed 1.7-h transit.

**(b)** Residuals.

**(c)** The phase-0.5 window, for INV and the raw box (REF is not plotted). No occultation is seen: the fitted occultation-shaped depths are +0.3 ± 0.8 ppt (INV), +0.3 ± 0.9 ppt (REF, not shown) and −1.7 ± 1.1 ppt (raw box). The dashed curve shows a 3-ppt occultation, about the 3σ limit.

The INV counts here (18 transits) are those of clean/; per_transit.png uses the referee's 19 INV transit times (referee3_clean q05).

Each ±0.30-d window has its own baseline removed: linear for INV and REF, quadratic plus pointing for the raw box. The baseline is fitted jointly with a free-amplitude transit template.

The raw-box curve was recomputed with the referee3_clean q02 method (`f0_rawap.py`). Its per-window flux losses match `q02_bigap_per_transit.csv` exactly.

Sources:
* `clean/lc_M3V_INV_clean.csv`, `clean/lc_M3V_REF_clean.csv`
* `referee3_clean/out/q02_bigap.json`, `referee3_clean/out/q08_global.json`, `referee3_clean/out/q12_secondary_ap.json`
* `clean/out/phase_scan.json`

## localisation.png
**(a)** The TESS Sector 39 median image near the target, with Gaia DR3 stars at the sector epoch. Labelled stars:
* the M3V target;
* the Mira, 17.8″ away (with its measured apparent image, shifted 7–8″ in S39 by lateral colour; 5–7.5″ across sectors);
* the 12″ pair A/B;
* the 3.1″, 5.5″ and 12″ faint stars.

All of these lie within about one 21″ pixel.

**(b)** Per-pixel dip S/N across the whole 21x21 cutout in S39. It is model-free: raw pixels, a common-mode correction and local baselines (referee2_search r07 method, recomputed for this sector). Only the pixels at the M3V/Mira position dim. The 3x3 box loses 13.1 ± 2.9 e-/s in S39 and 13.9 ± 1.4 e-/s over all six sectors. No pixel beyond 2.5 px dims in any sector.

**(c)** The combined six-sector localisation map (referee_loc). Contours show the injection-calibrated 68/95/99.7% regions, which reach 3.2″, 4.6″ and 8.1″ from the M3V. The 95% region includes the 3.1″ star. That star is disfavoured by its own hypothesis test (Δχ² +36 at baseline, 0/1000 injections; ~2.5–3.5σ across variants; panel d) and excluded by physics: it would need a ~50–54% flat-bottomed eclipse.

Panel b only confines the dip to the pixels holding the M3V and the Mira (0.9 px apart): the 3x3 box contains both. The separation between the M3V and the Mira comes from the injection-calibrated image fits of panels c and d.

**(d)** χ²(source on that star) − χ²(source on the M3V), observed, compared with the 2.5–97.5% range from 1,000 injections:
* with the M3V as the source (grey);
* with that star as the source (orange).

None of the injections on pair A, pair B, the Mira or the 3.1″ star reaches the observed values.

Sources:
* `referee_loc/out/inj_base.npz`, `locregion_base.json`, `fit_base.json`, `wcs_mira.json`, `data/gaia_400.csv`
* `referee2_search/out/r07_pixdepth.json`
* the S39 TESScut cube

## per_transit.png
**(a)** Per-transit mid-times minus the adopted ephemeris, for the cleaned INV and REF curves.
* Times and errors come from referee3_clean q05. The errors come from light-curve-level injections, one value per sector and pipeline. The times reproduce clean/ c06 to 0.1 min rms.
* The grey band is the adopted ephemeris ±1σ. The dashed line is the clean/ combined ephemeris.
* O−C rms is 8.3 min (INV) and 8.8 min (REF). The largest deviation is 2.0σ. Each pipeline's own linear fit gives χ² 8.9/17 (INV) and 12.6/19 (REF). There are no TTVs.

**(b)** Per-sector depths, all measured identically by referee3_clean q14 (fixed shape, linear baseline, fake-epoch errors), for:
* the two PRF photometries;
* two raw-pixel apertures.

The bands show the adopted depth, 15.8 ± 0.9 (stat) +1.5/−1.1 (sys) ppt.

The PRF photometries show a late-minus-early (S103+S104 vs S39+S66) excess: INV +4.1 ± 2.6 ppt; REF +6.1 ± 2.1 ppt with the scattered-light event S103 E58, or +3.4 ± 2.1 ppt without it (hollow square).

The PSF-independent raw apertures are consistent with a constant depth:
* 3x3 box: −1.0 ± 2.4 ppt, sector χ² 1.5/5;
* r ≤ 2.5 px: −4.9 ± 3.2 ppt, sector χ² 3.9/5.

So the late/early contrast is most likely a property of the PRF/scene-model photometry: paired with the apertures at the same epochs, INV is deeper by +6.6 to +7.4 ± 3.1–3.3 ppt late-minus-early (~2σ; split-independent test p 0.06–0.45), and the apertures alone cannot exclude a 3–4 ppt change. The split was chosen after looking at the data (look-elsewhere p 0.31–0.81).

Sources:
* `referee3_clean/out/q05_per_transit.csv`, `q05_timing.json`, `q14_contrast_all.json`, `q04_paired.json`, `q15_omnibus.json`, `q07_lee.json`
* `clean/out/ephemeris_final.json`

## false_positives.png
**(a)** Scenario posteriors from the sign-fixed Bayesian false-positive calculation, for the REF and INV light curves. Settings: Kipping planet eccentricity prior, short-period brown-dwarf eccentricity law, f_BD = 1.3e-4.
* FPP = P(not a planet on the M3V or its twin) = 0.0013 (REF) and 0.0031 (INV). With the Sagear & Ballard (2023) M-dwarf eccentricity prior, which referee3_bd recommends, 0.0016 and 0.0033.
* P(planet on the M3V or on an unresolved companion of it) = 0.9995 (0.9986–0.9996 over variants); including a planet on the bound pair star A, 0.9999.
* "Twin" means any companion within 0.75 mag of the M3V (a planet there would be 5.2–6.3 R_E).
* The largest false-positive terms are:
  * a planet on an unresolved ~0.2 Msun companion: 8.6e-4 (REF), 2.6e-3 (INV). The FPP code calls this PTP; the figure labels it "STP", its TRICERATOPS name, because TRICERATOPS uses PTP for a planet on the primary;
  * NTP on pair A: 3.7e-4 (REF), 4.2e-4 (INV).

Brown-dwarf channels are drawn twice:
* outline: the sign-fixed grid;
* filled: after referee3_bd's two calibrations, which together give a 0.6e-6 (REF) to 6.9e-6 (INV) total. The calibrations are:
  * the Gaia RVS factor (median 0.49, applied to BDs on the M3V and on a binary primary);
  * the injection-calibrated shape likelihood ratio (0.076 REF, 0.29 INV, applied to all BD channels).

Over the 64 variants the FPP ranges 0.0007–0.0037. With circular planet orbits it is 0.0007 (REF) and 0.0019 (INV). Over 1,500 joint draws of all priors, re-run after the sign fix in the final review (panel_fpp), the median is 0.0022 (5–95%: 0.0002–0.012).

Earlier FPP values are superseded: 0.011 (fpp/) used a mis-calibrated over-luminosity likelihood, and 0.003–0.007 (referee2_fpp) still carried the photo-eccentric sign error that bd/ found and fixed in both codes. TRICERATOPS (0.002–0.04) is unstable here and is not used.

**(b)** Expected FPP after follow-up, recomputed with the adopted composition by `../scripts/w05_fpp_ladder.py` (`../fpp_ladder_adopted.json`). It keeps the follow-up referee's rules and factors (speckle removes a fraction of the companion and bound-EB terms; a two-band transit applies the colour Bayes factor and removes pair A and the white dwarf) and swaps in the adopted brown-dwarf inputs (f_BD 1.3e-4, short-period eccentricity law; with and without the Gaia RV factor and the calibrated shape likelihood). The ranges cover those three brown-dwarf versions, the three Zorro contrast curves (speckle-only step; the combined steps use the conservative curve), the colour models and 5 or 10% depth-ratio precision. Open markers: the follow-up referee's published ladder, which used f_BD = 1e-3 and the uncalibrated shape grid after "now" (drawn only where it differs by more than 10%).

| Stage | REF | INV | Published ladder, REF / INV |
|---|---|---|---|
| Now | 1.3e-3 | 3.1e-3 | 1.4e-3 / 3.8e-3 |
| After an 8-m 832-nm speckle null | 4.9–6.4e-4 | 5.8e-4–1.0e-3 | 6.0–7.3e-4 / 1.0–1.5e-3 |
| After a two-band transit only | 3.2e-5–2.1e-4 | 4.0e-5–5.1e-4 | 1.6–3.2e-4 / 5.5e-4–1.0e-3 |
| After speckle + two-band transit | 3.1–6.1e-5 | 3.6e-5–1.2e-4 | 1.4–1.6e-4 / 4.6–5.2e-4 |
| … + brown-dwarf channel closed (by the ground transit shape or by RV) | 3.0–4.9e-5 | 3.0–8.3e-5 | 3.0–4.9e-5 / 3.0–8.6e-5 |

The floor of ~3e-5 is the background-star term, carried over from phase 2 and not reduced by any step.

The candidate is not statistically validated by community standards. There is no on-target ground-based detection and no real high-resolution imaging.

Sources:
* `referee3_bd/out/s02_fpp_fixed.csv`, `s12_summary.json`
* `referee3_followup/out/r05_hri_ladder.json`, `../fpp_ladder_adopted.json`
* `panel_fpp/out/p05_summary.json`, `p02_ranges.json`

## host.png
**(a)** Absolute M_Ks against BP−RP for 4,230 Gaia-clean nearby M dwarfs (referee_star comparison sample). The blue line is the single-star ridge (density peak). The dashed line is a ridge star's equal-mass twin. The target lies 0.84 mag above this ridge (0.75 mag with the local ridge of referee2_host), close to the twin line.

**(b)** The G−Ks colour residual at fixed BP−RP against the M_Ks offset, for the 2,126 clean stars within 0.1 mag of the target's colour. Single stars follow a metallicity trend (blue: running median and r02 fit), with scatter 0.044 mag. The dashed orange track shows unresolved binaries made of a normal primary and a companion of mass ratio q. It returns to r ≈ 0 for a near-equal twin.

The target (r = +0.097) sits on the metal-rich single-star trend (+0.070 ± 0.044 expected) and about 2σ from the twin locus. Hence P(twin) ≈ 0.3–1%, or 0.2–5% allowing photometric systematics (referee2_host). The colour test cannot exclude a companion that supplies up to ~1/3 of the Ks light.

The ridge, residual definitions and sample cuts reproduce `referee2_fpp/scripts/r02_cmd_colour.py`: slope 1.5309, target x = −0.838, r = +0.097.

Sources:
* `referee_star/data/mdw_cmd.csv`
* `referee2_fpp/out/r02_cmd_colour.json`
* `referee2_host/README.md`

## scenarios.png
**(a)** The cleaned REF light curve (12-min bins; per-window baselines of the planet fit) with four best-fit models, recomputed from saved parameters:
* a planet: referee3_bd s03 maximum-likelihood fit with free density (b 0.24). The adopted posterior median with the stellar-density prior is b ≈ 0.52; the shapes are nearly identical.
* a white dwarf totally occulted by the M3V: referee2_shape, 0.38 Msun, 3-min ingress;
* a brown dwarf grazing the M3V: referee3_bd s03, k 0.19, b 0.96;
* a V-shaped grazing body: referee2_shape.

**(b)** Residuals from the planet model, with each alternative model minus the planet. The planet's flat bottom and 11–14-min ingress fit best.

Δχ² against the planet:

| Model | Δχ² | Data | Notes |
|---|---|---|---|
| Grazing brown dwarf | +9.2 | cleaned REF, calibrated noise | Injection-calibrated likelihood ratio 0.004 (REF) / 0.03 (INV); needs e ≥ 0.5 |
| V shape | +5.7 | phase-1 REF | +6.7 on the cleaned curve; a V-shaped truth reproduces this in 0/25 REF injections |
| White dwarf | +1.5 | phase-1 REF | Light curve only; see below |

The white dwarf case is excluded by other data, not by the light-curve shape:
* light curve + secondary + Gaia RUWE give Δχ² ≥ 9–11;
* adding Gaia XP gives 10–12;
* the Gaia RVS error adds ≥ 10.

Sources:
* `referee3_bd/out/s03_fit_real.json`, `s11_calibration.json`
* `referee2_shape/out/scen_ref_fake.json`, `wd_combined.json`
* `clean/lc_M3V_REF_clean.csv`

## followup_2027.png
Transit windows in the 2027 and 2028 seasons (March–September), from the write-up's per-site tables (`../observability/obs_*.csv`, adopted ephemeris). A window is T14 + ±3σ timing + 1.5-h baselines, 5.0–5.2 h in total, with the target at airmass < 2 (altitude > 30°) throughout.

Symbols:
* filled circle: the whole window is in astronomical darkness (Sun < −18°);
* triangle: the whole window is dark only with Sun < −12°;
* tick: partial window (at least 50% of it dark).

Marker colour gives the Moon's illumination at mid-transit. The Chile row stands for CTIO/Cerro Pachón, Las Campanas and La Silla, which have identical window lists; Paranal differs at two epochs.

The first full windows come in April 2027. The 2027 season offers 10 distinct epochs with a full dark window from Chile, SAAO or Siding Spring, and 14 including Dome C. 2028 offers 9 (14).

Highlighted new-Moon events (mid-transit to the nearest minute):
* E126, 2027-06-01 12:56 UTC, Siding Spring (MuSCAT4);
* E127, 2027-06-07 05:32 UTC, Chile;
* E137, 2027-08-03 03:32 UTC, Chile.

1σ timing errors are 2.8–3.7 min in 2027 and 3.7–4.6 min in 2028; plan ±3σ margins (column `sigma3_min`).

Sources:
* `../observability/obs_*.csv`, `../scripts/w02_checks.json`
* `../ephemeris_adopted.json`
