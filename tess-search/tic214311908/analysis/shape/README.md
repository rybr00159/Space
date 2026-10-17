# Transit model, stellar-density test and blend physics for TIC 214311908

Folder: `tic214311908/analysis/shape/`
(AI-run analysis, 2026-10-17. Nothing was submitted, uploaded or posted anywhere, and no one was contacted.)

Inputs are the phase-1 deliverables: `../phot/lc_M3V_all.csv` (INV, the investigator's model-B deblended light curve,
rescaled from 935 to 979 e-/s), `../referee_phot/lc_M3V_referee_windows.csv` (REF, the referee's independent
local-window photometry), `../star/star_params.json` and the referee results. All significances quoted as dchi2
use per-point errors that include the red-noise factor beta (see Data). Where it mattered, they were calibrated
with injections (s13) or with the referee's fake-window errors.

## Bottom line

1. **Transit model (both photometries agree).** Fitted freely (rho*, b and LD with a Claret prior), the two
   photometries give the same answer:
   * Rp/R* = 0.119 +- 0.006 in both;
   * T14 = 1.69 (+0.17/-0.12) h (INV) and 1.65 (+0.13/-0.09) h (REF);
   * T23/T14 = 0.73 (INV) and 0.75 (REF);
   * b = 0.49 (+0.27/-0.35) and 0.40 (+0.27/-0.27);
   * rho_circ = 10.3 (+5.4/-5.9) and 12.3 (+4.0/-5.8) g/cm3.
   All INV-REF differences are below 0.8 sigma (T0 differs by 1.5 min).

   With the single-star density prior: **Rp = 4.94 +- 0.26 R_earth, b = 0.53 +- 0.10, a/R* = 25.3 +- 0.8,
   i = 88.81 +- 0.24 deg, T14 = 1.70 +- 0.07 h, Teq = 471 +- 10 K**. For an equal-twin host:
   **Rp = 5.23 +- 0.25 R_earth, b = 0.20 +- 0.14, a/R* = 29.7 +- 0.9, Teq = 420 +- 9 K.** The limb darkening is
   not constrained by the data; its posterior equals the prior.
2. **The stellar-density (photo-eccentric) test is not informative.**
   * rho_circ is consistent with both the single-star value (9.4) and the twin value (15.5). Nested-sampling
     evidences: ln B(single/twin) = -0.43 (INV) and -0.53 (REF), i.e. no preference.
   * Eccentricity is unconstrained: e < 0.74 at 95% for a uniform prior, and P(e < 0.1) = 0.20 against a prior
     of 0.11.
   * Injections confirm it: an injected rho = 9.5 planet is recovered at 4.9-15.9 g/cm3.
3. **The maximum-dilution test from shape alone is weak.**
   * If the occulter is smaller than the eclipsed source (transit geometry), dchi2 reaches 9 only at
     f_src ~ 0.031-0.04, i.e. dT ~ 3.5-3.8 mag.
   * If the occulter is larger (occultation geometry), partial occultations fit within dchi2 +1.8 to +7.8 for any
     f_src >= 0.018. A box-shaped total occultation at f_src = 0.016 (dT = 4.5) fits within +0.7/+2.9.
   * So the light-curve shape alone excludes neither the 3.1" star nor companions with dT 3.5-4.5. They are
     excluded by Kepler's law combined with stellar physics, photometry, the secondary eclipse and astrometry
     (items 4-5).
4. **Explicit false-positive (FP) fits.** Each scenario's dchi2 is relative to the planet, and includes the
   non-light-curve terms (secondary eclipse, photometry, parallax, astrometry) that the scenario must satisfy.
   * **White dwarf (WD) totally occulted by the M3V: the most dangerous scenario.** On the light curve alone it
     fits as well as the planet: dchi2_LC = +0.7 (INV) and +2.9 (REF); injections give a mean of +0.8 for planets
     and -2.9 for WDs.
     * It is excluded by two other measurements together. For M_WD <= 0.35 Msun the WD is large, so its
       transit across the M3V at phase 0.5 would be 1.6-3.3 ppt deep, which is not seen. For M_WD >= 0.3 Msun
       the 5.69-d photocentre orbit (0.35-0.8 mas) would raise RUWE above the observed 1.076.
     * Combined minimum: dchi2 = 13.2 (sigma_CCD 0.35 mas, at 0.3 Msun) or 15.2 (sigma_CCD 0.22 mas). That is a
       likelihood ratio of ~700:1 to 2000:1 against, before any prior, and the models have the same number of
       free parameters.
     * The Gaia RV error of 2.2 km/s sits at the 38th percentile, against K1 = 35-106 km/s for such a WD; this
       argument is qualitative.
     * A WD around one twin component gives dchi2 +11 to +18 (the transit is too short); a WD occulted by a
       late-M companion with dT 3.5-4.5 gives dchi2 +114 (T14 = 0.8 h).
   * **Bound hierarchical EBs.** These use the main-sequence Ks-band luminosity budget, with masses and radii
     from Mann+15/19.
     * Best: two 0.11 Msun M dwarfs at P = 11.38 d, diluted by a 0.34 Msun primary: dchi2 = +27.5 (INV) and
       +33.7 (REF). The eclipse is V-shaped.
     * A brown dwarf grazing the M3V or a twin component: dchi2 >= +84 (REF), because the grazing chord at the
       stellar density is too short.
     * An EB that supplies all the light, with no third star: dchi2 >= +46 (0.28 + 0.26 Msun, grazing, 11.38 d).
   * **Background star 3.1" NE (Gaia DR3 ...260544).** Its photometry is a single K3-4 dwarf at 1.6 kpc:
     chi2 = 1.9 for g, r, i, z and the parallax.
     * Best physical EB: a 0.3 Msun pre-ELM white dwarf (R 0.18, 8000 K) totally or partially occulted by a K7
       dwarf, dchi2 = +19.6 (LC +10.1).
     * Main-sequence pairs: +66 (P) and +93 (2P). WD: +58. Hot subdwarf: +84. Brown dwarf: +143.
     * The pixel localisation (phase 1) disfavours this star independently, at about 3-4 sigma.
   * **Background star 5.5" (...262208).** Its photometry is a single G0 dwarf at 3.8 kpc. Best physical EB: a
     pre-ELM white dwarf grazing behind an M dwarf, dchi2 = +27.4. sdB: +40. WD: +55. Main-sequence pairs: +76 and
     +97. Pixels (phase 1): p ~ 0.001.
5. **Secondary eclipse and phase curve.**
   * Phase 0.5: -0.16 +- 0.79 ppt (INV, empirical error) and +0.29 +- 0.77 ppt (REF). The 3-sigma upper limit is
     2.2 ppt, and 4.0 ppt at any phase.
   * The occulter must therefore have a TESS-band brightness temperature below 2460 K. This excludes any star of
     spectral type M8V (2570 K) or earlier as the occulter.
   * The geometric-albedo limit (A_g < 77) and the expected planetary signals (< 30 ppm) are meaningless at this
     precision.
   * Ellipsoidal and beaming signals of any companion up to 0.4 Msun are below 1 ppt, while the photometry
     limits P and P/2 sinusoids only to 4.2 and 2.1 ppt (95%). So the phase curve excludes nothing.
   * Brown dwarfs and stars as the transiting body are excluded by size and duration instead:
     * with the density prior, k <= 0.133 (3 sigma);
     * the BD minimum radius needs k >= 0.197;
     * forcing k = 0.20 costs dchi2 = +45 (INV k-profile).
6. **Per sector, the transit is consistent.**
   * Depth chi2 = 3.4/5 (REF), 8.4/5 (INV) and 7.1/5 (INV model A).
   * Duration chi2 = 2.2/5 (REF), 7.7/5 (INV) and 1.5/5 (INV-A).
   * Depths rise from 13-15 ppt in S39/S66 to 17-20 ppt in S103/S104 in both photometries (about 1.5-2 sigma).
   * The INV odd-even difference (+5.0 +- 2.1 ppt) is again absent in REF (+1.3 +- 1.9 ppt).

**Evidence against the planet, or still open:**
* The shape does not discriminate a planet from a box-shaped occultation of a hot compact star. That class is
  excluded only by the secondary eclipse plus Gaia astrometry, at a combined dchi2 = 13-15. The astrometric
  likelihood comes from a simplified RUWE simulation, and for M_WD >= 0.45 it rests on 0 of 400 simulations
  passing.
* "Flat-bottomed" is only moderately supported by physical models. A V-shaped grazing transit is worse by
  dchi2 = 3.0 (INV) and 5.4 (REF); injected planets give 5.4 +- 3.2. In the free fit, grazing solutions hold 4.9%
  (INV) and 2.2% (REF) of the posterior. The phase-1 dchi2 of 18-38 compared trapezoids, not physical models.
* The density test cannot distinguish single from twin hosts, nor constrain e.
* The light curve allows any system in which a hot compact star (WD/ELM, ~1.6-3% of the aperture TESS flux) is
  totally occulted by an M dwarf of R ~ 0.4 Rsun (s12 `_lconly`: dchi2 0.0). Bound versions fail on luminosity
  (such an occulter would outshine the system), duration (twin) or secondary + RUWE (the M3V itself). Catalogued
  background versions fail on photometry.
  * An **uncatalogued** background M dwarf + WD within ~1.5" with T <~ 17.5 is not excluded by anything here.
    Its prior is small (0.0037 stars/arcsec^2 brighter than T = 17.47, times the fraction that are eclipsing
    M + WD pairs at P = 5.7 d), but it belongs in the false-positive probability calculation.
* Phase scans contain excursions of 3.3-3.6 sigma that have opposite signs in the two photometries (phase 0.57).
  At that level the noise is not Gaussian. A 1.3-1.9 ppt dip at phase 0.490-0.495 appears in all three
  photometries of the same pixels (1.8-2.4 sigma, global p = 0.16).

## 1. Data, noise and model (s00-s03)

**Data.**
* INV: quality_ok == 1, flare_flag == 0, transit windows +-0.30 d (3696 points, 22 windows). Fractional flux is
  rescaled to the 979 e-/s (PRF-unit) normalisation, f = 1 + (flux_norm - 1) x 935/979.
* REF: all 4663 points (21 windows), already 1 + dflux/979.
* Both data sets get a per-window constant offset, profiled analytically.

**Noise** (`out/noise.json`).
* Per-sector white scale (OOT robust std / median error) is 0.95-1.07.
* Red-noise factor beta:
  * INV: from the time-averaging method on the full out-of-transit series, giving 1.00, 1.08, 1.24, 1.21, 1.27, 1.13
    for S12-S104. Phase 1 needed a ~1.3x inflation per transit, consistent with these.
  * REF: the in-window beta is biased low, because the per-pixel quadratic baselines absorb slow noise. The
    referee's fake-window calibration (err_fake/err_white) is used instead: 1.11, 1.47, 1.36, 2.01, 1.02, 1.58.

**Model** (`scripts/shapelib.py`).
* batman with quadratic limb darkening and exposure integration, using 31/11/5 sub-samples for the
  1800/600/200-s cadences.
* batman was validated against brute-force 2-D integration for radius ratios k = 0.1-5: maximum error 3e-5
  (`out/batman_check.json`). It is used for k > 1 (occultations of a small star by a larger body) as well.

**Limb darkening** (`out/ld_claret17.json`).
* Claret (2017) PHOENIX TESS tables (VizieR J/A+A/600/A30) at 3350 K and logg 4.84: u1 = 0.125, u2 = 0.533.
* The spread over LSM/PCM, quasi-spherical/r-method, Teff +-2 sigma and logg 4.5-5.0 is u1 0.11-0.18 and
  u2 0.43-0.61.
* Prior: N(0.125, 0.1) and N(0.533, 0.1). For the twin host: 0.125 and 0.541.

**Sampler.** dynesty, 500 live points, rwalk, dlogz 0.05 (s03). The free-LD (Kipping q1, q2) and eccentric runs
were started but stopped to free the shared CPUs. Eccentricity is handled in s04 by importance sampling instead.

| Fit (s03) | Data | Rp/R* | b | rho* (g/cm3) | T14 (h) | T23/T14 | ln Z | chi2_min / N |
|---|---|---|---|---|---|---|---|---|
| free (rho_circ) | INV | 0.119 +0.008/-0.005 | 0.49 +0.27/-0.35 | 10.3 +5.4/-5.9 | 1.69 +0.17/-0.12 | 0.73 +0.06/-0.20 | -1100.37 | 2177.3 / 3696 |
| free (rho_circ) | REF | 0.119 +0.005/-0.005 | 0.40 +0.27/-0.27 | 12.3 +4.0/-5.8 | 1.65 +0.13/-0.09 | 0.75 +0.03/-0.12 | -1274.04 | 2523.1 / 4663 |
| single (rho prior 9.41 +- 0.88) | INV | 0.118 +- 0.005 | 0.52 +0.08/-0.11 | 9.5 +- 0.9 | 1.71 +- 0.07 | 0.72 +- 0.04 | -1099.05 | 2177.4 |
| single | REF | 0.120 +- 0.005 | 0.54 +0.07/-0.08 | 9.5 +- 0.9 | 1.69 +- 0.06 | 0.71 +- 0.03 | -1272.55 | 2523.2 |
| twin (rho 15.5 +- 1.45, dilution 2) | INV | 0.165 +- 0.006 (undiluted) | 0.20 +0.15/-0.13 | 15.2 +- 1.3 | 1.67 +- 0.05 | 0.71 +0.01/-0.03 | -1098.62 | 2177.3 |
| twin | REF | 0.168 +- 0.006 | 0.20 +- 0.14 | 15.4 +- 1.3 | 1.67 +- 0.05 | 0.70 +0.01/-0.03 | -1272.02 | 2523.1 |

**Ephemeris (free fits).**
* INV: T0 = 3840.91496 +0.00140/-0.00149 BTJD_TDB, P = 5.6915098 +- 0.0000092 d.
* REF: T0 = 3840.91603 +- 0.00129, P = 5.6915057 +- 0.0000079 d.
* The two T0 differ by 1.5 min (0.6 sigma) and the two P by 0.3 sigma. Both agree with phase 1.

**Photometry comparison.** INV minus REF:
* Rp/R*: 0.000 (0.0 sigma).
* b: +0.09 (0.2 sigma).
* rho_circ: -2.0 g/cm3 (0.3 sigma).
* T14: +0.04 h (0.3 sigma).

So the physical model removes the T14 difference of the trapezoid fits (phase 1: 1.83 h against 1.57 h). The
free-fit posteriors are broad because b, rho and LD trade off.

**Limb darkening and tails.**
* The u1, u2 posteriors equal the prior: the data do not constrain limb darkening.
* Grazing solutions (b > 1 - k) hold 4.9% (INV) and 2.2% (REF) of the free posteriors. They extend to k = 0.40
  and rho_circ ~ 1-5 g/cm3. They disappear once the single-star density prior is applied: 3-sigma b = 0.04-0.72,
  k <= 0.133.

## 2. Stellar-density (photo-eccentric) test (s04, `out/photoecc.json`)

**rho_circ** (free fits):
* INV: median 10.3, 68% 4.5-15.7, 95% 1.7-19.1 g/cm3.
* REF: median 12.3, 68% 6.5-16.4, 95% 2.8-19.5 g/cm3.
* Single-star value 9.4 (8.6-10.3); twin value 15.5. Both lie inside the 68% interval of both photometries.

**Host comparison.**
* Nested-sampling evidences, including the twin's dilution and the larger undiluted k: ln B(single/twin) = -0.43
  (INV) and -0.53 (REF), so B = 0.65 and 0.59. There is no preference.
* The rho-only marginal-likelihood ratio gives the same, 0.68 and 0.50.

**Eccentricities implied** (uniform e < 0.9, uniform omega, Dawson & Johnson 2012):

| Host | Data | median e | 68% | 95% upper | P(e < 0.1) | e_min at median rho_circ |
|---|---|---|---|---|---|---|
| single | INV | 0.24 | 0.08-0.53 | 0.74 | 0.20 | 0.03 |
| single | REF | 0.23 | 0.08-0.53 | 0.74 | 0.20 | 0.09 |
| twin | INV | 0.21 | 0.05-0.54 | 0.74 | 0.29 | 0.13 |
| twin | REF | 0.18 | 0.04-0.50 | 0.71 | 0.33 | 0.08 |

* The prior gives P(e < 0.1) = 0.11, so the data mildly favour low e for both hosts.
* A circular orbit is consistent with both hosts; no eccentricity is required.
* The single host at e = 0 sits at the 44th (INV) and 30th (REF) percentile of the rho_circ posterior. The twin
  sits at the 83rd and 78th percentile.

**Injection check** (s13, 8 fake epochs; planet injected with rho = 9.47, b = 0.52).
* Recovered rho = 10.7 +- 4.9 (4.9-15.9) and b = 0.40 +- 0.34. This matches the posterior widths, so the
  density test is calibrated and genuinely weak.
* A box-shaped WD occultation injected into the same epochs is fitted by the planet model with rho = 17.2 +- 1.2
  and b ~ 0.1.

## 3. Maximum dilution from the shape (s05, `out/dilution_profile_{inv,ref}.json`, `figs/dilution_profile.png`)

The model is any opaque occulter of radius ratio k in front of a source that supplies a fraction f_src of the
aperture flux. Everything else is free: b, a/R_src and any LD (uniform q1, q2). The occulter's own light only
adds dilution, so this family contains every spherical blend at the primary eclipse. Minimum chi2 over a grid of
f_src (multi-start Powell, with a DE check):

| f_src | dT | INV transit (k<1): dchi2, k, T14 | INV occultation (k>1): dchi2 | REF transit (k<1): dchi2, k, T14 | REF occultation (k>1): dchi2 |
|---|---|---|---|---|---|
| 1.0000 | -0.00 | -0.4, 0.13, 1.81 h | +4.9 | -0.2, 0.12, 1.60 h | +7.8 |
| 0.5000 | 0.75 | -0.2, 0.18, 1.78 h | +5.1 | -0.2, 0.17, 1.61 h | +7.8 |
| 0.3000 | 1.31 | -0.3, 0.21, 1.83 h | +4.8 | +0.1, 0.23, 1.67 h | +7.6 |
| 0.2000 | 1.75 | -0.5, 0.28, 1.86 h | +4.8 | +0.7, 0.28, 1.74 h | +7.5 |
| 0.1300 | 2.22 | -0.4, 0.35, 1.87 h | +4.7 | +1.8, 0.34, 1.88 h | +7.7 |
| 0.0900 | 2.61 | -0.1, 0.42, 1.92 h | +4.7 | +3.0, 0.42, 1.95 h | +7.4 |
| 0.0600 | 3.05 | +1.8, 0.52, 2.02 h | +4.6 | +5.1, 0.52, 2.07 h | +7.4 |
| 0.0450 | 3.37 | +4.3, 0.60, 2.14 h | +4.2 | +7.4, 0.61, 2.14 h | +7.2 |
| 0.0350 | 3.64 | +7.5, 0.69, 2.22 h | +4.0 | +10.8, 0.70, 2.23 h | +7.0 |
| 0.0310 | 3.77 | +9.1, 1.00, 2.49 h | +3.8 | +12.4, 1.00, 2.29 h | +6.8 |
| 0.0225 | 4.12 | +9.6, 1.00, 2.52 h | +3.4 | +13.0, 1.00, 2.52 h | +6.4 |
| 0.0180 | 4.36 | +10.4, 1.00, 2.58 h | +1.8 | +14.0, 1.00, 2.59 h | +4.0 |

* The best occultation-regime solutions on this grid are grazing partial occultations by a much larger body
  (T14 ~ 2 h). The box-shaped total occultation of a source with f_src = 0.0157 (dT = 4.5) is the WD scenario of
  section 4a: dchi2_LC = +0.7 (INV) and +2.9 (REF).

* With k < 1 (planet-like or annular), the 3-sigma (dchi2 = 9) limit on the undiluted depth is reached at
  f_src ~ 0.031-0.04: dT_max ~ 3.5-3.8, with an undiluted depth of ~40-50% and k -> 1 (V-shaped).
* With k > 1, grazing partial occultations of a small source by a much larger body fit within dchi2 +1.8 to +7.8
  for every f_src from 1 down to 0.018. The box-shaped total occultation at f_src ~ 0.016 fits within +0.7/+2.9.
* **The shape therefore does not exclude a source as faint as dT ~ 4.5, the limit set by a 100% deep eclipse.** This covers the 3.1" star
  (dT = 3.77, needs 50-60%), the 5.5" star (dT = 4.08, needs 68-76%) and bound companions with dT 3.5-4.5. What
  remains is physics: which real binary can make a 50-100% deep, 1.5-1.7 h, flat eclipse at P = 5.69 d with no
  secondary? Item 4 answers this.
* With the rho prior on the M3V, the k-profile gives k = 0.10 (+1.4), 0.14 (+12.2), 0.17 (+38.9), 0.20 (+45.3)
  and 0.23 (+47.6). The twin host (rho 15.5, dilution 2) is best at k = 0.17 and +44.9 at k = 0.23 (INV).

## 4. Explicit eclipsing-binary scenarios (s06, s06b, s08, s12, s16; summary `out/scenario_table.csv`, `out/summary.json`)

dchi2 is relative to the free planet fit of the same light curve. The total column also includes:
* the measured phase-0.5 depth;
* stellar radius and mass priors;
* for background stars, the DECaPS2 g, r, i, z photometry (0.05 mag floor) and the Gaia parallax, relative to a
  single-star fit of the same photometry (s12b: bg3 chi2 1.87, a K3-4V at 1.57 kpc with E(B-V) 0.03; bg5 chi2 1.47,
  a G0V at 3.8 kpc with E(B-V) 0.15);
* for the bound WD, -2 ln L of the Gaia RUWE.

All scenario fits use differential evolution or a dense grid, with every parameter free within physical bounds
and generous nuisance ranges: E(B-V) up to 0.30 (SFD gives 0.20), and a +-10% source flux.

| Group | Scenario | Data | dchi2 light curve | dchi2 total | Best-fit configuration |
|---|---|---|---|---|---|
| background EB | bg3_celm | REF | +10.1 | +19.6 | A 0.30 Msun, 0.175 Rsun, 8000 K; B 0.66 Msun, 0.64 Rsun; T14 2.21 h; sec 1.45 ppt; phot chi2 4.1, plx 3.4 |
| background EB | bg3_celm_lconly | REF | -0.0 | -0.1 | A 0.27 Msun, 0.060 Rsun, 11681 K; B 0.41 Msun, 0.42 Rsun; T14 1.60 h; sec 0.34 ppt; phot chi2 721.8, plx 21.8 |
| background EB | bg3_cms_unphysical | REF | +3.9 | +11.1 | A 1.40 Msun, 0.204 Rsun, 7683 K; B 0.67 Msun, 0.65 Rsun; T14 1.99 h; sec 1.91 ppt; phot chi2 2.6, plx 1.7 |
| background EB | bg3_csdb | REF | +16.0 | +84.2 | A 0.55 Msun, 0.100 Rsun, 20000 K; B 0.78 Msun, 0.75 Rsun; T14 1.89 h; sec 0.50 ppt; phot chi2 60.2, plx 0.9 |
| background EB | bg3_cwd | REF | +21.1 | +58.1 | A 0.17 Msun, 0.055 Rsun, 21605 K; B 0.66 Msun, 0.64 Rsun; T14 1.73 h; sec 0.21 ppt; phot chi2 23.2, plx 6.7 |
| background EB | bg3_cwd_lconly | REF | +0.0 | -0.1 | A 0.17 Msun, 0.055 Rsun, 11776 K; B 0.39 Msun, 0.40 Rsun; T14 1.60 h; sec 0.33 ppt; phot chi2 15853.6, plx 5.3 |
| background EB | bg3_msbd | REF | +133.2 | +143.2 | A 0.73 Msun, 0.713 Rsun, 4600 K; B 0.08 Msun, 0.25 Rsun; T14 3.23 h; sec 0.00 ppt; phot chi2 3.0, plx 0.0 |
| background EB | bg3_msms11 | REF | +86.8 | +92.5 | A 0.73 Msun, 0.712 Rsun, 4590 K; B 0.74 Msun, 0.72 Rsun; T14 5.15 h; sec 0.00 ppt; phot chi2 2.6, plx 3.4 |
| background EB | bg3_msms5 | REF | +56.6 | +65.5 | A 0.81 Msun, 0.776 Rsun, 5034 K; B 0.46 Msun, 0.47 Rsun; T14 3.63 h; sec 2.41 ppt; phot chi2 0.6, plx 1.8 |
| background EB | bg5_celm | REF | +17.7 | +27.4 | A 0.30 Msun, 0.351 Rsun, 7087 K; B 0.37 Msun, 0.36 Rsun; T14 2.58 h; sec 0.94 ppt; phot chi2 0.8, plx 9.6 |
| background EB | bg5_csdb | REF | +12.8 | +40.1 | A 0.55 Msun, 0.100 Rsun, 20000 K; B 0.70 Msun, 0.70 Rsun; T14 1.86 h; sec 0.36 ppt; phot chi2 13.2, plx 6.7 |
| background EB | bg5_cwd | REF | +19.1 | +55.2 | A 0.17 Msun, 0.055 Rsun, 32662 K; B 0.69 Msun, 0.67 Rsun; T14 1.72 h; sec 0.12 ppt; phot chi2 17.0, plx 11.6 |
| background EB | bg5_msbd | REF | +186.7 | +201.6 | A 0.83 Msun, 0.787 Rsun, 5122 K; B 0.08 Msun, 0.25 Rsun; T14 3.47 h; sec 0.00 ppt; phot chi2 6.5, plx 1.0 |
| background EB | bg5_msms11 | REF | +94.3 | +96.5 | A 0.86 Msun, 0.797 Rsun, 5170 K; B 0.88 Msun, 0.81 Rsun; T14 5.47 h; sec 0.00 ppt; phot chi2 3.5, plx 0.1 |
| background EB | bg5_msms5 | REF | +64.2 | +76.3 | A 0.88 Msun, 0.822 Rsun, 5295 K; B 0.55 Msun, 0.56 Rsun; T14 3.95 h; sec 2.42 ppt; phot chi2 4.5, plx 0.2 |
| compact star occulted (bound) | lateM | INV | +92.3 | +94.7 | M_A 0.15, R_A 0.0232, R_B 0.130, b 0.00, T14 0.81 h, ingress 7.4 min, sec 1.06 ppt, e 0.00 |
| compact star occulted (bound) | single | INV | +0.7 | +3.1 | M_A 0.47, R_A 0.0143, R_B 0.390, b 0.00, T14 1.39 h, ingress 3.0 min, sec 0.88 ppt, e 0.00 |
| compact star occulted (bound) | single_Rfree | INV | +0.2 | +0.2 | M_A 0.32, R_A 0.0050, R_B 0.381, b 0.24, T14 1.38 h, ingress 1.1 min, sec -0.24 ppt, e 0.00 |
| compact star occulted (bound) | twin | INV | +7.7 | +17.9 | M_A 0.25, R_A 0.0190, R_B 0.304, b 0.00, T14 1.31 h, ingress 4.6 min, sec 1.97 ppt, e 0.00 |
| compact star occulted (bound) | twin_Rfree | INV | +2.3 | +4.5 | M_A 0.15, R_A 0.0050, R_B 0.302, b 0.00, T14 1.34 h, ingress 1.3 min, sec 0.04 ppt, e 0.00 |
| compact star occulted (bound) | any | REF | +2.9 | +2.9 | M_A 0.60, R_A 0.0124, R_B 0.382, b 0.00, T14 1.28 h, ingress 2.4 min, sec 0.27 ppt, e 0.00 |
| compact star occulted (bound) | lateM | REF | +112.9 | +113.9 | M_A 0.15, R_A 0.0232, R_B 0.130, b 0.00, T14 0.81 h, ingress 7.4 min, sec 1.06 ppt, e 0.00 |
| compact star occulted (bound) | single | REF | +2.9 | +2.9 | M_A 0.61, R_A 0.0124, R_B 0.381, b 0.02, T14 1.28 h, ingress 2.4 min, sec 0.27 ppt, e 0.00 |
| compact star occulted (bound) | single_ecc | REF | +1.8 | +1.8 | M_A 0.78, R_A 0.0109, R_B 0.381, b 0.34, T14 1.40 h, ingress 2.6 min, sec -0.34 ppt, e 0.20 |
| compact star occulted (bound) | twin | REF | +4.1 | +11.1 | M_A 0.25, R_A 0.0189, R_B 0.302, b 0.00, T14 1.30 h, ingress 4.6 min, sec 1.98 ppt, e 0.00 |
| bound hierarchical EB | best P5.69 | INV | +47.6 | +47.6 | m1 0.15, BD m2 0.07, k 0.42, m0 0.3405486778374248, f_src 0.107, P 5.69, b 0.00, T14 1.37 h |
| bound hierarchical EB | best P11.38 | INV | +27.5 | +27.5 | m1 0.11, MS m2 0.11, k 1.00, m0 0.34145592097617955, f_src 0.075, P 11.38, b 0.77, T14 1.70 h |
| bound hierarchical EB | best occulter_BD | INV | +47.6 | +47.6 | m1 0.15, BD m2 0.07, k 0.42, m0 0.3405486778374248, f_src 0.107, P 5.69, b 0.00, T14 1.37 h |
| bound hierarchical EB | best no_third_light | INV | +46.3 | +46.3 | m1 0.28, MS m2 0.26, k 0.94, m0 nan, f_src 1.000, P 11.38, b 1.70, T14 1.42 h |
| bound hierarchical EB | best P5.69 | REF | +58.4 | +58.5 | m1 0.15, BD m2 0.07, k 0.42, m0 0.3405486778374248, f_src 0.107, P 5.69, b 0.00, T14 1.37 h |
| bound hierarchical EB | best P11.38 | REF | +33.7 | +33.7 | m1 0.11, MS m2 0.11, k 1.00, m0 0.34145592097617955, f_src 0.075, P 11.38, b 0.74, T14 1.71 h |
| bound hierarchical EB | best occulter_BD | REF | +58.4 | +58.5 | m1 0.15, BD m2 0.07, k 0.42, m0 0.3405486778374248, f_src 0.107, P 5.69, b 0.00, T14 1.37 h |
| bound hierarchical EB | best no_third_light | REF | +54.8 | +54.8 | m1 0.28, MS m2 0.26, k 0.94, m0 nan, f_src 1.000, P 11.38, b 1.70, T14 1.44 h |

Scenario labels:
* bg3 / bg5 = the 3.1" and 5.5" background stars.
* cwd / csdb / celm = white dwarf / hot subdwarf / pre-ELM white dwarf, totally or partially occulted by a
  main-sequence star.
* msms5 / msms11 = two main-sequence stars at P / 2P.
* msbd = main-sequence star + brown dwarf.
* cms_unphysical = compact star with free M, R and Teff (unphysical optimum).
* single / twin / lateM / any = WD occulted by the M3V / by one twin component / by a bound late-M dwarf / by an
  M dwarf of free size.
* _Rfree = WD radius free down to 0.005 Rsun (unphysical).
* _ecc = eccentric.
* _lconly = photometry and parallax terms dropped; the 'total' column then compares light curve + secondary only.

For the WD scenarios the 'total' column does not yet include the astrometric term; see 4a.

### 4a. Hot compact star totally occulted by the M3V (bound; s06, s06b, s16)

**The light curve cannot distinguish this scenario from a planet.**
* The occultation is box-shaped, with 1-3 min ingress and T14 = 1.28-1.45 h; the planet has 13-15 min ingress
  and T14 = 1.7 h. Both have the same width at half depth.
* Real data: dchi2(WD - planet) = +0.79 (INV) and -0.33 (REF).
* Injections (s13):
  * planet truth: +0.8 +- 1.4 (-1.4 to +2.8);
  * WD truth (0.47 Msun box): -2.9 +- 1.6 (-5.3 to -1.1).
* So the INV value sits 2.3 sigma above the WD-truth mean and at the planet-truth mean, a shape likelihood ratio
  of roughly 15:1 for the planet in INV only. REF has no injection calibration.

**Discriminants that work** (REF, WD radius tied to its mass; `out/wd_mass_profile_ref_cool.json`, `out/wd_checks.json`):

| M_WD (Msun) | dchi2 LC + secondary + priors | secondary (ppt) | WD Teff needed (K) | -2 ln L RUWE (sigma_CCD 0.35 / 0.22) | total (0.35 / 0.22) |
|---|---|---|---|---|---|
| 0.17 | +17.9 | 3.25 | 4662 | 1.6 / 6.2 | +19.5 / +24.1 |
| 0.20 | +13.9 | 2.88 | 4776 | 2.0 / 8.0 | +16.0 / +21.9 |
| 0.25 | +9.5 | 2.38 | 4957 | 4.2 / 10.1 | +13.8 / +19.7 |
| 0.30 | +6.7 | 1.97 | 5129 | 6.4 / 12.3 | +13.2 / +19.0 |
| 0.35 | +4.9 | 1.62 | 5297 | 8.4 / 12.3 | +13.3 / +17.2 |
| 0.45 | +3.1 | 0.99 | 5693 | 12.3 / 12.3 | +15.4 / +15.4 |
| 0.60 | +2.9 | 0.29 | 6291 | 12.3 / 12.3 | +15.2 / +15.2 |
| 0.80 | +5.7 | 0.03 | 6271 | 12.3 / 12.3 | +17.9 / +17.9 |
| 1.00 | +13.6 | -0.82 | 7499 | 12.3 / 12.3 | +25.9 / +25.9 |
| 1.20 | +28.3 | -1.59 | 10324 | 12.3 / 12.3 | +40.6 / +40.6 |

* Low-mass (He) WDs are large, so their transit across the M3V at phase 0.5 (1.6-3.3 ppt) is excluded by the
  measured -0.16 +- 0.79 and +0.29 +- 0.77 ppt.
* CO WDs (>= 0.45 Msun) would move the photocentre by 0.45-0.8 mas on a 5.69-d orbit. Simulated RUWE has median
  1.3-2.2; observed is 1.076, and 0 of 400 simulations reach RUWE <= 1.10.
* Minimum combined dchi2 = **13.2** (0.30 Msun, conservative sigma_CCD) or **15.2** (0.6 Msun).
* The simulation follows referee_star r10: 41 FoV transits x 9 CCDs, a 5-parameter fit and an empirical RUWE0
  distribution. It is not the Gaia pipeline.

**Supporting checks** (not used in the dchi2):
* K1 = 35-106 km/s for WD masses 0.2-1.2. The Gaia DR3 combined-spectrum RV error (2.22 km/s, rv_method_used = 2)
  is at the 38th percentile of 1642 M dwarfs of the same grvs_mag. Smearing by K1 >~ 35 km/s should have broadened
  the combined CCF; this is qualitative.
* SkyMapper v - g is consistent with the empirical M-dwarf locus (n = 411, BP-RP 2.66-2.86): residual -0.08 with a
  robust sd of 0.13. This rules out the blue excess of a >= 0.8 Msun WD (2.5-6 sigma) but not lighter ones.
* SkyMapper g - i is 0.12 mag bluer than the locus (robust sd 0.019, 1.7th percentile). A 0.6 Msun WD would give
  -0.11, but DECaPS2 g shows a deficit relative to the model, not an excess. This is recorded as unexplained,
  not as evidence.

**Other occulters.** A WD around one twin component gives +11.1 (REF) and +17.9 (INV): T14 = 1.3 h is too short
and the secondary is 2 ppt. A WD occulted by a bound late-M dwarf gives +114 (T14 = 0.81 h). Letting the occulter
radius and mass float ('any') returns the M3V itself (+2.9). An eccentric orbit (REF; e = 0.20, omega = 261 deg, M_WD = 0.78) improves the light-curve fit only to +1.8, and the 0.78 Msun WD it needs is excluded by RUWE (-2 ln L >= 12.3). Post-common-envelope WD + M binaries at 5.7 d are circular in any case.

### 4b. Bound hierarchical eclipsing binaries (s08, `out/bound_eb_{ref,inv}.csv`)

**Grid.**
* Eclipsed star m1 (0.08-0.40 Msun) and occulter m2: main sequence up to m1, or a brown dwarf with
  R = 0.075-0.105 and L = 0.
* The third light m0 takes whatever Ks flux is left of the observed MKs = 6.44.
* Mann+19 masses, Mann+15 radii; T - Ks and Teff from Mamajek. Kepler fixes a/R.
* P = 5.69 d (secondary at phase 0.5, constrained) or 11.38 d (alternating eclipses, both fitted). b is free.
* 576 combinations.

**Best per class:**

| Class | Best configuration | dchi2 INV | dchi2 REF |
|---|---|---|---|
| any MS + MS at 11.38 d | 0.11 + 0.11 Msun, diluted by 0.34 Msun, b 0.74, V-shaped | +27.5 | +33.7 |
| BD occulter (any host) | 0.15 Msun + 0.07 Msun BD, diluted by 0.34 Msun | +47.6 | +58.5 |
| BD grazing the M3V (m1 0.36-0.38) | grazing, T14 1.2 h | | +84 to +90 |
| BD grazing a twin component (m1 0.26, m0 0.27) | | | +86 |
| no third light | 0.28 + 0.26 Msun at 11.38 d, grazing | +46.3 | +54.8 |
| MS + MS at 5.69 d | small M dwarf grazing the M3V, secondary 2.2 ppt | | >= +98 |

These configurations fail on three counts:
* At the stellar densities, a 1.6-1.7 h eclipse needs a central chord.
* A central chord with k ~ 1 is V-shaped, or with k ~ 0.3 is too deep.
* Stellar occulters bring a secondary eclipse.

A first version that applied an optical (Mamajek M_G, BP-RP) luminosity budget rejected the over-luminous single
M3V itself. It is kept as `out/bound_eb_ref_mamajekG.json` for the record (best +27.7); only the Ks-budget version
is used.

### 4c. Background stars (s12, s12b; `out/bg_eb_*.json`)

**bg3 (Gaia DR3 5925373272827260544, 3.1" NE; F = 30.3 e-/s, f = 0.031; needs >= 50%).**
* Its DECaPS2 colours (g-r 0.81, r-i 0.25, i-z 0.20) are those of a single K3-4 dwarf. A hot compact component that
  supplies >= 50% of the TESS light would make it much bluer.
* Every physical class has a large total dchi2:
  * pre-ELM (0.3 Msun, 0.18 Rsun, 8000 K, grazing behind a K7V): +19.6;
  * WD: +58 (the light-curve fit forces an ELM with 0.055 Rsun at 21,600 K, which the photometry and parallax
    reject);
  * sdB: +84;
  * MS + MS: +66 (5.69 d; T14 3.6 h) and +93 (11.38 d; T14 5.2 h);
  * MS + BD: +143.
* Main-sequence K-dwarf eclipses at P = 5.69 d last >= 2.3 h unless grazing; the eclipse would have to be
  50-60% deep; and an occulting K dwarf brings a secondary.
* An unconstrained 'compact' class finds a non-existent object (1.4 Msun, 0.2 Rsun, 7700 K) at +11.1.
* **The exclusion of bg3 is photometric, not light-curve based.** With the photometry and parallax dropped
  (`_lconly`), an ELM/He WD (0.17-0.27 Msun, 0.055-0.06 Rsun, ~11,700 K) occulted by a 0.39-0.41 Msun M dwarf fits
  the light curve exactly as well as the planet (dchi2 0.0). That is the same geometry as the bound WD scenario,
  since only the occulter's density matters. Such a system cannot be bg3: it would be an M dwarf + WD with
  r - i ~ 1, while bg3 has r - i = 0.25 and K-dwarf luminosity at its parallax (photometric chi2 720-15,850).
* Phase 1 pixel localisation disfavours bg3 at about 3-4 sigma (injection Mahalanobis^2 14.9 against 2.4; host
  test dchi2 9.9). That evidence is independent of this.

**bg5 (Gaia DR3 5925373272827262208, 5.5"; F = 22.1 e-/s, f = 0.0225; needs >= 68%).** Its photometry is a single G0 dwarf at ~3.8 kpc (chi2 1.5).
* Best physical EB: a 0.3 Msun pre-ELM (0.35 Rsun, 7100 K) grazing behind an M3 dwarf at 1.5 kpc, dchi2 = +27.4
  (LC +17.7; the parallax disagrees, chi2 9.6).
* sdB: +40.1. WD: +55.2. Main-sequence pairs: +76.3 (P) and +96.5 (2P). BD: +202.
* Phase-1 pixels disfavour it at p ~ 0.001.

## 5. Secondary eclipse and phase curve (s10, `out/secondary_phase.json`, `out/secondary_scan.csv`, `figs/secondary_scan.png`)

**Occultation depth.**
* Measured with the planet's own shape (uniform disk) at phase 0.5: INV -0.16 +- 0.79 ppt. The error is
  empirical, from 400 random phases; the formal error is 0.95. REF (referee r18): +0.29 +- 0.77 ppt.
* 3-sigma upper limit: 2.2 ppt (INV) and 2.6 ppt (REF).
* Scanning all phases for eccentric orbits:
  * The largest positive excursion is +1.67 ppt at phase 0.492 (2.1 sigma, global p = 0.16 for 71 independent
    trials). The same bump is seen in the investigator's scan (1.9 ppt at 0.490) and the referee's (1.3 ppt at
    0.495).
  * The largest negative excursion is -2.9 ppt at phase 0.576 (3.6 sigma). The referee scan shows a +2.4 ppt dip
    there instead.
  * The any-phase upper limit is 4.0 ppt.

**Planet limits** (Rp/a = 0.0054).
* A_g < 77 (3 sigma): not constraining. A_g = 1 would give 29 ppm, and thermal emission at Teq ~ 470-500 K is
  < 1 ppm in TESS.
* TESS-band brightness temperature < 2460 K (blackbody vs a 3350 K blackbody star).

**Companions of the M3V excluded.**
* By the secondary: any self-luminous transiting body hotter than ~2460 K, i.e. M8V (2570 K) and earlier. This
  assumes the occulted areas at the two conjunctions are equal (circular orbit, non-grazing).
* By size: with the density prior, k <= 0.133 (3 sigma). An old brown dwarf (R >= 0.075 Rsun) needs
  k >= 0.197 on the single star and >= 0.26 on a twin component. Forcing such k costs dchi2 >= 45 (INV k-profile)
  and >= 84 in the physical grid (s08, REF).
* Low-mass stars additionally fail on the secondary.

**Phase curve** (raw INV flux; per-orbit quadratic baselines; empirical errors from 300 random periods of 2-12 d).
* Amplitude at P: 2.9 ppt (21% of random periods are larger). At P/2: 0.4 ppt.
* 95% / 99.7% limits: 4.2 / 10.3 ppt at P, and 2.1 / 2.6 ppt at P/2.
* Expected semi-amplitudes (beaming / ellipsoidal):
  * 70 M_J brown dwarf: 0.25 / 0.014 ppt;
  * 0.1 Msun star: 0.35 / 0.02 ppt;
  * 0.37 Msun twin: 0.97 / 0.05 ppt;
  * 0.6 Msun WD: 1.31 / 0.06 ppt;
  * K1 = 13.8, 19.7, 53.9 and 73 km/s respectively.
* **No companion class is excluded by the phase curve.** With the Mira 18" away, slow variability below ~1-2%
  cannot be measured, and the detrended curve removes P/2 signals by construction.

## 6. Per-sector consistency and adopted parameters (s09, s15)

**Per sector** (`out/per_sector.csv`, `out/per_sector.json`, `figs/per_sector.png`).
* The global shape is fixed and each sector gets k alone (depth) or k plus a duration scale. Errors are
  profile-likelihood errors with the beta-inflated noise.

| Sector | INV depth (ppt) | INV-A depth | REF depth | REF duration scale |
|---|---|---|---|---|
| 12 | 13.0 +- 4.4 | 6.9 (no error: profile failed) | 14.1 +- 4.1 | 0.64 +- 0.27 |
| 39 | 13.0 +- 2.0 | 15.5 +- 3.0 | 14.7 +- 1.9 | 0.96 +- 0.09 |
| 66 | 13.3 +- 2.8 | 8.7 +- 3.7 | 13.8 +- 2.4 | 1.04 +- 0.09 |
| 93 | 16.1 +- 2.2 | 22.7 +- 5.1 | 14.8 +- 2.7 | 0.98 +- 0.11 |
| 103 | 19.4 +- 3.5 | 20.6 +- 5.0 | 17.2 +- 1.9 | 1.01 +- 0.09 |
| 104 | 20.5 +- 2.3 | 18.9 +- 4.1 | 18.4 +- 2.1 | 1.02 +- 0.07 |
| all | 15.9 +- 1.0 | 15.6 +- 1.7 | 15.8 +- 0.9 | 1.00 +- 0.04 |
| odd / even | 18.3 +- 1.4 / 13.3 +- 1.5 | 19.0 +- 2.3 / 14.1 +- 2.5 | 16.4 +- 1.3 / 15.1 +- 1.4 | 1.06 / 0.97 |

* Depth consistency chi2 = 8.4/5 (INV), 7.1/5 (INV-A), 3.4/5 (REF). Duration consistency 7.7/5, 1.5/5, 2.2/5.
* S12 (30-min cadence, 2-3 points per transit) gives a short duration scale (0.59 +- 0.19 INV, 0.64 +- 0.27 REF).
  That is the only excursion near 2 sigma, and it is consistent with the coarse sampling.
* Depths in S103/S104 (17-20 ppt) exceed those in S39/S66 (13-15 ppt) in all three photometries, at about
  1.5-2 sigma.
* Odd minus even: INV +5.0 +- 2.1 ppt, INV-A +4.9 +- 3.4, REF +1.3 +- 1.9 ppt. This reproduces the phase-1 finding
  that only model B shows it.
* By cadence: 30 min (S12) 13-14 ppt, 10 min (S39) 13-15.5 ppt, 200 s (S66-S104) 16.3-17.3 ppt. These agree within
  ~1.5 sigma, and the difference is the same time trend. Model A (all-free deblending, different template
  correlations) gives the same mean depth (15.6 +- 1.7 against 15.9 +- 1.0 ppt), so there is no sign of an
  aperture- or template-dependent dilution.

**Adopted parameters** (`out/final_params.json`).
* Mean of the INV and REF 'single' and 'twin' posteriors.
* Uncertainty = the larger of the two 68% half-widths, added in quadrature to half their difference.
* R* = 0.381 +- 0.012 (single) or 0.289 +- 0.009 (twin); Teq for A_B = 0 and full redistribution.

| Parameter | Single M3V host | One component of an equal twin |
|---|---|---|
| Rp/R* (undiluted) | 0.1189 +- 0.0050 | 0.1664 +- 0.0062 |
| **Rp (R_earth)** | **4.94 +- 0.26** | **5.23 +- 0.25** |
| b | 0.53 +- 0.10 | 0.20 +- 0.14 |
| a/R* | 25.3 +- 0.8 | 29.7 +- 0.9 |
| i (deg) | 88.81 +- 0.24 | 89.61 +- 0.28 |
| T14 (h) | 1.70 +- 0.07 | 1.67 +- 0.05 |
| T23 (h) | 1.21 +- 0.11 | 1.18 +- 0.06 |
| a (au) | 0.0449 +- 0.0019 | 0.0398 +- 0.0017 |
| S (S_earth) | 8.2 +- 0.7 | 5.2 +- 0.5 |
| **Teq (K)** | **471 +- 10** | **420 +- 9** |
| rho* (prior-dominated) | 9.5 +- 0.9 | 15.3 +- 1.3 |
| P (d), T0 (BTJD_TDB) | 5.691508 +- 0.000009, 3840.9155 +- 0.0015 (INV/REF mean; REF alone 3840.91594 +- 0.00132) | same |

**Caveats on the parameters.**
* Rp scales as sqrt(depth). The depth normalisation (979 e-/s) carries +-1.0 ppt of systematic error, i.e.
  +-3% in Rp, which the table does not include.
* The single-star Rp is lower than the referee_star 5.12 R_earth because the fitted depth is lower: 15.8 ppt
  with LD, not 16.4.
* Eccentric orbits would change b and a/R*. The data do not constrain them.

## Key numbers

| Quantity | Value | Method | Uncertainty | File |
|---|---|---|---|---|
| Rp/R* (free) | 0.119 (INV), 0.119 (REF) | dynesty, batman, Claret LD prior | +0.008/-0.005, +-0.005 | out/fit_{inv,ref}_free.json |
| rho_circ | 10.3 (INV), 12.3 (REF) g/cm3 | same | +5.4/-5.9, +4.0/-5.8; injection scatter 4.9 | same; out/inject_planet_wd.json |
| ln B single/twin | -0.43 (INV), -0.53 (REF) | nested-sampling evidences | +-0.2 | out/fit_*_{single,twin}.json, out/photoecc.json |
| e (single host) | < 0.74 (95%), P(e<0.1) = 0.20 | photo-eccentric importance sampling | | out/photoecc.json |
| Max dT of a blended source, transit geometry | ~3.5-3.8 (3 sigma) | chi2 profile in f_src | dchi2 = 9 threshold | out/dilution_profile_*.json |
| Max dT, occultation geometry | 4.5 (no shape limit; set by a 100% deep eclipse) | same | box occultation dchi2 +0.7/+2.9 | same; out/compact_*_single.json |
| WD occulted by the M3V | dchi2_LC +0.7/+2.9; combined with secondary + RUWE >= 13.2 | s06/s06b/s16 | RUWE simulation simplified | out/summary.json |
| Bound EB (any) | dchi2 >= +27.5 (INV), +33.7 (REF) | physical grid, Mann relations | | out/bound_eb_*.json |
| bg3 EB (any physical class) | dchi2 >= +19.6 (REF) | joint LC + photometry + parallax | | out/bg_eb_ref_bg3_*.json |
| bg5 EB | dchi2 >= +27.4 (REF) | same | | out/bg_eb_ref_bg5_*.json |
| Secondary (phase 0.5) | -0.16 +- 0.79 ppt (INV); +0.29 +- 0.77 (REF) | template fit; 400 random phases | empirical | out/secondary_phase.json |
| T_b of occulter | < 2460 K (3 sigma) | blackbody ratio | | same |
| Phase-curve limits | A(P) < 4.2 ppt, A(P/2) < 2.1 ppt (95%) | sinusoid regression, random-period calibration | | same |
| Per-sector depth chi2 | 3.4/5 (REF), 8.4/5 (INV) | per-sector k | profile errors | out/per_sector.json |
| Rp (single / twin) | 4.94 +- 0.26 / 5.23 +- 0.25 R_earth | posterior + R* | +3% depth-normalisation systematic | out/final_params.json |
| Teq (single / twin) | 471 +- 10 / 420 +- 9 K | | | same |

## Scripts (`scripts/`, run with /home/user/work/venv/bin/python, OMP_NUM_THREADS=1-2)

| Script | Purpose | Output |
|---|---|---|
| `shapelib.py` | LC class (data, batman per cadence, per-window offsets), transit and EB flux models, durations, density conversions | (library) |
| `eem.py` | Mamajek sequence interpolation; Stassun G - T; Mann+15/19 MKs-mass-radius relations | (library) |
| `specs.py` | BT-Settl (< 4000 K, ../star/data/models) and ATLAS9 ck04 (data/ck04/, downloaded from STScI CDBS) surface fluxes; DECam/SkyMapper/TESS band fluxes (SVO filters in data/filt_*.dat) | (library) |
| `s00_ld.py` | Claret 2017 TESS LD (VizieR) interpolation | out/ld_claret17.json, data/claret17_tableab.csv |
| `s01_data.py` | data sets, white and red noise | data/{inv_win,inv_full,ref_win}.npz, out/noise.json |
| `s02_batman_check.py` | batman vs brute force, k = 0.1-5 | out/batman_check.json |
| `s03_planet_fit.py` (+ `run_s03.sh`) | dynesty fits: free, single, twin (also freeLD, single_ecc, twin_ecc options, not completed) | out/fit_{inv,ref}_{model}.{json,npz}, out/raw_*.pkl |
| `s04_photoecc.py` | density test, eccentricity posteriors, evidence table | out/photoecc.json |
| `s05_dilution_profile.py` | chi2 vs source flux fraction (both geometries); k-profiles | out/dilution_profile_{inv,ref}.json, out/s05_*.log |
| `s06_compact.py` (+ `run_s06.sh`, `run_s06b.sh`) | compact star occulted by the M3V / twin / late-M / any; WD M-R relation; secondary with self-lensing | out/compact_{ds}_{scenario}.json, out/s06.log |
| `s06b_wd_mass_profile.py` | chi2 vs WD mass (cool and hot radius branches) | out/wd_mass_profile_ref_{cool,hot}.json |
| `s07_gaia_bg.py` | Gaia DR3 parameters of bg3/bg5/M3V; SFD extinction | out/gaia_bg.json, data/gaia_bg.csv |
| `s08_bound_eb.py` | bound hierarchical EB grid (Ks luminosity budget) | out/bound_eb_{ds}.{json,csv} |
| `s09_per_sector.py` | per-sector / cadence / odd-even depth and duration | out/per_sector.{csv,json} |
| `s10_secondary_phase.py` | secondary depth, phase scan, albedo/T_b limits, phase-curve limits, predicted beaming/ellipsoidal | out/secondary_phase.json, out/secondary_scan.csv |
| `s11_specs_check.py` | validates synthetic photometry: reproduces the star agent's SED model magnitudes of the M3V to < 0.01 mag | out/specs_check.json |
| `s12_bg_eb.py` (+ `run_s12.sh`) | background-star EB fits joint with photometry and parallax | out/bg_eb_{ds}_{src}_{class}.json |
| `s12b_bg_single.py` | single-star fits of the background stars' photometry (null) | out/bg_single.json |
| `s13_inject.py`, `s13b_real.py` | injection-recovery at 8 fake epochs (planet and WD truths); the same fits on real data | out/inject_planet_wd.json, out/inject_rows_planet_wd.csv, out/inject_real.json |
| `s14_figs.py` | figures | figs/*.png |
| `s15_final.py` | adopted parameters | out/final_params.json |
| `s16_wd_checks.py` | WD scenario: RUWE simulation, SkyMapper colour locus (CDS XMatch with II/379/smssdr4), RV percentile | out/wd_checks.json, data/sm4_xmatch.csv |
| `s17_summary.py` | scenario table | out/scenario_table.csv, out/summary.json |

## Figures (`figs/`)

* `fold_models.png`: both photometries folded, with the planet, WD-occultation and best geometric-blend models.
* `rho_posteriors.png`: rho_circ posteriors against the single and twin densities, and b against rho.
* `dilution_profile.png`: chi2 against source dT for both geometries, with bg3 and bg5 marked.
* `per_sector.png`: per-sector depths and duration scales.
* `secondary_scan.png`: occultation depth against phase.

## Caveats

* All chi2 use per-point errors inflated by per-sector beta. REF's beta comes from the referee's fake-window
  ratios, which are noisy (S93: 2.0). The dchi2 values are therefore conservative but approximate. Only the
  planet/WD shape comparison and the density recovery are injection-calibrated (INV, 8 epochs per truth).
* Scenario fits are maximum-likelihood fits over bounded physical parameter spaces (differential evolution,
  multi-start Powell, or grids). A better optimum could exist in principle. The bounds were set generously:
  E(B-V) up to 0.30, source flux +-30%, WD radius up to 2.5x Nauenberg. Several optimisers hit those bounds, which
  shows the scenarios were pushed to their most favourable corners.
* Model stellar colours (ATLAS9 logg 5 for compact stars; BT-Settl for M dwarfs) carry ~0.05-0.1 mag
  systematics. WD and sdB colours at logg 7-8 differ slightly from logg 5.
* The RUWE simulation is the simplified referee_star r10 model (scan law, parallax factors and the epoch count of
  this source are not reproduced). A dedicated Gaia-epoch simulation, or Gaia DR4 epoch astrometry, would make
  the WD exclusion firmer.
* The background-star EB fits were run on REF only. The planned INV cross-checks and the bg5 LC-only fits
  (`run_s12.sh`, last lines) were stopped to free the shared CPUs; the bound scenarios were run on both.
* Eccentric false-positive orbits were tested only for the WD scenario (the REF eccentric WD fit: dchi2 +1.8 before astrometry). Post-interaction compact
  binaries at P = 5.7 d are expected to be circular.
* The light curves come from the same pixels. INV and REF are two analyses, not independent measurements.

