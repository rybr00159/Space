# Referee 3: adversarial check of `bd/` (brown-dwarf and companion-hosted channels) for TIC 214311908

Folder: `tic214311908/analysis/referee3_bd/`

AI-run referee analysis, 2026-10-17. Nothing was posted, uploaded or submitted, and no one was contacted. Network use
was read-only: arXiv API abstracts, the NASA Exoplanet Archive TAP service, the Gaia DR3 archive and one web search
(Gaia DR4 date).

All code is my own (`scripts/`). I read the analyst's scripts in `../bd/scripts/` to look for errors but did not import
or run them. `../referee2_fpp/scripts/r05_fpp.py` is imported READ-ONLY (as the analyst's b10 did) to re-evaluate that
framework with the sign corrected by my own patch. Nothing is written outside this folder. Light curves: the CLEANED
curves of `../clean/` (`lc_M3V_REF_clean.csv`, `lc_M3V_INV_clean.csv`), which the analyst did not have; REF fake
windows from `../referee_phot/data/win_s*.npz`.

## Verdict

1. **The analyst's main finding is correct.** Both FPP codes (`../fpp/scripts/f04_bayes.py`, `../referee2_fpp/scripts/r05_fpp.py`)
   evaluate the light curve at rho_c = rho_host g^-3, with g = (1 + e sin w)/sqrt(1 - e^2). The correct relation is g^+3.
   * I checked it with my own Kepler-orbit sky-projection code, not batman: rho_c/rho = g^3 to < 0.2% for 7 (e, w) cases
     (`out/s01_sign_tides.json`).
   * With the wrong sign, the transit-weighted prior probability of rho_c/rho* <= 0.2 is 8.0% instead of 0.89% (Kipping
     e). That is a factor 9 in favour of a grazing eccentric BD.
2. **The sign-fixed numbers of b10 are reproduced** within Monte Carlo noise. r05 uses only 2,500 e draws; I used 10,000-20,000
   (`out/s02_fpp_fixed.csv`).
   * REF: FPP 0.00151 (b10 0.00144).
   * INV: FPP 0.00359 (b10 0.00345).
   * BD on the M3V, f_BD 1e-3: 1.2e-4 (b10 0.98e-4).
3. **The analyst's FPP is upheld, now as a full computation rather than a rescaling.** Settings: sign fixed, f_BD 1.3e-4,
   Ferreira+26 short-period e law for BDs, beta noise.

   | Planet e prior | FPP, REF | FPP, INV |
   |---|---|---|
   | Kipping | 0.0013 | 0.0031 |
   | Sagear & Ballard (2023) M-dwarf Rayleigh 0.19 | 0.0016 | 0.0033 |
   | Circular | 0.0007 | 0.0019 |

   * The analyst's 0.0013-0.0030 stands.
   * P(planet on some star of the bound system) = 1 - FPP_system = 0.9995 (REF) / 0.9995 (INV) in the adopted settings, and
     0.9986-0.9996 over all 64 variants (analyst >= 0.997-0.999).
4. **The BD channels are smaller than the analyst's "most likely 1e-5 to 5e-5".** The new evidence is the BD-truth
   calibration that bd/ did not run, plus an independent evidence calculation on the cleaned curves.
   * BD-truth injections never reach the observed shape statistic: max Lambda_G 3.2 (REF, n = 60) and 3.4 (INV, n = 28),
     against 9.2 and 6.6 observed.
   * Planet-truth injections reproduce the observed values: REF median 9.5, observed at quantile 0.48; INV median 5.6,
     quantile 0.61.
   * Calibrated LR(BD/planet) for the shape: 0.004 (REF) and 0.03 (INV). The r03 grids behind r05 carry 0.05 and 0.11,
     so they are 13x (REF) and 3.4x (INV) too lenient.
   * BD total after the Gaia RVS factor (median 0.5):
     * 8e-6 (REF) / 2.4e-5 (INV) in the sign-fixed r05;
     * 0.6e-6 / 7e-6 with my calibrated shape LR.
   * With the pessimistic long-period e law: <= 8.5e-5.
   * Most likely range: **~1e-6 to 2e-5**. Upper envelope: **<= 1.5e-4**, reached only for INV with the long-period e law
     and no Gaia factor.
   * After the fix, the BD channel no longer sets the floor after imaging. Speckle + two-band + the fixed sign gives
     FPP ~1.1e-4, dominated by the WD-occultation term.
5. **PTP (a planet on an unresolved ~0.20 Msun bound companion) is the only significant residual.**
   * Kipping e: 8.6e-4 (REF), 2.6e-3 (INV).
   * Rayleigh 0.19: 1.2e-3 / 2.85e-3.
   * Circular orbits: 3e-5 / 1.15e-3.
   * In REF it is mainly an eccentric, apoastron-transit channel. A circular planet on a rho ~26 g/cm3 companion transits
     too briefly. My own representative-companion evidence (s06) gives the same 15-30x circular/eccentric contrast.
   * An 8-m 832-nm speckle null removes 93% of it, giving FPP 0.0005 (REF) / 0.0007 (INV).
6. **Errors and overstatements found in bd/**, each detailed in sections 1-2:
   * a. **REF Lambda_G = 13.6 depends on the noise model.**
     * On the cleaned REF curve with clean/'s fake-epoch beta it is 9.2; with white noise 21.0.
     * INV agrees: 6.6 against 6.8.
   * b. **"With white noise every BD term is below 1e-6" is false for INV.**
     * b10's own output gives 1.0e-4 (f_BD 1e-3, INV white, sign fixed).
     * My grid gives 1.1e-5 at f_BD 1.3e-4.
   * c. **The tidal conclusion is overstated.**
     * The small-e Jackson+08 rate underestimates damping at e = 0.5 by 13-17x (Hut 1981 eccentricity functions).
     * For Q'* = 1e6, tau_e at e = 0.5 is 1.1-1.5 Gyr, not 18 Gyr.
     * An initial e of 0.55-0.7 falls to e ~0.37 within 1 Gyr.
     * Tides therefore disfavour a present-day e >= 0.5 unless Q'* >~ 1e7 or the system is young.
   * d. **The companion-host solution k = 0.40, D ~0.09 is not "self-consistent with a ~0.9 R_J body on a ~0.2 Msun
     companion".**
     * A 0.2 Msun companion supplies D ~0.25 of the TESS flux (r04, Mann relations).
     * D = 0.09 means a ~0.14 Msun companion of R ~0.17 Rsun, so k = 0.40 is a 0.66 R_J body: a giant planet, which PTP
       already counts, not a BD.
     * A self-consistent BD (R >= 0.75 R_J) needs k >~ 0.45 there: +7.4 (REF) / +2.4 (INV) by the analyst's own profile.
     * A full transit of a 0.85 R_J BD gives the observed ~15 ppt only for m2 ~0.10-0.12 Msun, with k ~0.6-0.7, which the
       flat-bottomed shape disfavours.
   * e. **The required eccentricity is e >= 0.50 at exact apoastron, not "e >~ 0.55".** rho_c/rho* = 0.19 at the
     profile minimum. Over the Delta-chi2 <= 4 range of rho_c (~1.15-2.55 g/cm3 REF, ~0.95-2.5 INV) the bound is e >~ 0.42-0.62.
   * f. **f_BD uses Bryant+23's full-sample giant rate (0.194%).** Their 0.26-0.42 Msun bin, which contains this 0.37 Msun
     host, gives 0.108 +- 0.083%, so f_BD ~7e-5. This sits inside the quoted 68% range.
   * g. **"2-3 epochs at 0.3-1 km/s" is too loose at the low-mass end.**
     * At 1 km/s, only 39-63% of 13 M_J BDs are detected (FAP 1e-3), even at fixed phases.
     * Closing the channel for 13 M_J needs <= 0.5 km/s (grazing configuration, 2 quadrature epochs) or <= 0.3 km/s (3
       epochs at phases 0.2/0.5/0.8, any BD e).
     * For the empirical BD mass distribution, <= 0.7-0.9 km/s suffices.
     * Random epochs fail: 2 random epochs miss 35-55% of 13 M_J BDs even at 0.3-0.5 km/s.
7. **Not re-derived here; carried from referee2_fpp:** the NTP pair-A term (3-7e-4, now the second-largest FP), the
   companion configuration priors, f_N = 0.003 and the WD/HEB terms. The FPP numbers inherit them.

## 1. Key numbers: analyst vs referee

| Quantity | bd/ (analyst) | Referee 3 | Method / file | Verdict |
|---|---|---|---|---|
| Photo-eccentric sign in f04/r05 | wrong (g^-3); correct is g^+3 | rho_c/rho = g^+3 to < 0.2% (7 cases) | own Kepler-orbit code, `out/s01_sign_tides.json` | **upheld** |
| P_tw(rho_c/rho* <= 0.2) | 0.7-0.9% (2.4% long-period law) | 0.67% (BD short), 0.89% (Kipping), 1.0% (Van Eylen), 2.8% (BD long); wrong sign: 8.0% (Kipping) | 2e6 draws, `s01` | upheld |
| b10 sign-fixed r05, REF / INV | FPP 0.0014 / 0.0034; BD 1.3e-4 / 4.3e-4 (f_BD 1e-3) | 0.0015 / 0.0036; BD 1.6e-4 / 4.8e-4 (10-20k draws; r05's 2.5k draws give ~20% MC noise) | `out/s02_fpp_fixed.csv` | upheld |
| FPP (task definition; sign fixed, f_BD 1.3e-4, BD short e law) | ~0.0013 (REF) - 0.0030 (INV), by rescaling | REF 0.0013 (Kipping) / 0.0016 (Rayleigh 0.19) / 0.0007 (circ); INV 0.0031 / 0.0033 / 0.0019; white 0.0009-0.0036; all 64 variants 0.0007-0.0036 | full r05 evaluation, `out/s02_fpp_fixed.csv`, `out/s12_summary.json` | **upheld** (now a full computation) |
| P(planet on some star of the bound system) | >= 0.997-0.999 | 0.9995 adopted (REF and INV); 0.9986-0.9996 over all 64 variants | 1 - FPP_system, `s02` | upheld (conservative) |
| PTP | 8.4e-4 (REF) - 2.5e-3 (INV) | Kipping 8.6e-4 / 2.6e-3; Rayleigh 0.19 1.2e-3 / 2.85e-3; circular 3.0e-5 / 1.15e-3; white REF 3.8e-4 | `s02`; representative companion `out/s06_evidence.json` | upheld, range widened; strongly e-prior dependent |
| Lambda_G (grazing BD - planet), beta noise | 13.6 (REF) / 6.8 (INV) | 9.2 (REF, cleaned, clean/ beta) / 6.6 (INV); white 21.0 / 7.8; k >= 0.20 bound: 9.3 / 6.7 | own fits, `out/s03_fit_real.json` | **revised** (REF noise-model dependent) |
| Grazing-BD solution | k 0.19 (bound), b 0.97, rho_c 1.8, e >~ 0.55 | k 0.19, b 0.964-0.975, rho_c 1.71-1.81, T14 2.08 h; Delta-chi2 <= 4 for rho_c ~0.95-2.55 -> e >= 0.50 at the minimum (apoastron), >= 0.42-0.62 over the range | `s03` profile | revised (e bound) |
| Planet-truth Lambda_G | 11.6 +- 5.7 (REF, n 40), 8.0 +- 5.1 (INV, n 46 strongly overlapping fake sets) | 9.9 +- 4.2 (REF, n 50; obs at quantile 0.48), 6.3 +- 3.7 (INV, n 28; quantile 0.61) | `out/s04_inject_*_planet.csv` | upheld |
| BD-truth Lambda_G | not run | REF n 60: mean 0.40, p95 2.2, max 3.2; INV n 28: mean 0.55, max 3.4; 0 trials >= observed | `out/s04_inject_*_grazBD.csv` | **new** |
| Calibrated shape LR (BD/planet) | not available | REF 0.0038 (0.004-0.015 with chi2_1 tail), INV 0.032 (0.03-0.06); formal exp(-Lambda/2) 0.010 / 0.036; r03 grids 0.05 / 0.11 | `out/s11_calibration.json` | **new**: r05 BD terms 3-13x too lenient |
| Independent evidence ratio Z(BD on M3V)/Z(Neptune on M3V) | - | REF 5.4e-5 (BD short), 2.0e-4 (long); INV 1.65e-4 / 6.3e-4; 99.4% of the BD evidence from rho_c < 3 | own grid, `s06` | new |
| BD on M3V posterior, f_BD 1.3e-4, before Gaia | ~1e-5 (REF) - 3e-5 (INV) after Gaia factor | 2.3e-6 (REF) / 6.9e-6 (INV) short law; 8.5e-6 / 2.6e-5 long law (independent); sign-fixed r05: 1.1e-5 / 3.3e-5 | `out/s12_summary.json` | revised down ~4.5x |
| All BD channels | <~ 1e-4, most likely 1e-5 to 5e-5 | after Gaia: 0.6e-6-8e-6 (REF), 7e-6-2.4e-5 (INV), short law; <= 8.5e-5 long law; envelope 1.5e-4 | `s12` | upper envelope upheld; likely range revised to ~1e-6-2e-5 |
| f_BD (3-10 d, per M dwarf) | 1.3e-4 (68% 5e-5 to 2.5e-4) | inputs verified: Bryant+23 0.194 +- 0.072% (all), **0.108 +- 0.083% for 0.26-0.42 Msun** -> ~7e-5; Csizmadia+15 FGK ~0.2% at P < 10 d, ~6x below hot Jupiters | arXiv 2303.00659, 1508.05763 (`lit/`) | upheld (possibly 2x high) |
| e of known M-dwarf BDs at P < 16 d | none > 0.36; 2 at 3-10 d (e 0.337, 0.354) | TOI-2119 b e = 0.337 (P 7.20 d, 64.4 MJ, host 0.53 Msun) and TOI-5610 b e = 0.354 (P 7.95 d, 40.4 MJ) verified; NEA: giants/BDs around M dwarfs at P < 12 d all e <= 0.38 | arXiv 2202.08842, 2501.16554; `data/nea/` | upheld |
| P(e > 0.5) | 0.08-0.09 (short), 0.41 (long) | 0.090 (short), 0.41 (long), 0.10 (Kipping) | `s12` | upheld |
| Tidal e damping, 65 MJ, Q'* 1e6 | 18 Gyr; "tides cannot exclude e ~ 0.5" | small-e 19.4 Gyr, but at e = 0.5 1.1-1.5 Gyr (Hut factor 13-17); e 0.55-0.7 -> 0.37 after 1 Gyr, 0.23 after 5 Gyr; Q'* 1e7: 0.50 after 1 Gyr | `s01` | **revised** |
| BD-tide importance | only if Q'_BD <~ 1e5 | true for 65 MJ (Q'_BD < ~3e4 Q'*); for 13 MJ, 1 R_J the BD tide dominates for Q'_BD <~ 2e6 | `s01` | minor revision |
| Gaia RVS error percentile | 47th (3,778 stars) | 46th (6,555 stars, own sample, rv_method_used = 2, grvs +-0.3, S/N +-20%); 37th of 10,455 | `out/s05_gaia_rv.json` | upheld |
| Gaia RV LR(K) | ~1 for K <= 6; 0.05-0.92 at 12; 0.003-0.46 at 16 | ~1 for K <= 6; 0.10-0.95 at 12; 0.003-0.66 at 16 (CCF width 13-20 km/s, peak, e, paired sampling) | `s05` | upheld |
| Mass-averaged Gaia factor (BD on M3V) | 0.3-0.6 | empirical BD masses 0.25-0.80 (0.41-0.57 at CCF width 16 km/s); uniform 13-80 MJ 0.47-0.90 | `s05` | upheld |
| K1 of a BD on the M3V | 2.9-15.5 km/s | 2.81 (13 MJ) - 15.6 (80 MJ) km/s circular; x1.2 at e = 0.55 | `s05`, `s10` | upheld |
| Photocentre semi-axis | 0.023-0.117 mas | 0.023 (13 MJ) - 0.126 (80 MJ) mas | `s05` | upheld |
| Companion-host k = 0.40 solution | D ~0.09, "self-consistent with ~0.9 R_J on ~0.2 Msun" | a 0.2 Msun companion has D ~0.25; D 0.09 -> ~0.14 Msun, occulter 0.66 R_J (giant planet, PTP); a BD needs k >~ 0.45 (+7.4 / +2.4) | r04 draws; analyst's own profile | **refuted** (the physics statement) |
| White-noise BD terms | "below 1e-6" | REF yes (6.6e-8); **INV 1.0e-4** (b10's own output, f_BD 1e-3) | `../bd/out/b10_r05_eccsign.json`, `s02` | **refuted** for INV |
| RV design | "2-3 epochs at 0.3-1 km/s at fixed phases" | 13 MJ grazing configuration: 2 quadrature epochs need sigma <= 0.5 km/s (99%); 3 epochs (0.2/0.5/0.8) <= 0.38; BD e law, 13 MJ: 3 epochs <= 0.24-0.3 km/s; empirical masses <= 0.7-0.9 km/s; at 1 km/s only 39-63% of 13 MJ detected | `out/s10_rv_design.json` | revised |
| Imaging vs PTP | 8-m speckle null + two-band | PTP companions: sep 0.047-0.48" (median 0.10"; 18% < 0.05"), dI 0.90-1.41, dKs 0.71-1.10, m2 0.19-0.22 Msun, a 6.7-44 au; detected by 8-m 832-nm speckle 93%, 8-m Ks AO 81%, 4-m I speckle 55%; FPP after speckle null 0.0005 (REF) / 0.0007 (INV); + two-band 5%: 1.1e-4 | `out/s07_ptp_followup.json` | upheld, quantified |

## 2. Details

### 2.1 Sign error (s01, s02)

**Sign check.** My own code integrates the Kepler orbit and projects it on the sky. The transit is at true anomaly
f = 90 deg - w, and the impact parameter is fixed at conjunction.

| (e, w) | rho_c/rho measured | g^+3 | g^-3 |
|---|---|---|---|
| (0.5, 90 deg) | 5.188 | 5.196 | 0.192 |
| (0.5, 270 deg) | 0.1925 | 0.1925 | 5.196 |

**Effect of the error.** The r05/f04 relation pairs the periastron boost in transit probability, (1 + e sin w)/(1 - e^2), with
long durations. For the Kipping e law this raises the transit-weighted probability of rho_c/rho* <= 0.2 from 0.89% to
8.0%. For the uniform 0-0.9 law it raises it from 3.2% to 41%.

**Reproduction (s02).** I imported r05 read-only and patched ecc() to return 1/g, so its rho_h*ge^-3 becomes rho_h*g^+3.

| Run | FPP | Published / b10 |
|---|---|---|
| Wrong sign, 2,500 draws, new seed (REF) | 0.00315 | 0.00313 |
| Wrong sign, 20,000 draws (REF) | 0.00299 | 0.00313 |
| Sign fixed (REF) | 0.00151 | b10 0.00144 |
| Sign fixed (INV) | 0.00359 | b10 0.00345 |

Between seeds, BD_on_M3V varies by ~20% at 2,500 draws.

**Full grid.** 4 light curves x 4 planet e priors x 3 BD e laws x 4 f_BD values (`out/s02_fpp_fixed.csv`). The FPP is
dominated by PTP and by the NTP pair-A term (3-7e-4). f_BD between 5e-5 and 1e-3 changes the FPP by only 0.00014 (REF) / 0.00043 (INV)
once the sign is fixed.

### 2.2 Light-curve shape on the cleaned curves (s03, s04, s11)

**Data.**
* REF: 21 windows, N = 4,663.
* INV: 18 windows, N = 3,508. Windows need >= 2 in-transit points and >= 10 points on each side.
* Errors: flux_err x clean/'s beta_transit.
* Baselines: per-window linear.

**Template check.** The templates match batman's supersampled exposures to <= 6 ppm.

**Best fits.**

| Fit | k | b | rho_c (g/cm3) | T14 (h) | Flat fraction |
|---|---|---|---|---|---|
| Planet, REF | 0.118 | 0.24 | 14.9 | 1.61 | 0.78 |
| Planet, INV | 0.119 | 0.69 | 5.3 | 1.83 | 0.62 |
| Grazing BD (both) | 0.19, at the bound | 0.96-0.98 | 1.7-1.8 | 2.08 | - |

Lambda_G is 9.2 (REF) and 6.6 (INV). With white noise it is 21.0 and 7.8.

**Why the analyst's REF value is higher.** bd/ quotes 13.6 for REF. It used its own fake-window beta (1.11-1.37). clean/'s
calibrated beta is 1.02-1.62. The INV values agree.

**Diluted large body.** Lambda at the k = 0.25 bound with D ~0.21-0.23 is 1.3 (REF) and -0.2 (INV). A diluted body fits
as well as a planet. Physically that is a Saturn-size planet on a ~0.2 Msun companion, which is PTP, not a BD.

**Injections.**
* Set-up:
  * truths are the best planet and the best grazing-BD fits of the real data;
  * REF uses referee_phot fake windows with the real sector mix;
  * INV uses the cleaned curve at 28 fake ephemerides spaced 0.02 in phase, so the in-transit parts are independent.
* Under BD truth Lambda_G is chi2-like with few effective dof. A scaled-chi2 fit gives s ~1.0-1.3 and nu ~0.4.
* Under planet truth it is broad: REF 9.9 +- 4.2; INV 6.3 +- 3.7, with 29% of INV trials below 4.
* The observed values sit at the planet-truth median, and above every BD-truth trial.

**Calibrated LR.** LR = f_BD(Lambda_obs)/f_planet(Lambda_obs).
* REF: 0.0038. With a chi2_1 tail instead of the fitted nu: 0.015.
* INV: 0.032 (0.056).
* For comparison, the formal exp(-Lambda/2) is 0.010 / 0.036.

**Caveat on the analyst's INV injections.** They used 48 fake ephemerides with phase spacing 0.011 (0.06 d) and ±0.3-d
windows. Neighbouring sets share most of their data, so the effective number of independent trials is far below 46.

### 2.3 Independent evidence for the grazing BD (s06)

**Grid.** I built a likelihood grid on the cleaned curves:
* k: planet 0.095-0.175, BD 0.200-0.275, diluted host 0.15-0.45;
* b: 64 steps;
* rho_c: 46 log steps over 0.5-80 g/cm3.

**Integration.** I integrated over the correct-sign photo-eccentric mapping with transit-probability weights.
* Priors: Neptunes flat in ln R over 4-7 R_E; R_BD ~ N(0.88, 0.08) R_J; rho* ~ N(9.42, 0.88).

**Result.** Z(BD on M3V)/Z(Neptune on M3V):

| e law | REF | INV |
|---|---|---|
| BD short-period | 5.4e-5 | 1.65e-4 |
| Kipping | 7.4e-5 | 2.1e-4 |
| BD long-period | 2.0e-4 | 6.3e-4 |

* 99.4-99.9% of the BD evidence comes from rho_c < 3 g/cm3, i.e. from e >~ 0.4 apoastron transits.
* Posterior with f_BD/f_N = 1.3e-4/0.003: 2.3e-6 (REF) and 6.9e-6 (INV) for the short law.
* That is ~4.5x below the sign-fixed r05 values (1.1e-5 / 3.3e-5). Most of the gap is the r03 grids' noise model, whose
  best grazing dchi2 is 5.8 / 4.3 against my 9.6 / 7.0.

**PTP check.** For a representative PTP companion (0.20 Msun, R 0.22 Rsun, rho 26.5 g/cm3, D 0.25), the evidence of a
4-10 R_E planet relative to the M3V Neptune, per unit configuration prior:

| Planet e prior | REF | INV |
|---|---|---|
| Kipping | 0.057 | 0.164 |
| Rayleigh 0.19 | 0.074 | 0.235 |
| Circular | 0.0039 | 0.0062 |

This confirms that PTP is driven by eccentric (apoastron) orbits and is 15-30x weaker for circular orbits.

### 2.4 Gaia (s05)

**What DR3 provides.** The target's rv_amplitude_robust and rv_chisq_pvalue are NaN, as expected for rv_method_used = 2.
The pair members A and B have no DR3 RV, so the systemic velocities cannot be compared.

**Smearing model.**
* Summed Gaussian CCFs over 21 transits, with the Zucker curvature error.
* Transits drawn either as independent random phases, or as Gaia-like pairs 106.5 min apart. The pairs make the smearing
  more variable and the LR weaker.

**CCF width.** For RVS (R ~11,500, sigma_LSF ~11 km/s), a template-matched CCF should have sigma ~15-16 km/s.

**Result.** At 16 km/s the mass-averaged factor for the empirical BD masses is 0.41-0.57, which I adopt as 0.5. Low-mass
BDs (K < 8 km/s, M <~ 38 MJ) are unconstrained.

### 2.5 Tides (s01)

**Method.** Jackson+08 eq. 2 in the small-e limit, multiplied by the Hut (1981) constant-time-lag eccentricity function
(1 - e^2)^-6.5 [f3(e) - 11/18 (1 - e^2)^1.5 f4(e) Omega/n] / (1 - 11/18 Omega/n). The function equals 13.3 (Omega = 0) or
17.4 (P_rot = 10 d) at e = 0.5, and 47-65 at e = 0.6.

**e-folding time for 65 MJ at e = 0.5:**

| Q'* | tau_e |
|---|---|
| 1e5 | 0.11-0.15 Gyr |
| 1e6 | 1.1-1.5 Gyr |
| 1e7 | 11-14 Gyr |

**Integrated evolution, from e0 = 0.55-0.7:**

| Q'* | e after 1 Gyr | e after 5 Gyr |
|---|---|---|
| 1e6 | 0.37 | 0.23 |
| 1e7 | 0.50-0.54 | 0.42-0.43 |

**Caveat.** Ferreira & Rice's Q* = 10^6.0 +- 0.1 (abstract verified) is calibrated within their own tidal model. If that
model is the small-e one, the effective Q to use with Hut's functions differs. The FPP uses empirical e laws, so this
does not change any posterior. It changes the qualitative claim.

### 2.6 RV design (s10), for observers

**Target.** V ~15.6 (Gaia transformation), I_C ~12.9, J 11.36, Ks 10.48, Dec -51 deg.
* Season (from ../ctx/): April-August from Chile and South Africa.

**Ephemeris.** Use the final cleaned ephemeris (../clean/out/ephemeris_final.json): T0 = 3840.91488 BTJD_TDB, P = 5.6915020 d.
Its 2027 timing error of ~2 min is negligible for RV phase scheduling.

**Amplitudes.**
* A BD on the M3V: K1 = 2.8 (13 MJ) to 15.6 (80 MJ) km/s; at the grazing configuration (e 0.5-0.6, w 250-290 deg)
  multiply by 1.2.
* The planet itself: K ~7 (10 M_E) to 21 (30 M_E) m/s.

**Recommended.** 3 epochs at phases 0.2, 0.5 and 0.8 after mid-transit, each with sigma <= 0.3 km/s. This detects
>= 98.5% of 13-MJ BDs for any e from the BD law, and >= 99.99% of the grazing configurations, at FAP 1e-3.
* With sigma <= 0.7-0.9 km/s, 99% of BDs drawn from the empirical M-dwarf BD masses are detected.
* Two epochs at quadrature (phases 0.25 and 0.75) with sigma <= 0.5 km/s suffice for the light-curve-allowed grazing
  configuration at 13 MJ (99.4%).

**Precision is not the problem.** These precisions are modest for R >~ 20,000 spectrographs on 4-8 m telescopes at
J ~11.4 / V ~15.6, but I did not run instrument exposure-time calculators. The Mira 17.8" away contributes << 1e-3 of the
light in a 1-2" fibre in normal seeing (estimated, not measured).

**What RV does not do.** It does not address PTP at 5.69 d.
* A long-term RV trend does: 74% of the posterior-weighted PTP companions give |dv/dt| > 10 m/s/yr, and 56% give
  > 30 m/s/yr.
* So a planet-mass campaign at <= 5 m/s spanning a year tests PTP for free in most configurations.

## 3. Missed tests and open issues

1. **BD-truth calibration of the shape statistic.** This was the step the task asked for and the one that matters. It
   is now done (s04/s11). It lowers the BD channels 3-13x below the r05 grid values.
2. **Gaia DR4, 2 December 2026, 66 months, with epoch astrometry and epoch RVs/RVS spectra.** ESA announcement, via a web
   search on 2026-10-17; whether epoch RVs will cover grvs = 12.5 stars is unverified. It arrives six weeks from now, at no
   cost.
   * Photocentre orbit. A 60-80 MJ BD gives 0.10-0.13 mas. With ~2x more CCD transits than DR3 at sigma_CCD ~0.2 mas, a
     5.69-d orbit fit could reach S/N ~5-9 for massive BDs (rough estimate).
   * Epoch RVs would constrain K >~ 10 km/s directly.
   * The DR4 acceleration and RUWE also test the PTP companion: a ~0.2 Msun companion at 7-44 au.
3. **Physical self-consistency of companion-hosted BDs.** D and k must follow from the companion mass. In the analyst's
   shape profile this was not imposed. r05 does impose it.
4. **The planet e prior for PTP.** PTP changes by 15-30x between circular and eccentric priors. M-dwarf-specific laws
   exist and are verified:
   * Sagear & Ballard 2023, arXiv:2305.17157: singles Rayleigh sigma 0.19.
   * Sagear+25, arXiv:2507.07169: higher e for R > 3.5 R_E.
   * They should be the adopted prior, not Kipping (2013), which comes from RV giant planets. The FPP barely changes
     (0.0016 against 0.0013 REF).
5. **Planet occurrence on a ~0.2 Msun host.** f_N and beta_N for 6 R_E planets at 5.7 d around 0.2 Msun stars set PTP
   linearly. They are carried from referee2_fpp and not re-verified.
6. **NTP pair A** is now the second-largest FP term (3-7e-4). It was not re-derived by bd/ or here.
7. **Long-term RV trend as a PTP test** (s07). Not mentioned by bd/.
8. **The INV-REF discrepancy** (b 0.69 vs 0.24, rho_c 5.3 vs 14.9) still drives a ~2.3x FPP difference. bd/ and clean/
   both note it; a joint treatment of the two photometries would remove it.

## 4. Scripts and outputs (Python: /home/user/work/venv/bin/python, OMP/MKL/OPENBLAS threads = 1)

| Script | What it does | Outputs |
|---|---|---|
| `scripts/rlib3.py` | constants, final cleaned ephemeris (T0 3840.91488 BTJD_TDB, P 5.6915020 d), helpers | - |
| `scripts/s01_sign_tides.py` | (a) photo-eccentric sign check with my own Kepler-orbit sky-projection code (no batman); (b) transit-weighted P(rho_c/rho* <= x) for 5 e priors, correct and wrong sign; (c) tidal e-damping: Jackson+08 small-e rate x Hut (1981) eccentricity functions, and e(t) integration | `out/s01_sign_tides.json`, `logs/s01.log` |
| `scripts/s02_fpp_fixed.py` | imports r05 read-only; my sign patch (10,000 e draws); reproduction of r05 as published and of bd/b10; full grid 4 light curves x 4 planet e priors (incl. Sagear & Ballard 23 Rayleigh 0.19) x 3 BD e laws x 4 f_BD | `out/s02_fpp_fixed.{csv,json}`, `logs/s02.log` |
| `scripts/lc3.py` | own likelihood engine: batman on a fine grid + exact box exposure integration (validated < 6 ppm vs batman supersampling, `logs/t_lc3.log`), per-window linear baselines profiled, template bank + Nelder-Mead fits | - |
| `scripts/t_lc3.py` | validation of lc3 | `logs/t_lc3.log` |
| `scripts/s03_fit_real.py` | fits of the cleaned REF/INV curves (beta and white): planet, grazing BD (k >= 0.19 and >= 0.20), diluted large body; rho_c profile of the grazing BD | `out/s03_fit_real.json`, `logs/s03.log` |
| `scripts/s04_inject.py` | fake-epoch injections, truths = best planet and best grazing-BD fits of the real data; refits both scenarios | `out/s04_inject_{ref,inv}_{grazBD,planet}.csv`, `logs/s04_*.log` |
| `scripts/s05a_gaia_query.py` | Gaia DR3 queries: target + neighbours; independent RV-error comparison sample (10,455 stars) | `data/gaia_target.csv`, `data/gaia_rvcomp.csv` |
| `scripts/s05b_gaia_rv.py` | own RVS combined-CCF smearing likelihood ratio (21 independent or Gaia-like paired transits; CCF width 13/16/20 km/s, peak 0.3/0.6; e 0 and 0.55), mass-averaged; photocentre semi-axes | `out/s05_gaia_rv.json`, `logs/s05b.log` |
| `scripts/s06_evidence.py` | own evidence ratio BD-on-M3V / Neptune-on-M3V from a (k, b, rho_c) likelihood grid of the cleaned curves; representative PTP companion (0.20 Msun, D 0.25) per e prior | `out/s06_evidence.json`, `logs/s06.log` |
| `scripts/s07_ptp_followup.py` | sign-fixed r05: posterior-weighted PTP companion properties, fraction detectable by three imaging set-ups, RV trend, FPP what-ifs | `out/s07_ptp_followup.json`, `logs/s07.log` |
| `scripts/s09_fetch_abs.py` | arXiv abstract fetcher (verification of literature numbers) | `lit/abs_*.xml`, `logs/s09*.log` |
| `scripts/s10_rv_design.py` | RV epochs/precision needed to detect any BD on the M3V (grazing configuration and BD e law; 13 MJ and empirical masses) | `out/s10_rv_design.json`, `logs/s10.log` |
| `scripts/s11_calibrate.py` | calibrated shape likelihood ratio from the s04 injections | `out/s11_calibration.json` |
| `scripts/s12_summary.py` | combination: sign-fixed FPP rows, Gaia RV factor, shape calibration, independent BD estimate, P(e > 0.5) | `out/s12_summary.json` |
| `scripts/lane_a.sh`, `scripts/lane_b.sh` | job queues | `logs/lane_*.log` |

Other files: `data/nea/*.csv` (NASA Exoplanet Archive queries: transiting companions > 10 MJ and giants around M
dwarfs), `lit/` (arXiv abstracts used).
