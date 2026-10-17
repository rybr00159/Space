# False-positive probability of the 5.69-d signal on TIC 214311908 (`fpp`)

Folder: `tic214311908/analysis/fpp/`
(AI-run analysis, 2026-10-17. Nothing was submitted, uploaded or posted anywhere, and no one was contacted.)

Two independent calculations: (1) an explicit Bayesian scenario comparison written here (`f01`-`f04`), using the
phase-1/phase-2 results with their referee corrections; (2) TRICERATOPS (`f05`-`f07`) on the deblended light curve,
run 59 times in 10 configurations.

## Bottom line

* **Bayesian scenario comparison (baseline): FPP = 0.0111, NFPP = 6.3e-04.** P(planet on the M3V or
  its twin) = 0.9889: on a single M3V 0.671, on the primary of an unresolved binary
  0.186, on a near-equal twin companion 0.132.
* **Dominant residual false positive: PTP, a planet transiting a fainter unresolved bound companion**
  (posterior 0.0104, 93 % of the FPP). Its typical configuration is a
  0.21 M_sun companion (16-84 %: 0.20-0.22), Delta T ~ 0.75-1.5 mag, at
  0.05-0.64" (3-40 au; median 0.14"),
  hosting a 5.8 R_e planet. That is still a planet, on the wrong star: counting any
  planet-sized body in the bound system as "planet", FPP_system = 7.2e-04. Next: NTP on pair A
  (6.3e-04, set by an extrapolated pixel likelihood ratio), hierarchical EB (5.1e-05),
  BTP on an uncatalogued star (3.1e-05), WD occulted by the M3V (1.2e-05).
  The two catalogued background stars contribute < 1e-8; the Mira 0 (physically impossible host; pixels log10 LR = -26).
* **Sensitivity.** Over 3000 joint draws of every uncertain prior and likelihood choice: FPP median 0.0067,
  16-84 % 0.0024-0.0157, 5-95 % 0.0008-0.0245;
  **82 % of draws give FPP < 0.015**; NFPP median 1.7e-05 (95th percentile
  6.7e-04; 99 % < 1e-3). 16 of 96 named
  variants reach FPP >= 0.015, e.g. no archival imaging/astrometric constraints (0.045), a twin
  boundary of Delta T 0.5 instead of 0.75 (0.024), or every input at its plausible FP-favouring end
  (0.053). The FPP is most sensitive to the twin boundary,
  the CMD scatter sigma_s, the host-mass scaling of Neptune occurrence and the planet suppression in close binaries.
* **Definition matters.** If a planet on a near-equal twin also counted as a false positive (the TRICERATOPS convention for
  "STP"), FPP_strict = 0.143 (draws: median 0.13, 5-95 % 0.03-0.33),
  because the -0.75 mag over-luminosity makes an unresolved twin likely: P(binary) = 0.47,
  P(twin, Delta T <= 0.3) = 0.28.
* **TRICERATOPS** (59 runs in 10 configurations). Main configuration (12 seeds, N = 2e5): FPP = 0.0123 +- 0.0087 (median 0.0116, range 0.0017-0.0277), NFPP = 1e-05;
  67 % of seeds below 0.015. The run-to-run scatter is large (rare high-likelihood companion draws). Its FPP comes
  from STP (planet on an unresolved companion; the same scenario as our PTP) and from x2P eclipsing binaries whose second
  eclipse TRICERATOPS does not model; with circular EB orbits FPP = 0.0098 +- 0.0122 (median 0.0035, range 0.0002-0.0330), all STP. With every TIC neighbour kept
  (TRICERATOPS' aperture treatment) FPP = NFPP = 0.75 (pair A, two bright stars 33" and 53" away); these stars
  are excluded by the pixel localisation at Delta chi2 >= 126-592; weighting each by its pixel likelihood ratio gives NFPP =
  1.0e-03 (seeds 7e-05-4e-03) with the heavy-tailed pair-A ratio 10^-3.2, or
  5e-07 with the Gaussian one (10^-6.7): almost all of it is a planet on pair A.
* **Not "statistically validated" by community standards now.** The FPP sits at, not safely below, the Giacalone et al.
  (2021) threshold (Bayesian baseline 0.011; 18 % of plausible prior choices and
  33 % of TRICERATOPS seeds exceed 0.015); there is no high-resolution imaging (the IPD/DECaPS2 "contrast curve" is
  a proxy); the host may be an unresolved twin; and the depth comes from deblended photometry 18" from a Mira.
  NFPP < 1e-3 is met in the baseline (6.3e-04) but rests on an extrapolated pixel likelihood ratio for pair A
  (TRICERATOPS with all neighbours, pixel-weighted: 1.0e-03).
* **What would most reduce the FPP:** a ground-based, seeing-limited transit observed simultaneously in a blue and a red
  band (one 1-m class telescope, 2027 April-August; 17 observable events). It places the dip on the M3V (resolving pair A and
  the Mira) and tests achromaticity: a planet on a ~0.21 M_sun companion would give depth(V)/depth(I) = 0.80
  (0.72-0.83). With equal depths to 5 %, FPP 0.0111 -> 1.2e-04
  (realistic-pessimistic priors: 0.053 -> 0.0017). 8-m speckle imaging alone gives
  0.0017 (0.021) and is in any case expected for a validation;
  high-resolution spectroscopy (SB2/RV) is needed to settle the twin question (FPP_strict 0.14 -> 0.017 with all three).

## 1. Inputs (all read-only; provenance)

| Input | Value used | Source |
|---|---|---|
| Light curves | INV = phot model B with the 0.6485-d RRab leak removed (`search/data/m3v_clean.npz`, `det_clean`, 979 e-/s scale), flares removed, sigma x white scale x red-noise beta (`shape/out/noise.json`); REF = referee local-window photometry (`shape/data/ref_win.npz`) | phot, search, referee_phot, shape |
| Ephemeris | T0 = 3840.9155 BTJD_TDB, P = 5.6915087 d | shape/out/final_params.json, phot/out/transits.json |
| Depth | 15.3 e-/s of 979 e-/s = 1.56 % (LD fit 15.8-15.9 ppt) | phot, shape |
| Star | R 0.381, M 0.370, rho 9.42 +- 0.88 g/cm3, M_Ks 6.4406; twin 2 x (0.289 R_sun, 0.264 M_sun, rho 15.45) | star/star_params.json |
| Over-luminosity | Delta M_Ks = -0.753 mag (referee) | referee_star/out/r09_cmd.json |
| Pixel localisation | 1000 injections per source through the referee's image pipeline; real combined Delta chi2; real 0.5" sky map | referee_loc/out/inj_base.npz, fit_base.json |
| Gaia IPD/RUWE companion limits | LR per (separation, Delta I) cell | star/out/ipd_limits.json |
| Astrometric limits | L_combined(M2, a) (RUWE + DR3-DR2 PM) | star/out/astrometric_limits.json |
| DECaPS2 contrast | referee-corrected 5-sigma: i 4.24/5.07/5.80/5.58, z 3.96/4.53/4.96/5.24 at 1/1.25/1.5/2"; investigator beyond 2.5" | referee_star/out/r04_decaps_inject.json, star/out/decaps_contrast.json |
| DSS look-back | no unbound star T <~ 18.3 within 3" of the 2019-2026 position (except near catalogued stars) | referee_star/out/r03_dss_limits.json |
| Star counts | 0.00123 / 0.00209 / 0.00368 arcsec^-2 for T <= 16 / 16.72 / 17.47 | star/out/blend_density.json |
| Physical EB / blend fits | Delta chi2 vs the planet for every EB class (WD occultation, hierarchical EB, BD, background EBs), incl. secondary, photometry, parallax, RUWE | shape/out/summary.json, scenario_table.csv |
| Occurrence rates | f(4-7 R_e, P 3-10 d) ~ 0.1-1 % (best 0.3 %); giants 0.1-0.3 %; Winters+19 multiplicity 26.8 %, log-normal a (peak 20 au, sigma 1.16 dex); Shan+15 0.11 close companions per M dwarf | ctx/README.md, referee_ctx |

## 2. Bayesian scenario comparison (`f01`-`f04`)

### 2.1 Formula
For every scenario S, E_S = prior(S exists, per unit ln P at 5.69 d, per target) x likelihood of all data under S,
integrated over S's internal parameters. Every E_S uses the same period normalisation and the same light-curve reference
chi2 (the best planet fit on the grid). Posterior P(S) = E_S / sum E.

* **Planet scenarios** use a light-curve likelihood grid (`f02`): dark body of radius ratio k on a host supplying a fraction
  D of the deblended flux, with circular-equivalent density rho_circ and impact parameter b (15 D x 26 rho x 39 b; k
  optimised with a Laplace width; per-window offsets profiled; batman with exposure integration; both photometries).
  The evidence integrates over b (prior density (R/a)(1+e sin w)/(1-e^2)), over k with the occurrence prior per ln R
  (classes 2-4 / 4-7 / 7-10 / 10-22 R_e, host-mass scaling (M/0.37)^beta), and over e, w (Kipping 2013 beta prior;
  rho_circ = rho_host g^-3, photo-eccentric effect). Grid sanity checks reproduce the shape agent: best b 0.55 (0.53 +-
  0.10), k 0.118 (0.119), ln B(single/twin) = -0.17 (e prior) / -0.60 (circular) against their -0.43/-0.53.
* **Bound-companion configurations** (`f03`): 400,000 Monte-Carlo companions (Winters+19 separations, q uniform
  re-weighted in f04 to power law q^gamma plus a twin excess F_twin at q > 0.95; isotropic orbits; primary mass solved so that
  the total M_Ks = 6.4406). Each draw carries the likelihoods of the Gaia IPD cells, the astrometric grid (mapped by the
  photocentre wobble |B - beta|, which vanishes for twins; the phase-1 grid assumed a 0.37 Msun primary), the DECaPS2
  limits and the CMD over-luminosity: the companion explains dM_bin(q) of the -0.753 mag offset and the rest must be the
  primary's own offset eps ~ N(0, sigma_s), sigma_s = 0.30 mag. (The empirical offsets of ~1000 RUWE<1.1 M dwarfs at the
  same colour have a compact bright side, 16th percentile -0.27 to -0.34 mag, and a long metal-poor faint tail, so the
  faint side could not be mirrored; sigma_s is varied 0.25-0.50.) Planet suppression in binaries S(a) = S_close (0.25)
  for a <= 10 au rising log-linearly to 1 at 200 au; Holman & Wiegert (1999) S-type stability.
* **EB scenarios** use the shape agent's maximum-likelihood Delta chi2 (dataset-matched), times an Occam factor
  f_occ = 0.1 (the fraction of the class's parameter space that reproduces the depth and duration; generous compared with
  the planet's own ~0.05), times literature rates and eclipse probabilities.
* **Background / neighbour scenarios** add the injection-calibrated localisation likelihood ratios and the archival
  imaging limits.

### 2.2 Localisation likelihoods (`f01`, `out/loc_likelihoods.json`)
Statistic s = chi2(X) - chi2(M3V) from the referee's combined image fits; LR = p(s_obs | X)/p(s_obs | M3V) with both
densities fitted to the 1000 injections (Gaussian and Student-t; the t-fit is adopted where it is heavier-tailed).

| Source | s_obs | s under X true (2.5/50/97.5/99.85 %) | P(s <= s_obs | M3V) | log10 LR (Gauss / t) |
|---|---|---|---|---|---|
| pair A | +126.2 | -201 / -114 / -58 / -35 | 0.12 | -6.74 / **-3.22** |
| pair B | +145.5 | -212 / -122 / -66 / -35 | 0.15 | -10.4 / -10.4 |
| Mira | +284.8 | -300 / -201 / -127 / -87 | 0.14 | -26.3 / -26.3 |
| G17.5 star (3.1") | +35.8 | -83 / -44 / -7 / +6 | 0.31 | -3.88 / -3.87 |
| G17.7 star (5.5") | +86.6 | -157 / -100 / -47 / -11 | 0.37 | -10.2 / -10.2 |
| G17.0 star (12") | +189.3 | -265 / -170 / -106 / -74 | 0.13 | -17.2 / -17.2 |

0/1000 injections on any alternative reach the observed value, so all |log10 LR| beyond 3 are extrapolations; the
sensitivity analysis uses -1.4 to -4.2 for the 3.1" star (referee variants; loc AMP gives -1.63) and -2 to -6.7 for pair A.
Uncatalogued source: effective area A_eff = sum exp(-Delta chi2_map/2 kappa) dA = 18.2 arcsec^2 at kappa = 2.2 (the
referee's 95 % injection threshold / Gaussian), 17.6-26.7 for kappa 1.5-4.8: 2.1 (r < 1"), 2.5 (1-1.25"), 10.0
(1.25-3"), 3.4 (3-5"), 0.1 arcsec^2 (5-8").

### 2.3 Priors and likelihood terms of every scenario (baseline)

| Scenario | Prior ingredients (baseline; range in sensitivity) | Likelihood ingredients | Posterior |
|---|---|---|---|
| TP single | (1 - f_mult) x CMD(single) x f(4-7 R_e) = 0.3 % (0.1-1 %) per 3-10 d | LC grid, rho = 9.42 | 0.6708 |
| TP diluted (planet on the primary of an unresolved binary) | f_mult 0.27 (0.20-0.45) x Winters a, q x IPD x astrometry x DECaPS2 x CMD x S(a) x stability | LC grid at D = 1 - D_c, rho_1 | 0.1859 |
| TP on twin (planet on a companion with Delta T <= 0.75) | same, companion-hosted | LC grid at D = D_c, rho_2 | 0.1322 |
| **PTP: planet on a fainter bound companion (Delta T > 0.75)** | same; occurrence by planet radius and companion mass (beta_N 0, beta_J 1.5) | LC grid at D = D_c, rho_2 | **0.0104** |
| EB: WD totally occulted by the M3V | WD+M dwarf rate 2e-5 per ln P (2e-6-2e-4; stress 2e-3, 1e-2) x f(T_WD) 0.3 x p_ecl 0.029 x f_occ 0.1 | Delta chi2 13.2 (LC + secondary + RUWE; stress 6.7-8) | 1.20e-05 |
| EB: WD occulted by one twin component | same x 2 x P(twin config) | Delta chi2 17.9 (INV) / 11.1 (REF) | 1.02e-06 |
| EB: BD or M dwarf eclipsing the M3V (also 2P) | 0.024 close companions per ln P x 0.05 x f_occ | Delta chi2 >= 45 (INV) / 84 (REF) | 1.43e-08 |
| HEB: bound companion that is itself an EB | P(companion able to supply the flux) x 0.024 x 0.06 x f_occ | Delta chi2 27.5 (INV) / 33.7 (REF) | 5.06e-05 |
| BEB on the 3.1" star | field EB rate 0.0017 per ln P x (MS+MS 0.99, compact 0.01) x f_occ | Delta chi2 65.5 (MS) / 19.6 (pre-ELM) incl. photometry + parallax; pixels 10^-3.87 | 8.88e-09 |
| BEB on the 5.5" star | as above | Delta chi2 76.3 / 27.4; pixels 10^-10.2 | 8.86e-17 |
| BTP on an uncatalogued star | stars T <= 16 in A_eff x miss prob. (DSS look-back eps 0.02 within 3", DECaPS2, Gaia/IPD) x M-dwarf fraction 0.2 x hot-Jupiter rate | planet shape fits (k < 0.43) | 3.11e-05 |
| BEB uncatalogued (M dwarf + hot WD, fits the LC perfectly) | stars T <= 17.47 in A_eff x miss prob. x 0.2 x WD+M rate x 0.1 x 0.029 x f_occ | Delta chi2 0 | 1.79e-07 |
| BEB uncatalogued (MS + MS) | as above x 0.0017 | Delta chi2 >= 30 | 8.01e-09 |
| NTP pair A (5.2 R_e planet on the bound 0.30 M_sun pair star; S(53 au)) | as TP | LC grid at rho_A; pixels 10^-3.22 | 6.27e-04 |
| NEB pair A | WD rate x ... | Delta chi2 6.7; pixels 10^-3.22 | 2.61e-07 |
| NTP / NEB pair B (needs 12.5 % depth) | | pixels 10^-10.4 | 1.90e-12 / 1.70e-14 |
| Mira | a 5.69-d orbit (a ~ 15 R_sun) is inside the Mira (hundreds of R_sun) | pixels 10^-26 | 0 |

Configuration posterior (non-transit data only weighted by the transit evidences; baseline): single star
0.53; unresolved bound companion 0.47, of which near-equal (Delta T <= 0.3)
0.28, Delta T 0.3-2 0.17, Delta T 2-4.5 0.020.
The over-luminosity drives this: at sigma_s = 0.30 mag a single star at -0.753 mag is 23x less probable than an equal twin
(`out/companions_summary.json`: p_single/p_twin = 0.043). The archival constraints (IPD x astrometry x DECaPS2) keep 17 % of
the Winters+19 companion prior; with the CMD likelihood, companions at 6-65 au and < 0.15 au are the most favoured.

### 2.4 Results

| Quantity | Value | Method | Uncertainty | File |
|---|---|---|---|---|
| FPP (planet on M3V or twin, Delta T <= 0.75) | **0.0111** | explicit Bayesian evidences, INV light curve | sensitivity: 16-84 % 0.0024-0.0157 | out/fpp_baseline.json, out/fpp_random.csv |
| FPP with REF light curve | 0.0081 | same, referee photometry | REF draws 16-84 %: 0.0017-0.0132 | out/fpp_variants.csv |
| NFPP | **6.3e-04** | pair A NTP with t-fit pixel LR 10^-3.22 | draws 95th pct 6.7e-04; LR 10^-6.7 (Gaussian) -> 2e-07; LR 10^-2 (stress) -> 0.010 | same |
| FPP_strict (twin = FP) | 0.143 | same | draws 5-95 % 0.034-0.331 | same |
| FPP_system (any planet in the bound system) | 7.2e-04 | same | draws 95th pct 4.2e-03 | same |
| Dominant residual | PTP 0.0104: 5.8 R_e planet on a 0.21 M_sun companion, Delta T 1.04, 0.14" | posterior-weighted companion draws | by separation: < 4 au 2.3e-04, 4-6 au 1.2e-03, 6-65 au 7.0e-03, > 65 au 2.0e-03 | out/fpp_baseline.json, out/ptp_weights.npz, figs/ptp_map.png |

* The PTP probability sits where the archival data are weakest (`figs/ptp_map.png`): inside 0.08" (no IPD calibration,
  no DECaPS2 sensitivity, astrometry blind to near-equal companions) and at 0.1-0.6" with Delta I ~ 1, where the IPD cells
  only disfavour companions by factors 5-15 (LR 0.06-0.2).
* Planets on companions need Rp = 5.7-6.1 R_e, still in the Neptune class, and the companions' higher density (rho ~ 15-25
  g/cm3) is compatible with the measured duration once moderate eccentricity is allowed: with circular orbits the PTP share
  falls to 0.0033 (FPP 0.0042).
* The EB scenarios that fit the light curve as well as a planet (a WD totally occulted by an M dwarf) are suppressed by
  their rates and by the secondary-eclipse + RUWE terms (bound) or by photometry/parallax/pixels (catalogued background
  stars) and the DSS look-back (uncatalogued unbound stars within 3").

### 2.5 Sensitivity analysis (`out/fpp_variants.csv`, `out/fpp_random.csv`, `figs/variants.png`, `figs/fpp_distribution.png`)

**Named one-at-a-time variants** (96 rows incl. stress tests and what-ifs): FPP 7.6e-05 -
0.42. Selected:

| Variant | FPP | FPP_strict | NFPP |
|---|---|---|---|
| baseline | 0.0111 | 0.143 | 6.3e-04 |
| data REF | 0.0081 | 0.127 | 6.9e-04 |
| e prior Van Eylen | 0.0128 | 0.143 | 6.1e-04 |
| e circular | 0.0042 | 0.177 | 8.3e-04 |
| f_N=0.001 | 0.0114 | 0.143 | 6.3e-04 |
| f_N=0.01 | 0.0110 | 0.143 | 6.3e-04 |
| beta_N=-1.0 | 0.0171 | 0.189 | 6.9e-04 |
| beta_N=2.0 | 0.0046 | 0.078 | 4.9e-04 |
| f_mult=0.2 | 0.0085 | 0.109 | 6.0e-04 |
| f_mult=0.45 | 0.0173 | 0.226 | 7.1e-04 |
| F_twin=0.0 | 0.0130 | 0.119 | 6.1e-04 |
| F_twin=0.3 | 0.0080 | 0.184 | 6.5e-04 |
| S_close=0.05 | 0.0080 | 0.072 | 6.7e-04 |
| S_close=1.0 | 0.0167 | 0.270 | 5.5e-04 |
| sig_s=0.25 | 0.0176 | 0.291 | 7.5e-04 |
| sig_s=0.4 | 0.0051 | 0.050 | 5.5e-04 |
| IPD LR84 | 0.0134 | 0.156 | 6.2e-04 |
| IPD off | 0.0304 | 0.273 | 5.8e-04 |
| astrometry soft | 0.0128 | 0.159 | 6.5e-04 |
| astrometry off | 0.0249 | 0.192 | 7.3e-04 |
| DECaPS2 off | 0.0351 | 0.295 | 4.1e-04 |
| IPD+astrometry+DECaPS2 off | 0.0448 | 0.350 | 4.6e-04 |
| dT_twin=0.5 | 0.0237 | 0.143 | 6.3e-04 |
| dT_twin=1.0 | 0.0054 | 0.143 | 6.3e-04 |
| r_wd=0.0002 | 0.0112 | 0.143 | 6.3e-04 |
| r_wd=0.002 | 0.0124 | 0.144 | 6.5e-04 |
| r_wd=0.01 | 0.0177 | 0.149 | 7.5e-04 |
| dchi2_wd=6.7 | 0.0114 | 0.144 | 6.3e-04 |
| f_occ=1.0 | 0.0117 | 0.144 | 6.3e-04 |
| dchi2_heb=20.0 | 0.0132 | 0.145 | 6.3e-04 |
| lr_bg3=10^-1.40 | 0.0111 | 0.143 | 6.3e-04 |
| lr_pairA=10^-2.00 | 0.0206 | 0.152 | 1.0e-02 |
| lr_pairA=10^-6.74 | 0.0105 | 0.143 | 1.9e-07 |
| kappa=4.8 | 0.0111 | 0.143 | 6.3e-04 |
| eps_dss=0.3 | 0.0115 | 0.144 | 6.3e-04 |
| OPTIMISTIC combined | 0.0015 | 0.011 | 4.0e-04 |
| PESSIMISTIC realistic (each input at its plausible FP-favouring end) | 0.0527 | 0.421 | 2.9e-03 |
| EXTREME stress (all FP-favouring ends at once, incl. r_wd=2e-3, f_occ=1) | 0.4239 | 0.509 | 1.4e-02 |

The EXTREME stress row combines every FP-favouring end at once, including a WD+M dwarf rate 100x the PCEB-based estimate
and no Occam penalty for EBs; it is dominated by "WD occulted by the M3V" and is not a plausible prior.

**Joint random draws** (3000; `out/fpp_random.csv`, `figs/fpp_distribution.png`): every uncertain input drawn at once
(light curve INV/REF; e prior Kipping/Van Eylen/circular; f(4-7 R_e) 0.1-1 %; f(7-10) 0.01-0.3 %; f(10-22) 0.05-0.5 %;
host-mass slopes; f_mult 0.22-0.40; separation and q distributions; twin excess 0-0.3; S_close 0.05-1; sigma_s 0.25-0.40;
IPD LR or 84 % bound; astrometry nominal or softened; twin boundary 0.5-1.0 mag; WD+M rate 2e-6-2e-4; WD Delta chi2 8-15.2;
f_occ 0.03-0.3; HEB Delta chi2 20-33.7; pixel LRs of the 3.1" star 10^-4.2-10^-1.4 and pair A 10^-6.7-10^-3; kappa;
DSS miss probability 0.005-0.1). FPP percentiles 5/16/50/84/95: 0.0008, 0.0024, 0.0067, 0.0157, 0.0245; P(FPP < 0.015) = 0.82;
NFPP 5/50/95: 3.4e-07, 1.7e-05, 6.7e-04; P(NFPP < 1e-3) = 0.99. Dominant FP:
PTP_companion 2739, HEB_bound 140, NTP_pairA 43, BTP_uncatalogued 34, EB_WD_on_M3V 29, EB_WD_on_twin 15. Spearman correlation of FPP with the inputs: twin boundary
-0.45, sigma_s -0.38, beta_N -0.33, S_close
+0.20, f_mult +0.12, f_N -0.09; all others |rho| < 0.1.

### 2.6 What would reduce the FPP (`WHAT-IF` rows of `out/fpp_variants.csv`)

Null results assumed. Speckle/AO: typical 8-m 5-sigma curve at 832 nm (Delta I 2.0/3.0/4.0/4.5/5.0/5.5/6.0 at
0.03/0.05/0.1/0.2/0.5/1/2"). Multicolour: depth(V)/depth(I) = 1 to 5 % (companion-hosted planets give 0.72-0.83; Mamajek
V and Ic light fractions). Spectra: no SB2 or RV variation, excluding Delta T < 2 companions inside 3 au. Ground transit
on target: the multicolour test plus pair A/B resolved (11.7").

| Follow-up (null result) | FPP (baseline) | FPP_strict | NFPP | FPP (realistic-pessimistic priors) |
|---|---|---|---|---|
| none | 0.0111 | 0.143 | 6.3e-04 | 0.0527 |
| speckle/AO imaging null | 1.7e-03 | 0.043 | 6.1e-04 | 0.0213 |
| multicolour transit depths equal to 5% | 7.9e-04 | 0.142 | 6.7e-04 | 0.0053 |
| high-res spectra: no SB2/RV (dT<2, a<3 au) | 1.2e-02 | 0.123 | 5.5e-04 | 0.0632 |
| imaging + multicolour | 7.0e-04 | 0.042 | 6.2e-04 | 0.0064 |
| imaging + multicolour + spectra | 5.6e-04 | 0.017 | 5.2e-04 | 0.0051 |
| ground multicolour transit on target (chromatic test + pair resolved) | 1.2e-04 | 0.142 | 1.1e-08 | 0.0017 |
| ground multicolour transit on target + speckle imaging | 7.6e-05 | 0.042 | 1.0e-08 | 0.0019 |

* **Most effective single observation: a ground-based multi-band transit on the M3V.** It removes PTP through the
  chromatic test and the pair-A NTP by resolving the pair; the residual is then the hierarchical EB and uncatalogued-star
  scenarios (~1e-4). Precision needed: depth ratio to ~5 %, i.e. ~4-6 events in two bands with a 1-m telescope (per-event
  S/N 11-25 for a 1.5 % depth, ctx/out/followup.json).
* **Speckle/AO imaging** removes most PTP configurations (those at >= 0.03") and is the standard requirement for
  validation, but leaves the < 0.03" (< 2 au) twin/companion space.
* **High-resolution spectroscopy** barely changes FPP (the PTP mass lies at 3-40 au) but is the only way to settle the
  twin question and thus FPP_strict (with imaging + multicolour + spectra: 0.017); it also measures [Fe/H] (the
  alternative explanation of the over-luminosity) and an RV mass.

## 3. TRICERATOPS (`f05`-`f07`, `out/tri_summary.json`, `out/tri_runs.csv`, `figs/tri_runs.png`)

### 3.1 Setup
* Light curve: INV (phot model B, RRab removed), sectors 66/93/103/104 (200 s), folded on the ephemeris, binned at 200 s
  in phase (216 bins, |dt| < 0.25 d), exposure time 200 s with 7 sub-samples; flux error = rms of the out-of-transit bins
  (0.0043; the propagated bin errors overestimate the scatter, chi2/dof 0.6). Variants: REF (E1), all 600-s-or-shorter
  sectors binned at 600 s (E2), propagated errors (E4). S12 (30 min) is excluded: TRICERATOPS takes one exposure time.
* Depth 0.0157 (on-target; dilution_corrected = True). Target overrides: M 0.370, R 0.381, Teff 3350 K, plx 15.549 mas.
* Apertures: per sector the pixels whose centres lie within 3 px of the M3V (26-29 px), the M3V part of the phot fit
  region (`p05_phot.py`, RFIT = 3); E5 uses 1.5 px.
* Pixel data: the phase-1 TESScut cubes served through a stand-in for `lightkurve.search_tesscut` (no download).
* Neighbours: TIC 1515709367 is a TIC duplicate of pair A (same Gaia DR3 source) and is always removed. "loc" keeps only
  the 3.1" (TIC 1515709327) and 5.5" (TIC 1515709330) stars, the only catalogued stars inside the 99.7 % pixel-localisation
  region; "keep" (D, E5) keeps all 55 TIC stars that could host the dip.
* Contrast curve (`data/tri/cc_cons.txt`, filter TESS): **not a high-resolution imaging contrast curve.** It is the
  conservative (running-minimum) envelope of the Gaia IPD/RUWE calibration cells (contrast at which companions are
  disfavoured with LR_84 < 0.2: 0 inside 0.4", 1 at 0.4-0.6", 2 at 0.6-0.9", 3 at 0.9-1") and the referee-corrected
  DECaPS2 5-sigma limits (4.24 at 1", 5.07 at 1.25", 5.6-5.8 at 1.5-2", 6.2-8.6 at 2.5-6"). The optimistic version
  (`cc_opt.txt`, E3) uses the running maximum (2-3 mag from 0.08").
* N = 2e5 draws (3e4 for the runs with all neighbours), parallel = True; seeds 1-12.

### 3.2 Results (`out/tri_summary.json`)

| Group | Configuration | runs | N | FPP mean +- sd (median; min-max) | NFPP | FPP without x2P | FPP_system | NFPP with pixel LRs | runs < 0.015 |
|---|---|---|---|---|---|---|---|---|---|
| A | main: INV 200-s, conservative contrast curve, neighbours inside the localisation region | 12 | 2e+05 | 0.0123 +- 0.0087 (0.0116; 0.0017-0.0277) | 9.9e-06 | 0.0072 | 0.0053 | 1.3e-09 | 67 % |
| F | A with circular EB orbits | 10 | 2e+05 | 0.0098 +- 0.0122 (0.0035; 0.0002-0.0330) | 1.3e-11 | 0.0098 | 0.0000 | 1.7e-15 | 70 % |
| C | A without contrast curve | 10 | 2e+05 | 0.0208 +- 0.0129 (0.0216; 0.0031-0.0416) | 9.6e-06 | 0.0144 | 0.0067 | 1.3e-09 | 40 % |
| G | C with circular EB orbits | 3 | 2e+05 | 0.0029 +- 0.0036 (0.0009; 0.0008-0.0071) | 1.3e-11 | 0.0029 | 0.0006 | 1.7e-15 | 100 % |
| E1 | REF light curve | 3 | 2e+05 | 0.0075 +- 0.0040 (0.0092; 0.0029-0.0104) | 4.8e-08 | 0.0069 | 0.0006 | 6.4e-12 | 100 % |
| E2 | INV incl. S39, 600-s bins | 3 | 2e+05 | 0.0049 +- 0.0076 (0.0009; 0.0001-0.0136) | 1.0e-09 | 0.0048 | 0.0000 | 1.3e-13 | 100 % |
| E3 | optimistic contrast curve | 3 | 2e+05 | 0.0098 +- 0.0043 (0.0108; 0.0051-0.0136) | 8.7e-06 | 0.0024 | 0.0075 | 1.2e-09 | 100 % |
| E4 | propagated (larger) errors | 3 | 2e+05 | 0.0375 +- 0.0059 (0.0381; 0.0313-0.0430) | 2.9e-04 | 0.0146 | 0.0240 | 3.9e-08 | 0 % |
| D | all TIC neighbours (TRICERATOPS default) | 10 | 3e+04 | 0.7519 +- 0.1423 (0.7866; 0.4900-0.9222) | 7.5e-01 | 0.7504 | 0.7512 | 1.0e-03 | 0 % |
| E5 | all neighbours, 1.5-px aperture | 2 | 3e+04 | 0.6452 +- 0.2655 (0.6452; 0.4574-0.8329) | 6.4e-01 | 0.6437 | 0.6392 | 1.5e-03 | 0 % |

* Main configuration: FPP = 0.0123 +- 0.0087 (median 0.0116, range 0.0017-0.0277); mean scenario probabilities: TP 0.8660, DTP 0.0754, PTP 0.0463, STP 0.0070, EBx2P 0.0033, PEBx2P 0.0011, DEBx2P 0.0004.
* The FPP comes from **STP** (planet on the unresolved companion; best fits 6-7 R_e on 0.2-0.36 M_sun companions, e up to
  0.4) and from **x2P eclipsing binaries** (grazing near-equal M-dwarf pairs at 11.38 d with e = 0.2-0.46). TRICERATOPS'
  x2P likelihood fits only the primary-eclipse shape (`likelihoods.lnL_EB_twin`): it ignores that the second eclipse of an
  eccentric orbit is shifted by ~(2/pi) e cos w in phase and changes depth with e sin w, whereas our odd and even transits
  agree in depth (-0.6 +- 1.5 e-/s) and timing (no TTVs). M-dwarf binaries at 5.7-11.4 d are also tidally circularised.
  With circular EB orbits (F) the x2P terms vanish and FPP = 0.0098 +- 0.0122 (median 0.0035, range 0.0002-0.0330), entirely STP.
* **Run-to-run scatter** is large (0.0123 +- 0.0087 (median 0.0116, range 0.0017-0.0277) across 12 seeds) because STP/x2P evidences are dominated by a few rare
  high-likelihood draws. The mean over seeds is the relevant number; the standard error of the mean is 0.0025.
* Without the contrast curve (C) FPP rises to 0.021: more bright unresolved companions are allowed.
* With the propagated bin errors (0.0054 instead of the empirical 0.0043; E4) the light curve discriminates less and
  FPP = 0.037 +- 0.006 (0.015 without x2P). The REF light curve (E1) and the
  600-s binning with S39 (E2) give 0.0075 and 0.0049.
* With all TIC neighbours (D) TRICERATOPS returns FPP = 0.75, NFPP = 0.75: NTP on pair A and NTP/NEB on
  TIC 214311955 (T 13.85, 33") and TIC 214311945 (T 11.08, 53"). TRICERATOPS assigns every star a Gaussian PSF and an
  aperture share and has no access to the PSF-fit deblending or the localisation; the pixel analysis excludes these stars at
  Delta chi2 = 126 (pair A) to 563-592 (TIC 214311955, 214311941; referee fits), 0/1000 injections, and the whole-cutout
  scan (search/) rejects any source more than 2 px away by >= 144. Weighting each neighbour by its pixel likelihood ratio
  (injection-calibrated t-fit LRs for the pair and the faint stars; exp(-Delta chi2/9.6) or the search bound for the others,
  capped at pair A's) gives NFPP = 1.0e-03 (median 6e-04, seeds 7e-05-4e-03), almost all
  NTP on pair A; with the Gaussian-fit pair-A ratio (10^-6.7) it is 5e-07. TRICERATOPS weights pair A's planet
  scenario ~1.3x the on-target TP before localisation (its R^-7 radius prior favours the 4.0 R_e planet pair A would need
  over 5.0 R_e on the M3V), so the NFPP condition hinges on how strongly the pixels reject pair A.

### 3.3 TRICERATOPS limitations in this field
* It cannot use the deblended PSF-fit photometry: nearby-star scenarios are judged from TIC magnitudes, a Gaussian PSF
  (sigma 0.75 px) and an aperture. Next to a Mira that is 7-17x brighter, varies by 0.9 mag in T and whose PSF the SPOC
  PRF does not describe, the aperture shares (Mira 30 % of a 3-px aperture, M3V 6 %) and the required depths are not
  physical. The Mira (TIC T = 11.26, no TIC radius) is given a nominal 1 M_sun and R = 200 R_sun, so its NTP/NEB evidences
  are -inf - correct by accident.
* TIC positions are epoch 2000; the M3V (143 mas/yr) has moved 2.6-3.7" since, so its separations from neighbours are
  off by up to ~3.7" (the 3.1" star is listed at 4.3").
* TIC duplicates (pair A appears twice) double-count neighbours unless removed.
* One exposure time and one flux error: mixed 30-min/10-min/200-s data cannot be modelled together; red noise is not
  represented.
* x2P eclipsing binaries are judged on one eclipse only (above); EB eccentricities are free.
* No occurrence-rate priors: the planet-versus-binary rate priors are switched off in the code (`lnprior_Mstar_planet`,
  `lnprior_Mstar_binary` return 0), and the planet radius prior for M <= 0.45 M_sun falls as R^-7 between 3 and 6 R_e.
* Its companion prior does not know the -0.75 mag over-luminosity, so it does not weight the twin case.
* The background population comes from Gaia within 0.1 deg^2 with dwarf-sequence parameters (giants mis-described).
* The contrast curve is a proxy built from Gaia IPD and DECaPS2, not high-resolution imaging.

### 3.4 Comparison with the Bayesian calculation
Both methods find FPP ~ 0.01 with the same dominant scenario (a planet on an unresolved companion: STP in TRICERATOPS,
PTP here): TRICERATOPS main 0.012 (circular EBs 0.010), Bayesian baseline 0.011
(draws median 0.0067). The NFPP agrees too once the pixel localisation is used, and in both it is set by a planet on
pair A: Bayesian 6.3e-4 and TRICERATOPS (all neighbours, pixel-weighted) 1.0e-3 with the heavy-tailed pair-A ratio
10^-3.2; 2e-7 and 5e-7 with the Gaussian ratio 10^-6.7. (TRICERATOPS restricted to the stars inside the localisation
region gives 1e-5, but that mode removes pair A by construction.) NFPP < 1e-3 therefore holds only if the pixels reject
pair A at least as strongly as the t-fit extrapolation says.

## 4. Is the candidate "statistically validated"?

**No.** By the Giacalone et al. (2021) criteria (FPP < 0.015 and NFPP < 0.001), the candidate is at best marginal and
should be called a **likely planet / planet candidate, not a statistically validated planet**:
1. FPP is ~0.01, within a factor of 1.5 of the threshold: Bayesian baseline 0.0111, but 18 % of plausible
   prior/likelihood choices give FPP >= 0.015 (84th percentile 0.016); TRICERATOPS 0.012 +- 0.009
   with 33 % of seeds above 0.015.
2. Community validations rest on high-resolution imaging; none exists here. The contrast curve used is a Gaia-IPD/DECaPS2
   proxy, and without the archival constraints the Bayesian FPP is 0.045.
3. The host may be an unresolved binary: P(unresolved companion) = 0.47, P(near-equal twin, Delta T <= 0.3) = 0.28.
   If a planet on the twin counts as a false positive, FPP_strict = 0.14.
4. NFPP < 1e-3 is met only marginally and only through a pixel likelihood ratio for pair A that is extrapolated beyond
   1000 injections: Bayesian 6.3e-04 and TRICERATOPS (all neighbours, pixel-weighted) 1.0e-03 with the
   heavy-tailed ratio 10^-3.2; ~1e-6 with the Gaussian ratio 10^-6.7. A ground-based transit that resolves the M3V from pair
   A (11.7") would settle it.
5. The FPP framework assumes an astrophysical, strictly periodic signal. That rests on the phase-1 photometry next to a
   Mira (two independent pipelines, six sectors, constant depth), which no FPP code can test.

**Observation that most reduces the FPP:** a seeing-limited ground-based transit observed simultaneously in two widely
separated bands (e.g. g' and i'/z'), 4-6 events in the 2027 April-August window. On-target detection with equal depths
to ~5 % lowers the Bayesian FPP from 0.011 to ~1e-4 (realistic-pessimistic priors: 0.053 -> 0.0017).
For a validation claim that follows community practice it should be combined with 8-m speckle/AO imaging
(FPP -> 8e-05); high-resolution spectroscopy is needed to decide whether the host is a twin.

## 5. Evidence for and against, and caveats

**For the planet (from the FPP analysis):**
* Every EB class that can mimic the flat-bottomed 1.6-h dip fails on physics: bound WD occultations on the secondary +
  RUWE (Delta chi2 >= 13), hierarchical EBs on shape (>= 27.5), BD/M-dwarf eclipsers on size and duration (>= 45),
  background EBs on photometry, parallax and pixels (bg3: >= 19.6 plus 10^-3.9).
* The pixel localisation excludes the pair (10^-3.2 to 10^-6.7), the Mira (10^-26) and the 5.5" star (10^-10).
* The 3.1" star cannot host the dip physically (it would need a 54 % flat eclipse from a single K dwarf); its posterior is < 1e-8.
* The archival look-back (DSS) excludes unbound stars bright enough to host the dip within 3" except near catalogued stars.

**Against or open:**
* The largest open space is bound companions inside ~0.6" (< 40 au), especially < 0.08", with Delta T ~ 0.75-2. A planet on
  such a companion is still a ~6 R_e planet, but on another star.
* The over-luminosity (-0.75 mag) makes an unresolved twin likely; the planet's host, radius (4.9 vs 5.2 R_e) and
  insolation then change, and FPP_strict = 0.14.
* The NFPP depends on a pixel likelihood ratio for pair A extrapolated far beyond the 1000 injections (Gaussian 10^-6.7,
  Student-t 10^-3.2). The pair hosts would need only a 2.2-2.5 % (pair A) dip, so the pixels alone exclude them.
* TRICERATOPS FPP straddles 0.015 between seeds; its mean is close to the threshold.
* The analysis assumes the signal is astrophysical; the S103/S104 depth trend (1.5-2 sigma), the odd/even difference seen
  only in one photometry and the Mira's proximity are phase-1 issues outside the FPP.

**Caveats of the Bayesian model:**
* Occurrence of 4-7 R_e planets around mid-M dwarfs is poorly known (0.1-1 %; Gillis+2026 upper limit 3.4 % for
  3-6.5 R_e at P < 30 d); the FPP is insensitive to its level (it cancels between TP and PTP) but sensitive to how it
  scales with host mass (beta_N).
* The planet suppression in close binaries (S_close = 0.25, varied 0.05-1) and the 0.25 au-to-10 au shape are taken from
  memory of Kraus+2016 / Moe & Kratter 2021, not re-verified here (referee_ctx found the fetched Moe & Kratter abstract was
  the 2018 paper).
* The WD+M dwarf rate at P ~ 5.7 d (2e-5 per ln P) is an order-of-magnitude estimate from PCEB space densities
  (~1e-5 pc^-3) and period distributions, from memory; it is varied 2e-6-2e-4 and stress-tested at 2e-3 and 1e-2.
* EB evidences use maximum-likelihood Delta chi2 from the shape agent times an Occam factor (0.1; 0.03-1), not full
  marginal likelihoods; the planet evidences are marginalised on the grid. This favours the EBs.
* The Gaia IPD likelihood cells have few calibration pairs at Delta I > 3 and inside 0.15"; the astrometric grid is a
  simplified simulation (referee: too strict for 0.3 M_sun); both are varied.
* The CMD single-star scatter sigma_s (0.30 mag) is set from the bright-side width of nearby RUWE < 1.1 M dwarfs; it is
  the second most influential input.
* The light-curve grid uses fixed limb darkening for all hosts and fixed T0/P; companions of 0.15-0.3 M_sun have
  slightly different limb darkening (a small effect at this S/N).
* Bound-companion separations are drawn from the Winters+19 field distribution; companions in a system that already has a
  wide companion (pair A+B at 754 au) may be more frequent (f_mult varied to 0.45).

## 6. Key numbers

| Quantity | Value | Method | Uncertainty | File |
|---|---|---|---|---|
| FPP (Bayesian, task definition) | 0.0111 | explicit scenario evidences, INV | 16-84 % of 3000 prior draws 0.0024-0.0157; 5-95 % 0.0008-0.0245 | out/fpp_baseline.json, out/fpp_random.csv |
| P(FPP < 0.015) over prior draws | 0.82 | joint random draws | | out/fpp_summary.json |
| NFPP (Bayesian) | 6.3e-04 | pair A NTP x pixel LR | 95th pct 6.7e-04 | same |
| FPP_strict / FPP_system | 0.143 / 7.2e-04 | same, other definitions | strict 5-95 % 0.034-0.331 | same |
| Dominant residual | PTP 0.0104 (5.8 R_e planet on a 0.21 M_sun companion at 0.05-0.64") | posterior-weighted MC | dominant in 2739/3000 draws | out/ptp_weights.npz |
| P(unresolved companion) / P(twin) | 0.47 / 0.28 | CMD + IPD + astrometry + DECaPS2 + transit evidences | sigma_s 0.25-0.40 changes P(twin-hosted planet) 0.27-0.04 | out/fpp_baseline.json |
| Pixel log10 LR vs M3V | pair A -3.2 (t) / -6.7 (Gauss); pair B -10.4; Mira -26.3; 3.1" star -3.9; 5.5" star -10.2 | 1000 injections per source (referee pipeline) | extrapolated beyond 1e-3 | out/loc_likelihoods.json |
| Uncatalogued-source effective area | 18.2 arcsec^2 (kappa 2.2) | tempered Delta chi2 map | 17.6-26.7 (kappa 1.5-4.8) | out/loc_likelihoods.json |
| TRICERATOPS FPP (main) | 0.0123 +- 0.0087 | 12 seeds, N = 2e5 | range 0.0017-0.0277; SEM 0.0025 | out/tri_summary.json, out/tri_runs.csv |
| TRICERATOPS FPP, circular EBs | 0.0098 +- 0.0122 | 10 seeds | range 0.0002-0.0330 | same |
| TRICERATOPS FPP, no contrast curve | 0.0208 +- 0.0129 | 10 seeds | | same |
| TRICERATOPS NFPP, neighbours restricted to the localisation region (3.1" and 5.5" stars) | 1e-05 | main runs | pair A excluded by construction in this mode | same |
| TRICERATOPS FPP/NFPP, all TIC neighbours | 0.75 / 0.75 | 10 seeds, N = 3e4 | pixel-weighted NFPP 1.0e-03 (t) / 5e-07 (Gauss) | same |
| TRICERATOPS FPP with propagated errors | 0.0375 +- 0.0059 | 3 seeds, sigma 0.0054 instead of 0.0043 | without x2P 0.0146 | same |
| FPP after a ground multicolour transit on target | 1.2e-04 | what-if (depth ratio to 5 %, pair resolved) | realistic-pessimistic 0.0017 | out/fpp_variants.csv |
| FPP after 8-m speckle imaging | 0.0017 | what-if | realistic-pessimistic 0.021 | same |

## 7. Scripts and outputs

| Script | What it does | Outputs |
|---|---|---|
| `scripts/fpplib.py` | paths, constants, ephemeris, INV (RRab-cleaned) and REF window loaders | (library) |
| `scripts/f01_loc_likelihoods.py` | injection-calibrated localisation LRs of the catalogued alternatives; effective area for uncatalogued sources | `out/loc_likelihoods.json`, `logs/f01.log` |
| `scripts/f02_lc_grid.py <inv\|ref>` | light-curve chi2 grid over (D, rho_circ, b) with optimised k and Laplace widths | `out/lc_grid_inv.npz`, `out/lc_grid_ref.npz`, `logs/f02_*.log` |
| `scripts/f03_companions.py` | bound-companion Monte Carlo with IPD, astrometric (wobble-mapped), DECaPS2 and CMD likelihoods | `out/companions_mc.npz`, `out/companions_summary.json`, `logs/f03.log` |
| `scripts/f04_bayes.py` | scenario evidences, posteriors, FPP/NFPP (three definitions), 96 named variants, 3000 joint prior draws, what-if follow-up | `out/fpp_baseline.json`, `out/fpp_variants.csv`, `out/fpp_random.csv`, `out/fpp_summary.json`, `out/ptp_weights.npz`, `logs/f04.log` |
| `scripts/f05_tri_prep.py` | folded/binned light curves for TRICERATOPS (inv200, inv600, ref200); contrast curves from IPD + DECaPS2 (not high-resolution imaging) | `data/tri/lc_*.npz`, `data/tri/cc_cons.txt`, `data/tri/cc_opt.txt`, `data/tri/inputs.json` |
| `scripts/f06_tri_run.py` | one TRICERATOPS run (local TESScut cubes, apertures = pixels within r px of the M3V, stellar overrides, neighbour modes keep/loc/drop, optional circular EBs) | `out/tri/<run>.json`, `logs/tri_*.log` |
| `scripts/run_tri_batch.sh`, `run_tri_batch3.sh`, `run_tri_batch4.sh` | the run sets A, C, D, E1-E5 (stream 1), F, G (stream 2), E, D (stream 3; duplicates skipped) | `logs/tri_batch*.log` |
| `scripts/f07_tri_collect.py` | run-to-run statistics; FPP without x2P scenarios; FPP_system; localisation-weighted NFPP | `out/tri_runs.csv`, `out/tri_summary.json` |
| `scripts/f08_figs.py` | figures | `figs/*.png` |

Other files: `data/tri/214311908_gaia_background.csv` (TRICERATOPS' Gaia field population), `logs/` (console output of
every script and run). A first TRICERATOPS timing test (N = 20,000, all neighbours, pair A counted twice through a TIC
duplicate) gave FPP = NFPP = 0.64; it was deleted and replaced by the D runs.

Figures: `figs/scenario_posteriors.png` (baseline posteriors), `figs/fpp_distribution.png` (FPP and NFPP over 3000
prior draws, TRICERATOPS runs marked), `figs/variants.png` (all named variants), `figs/ptp_map.png` (where the residual
PTP probability lies in separation-contrast space, with the IPD cells, the DECaPS2 limit and a hypothetical speckle
curve), `figs/tri_runs.png` (TRICERATOPS FPP/NFPP per run and group).
