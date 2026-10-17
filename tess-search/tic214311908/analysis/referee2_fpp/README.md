# Referee 2: independent re-derivation of the false-positive probability (FPP) of TIC 214311908

Folder: `tic214311908/analysis/referee2_fpp/`
(AI-run referee analysis, 2026-10-17. Nothing was posted, uploaded or submitted, and no one was contacted. No network
access was used except what TRICERATOPS itself queries (MAST TIC / Gaia field population) in `r06b`.)

This is an adversarial check of `../fpp/`. All code here is my own (`scripts/`). I read the analyst's scripts to look
for errors but did not import or run them. Shared inputs are data products only: the phase-1 light curves, the
referee_loc injection sets, the Gaia IPD cells and astrometric grid of `../star/`, the comparison-star table
`../referee_star/data/mdw_cmd.csv`, the Mamajek table, the TESScut cubes and the noise calibrations of
`../referee_phot/` and `../referee2_shape/`.

## Verdict

1. **The analyst's scenario bookkeeping is correct.** When I feed the analyst's over-luminosity likelihood into my own
   framework, I reproduce the analyst's configuration and FPP numbers. My framework has its own light-curve grid,
   companion Monte Carlo and scenario code. With the analyst's Gaussian sigma_s = 0.30 at x = -0.753, I get
   P(binary) = 0.472 (analyst 0.47), P(dT <= 0.3) = 0.265 (0.28), FPP = 0.0099 (0.0111 INV / 0.0081 REF) and
   FPP_strict = 0.095 (0.143).
2. **The over-luminosity likelihood is wrong, and it drives the headline numbers.** The analyst fixed the single-star
   scatter at sigma_s = 0.30 mag and ignored colour. Two independent facts contradict this:
   * **sigma_s.** A mixture fit to the 2,126 Gaia-clean comparison M dwarfs (single stars plus an explicit
     unresolved-binary population) prefers a single-star bright-side width of 0.40 mag. sigma_s = 0.30 is worse by
     -2 dlnL = 21 (`out/r02_cmd_colour.json`, `figs/r07_cmd_colour.png`). The analyst's own percentiles (16th
     -0.34, median +0.05) already imply a width of ~0.39 about the median.
   * **The G-Ks colour residual at fixed BP-RP.** Single M dwarfs follow a tight metallicity trend: r = -0.078 x
     with scatter 0.044. The target (x = -0.84, r = +0.097) sits on the metal-rich extension of that trend, where
     +0.070 is expected. A near-equal twin keeps the colours of a normal star (r ~ 0.00-0.01), which puts it 2.0
     sigma off.
   * **Combined result.** P(binary) falls from 0.47 to 0.10 (0.08-0.13 over F_twin) and P(dT <= 0.3) from 0.28 to
     0.020. Only 4-14% of the remaining binaries have q >= 0.9.
   * **Independent confirmation.** The `../host/` analysis was finished before `fpp` but not used by it. It finds
     metal-rich colours and XP spectra and gets P(twin) ~ 1%.
3. **Corrected Bayesian FPP.**
   * Baseline: **0.0031 with the REF light curve and 0.0071 with INV** (beta-inflated noise in both).
   * 1,500 joint prior draws: median 0.0037, 16-84% 0.0007-0.011, 95th percentile 0.019; P(FPP < 0.015) = 0.92.
     4 of 70 named variants reach FPP >= 0.015, all of them extreme combinations.
   * FPP_strict (a planet on any companion counts as a false positive, as in TRICERATOPS) drops from 0.143 to
     **0.008 (REF) / 0.015 (INV)**. Over the draws its median is 0.018 (5-95% 0.003-0.12).
   * Dominant residual: still a planet on a fainter bound companion (PTP 0.0016). The companion is ~0.20 Msun,
     dT 1.2 (0.9-1.4), a = 17 au (7-45), separation 0.10" (0.05-0.49").
4. **A false-positive class was missed: brown dwarfs.** The analyst's EB list has no eccentric brown dwarf.
   * The physical HEB grid in `../shape` assumed circular orbits. Its BD-on-companion case failed on duration
     (T14 = 1.37 h). Eccentricity removes that objection, as it does for the analyst's own PTP planets.
   * Brown dwarfs give **~9e-4 (REF) to 2e-3 (INV)** of posterior, about a quarter to a third of the corrected
     FPP. This assumes f_BD = 1e-3 per star at 3-10 d and the Kipping e prior.
   * Most of it is a grazing BD on the M3V. That needs rho_circ ~ 1.7 g/cm3, so e >~ 0.5, plus a V shape costing
     dchi2_LC ~ 4.3-5.8 with beta-inflated noise. It vanishes for circular orbits.
   * **It is conditional on the noise model and e prior (r11, `out/r11_bd_noise.json`).** The best grazing-BD
     solution costs dchi2 5.8 (REF) / 4.3 (INV) with beta-inflated errors, but 21.1 (REF) / 7.3 (INV) with white
     errors. The BD total is then 9.5e-4 (REF beta), 2.2e-3 (INV beta), 6.5e-4 (INV white), 5e-7 (REF white), and
     < 1e-7 for circular orbits. The calibrated V-shape test of `../referee2_shape` (V-shaped truth reproduces the
     observed flat preference in 0/25 REF, 1/25 INV injections, so LR <~ 0.04) is about as strict as the REF-beta
     grid (0.054) and stricter than the INV-beta grid (0.12). The baseline BD share is therefore an upper-end value.
     Real BDs at P = 5.7 d around M dwarfs with e > 0.5 are rarer than the Kipping prior implies.
   * Where it survives, this class sets the floor after speckle and two-band follow-up (FPP ~ 1e-3, not 1e-4).
     Only RV (K ~ 13 km/s for 60 M_J, ~3 km/s for 13 M_J; a 5 R_e planet gives ~20 m/s) or a precise transit
     shape removes it.
5. **NFPP.** The analyst adopted pixel log10 LR(pair A) = -3.22, which is the **most lenient of the 19** referee_loc
   variants that have injections. The median over variants is -9.2 and the Gaussian baseline -6.7. Only a
   systematics-inflated Gaussian evaluated at the worst variant (s_obs = +75) gives -2.8. The resulting NFPP is
   4.5e-4 (t-fit, conservative), 1e-7 to 1e-9 typically, and 1.2e-3 only in that worst case. P(NFPP < 1e-3) = 0.97
   over the draws.
6. **TRICERATOPS.** Re-aggregating the analyst's 59 JSON outputs reproduces their statistics exactly. However:
   * the seed scatter is so large that configuration comparisons are not significant. Run F (contrast curve +
     circular EBs) gives 0.0098 against 0.0029 for run G (no contrast curve + circular EBs), which is the
     impossible ordering;
   * FPP changes by 3x with the error convention (0.012 to 0.037, E4);
   * my own re-run (REF light curve, own fold, 2.5-px aperture) gives 0.0016 +- 0.0020 with my red-star IPD
     contrast proxy, 0.0115 +- 0.0103 with the analyst's conservative proxy (reproducing the analyst's 0.012), and
     0.039 +- 0.012 with beta-inflated errors (section 6). TRICERATOPS therefore spans 0.002-0.04, set by the
     contrast proxy and the noise convention.
   The claim that the Bayesian and TRICERATOPS FPPs "agree" compares different definitions and different
   information. TRICERATOPS counts every companion-hosted planet, twins included, as a false positive, and it does
   not know about the over-luminosity. The two methods are not independent either: they share the light curve and
   the contrast proxy.
7. **"Statistically validated?"** Still **no** by community standards, but for different reasons than the analyst
   gives.
   * The Bayesian FPP is now ~0.003-0.007, comfortably below 0.015, and the twin question is largely settled
     (P ~ 1-2%).
   * What remains: (i) there is no high-resolution imaging, and the Giacalone et al. (2021) thresholds are
     calibrated for TRICERATOPS with such imaging; (ii) TRICERATOPS is unstable here (0.0005-0.04 depending on
     seed, noise convention and contrast proxy); (iii) the residual BD channel; (iv) the non-FPP photometric
     systematics next to the Mira.
   * Most useful observation: 8-m speckle/AO imaging (expected for any validation; FPP 0.0031 to 0.0015) together
     with a seeing-limited two-band ground transit (resolves pair A, tests achromaticity; 0.0010; jointly 0.0009).
     After both, the residual is the BD channel, which one or two RV epochs would close (FPP ~1e-4).

## 1. Scripts and outputs

| Script | What it does | Outputs |
|---|---|---|
| `scripts/rlib.py` | Paths, constants, Mann+19 mass-MKs and Mann+15 radius relations (own implementation), own parser of the Mamajek table (`../referee_star/data/EEM_dwarf.txt`), band offsets | (library) |
| `scripts/r01_loc_lr.py` | Pixel likelihood ratios of the catalogued alternatives (pair A/B, Mira, 3.1", 5.5", 12" stars) for every referee_loc variant with injections (19), with Gaussian, Student-t, generalised-Pareto-tail and systematics-inflated density models | `out/r01_loc_lr.json`, `logs/r01.log` |
| `scripts/r02_cmd_colour.py` | Mixture model of the Gaia-clean comparison M dwarfs in the plane of (M_Ks offset x, G-Ks residual r): single stars (split-normal + faint tail, metallicity colour trend) plus unresolved binaries (Mann/Mamajek shifts). Profile likelihood over the single-star bright-side width sigma_b, for F_twin = 0 / 0.1 / 0.25; P(binary | x) and P(binary | x, r) at the target | `out/r02_cmd_colour.json`, `logs/r02.log` |
| `scripts/r02_cmd_mixture.py` | First attempt (x only, free tail): degenerate, abandoned; kept for transparency | `logs/r02_first_attempt_degenerate.log` |
| `scripts/scratch_colour.py` | Exploratory table of the median G-Ks residual against the M_Ks offset, and the binary colour tracks | console only (numbers quoted in section 3) |
| `scripts/r03_lc_grid.py <ref\|inv> <beta\|white>` | Own light-curve grid. Dark body of radius ratio k on a host supplying a fraction D of the flux (12 D x 30 rho_circ x 31 b); k integrated numerically; per-window polynomial baselines profiled; midpoint exposure integration | `out/r03_grid_{ref,inv}_{beta,white}.npz`, `logs/r03_*.log` |
| `scripts/r04_companions.py` | Own bound-companion Monte Carlo (300,000 draws). Constraint likelihoods: IPD (all / 84% / red-star / beyond 0.25" only), astrometry (wobble-mapped), DECaPS2 (graded: 0.02 within 1.5 mag of the limit, 1e-3 beyond), Gaia-wide (> 2.5"), plus S-type stability and colour shifts | `out/r04_companions.npz`, `logs/r04.log` |
| `scripts/r05_fpp.py` | Scenario evidences: TP single / diluted / twin, PTP, brown dwarfs on every host (new), WD occultation, MS hierarchical EB, NTP pair A, background term. FPP, FPP_strict, FPP_system, NFPP; 70 named variants; 1,500 joint draws | `out/r05_baseline.json`, `out/r05_variants.csv`, `out/r05_random.csv`, `out/r05_summary.json`, `logs/r05.log` |
| `scripts/r06_tri_check.py` | Re-aggregates the analyst's TRICERATOPS JSONs and my re-runs from the stored scenario tables | `out/r06_tri_check.json`, `out/r06_tri_runs.csv` |
| `scripts/r06b_tri_run.py`, `run_tri_r2.sh`, `run_tri_r2b.sh` | My TRICERATOPS runs (venv_tri, NUMBA 2 threads): REF light curve folded and binned by my own code, the analyst's or my contrast proxy, 2.5-px aperture, localisation-region neighbours | `out/tri/r2_*.json`, `data/tri/lc_ref200_r2.npz`, `data/tri/cc_r2.txt`, `logs/tri_r2_*.log` |
| `scripts/r07_diagnostics.py` | Properties of the PTP residual (V/I depth ratio, Rp, separation); analyst-CMD vs referee CMD side by side; figures | `out/r07_diagnostics.json`, `figs/r07_cmd_colour.png`, `figs/r07_fpp_draws.png` |
| `scripts/r08_whatif_extra.py` | Follow-up what-ifs: two-band 10% vs 5% precision, speckle, combinations, with and without the BD class, for REF and INV | `out/r08_whatif.json` |
| `scripts/r09_aeff.py` | Re-computes the analyst's effective localisation area for uncatalogued sources from the referee_loc sky map | `out/r09_aeff.json` |
| `scripts/r11_bd_noise.py` | Robustness of the BD channel: r05 BD terms for REF/INV x beta/white noise x Kipping/Van Eylen/circular e, and the best grazing-BD dchi2 on the M3V in each grid | `out/r11_bd_noise.json`, `logs/r11.log` |
| `scripts/r10_t0_check.py` | FPP and ln B(single/twin) with the REF grid rebuilt at T0 = 3840.9160 (REF free-fit epoch; `r03` run with R03_T0, tag _T0ref) | `out/r10_t0_check.json`, `out/r03_grid_ref_beta_T0ref.npz` |

Other files: `data/tri/214311908_gaia_background.csv` (TRICERATOPS' Gaia field population, written by TRICERATOPS
into its working directory), `logs/` (console output of every script and run).

## 2. Key numbers (mine vs the analyst's)

| Quantity | Analyst (fpp) | Referee 2 | Method / file | Verdict |
|---|---|---|---|---|
| Scenario bookkeeping with the analyst's CMD form | FPP 0.0111 (INV) / 0.0081 (REF); P(bin) 0.47; P(dT<=0.3) 0.28 | FPP 0.0099, strict 0.095, P(bin) 0.472, P(dT<=0.3) 0.265 (REF) | own grid + MC + evidences, `out/r07_diagnostics.json` | upheld (no arithmetic error) |
| Single-star CMD width sigma_s | 0.30 (assumed) | 0.40 preferred; 0.30 worse by -2 dlnL = 21; 0.35-0.45 carry all profile weight | `out/r02_cmd_colour.json` | **revised** |
| G-Ks residual test (twin vs metal-rich) | not used | target r = +0.097; single trend at x_t: +0.070 +- 0.044; twin: ~0.00-0.01 | r02; `figs/r07_cmd_colour.png` | **missed test** |
| P(unresolved companion) / P(dT <= 0.3) | 0.47 / 0.28 | **0.10 / 0.020** (with all target-specific constraints; population-only 0.08-0.13) | r05 config; r02 | **revised** |
| FPP (task definition) | 0.0111 (INV) | **0.0031 (REF) / 0.0071 (INV)** | r05 baseline / variants | **revised** |
| FPP over joint draws | median 0.0067, 16-84% 0.0024-0.016, P(<0.015) 0.82 | median 0.0037, 16-84% 0.0007-0.011, 95% 0.019, P(<0.015) 0.92 | `out/r05_summary.json` | revised |
| FPP_strict (twin = FP) | 0.143 | **0.008 (REF) / 0.015 (INV)**; draws median 0.018 (0.003-0.12) | r05 | **revised** |
| FPP_system | 7.2e-4 | 1.5e-3 (incl. brown dwarfs) | r05 | revised |
| Dominant residual | PTP 0.0104; 0.21 Msun, dT 1.0, 0.05-0.64" | PTP 0.0016; 0.20 Msun, dT 1.2 (0.9-1.4), 0.05-0.49", a 7-45 au; V/I depth ratio 0.75 (0.68-0.81) | r07 | upheld (smaller) |
| Brown-dwarf channels | absent | 9.5e-4 (REF), 2.1e-3 (INV): grazing eccentric BD on M3V 7.7e-4; BD on a bound companion 1.1e-4. Beta-noise values; white noise 5e-7 (REF) / 6.5e-4 (INV); circular orbits < 1e-7 | r05, r08, r11 | **missed scenario** (upper-end) |
| NFPP | 6.3e-4 (t-fit LR 10^-3.22) | 4.5e-4 with the same LR; LR is the most lenient of 19 variants (median 10^-9.2); range 1e-9 to 1.2e-3 | r01, r05 | upheld as conservative |
| Pair A pixel log10 LR | -3.22 (t) / -6.74 (Gauss) | base reproduced; over 19 variants t: median -9.2 (-14.3 to -3.22); Gauss median -9.6; sys-inflated -6.1 (base) / -2.8 (worst variant) | `out/r01_loc_lr.json` | upheld, framing too pessimistic |
| 3.1" star pixel log10 LR | -3.87 | base -3.87; variants -4.2 to -0.5 (the -0.5 needs BP-RP 2.76, unphysical) | r01 | upheld (irrelevant: needs a 54% eclipse) |
| WD-occultation dchi2 | 13.2 | 10 (referee2_shape correction: 9-11; 6-8 at sigma_CCD 0.35): posterior 6.8e-5 | r05 | negligible either way |
| Radius prior inside 4-7 R_e | flat in ln R | R^-3 across 4-10 R_e: FPP 0.0031 to 0.0023 | r05 | sensitivity not explored by the analyst |
| IPD cells inside 0.25" | applied as calibrated | calibration stars are nearer, so at a given angle they have shorter periods and larger RUWE. Using cells only beyond 0.25" gives FPP 0.0062 | r05 | flag |
| TRICERATOPS main (A) | 0.0123 +- 0.0087 (12 seeds) | identical from the JSONs; bootstrap 68% of the mean 0.010-0.015 | `out/r06_tri_check.json` | upheld |
| TRICERATOPS re-run (own inputs) | - | see section 6 | r06b | |
| Effective area for uncatalogued sources | 18.2 arcsec^2 (kappa 2.2) | 18.2 (kappa 2.2), 20.5 (1.0), 26.7 (4.8); expected T <= 17.47 stars inside it before archival vetoes 0.067 | `out/r09_aeff.json` | upheld |
| Fixed T0 | 3840.9155 | baseline uses the phase-1 epoch 3840.9141; with the REF free-fit epoch 3840.9160: FPP 0.0035 (vs 0.0031), ln B(single/twin) -0.16 / -0.57 (Kipping / circular; analyst -0.09 / -0.50 REF); the analyst's CMD form on that grid gives FPP 0.0115, strict 0.126 (analyst REF: 0.0081, 0.127) | `out/r10_t0_check.json` | small effect; bookkeeping reproduced |
| What-if: speckle null | 0.0017 | 0.0015 (REF) / 0.0027 (INV); residual = BD | r08 | upheld |
| What-if: two-band ground transit (+ pair resolved) | 1.2e-4 | 0.0010 (5%) / 0.0011 (10%) REF; 0.0020-0.0023 INV; residual = BD; 1.1e-4 only once BDs are also excluded (speckle + two-band, BDs removed by RV) | r08 | **overstated** |

## 3. The over-luminosity: why P(twin) collapses (r02)

* **Comparison sample.** 2,126 Gaia-clean M dwarfs with |BP-RP - 2.756| < 0.1, RUWE < 1.15,
  ipd_frac_multi_peak = 0, gof_harmonic < 0.04 and non_single_star = 0, from `../referee_star/data/mdw_cmd.csv`.
* **Offsets.** The ridge comes from KDE modes in colour bins (slope 1.53), giving a target offset
  **x = -0.838** (the referee ridge gives -0.753; both are carried as variants).
* **Colour trend.** The median G-Ks residual climbs steadily with over-luminosity:
  * -0.041 at x 0.4..1.5;
  * +0.004 at -0.15..0.15;
  * +0.040 at -0.6..-0.4;
  * +0.057 at -0.75..-0.6;
  * +0.074 at -0.9..-0.75 (53 stars).
  Binaries cannot make this trend at the bright end. A q = 1 twin sits at (x, r) = (-0.753, 0.000), and q = 0.8 at
  (-0.66, +0.01). The over-luminous field stars at the target's x are therefore mostly metal-rich single stars.
* **Target.** It has r = +0.097, on the single-star trend. The `host/` agent shows that 2MASS Ks is sound and gets
  the same r from G-W1, which is independent of 2MASS.
* **Mixture fit.** The binary fraction among clean stars comes out at f = 0.030-0.040. My companion Monte Carlo
  predicts 0.04-0.05 for the same Gaia cuts, which is consistent. The bright-side width is sigma_b = 0.40, with
  profile weight 0.88 at 0.40 and 0.12 at 0.35. Population posteriors at the target:

  | F_twin | P(binary given x) | P(binary given x, r) | q >= 0.9 among binaries |
  |---|---|---|---|
  | 0 | 0.203 | 0.125 | 4% |
  | 0.1 | 0.185 | 0.102 | 8% |
  | 0.25 | 0.163 | 0.076 | 14% |

* **Caveats.** The split-normal plus tail model is a simplification: the faint tail hits its bounds, though it
  affects only the faint side. The profile preference for 0.40 over 0.30 is model-dependent, so it is varied
  (fixed sigma_b 0.2/0.3/0.4; x only; the analyst's Gaussian form).

## 4. Corrected scenario posteriors (r05 baseline: REF light curve, beta noise, Kipping e, colour CMD marginalised over sigma_b, F_twin 0.1)

| Scenario | Posterior |
|---|---|
| TP single | 0.957 |
| TP diluted (planet on the primary of a binary) | 0.035 |
| TP on twin (dT <= 0.75) | 0.0048 |
| **PTP (planet on a fainter bound companion)** | **0.0016** |
| BD grazing the M3V (eccentric) | 7.7e-4 |
| NTP pair A (t-fit LR) | 4.5e-4 |
| BD on a bound companion / on the primary of a binary | 1.1e-4 / 7.7e-5 |
| WD occulted by the M3V (dchi2 10) | 6.8e-5 |
| Background-star scenarios (carried at the analyst's level) | 3e-5 |
| MS hierarchical EB (dchi2 27.5) | 8e-6 |

Configuration (before the transit): P(single) 0.90, P(dT <= 0.3) 0.020, P(dT 0.3-0.75) 0.008, P(dT 0.75-2) 0.028,
P(dT > 2) 0.042.

**Named variants (70, `out/r05_variants.csv`)**: FPP range 7e-4 to 0.042. Values of 0.015 or more occur only for:
* the analyst's Gaussian CMD form with sigma_s 0.25 (0.022);
* all archival constraints off (0.027);
* realistic-pessimistic priors (0.028; 0.042 if the IPD is also dropped inside 0.25").

Selected variants:

| Variant | FPP |
|---|---|
| INV | 0.0071 |
| White noise | 0.0016 (REF) / 0.0052 (INV) |
| e circular | 0.0007 |
| beta_N = -1 / +2 | 0.0045 / 0.0019 |
| f_mult 0.45 | 0.0051 |
| S_close 1 | 0.0058 |
| CMD without colour term | 0.0051 |
| IPD LR84 | 0.0037 |
| IPD off | 0.011 |
| Astrometry off | 0.0076 |
| f_BD 0 / 3e-3 | 0.0022 / 0.0050 |
| dT_twin 0.5 | 0.0036 |

**Joint draws (1,500)**: the median by light curve is 0.0026 (REF, beta) and 0.0060 (INV, beta). The Spearman
correlations of FPP with the inputs are weak (|rho| <= 0.22; strongest beta_N, S_close, f_mult). Dominant false
positive: PTP in 1,006 draws, WD in 173, grazing BD in 161, MS-HEB in 68, NTP pair A in 67.

## 5. Other findings on the analyst's analysis

* **Correct.** The FPP definitions, the photo-eccentric treatment, the transit-probability factor, the q/a
  re-weighting, the S-type stability, the PTP localisation in separation and contrast, and the V/I chromatic
  argument for PTP.
* **Twin boundary dT <= 0.75.** This is generous. A planet on a dT = 0.75 companion would be ~5.8 R_e (grid k = 0.119)
  to ~6.3 R_e (depth-based k = 0.129), not 4.9-5.2. With the corrected CMD the choice hardly matters (dT_twin 0.5 gives 0.0036). With the analyst's CMD it
  moved the FPP across the threshold (0.024 at 0.5).
* **DECaPS2 floor.** Applying 0.02 to companions far brighter than the 5-sigma limit has a negligible effect
  (0.0031 either way).
* **IPD calibration distance mismatch.** This is not discussed by the analyst and is worth a factor of ~2 in FPP.
* **TRICERATOPS scatter.** A heavy-tailed estimator needs many more seeds, or a larger N, before configuration
  differences of 2-3x mean anything. The statements "C 0.021", "F 0.0098 entirely STP" and "G 0.0029" are within
  noise of each other.
* **Shared assumption.** The FPP framework assumes the signal is astrophysical and on-sky periodic. The analyst
  states this correctly; phase-1 issues sit outside it.

## 6. TRICERATOPS re-run (r06b)

Setup:
* light curve: REF, re-folded by my own code; flux_err = OOT rms = 0.00351 ('emp') or the propagated errors with
  fake-window beta ('prop');
* depth 0.01563;
* neighbours: the 3.1" and 5.5" stars only;
* N = 2e5;
* contrast curve: my red-star IPD + DECaPS2 proxy (2 mag inside 0.75", 3.5 mag at 0.8-0.9", DECaPS2 beyond 1"),
  or the analyst's conservative curve. Neither is a true high-resolution imaging curve.

Results, from `out/r06_tri_check.json` and `out/r06_tri_runs.csv` (FPP recomputed from the stored scenario tables):

| Group | Light curve / errors | Contrast proxy | Seeds | FPP mean +- sd (median; range) | Bootstrap 68% of mean | NFPP | Dominant |
|---|---|---|---|---|---|---|---|
| r2_emp | REF (mine), OOT rms 0.00351 | mine (red-star IPD cells + DECaPS2) | 6 | 0.0016 +- 0.0020 (0.0007; 3.5e-5 to 0.0048) | 0.0009-0.0024 | 4e-8 | STP, DTP |
| r2_emp_cccons | REF (mine), OOT rms 0.00351 | analyst's conservative cc_cons | 4 | 0.0115 +- 0.0103 (0.0105; 0.0019-0.023) | 0.0067-0.0164 | 4e-8 | STP |
| r2_prop | REF (mine), beta-inflated propagated 0.00642 | mine | 3 | 0.039 +- 0.012 (0.026-0.050) | 0.034-0.043 | 1.2e-3 to 2.1e-3 | EBx2P (single-eclipse x2P likelihood) |
| analyst A | INV, 0.0043 | cc_cons | 12 | 0.0123 +- 0.0087 | 0.0100-0.0148 | 1e-5 | STP, x2P |
| analyst E1 | REF (analyst fold), 0.0036 | cc_cons | 3 | 0.0075 +- 0.0040 | | 5e-8 | |
| analyst E4 | INV, propagated 0.0054 | cc_cons | 3 | 0.037 +- 0.006 | | 2.9e-4 | |

What this shows:
* **Contrast proxy.** The proxy curve alone changes TRICERATOPS' FPP by ~7x (0.0016 against 0.0115, same light
  curve and errors). The analyst's conservative curve is 0.01 mag inside 0.4". It comes from one Gaia-IPD cell
  (0.25-0.4", dI 0-1, LR_84 = 0.27) propagated inward by the running minimum.
* **Error convention.** A ~2x larger flux error raises the FPP to 0.03-0.05, through the x2P EB scenarios whose
  second eclipse TRICERATOPS does not test. The beta-inflated propagated errors overstate the scatter of the binned
  fold (OOT rms 0.0035 against 0.0064).
* **Net.** TRICERATOPS spans 0.0016-0.04 depending on two analyst choices that real high-resolution imaging and a
  properly calibrated noise model would remove. It cannot by itself support or refute validation here. The
  analyst's reported numbers are reproduced.
* **Bug of my own, fixed.** My first three 'prop' runs used a wrong bin-error formula (flux_err ~0.4, FPP 0.56).
  They are quarantined in `out/tri_buggy_prop/` and not used; the 'prop' rows above are the corrected re-runs.

## 7. Important missed tests and open issues (for the synthesis)

1. **Use the colour and metallicity information** (G-Ks and G-W1 residuals, XP shape; `../host/`) in the companion
   prior. This alone moves P(twin) from 0.28 to ~0.01-0.02 and FPP_strict from 0.14 to ~0.01.
2. **Fit, rather than assume, the single-star CMD width** with a binary-aware mixture (r02). sigma_s = 0.30 is
   disfavoured.
3. **Add brown dwarfs** (eccentric) as an FP class: on bound companions (HEB-BD) and grazing the M3V. Use
   BD-appropriate eccentricities. The Kipping prior used here is generous: the grazing case needs e >~ 0.5 at
   P = 5.7 d, which is likely rare for BDs. RV (K ~ 13 km/s for 60 M_J) excludes them cheaply. The channel is
   only significant with beta-inflated errors (r11); with white errors on REF it is 5e-7.
4. **IPD cells inside ~0.25".** The calibration stars are at smaller distances, so at a given angle they have
   shorter orbital periods and more RUWE excess. Re-calibrate in physical separation, or drop those cells
   (FPP x2).
5. **TRICERATOPS.** It needs a real contrast curve and either more seeds or larger N, since the mean over 12 seeds
   at N = 2e5 is only good to ~25%. The x2P scenarios should be filtered by the observed odd/even equality and
   timing.
6. **Light-curve choice.** INV and REF differ by ~2x in FPP (REF has lower per-transit noise; the INV curve still
   carries the RRab leak, here absorbed by quadratic baselines). Both should be quoted.
7. **T0 is fixed in both analyses.** It moves the FPP by ~10% here (r10). It should ideally be marginalised.
8. **Wording.** The "agreement" of the Bayesian FPP and TRICERATOPS should not be presented as independent
   confirmation.

## 8. Caveats of this re-derivation

* Several priors are carried over from the analyst and not re-verified:
  * planet occurrence classes and host-mass slopes;
  * S_close, the WD+M rate, the Occam factors and the field EB rates;
  * the BD rate f_BD = 1e-3 per star at 3-10 d (varied 3e-4 to 3e-3), which is my own order-of-magnitude choice.
* The CMD mixture is a parametric model. Its profile likelihood strongly prefers sigma_b ~ 0.40, but that preference
  is conditional on the split-normal + tail shape and on my binary q-distribution. The colour term is the more
  robust result, and it is independently confirmed by `../host/`.
* The astrometric grid and the IPD cells are phase-1 products, used here with the same wobble mapping as the analyst.
* My REF noise multipliers (fake-window err_fake/err_white, up to 2.65 in S93) may be conservative. The white-noise
  variant gives FPP 0.0016.
* The light-curve grid fixes the limb darkening at the M3V's values for every host.
* `out/r07_diagnostics.json` 'PTP_Rp_16_50_84' (6.2/6.5/6.8 R_e) uses k = sqrt(0.0157/Dc) x 1.03, which is ~8% larger
  than the grid's best k (0.119 at b ~ 0.5). With the grid k the PTP planet is ~5.7-6.2 R_e, as the analyst found.
