# Referee 2: the 'host' analysis (single metal-rich M3V or unresolved twin?)

Folder: `tic214311908/analysis/referee2_host/`

*AI-run adversarial referee check, 2026-10-17. Nothing was posted, uploaded or submitted anywhere, and nobody was
contacted. Inputs: public Gaia DR3 archive queries (r00, r06a), CDS XMatch 2MASS (r02), and the analyst's raw
downloads in `../host/data/` (`field_mdwarfs.csv` = raw Gaia-archive rows, `local_*.csv` = raw VizieR rows), used as raw
data only. Every derived quantity here is recomputed with my own code: local ridges, a label-free posterior, and my own
checks. Scripts r00-r04 were written earlier in this referee task, and r05-r09 in its final part.*

## Verdict in brief

* **Upheld:** the conclusion that the M3V is most likely a single, mildly metal-rich star, with P(twin) of order 1%.
  * My label-free estimate is **P(twin) ≈ 0.3–1%** if 2MASS H/Ks and unWISE W1 are right (r07).
  * It rises to **~1.5–2%** if Ks and W1 are both 0.03–0.04 mag too bright.
  * It rises to **~4–6%** if Ks alone is 0.04 mag too bright and W1 is ignored.
  * The analyst's 0.2–2.5% range is therefore roughly right in size, but its upper end should be ~5%.
* **Errors and overstatements in how the result was reached:**
  1. **The CMD-only twin probability is understated by about 2.7×.**
     * The analyst gives 6.3% (README: 6–10%). The correct value, with the analyst's own prior f_T = 2.35%, is
       **17%**; it is 11% for f_T = 1.5% (r07, r05).
     * The cause is h15's single-star CMD model, N(a + b[Fe/H], 0.29 mag), integrated over the kinematic [Fe/H]
       prior. The 0.29-mag scatter already contains [Fe/H] label noise, so the noise is counted twice.
     * As a result the model predicts **9.5%** of singles to be as bright as the M3V, whereas only **2.5–2.8% of all**
       clean field stars at its colour (singles and binaries together) are (r04, r05). The model's single-star density
       at ΔMKs = −0.75 is **2.0× the empirical density of all stars**.
  2. **The "several independent metallicity indicators" are essentially one measurement.**
     * The measurement is the M3V's G−IR colour excess at fixed BP−RP.
     * Once r(G−Ks) is known, the XP-shape [Fe/H] adds a factor of only **0.93–1.01** (16–84%: 0.82–1.16).
       Unconditionally it gives 0.21 (r07).
     * Joint and product likelihoods of r(G−Ks) and XP differ by about 3×: 0.072 vs 0.025 (r03).
     * r(G−W1) correlates with r(G−Ks) at r = 0.80 (r01b). Adding W1 changes the 2-D ratio only from 0.19 to 0.16
       (r07).
     * h15 multiplies XP, r(G−W1) and the wide-pair factor. This over-counts, which happens to offset error 1 roughly.
  3. **The population mixtures are not a robust independent confirmation.**
     * The analyst's h17b fit gives a global twin weight of 0.43–0.73%. That is below the analyst's own literature
       prior (1.9–2.4%), and P(twin) scales with it.
     * Window mixtures are not identifiable. In the Ks plane the twin fraction is 0–0.12, giving P(twin) 0–0.7%. In
       the W1 plane it is 0.28–0.46, giving P(twin) 2.3–4.2%, and 11–18% if W1 is 0.03 mag fainter (r01, r01c).
  4. **The wide-pair [Fe/H] (+0.13 ± 0.11) is biased high and weak.**
     * In h14, B's predicted colour is capped at BP−RP = 3.157 for [Fe/H] ≥ −0.15 (h14.log).
     * Component A has RUWE 13.7. If A is itself an unresolved binary, the inferred [Fe/H] is pushed up.
     * The calibrator outlier fraction is 35%.
     * Its likelihood factor, 0.36–0.45, should be treated as at most that, and possibly as 1.
  5. **Planet radii should come from the transit fit,** not from depth/1.04.
     * The shape analysis gives Rp = **4.94 ± 0.26 R_earth** (single, k = 0.1189) and **5.23 ± 0.25 R_earth** (twin,
       k = 0.1664).
     * The analyst gives 5.16 and 5.40–5.56. The relative change of +6% is upheld.
  6. **Minor overstatements.**
     * The RV error does not "exclude" twins inside 0.1 au. The analyst's own likelihood ratio at 0.1 au is 0.05–0.20
       (σ_CCF 11–20 km/s), and 0.25–0.70 at 0.4 au.
     * Astrometry excludes q 0.90–0.93 only near ~1 au, where L ≈ 0.02. At 0.3 and 3 au, L is 0.47–0.85 (h09.log).
       The summary's "q ≤ 0.93 at 0.4–2.4 au" is too strong.
* **New supporting evidence (missed tests):**
  * **Known unresolved near-equal pairs sit on the twin line.**
    * Three of the four Gaia DR3 SB2/SB2C M dwarfs at ΔMKs −0.67 to −1.10 have r(G−Ks) = +0.009 to +0.017. The
      fourth has −0.14. Clean stars at the same ΔMKs have +0.056 to +0.083 (r06b).
    * This validates the discriminant, although the sample is tiny.
    * The M3V (+0.095) is unlike them.
  * **The M3V's excess has the wavelength pattern of metallicity,** not of Mira-halo contamination.
    * The excess is J : H : Ks : W1 = 0.048 : 0.094 : 0.095 : 0.096.
    * The field metallicity pattern fits with χ² = 0.58 (4 points). A λ-scaled halo pattern gives χ² = 9.9, and a
      flat halo 3.0 (r09).
    * The 'd'-flagged 2MASS field stars show no bias: median −0.002 mag, 95% upper limit +0.026 (r02).
* **New caveats (missed tests):**
  * **The colour test cannot exclude a companion that supplies up to ~1/3 of the Ks light.**
    * The density ratio is ≥ 1 for a companion contribution D ≤ 0.45 mag (r09).
    * Such q ≈ 0.5–0.75 companions are limited only by Gaia astrometry, RV and IPD. A single star versus a
      full-luminosity twin is not the only choice.
  * **The kinematics lower the frequency of metal-rich singles.** The over-luminous fraction is 1.3% for |W_LSR| 25–70
    and 1.9% for v_tan 30–60 km/s, against 2.5% overall (r05). This raises P(twin) by ×1.3–1.9 for kinematic analogues
    (in r07, ×1.4–1.9).

## Key numbers (my values)

| Quantity | Value | Method | Uncertainty / calibration | File |
|---|---|---|---|---|
| M3V ΔMKs, ΔMG | −0.749, −0.685 (ZP-corrected parallax: −0.742, −0.678) | Local ridge (linear slope + KDE mode, \|ΔBP−RP\| < 0.12, clean) | ±0.03 | out/r01_population.json |
| Fraction of clean stars as bright | 2.5% (ΔMKs), 2.8% (ΔMG), n = 3808 | Empirical | Poisson ±0.3% | r01 |
| r(G−Ks), r(G−W1), r(G−H), r(G−J), r(Ks−W1), r(H−Ks) | +0.095, +0.096, +0.094, +0.048, +0.001, +0.001 | Local ridges | ±0.026–0.03 | r01 |
| DENIS−2MASS offset in image 1724791 | J +0.386 (n = 36), Ks +0.301 (n = 52); M3V +0.411 / +0.341 | Own cross-match | Robust sd 0.040 / 0.046; other images +0.08 / +0.05 | out/r04_checks.json |
| rv_error percentile | 61.6% (n = 979); 52% at matched colour; 59% among stars with > 40% deblended transits | Empirical | – | r04 |
| h15 single-star model vs data | Predicts 9.5% as bright; data 2.8% of all stars; density ratio 2.0 | Recomputed h15 model vs KDE | – | r04 |
| CMD-only P(twin) | 10.9% (f_T 1.5%), 17.1% (f_T 2.35%), 5.8–25% (f_T 0.8–3.5%) | f_T · n(x+0.753)/n(x), KDE of all clean stars | Ratio 7.27 (16–84%: 6.69–8.11) | out/r07_empirical_posterior.json |
| P(twin), (ΔMKs, r(G−Ks)) | 0.28% (f_T 1.5%), 0.44% (2.35%) | 2-D KDE ratio, dy = 0 (q = 1, twin-friendly) | Ratio 0.19 (0.14–0.24); y bandwidth 0.02 gives 0.44, i.e. 0.66–1.0% | r07 |
| P(twin), 3-D (+ r(G−W1)) | 0.25% / 0.39% | 3-D KDE | Ratio 0.16 (0.015–0.31) | r07 |
| P(twin), v_tan 30–60 subsample | 0.39–0.61% (2-D), 0.48–0.75% (3-D) | As above | Ratio 0.26 / 0.32 | r07, r05 |
| P(twin), photometric systematics | Ks +0.04: 4.0–6.3% (2-D); Ks +0.04 and W1 +0.03: 1.4–2.3% (3-D); Ks +0.03: 0.7–1.1% (3-D) | As above | – | r07 |
| XP factor given r(G−Ks) | 0.93–1.01 (0.82–1.16); unconditional 0.21 | Linear-Gaussian conditionals: normal (n = 219) vs window (n = 249) | Bootstrap 300 | r07, r03 |
| Joint vs product LR (r(G−Ks), XP) | 0.072 vs 0.025 (k = 1); 0.098 vs 0.035 (k = 1.5) | 2-D Gaussian mixture | Bootstrap | out/r03_joint.json |
| Window mixtures | Ks plane: f_T(window) 0–0.12, P(twin) 0–0.7%; W1 plane: 0.28–0.46, P 2.3–4.2% (W1 +0.03: 11–18%) | Fixed/free-width conditional mixtures | Not identifiable | out/r01c_fixedwidth.json, r01 |
| Gaia SB2s at ΔMKs −0.67 … −1.10 | r(G−Ks) = +0.016, +0.017, +0.009, −0.144; clean stars at the same ΔMKs +0.056–0.083 | Global ridges, Gaia DR3 NSS | n = 4 | out/r06b_sb2_colour.json |
| Astrometric-orbit binaries, ΔMKs −1.2 … −0.55 | Median r(G−Ks) +0.042 (n = 281) vs +0.063 clean | Same | – | r06b |
| Companion light budget | Density ratio ≥ 1 for D ≤ 0.45 mag (≤ 34% of Ks light); 0.19 at D = 0.75 | 2-D KDE vs D | – | out/r09_budget_pattern.json |
| IR-excess pattern | Metallicity χ² = 0.58; Mira halo 9.9 (λ-scaled), 3.0 (flat) | Pattern fit, 4 bands | Assumed Mira colours (shape only) | r09 |
| 2MASS 'd'-flag bias | −0.002 mag (95% upper limit +0.026), n = 12 | XMatch Cflg, clean field | Bootstrap | out/r02_tmass_dflag.json |
| Rp (from the shape analysis) | single 4.94 ± 0.26, twin 5.23 ± 0.25 R_earth | k from transit fits with density priors | – | ../shape/README.md |

## Scripts and outputs

* `scripts/r00_target_gaia.py` → `out/r00_target_gaia.json`: the M3V's Gaia DR3 row, fetched from the archive.
  * rv_nb_deblended_transits is 12 of 21.
  * rv_chisq_pvalue, rv_amplitude_robust and rv_renormalised_gof are not published, because it is a faint star.
* `scripts/r01_population.py` → `out/r01_population.json`, `figs/r01_population.png`:
  * local ridges, the target's offsets and the fraction as bright;
  * the r(G−Ks) profile against ΔMKs;
  * conditional window mixtures with a free single-star width (twin template = stars at x + 0.753), with variants;
  * a binary control: RUWE 1.4–8 and NSS stars.
* `scripts/r01b_planes.py` → `out/r01b_planes.json`: why the Ks and W1 planes disagree (window medians, correlation of
  r(G−Ks) and r(G−W1) = 0.80).
* `scripts/r01c_fixedwidth.py` → `out/r01c_fixedwidth.json`, `out/r01c.log`: window mixtures with the single-star width
  fixed to k = 1–2 × the template width, in both planes and with the photometric shifts.
* `scripts/r02_tmass_dflag.py` → `out/r02_tmass_dflag.json` (`data/r02_*.csv` are raw XMatch/VizieR rows): bias of
  2MASS 'd'-flagged photometry.
* `scripts/r03_joint_xp_colour.py` → `out/r03_joint.json`: correlation of XP-[Fe/H] with r(G−Ks), and the joint
  versus product likelihood ratio.
* `scripts/r04_rv_denis_checks.py` → `out/r04_checks.json`:
  * the rv_error percentile;
  * the DENIS frame offset, with my own cross-match;
  * the tail of h15's single-star model against the empirical tail.
* `scripts/r05_kinematics_overlum.py` → `out/r05_kinematics.json`:
  * target UVW (LSR: U 33, V 7, W 48 km/s; v_tot 59);
  * the over-luminous fraction, the component/target density ratio and the CMD-only P(twin) for kinematic subsamples.
* `scripts/r06a_sb2_query.py` → `data/sb2/nss_mdwarfs.csv`: Gaia DR3 NSS orbits of nearby M dwarfs with 2MASS
  (1018 rows; 12 SB2/SB2C).
* `scripts/r06b_sb2_colour.py` → `out/r06b_sb2_colour.json`: r(G−Ks) and ΔMKs of the SB2s and astrometric-orbit binaries
  relative to global ridges.
* `scripts/r07_empirical_posterior.py` → `out/r07_empirical_posterior.json`: the label-free posterior,
  P(twin) = f_T · n(x_t + Δx, y_t + Δy)/n(x_t, y_t).
  * It is computed in 1-D, 2-D and 3-D, with bandwidth, kinematic and photometric variants.
  * It includes the XP conditional factor and combined tables.
* `scripts/r08_binary_validation.py` → `out/r08_binary_validation.json`: colour residuals of IPD-multi-peak, RUWE-binary
  and NSS stars against ΔMKs.
  * The multi-peak sample is too small (n < 10).
  * The RUWE binaries' slope is −0.057 against −0.089 for clean stars.
* `scripts/r09_budget_pattern.py` → `out/r09_budget_pattern.json`: the companion light budget, and the wavelength pattern
  of the IR excess (metallicity versus Mira-halo contamination).

Run all scripts with `/home/user/work/venv/bin/python`, setting OMP_NUM_THREADS=2. No background jobs are left
running.

## Limitations of this referee check

* **Not re-derived:**
  * the XP spectra or the XP regression. I used the analyst's field predictions `h12b_xp_pred_field.csv` only, as data
    for the conditional test;
  * the Mira contamination of the XP spectrum;
  * the h09 astrometric simulations;
  * the h07 CCF simulation. I only checked it analytically: at 1 au R ≈ 1.2, at 0.1 au R ≈ 2.6, consistent with h07.
* **f_T (prior twin fraction surviving Gaia tests) is taken from the literature and the analyst's h15**: 0.8–3.5%.
  P(twin) scales linearly with it.
* **KDE tail densities rest on a handful of stars.** Only 0.4% of normal stars (3 of 749) are as red in r(G−Ks) as the
  M3V, so the 2-D/3-D ratios carry bootstrap ranges of a factor of 2 or more.
* **The SB2 validation has only four stars.** A larger set of literature SB2 or AO near-equal M-dwarf binaries is the
  obvious next test.
* **A high-resolution spectrum (SB2 and spectroscopic [Fe/H]), or Gaia DR4 epoch RVs and astrometry, would settle the
  question.**
