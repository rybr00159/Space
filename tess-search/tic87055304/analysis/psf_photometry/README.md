# Which star dims? Simultaneous PRF photometry of TIC 87055304 and its 9.3" neighbour

Investigator area: `/home/user/work/tic87/psf_photometry/`. All scripts are in `scripts/`. They run with
`/home/user/work/venv/bin/python`, OMP_NUM_THREADS=1, niced. Numbers are collected in `results.json`.
Ephemeris used throughout: P = 4.136195 d, T0 = 1658.41325 BTJD. Depth template: batman shape with
rp 0.0668, a/R* 4.30, b 0.85, normalised to unit depth and integrated over the exposure time.

## Verdict

**The transit signal comes from the neighbour TIC 87055315 (Gaia DR3 6720986208895208960), or from an unresolved
source within about 1" of it. It does not come from the target.** For the planet hypothesis on TIC 87055304
this result **supports a false positive**. Confidence is high. Four largely independent estimators agree in all five
sectors:

| test | target hosts | neighbour hosts |
|---|---|---|
| PSF photometry: fraction of the flux decrement on the target, f_target (expected 1 or 0) | excluded: f = 0.16 ± 0.09 (stat) ± 0.07 (sys), 7.3σ from 1 | consistent (1.4σ) |
| PSF photometry χ² (5 sectors, 4 dof) | 134.0 | 4.6 (Δχ² = 129) |
| target's own PSF light curve, depth (ppt; 4.1 expected if target hosts) | 0.64 ± 0.36 (9.6σ short) | – |
| neighbour's PSF light curve, depth relative to Gaia-predicted neighbour flux (ppt; 0 expected if target hosts) | 17.9 ± 1.6 (11σ) | 21.6 expected if all of it is on the neighbour |
| pixel-map (1-px apertures) χ², eps = 5 % PRF error (96 dof) | 188.3 | 53.2 |
| named apertures (target / neighbour / opposite pixel, 1 px and 2x2) χ² (15 dof) | 63.1 | 3.5 |
| free-position source of the pixel-map decrement (5 sectors combined) | Δχ² = 141.6 at target | 0.9" ± 0.9" from neighbour, Δχ² = 1.3 |
| pair-image motion relative to the field during transit (centroid test) χ² (10 dof) | 297.6 | 15.6 |

## 1. Method

### 1a. Static scene (`scripts/fit_static.py` → `static/`)
- **Inputs.** Median TESScut image per sector. Gaia DR3 stars with G < 19 within 330" (`gaia_330as_G19.csv`, VizieR
  I/355), propagated with proper motions to the sector epoch and placed with the cutout WCS.
- **PRF.** TESS_PRF (start_s0004 models, downloaded to `prf_files/`), interpolated to the target's CCD position.
- **Free parameters.**
  - Fluxes of all stars with T < 15.5 (30 % ridge prior to the Gaia prediction, but no prior for the target or the
    neighbour).
  - One zero-point for the fainter stars, and a planar background.
  - A global registration shift (dx0, dy0) and a Gaussian smear σ of the PRF, both fitted **only on pixels more than
    3 px from the target**, so that the registration does not depend on the target/neighbour pair.
- **Why the TESS PRF.**
  - The fitted smear is σ ≈ 0 px (0.05 px in S104).
  - The residuals within 2.5 px of the pair are 1.9–2.9 % of the peak.
  - No empirical PSF was needed. A Gaussian PSF, smeared PRFs and a PRF from a distant CCD position were run as
    robustness variants (section 4).
- **Registration accuracy.**
  - Global offsets are at most 0.11 px, with formal errors of 0.003 px.
  - Bright isolated field stars fitted with free positions scatter by 0.010–0.035 px about the registration. Projected
    on the target→neighbour axis this is 0.015–0.027 px.
  - The neighbour lies 0.45–0.47 px from the target, so a registration error ε shifts f_target by about ε/0.45.
    Measured: ±0.05 px gives ±0.11 to ±0.12 (variants `off_p05`/`off_m05`).
- **Static flux ratio.**
  - The fitted F_nb/F_t is 0.34, 0.11, 0.12, 0.15 and 0.16 (S13, S39, S66, S93, S104); Gaia/TIC predicts 0.190.
  - Fitting the pair alone with the ratio fixed to Gaia needs a pair shift of only +0.047, −0.028, −0.028, −0.021 and
    −0.013 px along the neighbour axis. This is consistent with the field-star scatter.
  - Because the PSF-normalised depth of the neighbour inherits this bias, I also quote depths **normalised to the
    Gaia-predicted flux** ("cat"), i.e. ΔF divided by F_pair·r/(1+r).

### 1b. Per-cadence photometry (`scripts/psfphot.py` → `lc/<variant>/psf_s*.npz`)
1. **Per-cadence registration** d(t), from the field pixels (more than 3 px from the target) using a second-order
   Taylor model of the static scene: I = a(M + d·∇M + ½ dᵀ∇²M d) + plane.
   - Pointing jitter has rms 0.013–0.023 px in x and 0.032–0.060 px in y.
   - The pair region (r < 2.5 px) is fitted the same way to give the pair's own shift, which is used in the centroid
     test.
2. **Flux fit.** In pixels within 3.5 px of the target, the per-cadence weighted least-squares fit (1/FLUX_ERR²) has:
   - free fluxes for the target, the neighbour, stars with T < 15.5 within 3.8 px and stars with T < 13.5 within 4.5 px
     (12–14 free stars in total);
   - all other stars fixed at their static fluxes and moved with d(t);
   - a planar background.
3. **Quality cuts.** Only QUALITY == 0 cadences are used. Cadences are masked when the background exceeds the median
   + 6 MAD (this removes scattered light, including the start of S104) or when χ² exceeds 5× the median.

### 1c. Depths (`scripts/depthfit.py`, `scripts/analyze_lc.py`)
- **Estimator.** A joint fit within ±0.45 d of each transit with a per-transit quadratic baseline and a shared
  amplitude A of the unit-depth template. A transit is used only if at least 50 % of its in-transit and 25 % of its
  out-of-transit cadences are present. The fit is run on:
  - absolute fluxes, which gives ΔF in e-/s;
  - PSF-normalised fluxes;
  - the pair sum.
- **Red-noise-aware errors.** The adopted error is the maximum of four estimates:
  - the white-noise error;
  - the per-transit scatter divided by √N;
  - the scatter of the same estimator at 65 control phases that avoid the transit and phase 0.5;
  - the scatter of single-epoch amplitudes at the control phases divided by √N.

  The control phases also give the target–neighbour covariance, ρ = −0.88 to −0.93.
- **Combining sectors.** f_target is obtained from a GLS profile likelihood in f, with ΔF_t = f·S_s and
  ΔF_n = (1−f)·S_s, where S_s is free for each sector.

## 2. Results: PSF light curves (nominal; `lc/nominal/depths.json`, figure `fold_psf_nominal.png`)

| sector | N_tr | ΔF_t (e-/s) | ΔF_n (e-/s) | f_target | target depth, PSF-norm (ppt) | neighbour depth, cat-norm (ppt) | pair depth (ppt) |
|---|---|---|---|---|---|---|---|
| 13 | 6 | −1.6 ± 3.3 | 14.3 ± 3.9 | −0.12 ± 0.26 | −0.51 ± 1.08 | 21.7 ± 5.9 | 3.07 ± 0.43 |
| 39 | 5 | 3.7 ± 3.0 | 12.6 ± 2.3 | 0.23 ± 0.17 | 0.99 ± 0.78 | 19.5 ± 3.6 | 4.04 ± 0.35 |
| 66 | 4 | 2.7 ± 3.6 | 12.1 ± 2.8 | 0.18 ± 0.23 | 0.72 ± 0.96 | 18.6 ± 4.3 | 3.62 ± 0.40 |
| 93 | 4 | 2.1 ± 3.0 | 10.9 ± 2.2 | 0.16 ± 0.22 | 0.61 ± 0.85 | 16.9 ± 3.3 | 3.22 ± 0.31 |
| 104 | 5 | 2.6 ± 2.0 | 11.0 ± 1.6 | 0.19 ± 0.14 | 0.75 ± 0.58 | 17.0 ± 2.5 | 3.35 ± 0.23 |
| all | 24 | 2.2 ± 1.3 | 11.7 ± 1.0 | **0.16 [0.07, 0.25]** | **0.64 ± 0.36** | **17.9 ± 1.6** | **3.44 ± 0.14** |

- **Expectations.**
  - If the target hosts, the target's PSF depth should be about 4.1 ppt, i.e. 3.44 × (1 + 0.19), and the neighbour's
    should be zero.
  - If the neighbour hosts, the neighbour's depth should be about 21.6 ppt, i.e. 3.44 × (1 + 1/0.19) relative to its
    Gaia flux, and the target's should be zero.
- **PSF-normalised neighbour depths.** These are 13.6, 45.5, 31.6, 19.5 and 18.9 ppt. They scatter between sectors
  (χ² = 14.7 for 4 dof) only because the PSF flux normalisation is biased (section 1a); the cat-normalised depths
  agree (χ² = 0.9 for 4 dof).
- **Per-transit amplitudes.** Shown in `per_transit_dF.png`. The 24 single transits cluster near (ΔF_n, ΔF_t) ≈
  (12, 2) e-/s along the anti-diagonal of flux-transfer noise. One S13 transit falls at (−4, 14).

**Out-of-transit validation** (`scripts/validate_oot.py`, `lc/nominal/oot_validation.json`):
- The point-to-point scatter of F_t and F_n is 1.0–1.3 times the formal LSQ errors (e.g. S104: 20.1 versus 15.0 e-/s).
- The empirical correlation of F_t and F_n is −0.85 to −0.88, against −0.80 to −0.91 formal.
- The rms of 4.5-h bins exceeds the white-noise expectation by a factor β = 1.3–3.0. This red noise is covered by the
  control-phase errors.
- Slow anti-correlated trends between F_t and F_n are flux transfer from residual registration error. The per-transit
  quadratic baselines remove them.

**Robustness variants** (`scripts/run_variants.sh`, `variants_f_target.png`): combined f_target with its 1σ range and
Δχ²(target − neighbour).

| variant | f_target | Δχ² |
|---|---|---|
| nominal | 0.16 [0.07, 0.25] | +129 |
| plus regressors dx, dy, bkg | 0.12 [0.02, 0.21] | +134 |
| no per-cadence registration | 0.17 [0.04, 0.29] | +53 |
| free per-cadence registration (degenerate; see note) | 0.20 [0.06, 0.33] | +39 |
| pair registered with Gaia ratio | 0.13 [0.04, 0.21] | +148 |
| Gaussian PSF (σ fitted) | 0.15 [0.05, 0.24] | +124 |
| PRF of CCD position (300, 300) | 0.03 [−0.08, 0.13] | +114 |
| PRF smeared by 0.25 px | 0.13 [0.03, 0.23] | +103 |
| PRF smeared by 0.50 px | −0.07 [−0.22, 0.07] | +72 |
| R_fit = 2.5 px | 0.10 [0.00, 0.19] | +103 |
| R_fit = 5 px | 0.17 [0.08, 0.25] | +127 |
| T_free = 14 | 0.08 [0.02, 0.14] | +265 |
| T_free = 16.5 | 0.18 [0.09, 0.27] | +71 |
| registration +0.05 px toward the neighbour | 0.28 | +86 |
| registration −0.05 px | 0.04 | +165 |
| registration +0.10 px | 0.39 | +51 |
| registration −0.10 px | −0.06 | +217 |
| background mask at 30 MAD | 0.13 [0.03, 0.22] | +131 |

- Note on free registration: free per-cadence shifts of the pair are degenerate with the target/neighbour flux split,
  which is why its errors are large.
- **Systematic error.** Across the non-deliberate variants f_target has mean 0.115 and standard deviation 0.075, so I
  adopt 0.07.
- **What H_target would need.** Making the target the host requires a registration error of about +0.37 px (about 8")
  toward the neighbour, against measured registration uncertainties of 0.02–0.03 px.

**Injection-recovery and leakage** (`scripts/injection.py`, `injection/`, `injection_leakage.png`):
- **Set-up.** The template transit was injected into the pixel data of ONE star at phase 0.27 or 0.73, with
  ΔF = 0.0036·(F_t + F_n), about 15 e-/s. The full pipeline was then re-run. In total 5 sectors × 2 phases × 2 hosts
  × 6 injection models = 120 runs.
- **Same PRF.** Recovery is 1.000 ± 0.003 and leakage is at most 0.005, as expected for a linear fit.
- **Mismatched models.** The injection used a 0.4-px smeared PRF, a Gaussian PSF, the PRF of another CCD position, or
  the star's position offset by ±0.05 px:
  - target-injected signals give f_target between 0.73 and 1.12;
  - neighbour-injected signals give f_target between −0.18 and +0.21.
- **Comparison with the data.** The real value, 0.16 ± 0.09, lies inside the neighbour-injected range and well away
  from the target-injected range.
- **Noisy recovery** with the same PRF: injection into the target gives f = 0.98; into the neighbour, f = −0.04.

## 3. Results: aperture photometry versus aperture position (`scripts/aperture.py`, `aperture/`)

- **Apertures.**
  - Single pixels and 2x2 boxes centred on (or nearest to) the target, the neighbour and the point opposite the
    neighbour, x_opp = 2x_t − x_n.
  - A 3x3 box on the target.
  - Every single pixel within 2.5 px of the target (the pixel map).
- **Photometry.** Aperture flux minus the per-cadence PSF-fit background plane. The depth is the ΔF from the same
  estimator; the relative depth is ΔF divided by the median aperture flux.
- **Predictions.** ΔF_ap = S·frac_host(ap), with frac from the PRF at the in-transit registration and S free.
- **Errors.** The control-phase covariance across apertures is shrunk 50 % towards the diagonal and scaled to the
  adopted errors. A PRF-model error term of eps × prediction is added, with eps = 0, 5 or 10 %.
- **Duplicate apertures.** Where two apertures are identical, for example when the target and the neighbour fall in the
  same pixel in S13, S39 and S66, only one is kept.

| sector | named apertures χ²: target / neighbour (dof) | pixel map χ²: target / neighbour (dof) |
|---|---|---|
| 13 | 7.0 / 0.7 (3) | 16.0 / 5.9 (20) |
| 39 | 24.6 / 1.0 (3) | 79.8 / 15.4 (20) |
| 66 | 3.3 / 0.2 (2) | 34.0 / 14.2 (20) |
| 93 | 13.6 / 0.7 (4) | 34.9 / 12.6 (18) |
| 104 | 14.5 / 0.9 (3) | 23.6 / 5.1 (18) |
| total | **63.1 / 3.5 (15)** | **188.3 / 53.2 (96)** |

All values use eps = 5 %. With eps = 0 % the totals are 78.5 / 4.1 and 194.3 / 54.2; with eps = 10 % they are
45.0 / 2.6 and 174.0 / 50.8. The neighbour hypothesis fits at χ² below the degrees of freedom, so the adopted errors are
conservative.

**Example depths** (ppt):
- S93: 1-px target 1.83 ± 0.67, 1-px neighbour 4.20 ± 0.71, 1-px opposite 1.15 ± 0.87.
- S104: 2x2 opposite 2.03 ± 0.29, 2x2 target 2.89 ± 0.22, 1-px neighbour 4.07 ± 0.81.

The depth rises toward the neighbour, as the neighbour hypothesis predicts.

**Figures.** `depth_vs_aperture.png` and `pixel_depth_vs_prediction.png`. In the second, measured pixel depths follow
the 1:1 line for the neighbour hypothesis and scatter about it for the target hypothesis.

**Alternative aperture runs.**
- `aperture_regs/` adds a quadratic jitter decorrelation in dx, dy. It gives the same conclusion: named 75.8/3.8,
  pixel map 217.3/51.0 (target/neighbour).
- `aperture_nobkg/` omits per-cadence background subtraction. Errors are 2–10× larger, and the conclusion is the same
  but weaker: named 16.8/1.2, pixel map 49.3/11.5 (target/neighbour). Its combined localisation is in
  `aperture_nobkg/localisation_combined.json`: 0.9" from the neighbour, Δχ² at the target = 38.6.

**Other stars as hosts** (pixel map, summed χ² over the 5 sectors):
- The neighbour is best, at 53.2.
- The next is Gaia DR3 6720986208880830336 (T = 17.35, 11.6"), at 91.4. It would need a 56–83 % deep eclipse.
- Every other star within 60" scores 164 or worse, or would need an eclipse deeper than 100 %. Examples:
  - the G 18.1 star at 5.5" (χ² 192);
  - the target itself (188).
- Any star fainter than T ≈ 17.6 cannot supply ΔF ≈ 12–14 e-/s even with a 100 % eclipse.

**Free-position localisation** (`scripts/localise_combined.py`, `localisation_combined.png`, `pixelmap_localisation.png`).
A single PRF at a free position fitted to the pixel decrements, with the five sectors combined in sky coordinates:
- best source position (E, N) = (5.8, 6.4)", with 1σ ranges E 5.0–6.8" and N 5.6–7.4";
- the neighbour lies at (6.65, 6.71)": 0.9" from the best position, Δχ² = 1.3;
- the target lies at Δχ² = 141.6, i.e. 11.6σ for 2 dof.

**Centroid test** (`scripts/centroid_test.py`, `centroid/`). This test does not use the flux split.
- **Observed shift.** During transit the pair image moves away from the neighbour by 0.80–1.19 mpx (±0.17–0.23)
  relative to the field, in every sector.
- **Predictions.** Applying the same linear fit to model difference images predicts 0.63–1.46 mpx away from the
  neighbour if the neighbour hosts. If the target hosts, it predicts 0.17–0.84 mpx toward the neighbour.
- **Fit.** Total χ² = 15.6 for the neighbour and 297.6 for the target (10 dof). The fraction of the way from the
  target to the neighbour is 1.22, 0.90, 0.85, 0.81 and 0.84.

## 4. The neighbour's own S13 QLP and eleanor-lite light curves (`scripts/neighbour_lcs.py`, `neighbour_lcs/`, `neighbour_own_lcs_s13.png`)

- **eleanor-lite.** The files for TIC 87055304 and TIC 87055315 contain **bit-identical** RAW_FLUX, CORR_FLUX and
  PCA_FLUX, with the same aperture ('1p25_circle_center') on the same postcard pixels. Both stars sit in one pixel, so
  the depth is identical: 2.93 ± 0.51 ppt raw, 3.09 ± 0.44 ppt corrected. **These light curves carry no information
  on the host.**
- **QLP depths.** SAP depth is 3.92 ± 0.27 ppt for the target and **21.65 ± 1.44 ppt** for the neighbour, a ratio of
  5.53 ± 0.52. The two QLP aperture centres differ by only (0.07, 0.15) px.
- **Why the neighbour's QLP curve shows a 2 % dip.** The depth ratio equals three other numbers:
  - the ratio of QLP's flux errors, 5.39;
  - the TIC flux ratio of the two stars, 10^(0.4·1.82) = 5.36;
  - the ratio found in transit/nb_lc_check.json, 5.53 ± 0.47.

  Each star's light curve is therefore built from an almost identical aperture flux but normalised to (or
  dilution-corrected by) that star's own catalogue flux. The ~2 % dip of the neighbour is the shared aperture dip scaled
  by about F_t/F_n. QLP would produce it whichever star hosts the signal, so **it is not evidence for either host.**
- **QLP target depth.** The same scaling explains why QLP's target SAP depth (3.9 ppt in S13) exceeds the eleanor raw
  depth (2.9 ppt) and the PSF pair depth (3.1 ppt in S13). QLP SAP depths, including the 3,580 ppm combined value in
  CONTEXT, are therefore dilution-scaled numbers.
- **Coincidence with the true depth.** The neighbour's true intrinsic depth from the PSF photometry, about 18–22 ppt,
  happens to be similar to QLP's 21.7 ppt, but for a different reason.

## 5. Secondary eclipse and odd/even in the separated light curves (`scripts/secondary_oddeven.py`, `secondary_scan_psf.png`)

Depths below are relative to the Gaia-predicted flux of each star, combined over 5 sectors.
- **Secondary at phase 0.5.**
  - Neighbour: 0.30 ± 1.87 ppt, against a primary of 17.9 ± 1.6 ppt. The 2σ upper limit is 4.0 ppt, i.e.
    secondary/primary < 0.23.
  - Target: 0.18 ± 0.42 ppt.
- **Eccentric secondary.** A scan over phases 0.30–0.70 (81 phases, primary duration) finds no feature: the maximum
  |S/N| is 1.55 for the neighbour (at phase 0.37) and 1.28 for the target. The 2σ upper limit anywhere is 7.0 ppt for
  the neighbour and 1.5 ppt for the target.
- **Odd/even.**
  - Neighbour: odd 16.2 ± 2.4, even 20.2 ± 2.2 ppt, difference −4.0 ± 3.2 ppt (1.3σ).
  - Target: odd 0.98 ± 0.54, even 0.30 ± 0.51 ppt.
- **Interpretation.** The signal on the neighbour is not an equal-mass eclipsing binary at P = 4.136 d, which would
  give a secondary of similar depth. Other stellar companions are not excluded:
  - a low-luminosity eclipsing companion (for example an M dwarf, whose secondary would be shallower than 4 ppt);
  - an EB at 2P with nearly equal components, which would show odd ≈ even; this is only weakly constrained here.
- **Size of the companion (outside this task's scope).** A depth of about 2 % on a star with R = 1.9 R_sun (TIC)
  implies a companion of about 0.27 R_sun. That is too large for a planet, so a stellar or brown-dwarf eclipser of the
  neighbour is the natural reading.

## 6. Comparison with transit/ (read only after all of the above was finished; `comparison_transit.json`)
- **transit/localize_twosource.json.** It gives f_target = 0.02 ± 0.04, against 0.16 ± 0.09 (stat) here. The
  conclusion is the same.
  - Its per-sector values scatter by χ² = 11.8 for 4 dof, so its errors look about 1.7× too small.
  - My `tfree14` variant, which fixes the 1.2-px T = 14.7 star as that analysis does, gives 0.08 ± 0.06.
- **transit/aperture_depth.json** (concentric apertures on the target). Its total χ² is 191.0 for the target and 64.1
  for the neighbour (40 apertures), the same direction as my 63.1 / 3.5 (named) and 188 / 53 (pixel map).
- **transit/nb_lc_check.json.** It agrees with section 4: QLP ratio 5.53 ± 0.47, eleanor ratio 1.0000.
- **transit/prf_localize.json** places the source *between* the stars, at (E, N) = (4.2, 3.3)" and about 7σ from the
  neighbour.
  - Its PRF/WCS registration was fitted on the out-of-transit image with free target and neighbour fluxes. That fit is
    degenerate with the pair centroid along the target–neighbour axis.
  - In S104 its x-registration (along the separation there) differs from my field-star registration by about 0.24 px,
    while its y-registration (perpendicular) agrees: 0.112 against 0.105 px.
  - With field-star registration, which excludes the pair, the source is 0.9 ± 0.9" from the neighbour. I regard the
    prf_localize position as biased by its registration.

## 7. Corrections to earlier claims
1. **"The 9" neighbour is excluded at >6σ"** (src/centroid.py, earlier README) is wrong. The neighbour is the
   preferred host, and the target is excluded:
   - PSF Δχ² = 129;
   - pixel-map Δχ² = 142;
   - centroid χ² 298 against 16.
2. **Separation.** It is 9.33" at PA 45° (neighbour to the NE) between the Gaia DR3 epoch-2016 positions, and
   9.37–9.44" at the TESS epochs. The 8.94" quoted in CONTEXT mixes the TIC (DR2-based, J2000) target position with
   the DR3 neighbour position. In TESS pixels the separation is 0.45–0.47 px.
3. **The ~2 % dip in QLP's neighbour light curve** is a normalisation effect, not host evidence (section 4). The
   eleanor-lite light curves of the two TICs are identical.
4. **QLP SAP depths** (3,580 ppm combined) are dilution-scaled. The depth of the blended target+neighbour pair from PSF
   photometry is 3.44 ± 0.14 ppt; the eleanor raw depth in S13 is 2.9 ppt.

## 8. Open issues / limits
- f_target = 0.16 ± 0.09 is 1.8σ above zero. This is within the PRF and registration systematics: variants range from
  −0.07 to 0.20, and neighbour-injected signals with mismatched PRFs give up to 0.21. With all errors, at most about
  40 % of the decrement can be on the target (2σ).
- The PSF fluxes of the neighbour, and so its PSF-normalised depths, depend on registration at the 0.03-px level. Use
  the cat-normalised depth, 17.9 ± 1.6 ppt for the neighbour's share, or about 21.6 ppt if all of the signal is on the
  neighbour.
- The source position is known to about ±0.9". A companion of the neighbour within about 1" cannot be distinguished
  from the neighbour itself. Anything at the target's position is excluded.
- **Decisive and cheap follow-up.** Seeing-limited ground-based photometry at a predicted transit resolves the 9.3"
  pair directly. The neighbour should dim by about 2 % and the target should not.

## Files
- **Scripts** (`scripts/`):
  - `common.py` — data, Gaia, PRF and template helpers;
  - `get_gaia.py`;
  - `fit_static.py`;
  - `psfphot.py`, run through `run_variants.sh`;
  - `depthfit.py`, `analyze_lc.py`;
  - `injection.py`;
  - `aperture.py` (environment variables AP_REGS, AP_BKG), `localise_combined.py`;
  - `centroid_test.py`, `validate_oot.py`;
  - `neighbour_lcs.py`, `secondary_oddeven.py`;
  - `figures.py`, `final_summary.py`.
- **Figures:**
  - `fold_psf_nominal.png` — the two PSF light curves folded per sector. There is also a `fold_psf_<variant>.png` for
    each variant.
  - `injection_leakage.png`, `variants_f_target.png`, `per_transit_dF.png`;
  - `depth_vs_aperture.png`, `pixel_depth_vs_prediction.png`;
  - `pixelmap_localisation.png`, `localisation_combined.png`;
  - `neighbour_own_lcs_s13.png`, `secondary_scan_psf.png`;
  - `static/static_s*.png`, `diag/s66_timeseries.png`.
- **Data products:**
  - `lc/<variant>/psf_s*.npz` — per-cadence fluxes, covariances, background and registration;
  - `lc/<variant>/depths.json`;
  - `results.json` — all numbers.
