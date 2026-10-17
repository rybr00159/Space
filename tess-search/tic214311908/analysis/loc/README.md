# TIC 214311908: which star dims? Independent pixel-level attribution

Folder: `tic214311908/analysis/loc/`
(AI-run analysis, 2026-10-16. Nothing was submitted or posted anywhere.)

## Bottom line

- **The pair and the Mira are excluded as the source.** Both of my image methods say so when the sectors are
  combined, and individually in S39, S93, S103 and S104. S12 and S66 are ambiguous on their own.
  - Combined difference-image χ² penalties relative to the M3V: pair A +133, pair B +150, Mira +275. The
    amplitude images give +65, +73 and +165.
  - In 1,000 injections of the observed signal on pair A, pair B or the Mira, none reaches these values in either
    method. Under those hypotheses the medians are −83, −91 and −153 (DIFF).
  - Injections on the M3V reproduce the real values: P(≤ real) = 0.31, 0.41 and 0.42 (DIFF) and 0.17, 0.17 and
    0.31 (AMP).
  - With the Mira's own empirical image as the PSF model, the Mira fits *worse*, not better (+291 instead of
    +275). PRF mismatch of the Mira does not explain its rejection.
- **The source is on the M3V, or within about 3″ of it.**
  - Combined free-fit position (difference images, M-dwarf colour): (+0.6, +0.4)″ from the M3V. The injection
    scatter for a source on the M3V is 1.0″ per axis.
  - 99.7% region: within 5.1″ of the M3V (difference images) or 9.2″ (amplitude images).
  - The M3V is the only catalogued star inside the DIFF region. The AMP 99.7% region also contains the G = 17.5
    star described below.
  - Posterior with equal priors over the 93 catalogued stars, injection-calibrated:
    - Difference images: P(M3V) = 0.99999, P(pair A+B) ≈ 6×10⁻²⁰, P(Mira) ≈ 2×10⁻⁴⁰, P(all others) ≈ 8×10⁻⁶,
      almost all of it on one G = 17.5 star.
    - Amplitude images: P(M3V) = 0.977, P(that G = 17.5 star) = 0.023, P(pair) ≈ 3×10⁻⁶, P(Mira) ≈ 5×10⁻¹⁵.
    - At such extremes the injection reliability test shows the posteriors are somewhat over-confident; see the
      Calibration section.
- **The one catalogued alternative not cleanly excluded is Gaia DR3 5925373272827260544.** It is a G = 17.5,
  BP−RP = 1.32 star 3.1″ from the M3V, with parallax 0.65 ± 0.11 mas (about 1.5 kpc, so a K dwarf). It was not
  mentioned in the earlier vetting.
  - It is disfavoured by Δχ² = +35 (difference images) and +19 (amplitude images). Injections on it reach the real
    value in 0/1000 and 10/1000 trials.
  - That discrimination comes mostly from Sector 39 and from the lateral-colour shift:
    - If the star were as red as the M3V, Δχ² drops to +8 and +7.
    - Without Sector 39 it drops to +17 and +5.
  - It is mainly excluded by **physics**: it would need a 51% deep, flat-bottomed, 1.5-h dip. No dark transiting
    body can produce that. Only a total eclipse of a hot compact companion supplying more than half of its light
    could.
  - A second star, G = 17.7 at 5.7″ (Gaia DR3 5925373272827262208), is excluded by both pixels (Δχ² +80 / +51)
    and physics (71% depth).
- **No evidence for an uncatalogued source.**
  - The best catalogued star (the M3V) fits as well as a free-position source: Δχ² 0.3 (difference images) and 1.2
    (amplitude images). For a source truly on the M3V, P(≥ these) = 0.78 and 0.63.
  - An uncatalogued source would have to lie within about 2–3″ of the M3V to be indistinguishable in TESS.
    - At 4.3″ the M3V would be excluded at 99.73% in half of the injection trials; at 6″, in 92% (difference images).
    - With the noisier amplitude images these distances are 7.1″ and more than 8″.
- **Only the M3V and the M dwarfs of the pair can physically host the dip.**
  - The required depth is 1.7% on the M3V (5.4 R⊕), 2.3% on pair A (5.2 R⊕) and 5.8% on pair B (6.1 R⊕).
  - Every catalogued star within 12″ of the M3V other than these needs ≥ 36% (51% and 71% for the two nearest).
    That is impossible for a dark body making a flat-bottomed dip whose flat part is ≥ 40% of the duration.
  - The brighter field stars, which would need only 0.3–30%, lie 23–95″ from the M3V. The pixels exclude them by
    Δχ² ≥ 451 (DIFF) and ≥ 267 (AMP). They also fail the circular-orbit test: a flat-bottomed 1.5-h transit at
    5.69 d needs a small, dense star.
  - The Mira (hundreds of R☉) cannot host a 5.69-d orbit at all.
- **Corrections to earlier numbers.**
  1. The Mira's colour shift is **5–7″, not 16–17″.**
     - Measured against a WCS calibrated on normal-colour stars, the Mira's image is displaced 5.0–7.4″ radially
       away from the camera boresight, with tangential components ≤ 0.6″, in all six sectors.
     - This matches the extrapolated lateral-colour curve (6.1–7.2″).
     - So the 16–17″ calibration error described in the vetting README was not caused by the colour shift alone.
     - The conclusion drawn there (calibrate on normal-colour stars) is still right.
  2. The Mira-flux-scaling argument is weak. Under the Mira hypothesis, the per-sector amplitudes compared with
     the Mira's TESS flux (6,100–13,900 e⁻/s in my fits) give χ² 11.2 for proportional against 8.6 for constant
     (5 degrees of freedom, difference images), and 5.7 against 6.8 (amplitude images). It does not discriminate;
     the pixel positions do.

## Data and calibration (steps 0–2)

### Data preparation
- **Data.** TESScut 21×21 px cutouts, Sectors 12, 39, 66, 93, 103 and 104.
  - Times: TESScut TIME with the FFI-direction barycentric correction replaced by the Roemer delay towards the target
    (Earth ephemeris, error < 1.3 s). The corrections are −19 to +31 s, negligible here.
  - Ephemeris: P = 5.6914968 d, T0 = BTJD 3840.9153, T14 = 1.63 h, flat fraction 0.7.
  - Transits used: 3, 4, 3–4, 3, 3–4 and 4 per sector (others fall in gaps).
- **Background.** A per-cadence background plane is fitted to the faintest 30% of pixels and subtracted.
  - This is essential. Without it, fake-window noise is 4–37 e⁻/s per pixel, 93–100% of it a uniform offset
    (`figs/noise_pca.png`).
  - With it, the pixel noise in fake windows is 1.02–1.40× the formal error and pixel-to-pixel correlations are
    ~0.03–0.08 (median).
- **WCS** (`c01_wcs.py` → `out/wcs_cal.json`).
  - Every Gaia star with T < 15.5 in the cutout is fitted with the TESS PRF ('PRF' package) at its Gaia position,
    propagated with proper motion and shifted by the lateral-colour curve. A common offset is fitted; the Mira's
    core is masked.
  - Offsets are −0.11 to +0.05 px in x and +0.07 to +0.18 px in y, with formal errors ±0.01 px.
  - Isolated calibrators agree to ≤ 0.02 px.
  - Red calibrators (BP−RP 1.83, 1.91, 2.41) show radial residuals of −0.14″ to +0.35″ after the colour curve, so
    the curve works for these cameras.
- **Lateral-colour curve.**
  - Taken from the TIC 71316387 campaign's whole-CCD field-star fits: −1.55″/mag up to BP−RP ≈ 2, steepening to
    −3.4″ at 2.48 and −3.95″ at 2.94.
  - It is scaled with field angle as θ^0.35. Field angles here are 8.6–13.7°.
  - Shifts: M3V (2.76) ≈ −3.5 to −4.1″; pair A (2.88) ≈ −3.6 to −4.3″; pair B assumed BP−RP 3.3 (no Gaia colour).
  - Uncertainty: about 0.5″ for C > 2.3. In the injections this is propagated as a random position perturbation;
    variants scale the curve by ×0.7 and ×1.3.
- **The Mira's apparent position** (`c02_mira.py` → `out/mira_pos.json`, `figs/mira_direct_s*.png`).
  - The direct image within 3.5 px is fitted with the Mira's position and flux free, all other stars fixed at
    catalogue flux, and the reference star and background free.
  - Offsets from Gaia (″), with radial component (negative = away from the boresight) and tangential component:

    | Sector | Offset (E, N) | Radial | Tangential |
    |---|---|---|---|
    | 12 | (+3.0, −4.3) | −5.2 | −0.6 |
    | 39 | (+2.9, −6.9) | −7.4 | +0.5 |
    | 66 | (−1.2, −5.2) | −5.4 | +0.2 |
    | 93 | (−1.7, −5.9) | −6.1 | +0.2 |
    | 103 | (−2.0, +4.7) | −5.0 | +0.6 |
    | 104 | (−6.7, −1.1) | −6.8 | −0.3 |

  - Formal errors are ±0.4–0.7″. Across 18 fit variants the offsets move by 1–2″, and 2-d chunks scatter by 0.3–1.2″.
  - The Mira's fitted TESS flux is 6,100–13,900 e⁻/s (T ≈ 10.1–11.0).
  - Fit residuals are 4.4–5.8% rms of the pixel flux in the core. That is the Mira's PRF mismatch: ±8% in single
    pixels, or about 800 e⁻/s.
- **PRF match for normal stars** (`c03_prfcheck.py`): residuals are 1–3% of the peak (S103: up to 7%).

### Image methods (`l1img.py`, `l2fit.py`)
- **DIFF: difference images.**
  - In-transit cadences: exposure-integrated template ≥ 0.75 (flat bottom).
  - Out-of-transit windows on both sides, from 0.5·T14 + 0.5·exposure + 12 min to that plus 0.25 d.
  - Per pixel, a straight line through the out-of-transit cadences (removes the Mira's slow variability) and
    D = (line − flux_in) / ⟨template⟩, so D is the full depth in e⁻/s.
  - Transits are combined with inverse-variance weights.
- **AMP: amplitude images.**
  - Per pixel, a regression on a cubic B-spline (knots every 0.5 d, per orbit segment) plus the exposure-integrated
    trapezoid template, solved as α = τᵀMy / τᵀMτ.
  - This is a global baseline, unlike DIFF's local lines.
- **Linearity.** Both methods are linear, so injection = noise image + A·PSF.
  - Checked by injecting into the full cube and running the whole pipeline (`r03_cubecheck.py`): response = A·PRF
    to < 1.1% of the peak per pixel.
  - Summed over all pixels there is a −2 to −3 e⁻/s uniform offset from background-plane leakage, which the free
    constant absorbs.
- **Noise maps.**
  - Formal per-pixel errors are inflated by the empirical ratio from single fake windows: 1.07–1.20 (DIFF) and
    1.2–1.3 (AMP), with mild per-pixel smoothing.
  - Fit region: pixels within 5.5 px of the Mira.
- **Model.** d = a·PSF(apparent position) + b for each hypothesis, with a weighted closed-form fit.
  - Hypotheses: 93 Gaia stars within 80″ of the Mira with T < 17.6 (any fainter star cannot lose 15 e⁻/s). The Mira
    is placed at its measured apparent position.
  - Free fit: a 1″ sky grid (±100″) co-moving with the M3V, for a source of BP−RP 2.76 (variant 1.3), refined
    quadratically.

## Results

### 1. The signal in every sector (`e02_noise.py`, `e05_perwindow.py` → `out/per_window.csv`, `figs/diff_images.png`)
- **Flux lost**, summed within 2.5 px of the M3V (DIFF, e⁻/s):

  | Sector | Flux lost |
  |---|---|
  | 12 | 16.8 ± 4.6 |
  | 39 | 13.6 ± 1.4 |
  | 66 | 17.0 ± 2.1 |
  | 93 | 15.1 ± 1.8 |
  | 103 | 8.9 ± 2.8 |
  | 104 | 11.5 ± 1.7 |

  AMP agrees: 17.4, 13.1, 17.9, 15.2, 6.9 and 12.3.
- **Individual transits.** The fitted amplitude on the M3V is 9–24 e⁻/s in all transits but one. The exception,
  S103 epoch 56, gives 5 ± 3 e⁻/s in both methods. In DIFF, the pair fits worse than the M3V in 19 of the 21
  transits and the Mira in all 21.
- **Fit quality.** A single PRF on the M3V leaves noise-like residuals in every sector: χ² = 86/90, 90/94, 86/93,
  74/94, 79/91 and 98/94. Pair and Mira fits leave dipoles (`figs/diff_images.png`).

### 2. Hypothesis χ² per sector and combined (`r01_sim.py`, `r02_post.py` → `out/post_{diff,amp}_prf_all.json`)
χ² minus that of the best catalogued star (DIFF / AMP):

| Sector | M3V | pair A | pair B | Mira | G17.5 star at 3.1″ | G17.7 star at 5.7″ |
|---|---|---|---|---|---|---|
| 12 | 4.1 / 5.8 | 10.4 / 8.5 | 14.7 / 12.5 | 24.6 / 23.8 | 1.5 / 4.7 | 2.4 / 6.7 |
| 39 | 0 / 0 | 45.0 / 32.3 | 49.1 / 35.8 | 103.8 / 77.2 | 18.3 / 13.7 | 38.3 / 29.2 |
| 66 | 0.3 / 0 | 1.0 / 5.0 | 2.4 / 5.5 | 10.8 / 19.6 | 0 / 1.6 | 4.7 / 6.2 |
| 93 | 0 / 0 | 21.7 / 6.5 | 25.3 / 7.3 | 51.2 / 15.6 | 2.6 / 1.0 | 11.8 / 3.5 |
| 103 | 0 / 1.9 | 17.4 / 5.6 | 20.1 / 4.7 | 30.1 / 14.7 | 2.2 / 0 | 2.3 / 3.7 |
| 104 | 0 / 0 | 41.8 / 15.0 | 42.9 / 14.9 | 58.8 / 21.6 | 15.1 / 5.5 | 24.5 / 9.4 |
| **combined** | **0 / 0** | **132.8 / 65.1** | **150.1 / 72.9** | **274.7 / 164.6** | **35.3 / 18.8** | **79.6 / 51.0** |

S12 (30-min cadence, 3 transits) and S66 are individually ambiguous. In S12 the best-fitting star is a G = 17.0 star
7″ from the Mira; the M3V is 4.1 (DIFF) and 5.8 (AMP) worse.

### 3. Injection–recovery calibration (item 3)
- **Design.**
  - Noise comes from fake transit windows at random times away from the real ones. Each draw has the same number
    of windows as the real data, at least 0.7 d apart, from 259–350 candidate windows per sector.
  - Each of the 93 stars is injected 1,000× with the observed total flux loss (14.3 e⁻/s DIFF, 14.5 AMP).
  - Each injection's position is perturbed by the systematic uncertainty: WCS 0.02 px; colour curve 0.15–0.5″
    radial; Mira 0.8″; pair B 1″.
  - Sector draws are paired by index, so the sector sum is a proper 6-sector realisation.
- **Distributions of the combined χ² differences** (2.5 / 50 / 97.5%):

  | True source | pairA − M3V | Mira − M3V | G17.5 − M3V |
  |---|---|---|---|
  | M3V (DIFF) | 91 / 146 / 199 | 202 / 283 / 359 | 8 / 43 / 77 |
  | M3V (AMP) | 44 / 84 / 133 | 121 / 180 / 248 | −10 / 24 / 57 |
  | pair A (DIFF) | −120 / −83 / −46 | | |
  | Mira (DIFF) | | −206 / −153 / −108 | |
  | G17.5 star (DIFF) | | | −65 / −32 / +0.6 (99.85%: +18) |
  | G17.5 star (AMP) | | | −54 / −19 / +12 (99.85%: +32) |
  | **real** | **132.8 (DIFF), 65.1 (AMP)** | **274.7, 164.6** | **35.3, 18.8** |

  The real data look exactly like an M3V source.
- **Posterior temperature.**
  - The posterior P_j ∝ exp(−χ²_j / 2κ) with κ = 1 is over-confident.
  - κ fitted by maximum log-score of the true source over all 93,000 injections: DIFF κ = 1.5 combined (per sector
    1.0–2.0); AMP κ = 2.5 (1.2–2.5).
  - Reliability at calibrated κ, for the 16 sources within 25″: predicted 0.97 → true 0.94 (DIFF) / 0.95 (AMP);
    predicted ≥ 0.99 → 0.998 / 0.977.
  - So extreme posteriors (1 − 10⁻⁵) should be read as "≥ 0.98–0.998".
- **Posteriors with equal priors over the 93 stars, calibrated κ.** Pair = A + B; others = all other catalogued
  stars.

  | | M3V | pair | Mira | G17.5 star (3.1″) | G17.7 star (5.7″) | others |
  |---|---|---|---|---|---|---|
  | DIFF S12 | 0.05 | 0.002 | 2e-6 | 0.17 | 0.11 | 0.66 |
  | DIFF S39 | 0.998 | 4e-7 | 9e-16 | 0.002 | 3e-6 | 7e-7 |
  | DIFF S66 | 0.21 | 0.30 | 0.015 | 0.22 | 0.07 | 0.19 |
  | DIFF S93 | 0.69 | 7e-4 | 3e-8 | 0.29 | 0.013 | 2e-4 |
  | DIFF S103 | 0.66 | 3e-5 | 1e-8 | 0.18 | 0.17 | 5e-5 |
  | DIFF S104 | 0.993 | 2e-6 | 3e-9 | 0.007 | 3e-4 | 2e-7 |
  | **DIFF combined** | **0.999992** | **6e-20** | **2e-40** | **8e-6** | **3e-12** | **4e-20** |
  | AMP S12 | 0.04 | 0.015 | 2e-5 | 0.06 | 0.03 | 0.86 |
  | AMP S39 | 0.99 | 3e-5 | 7e-12 | 0.010 | 6e-5 | 2e-5 |
  | AMP S66 | 0.32 | 0.22 | 0.006 | 0.23 | 0.09 | 0.13 |
  | AMP S93 | 0.26 | 0.13 | 0.012 | 0.21 | 0.13 | 0.25 |
  | AMP S103 | 0.13 | 0.11 | 0.005 | 0.20 | 0.08 | 0.47 |
  | AMP S104 | 0.67 | 0.03 | 0.003 | 0.17 | 0.06 | 0.07 |
  | **AMP combined** | **0.977** | **3e-6** | **5e-15** | **0.023** | **4e-5** | **1e-6** |

- **Independent second likelihood.**
  - Per-sector free-fit positions are compared with Gaussian-plus-outlier models of the injected positions for each
    star. The models are built on odd-numbered injections and tested on even-numbered ones, where the M3V is
    recovered as the most likely star (MAP) in 99% (DIFF) and 94% (AMP) of trials.
  - Results, combined: DIFF P(M3V) = 0.99997, G17.5 star 3×10⁻⁵, pair 1×10⁻¹⁷, Mira 5×10⁻²⁴. AMP P(M3V) =
    0.995, G17.5 star 0.005, pair 8×10⁻¹³, Mira 1×10⁻¹⁶.
- **Localisation region** (combined Δχ² map, M-dwarf colour; thresholds = quantiles of Δχ² at the true position for
  M3V injections).
  - Thresholds: DIFF 2.9 / 9.8 / 24.6 and AMP 3.8 / 14.9 / 39.7 (Gaussian 2.3 / 6.2 / 11.8; the noise is
    heavier-tailed).
  - Ring injections confirm these at other positions: 95% coverage 0.95–0.97, 99.7% coverage 0.994–0.998.

    | Method | 68% | 95% | 99.7% | Area of 99.7% region |
    |---|---|---|---|---|
    | DIFF | 2.2″ | 3.2″ | 5.1″ | 59 arcsec² |
    | AMP | 3.6″ | 6.1″ | 9.2″ | 174 arcsec² |

    Distances are the largest distance of each region from the M3V.
  - DIFF: free fit at (+0.6, +0.4)″ from the M3V; injections on the M3V give bias (−0.3, −0.1)″ and scatter 1.0″
    per axis.
  - For a normal-colour source (BP−RP 1.3) the DIFF region centres at (0, −1.5)″ and reaches 2.5 / 3.8 / 5.4″ from
    the M3V (`v02_colourmap.py`).
  - Only the M3V lies inside the DIFF 99.7% region. The AMP 99.7% region also contains the G17.5 star, by its own
    χ² (`figs/localisation.png`).

### 4. Free-amplitude multi-source fits (`p01_joint.py` → `out/joint_fits.json`)
- **Joint fits of the six sectors** (shared amplitudes, e⁻/s):
  - DIFF: M3V 14.7 ± 1.8, pair −2.0 ± 5.1, Mira +2.0 ± 3.8.
  - AMP: M3V 11.9 ± 2.3, pair +7.4 ± 6.5, Mira −4.6 ± 4.9.
  - The pair and Mira amplitudes are strongly anti-correlated (ρ = −0.95), so these errors are large. The M3V
    amplitude is still significant at 8σ (DIFF) and 5σ (AMP).
- **Non-negative fit with all 93 stars at once** (shared across sectors):
  - DIFF puts 13.7 e⁻/s on the M3V and < 0.7 on any other star.
  - AMP puts 9.8 on the M3V, 2.7 on the G17.5 star and ≤ 1.8 on the others.

### 5. Is the depth physically possible? (`p02_shape.py`, `p03_physics.py` → `out/shape_profile.json`, `out/physics_table.csv`)
- **Shape** (matched-filter light curve at the M3V, 21 transits, red-noise-inflated errors).
  - Depth 15.0 ± 0.8 e⁻/s, T14 = 1.48 ± 0.08 h, T23/T14 = 0.73 ± 0.09.
  - A V shape is worse by Δχ² = 31. T23/T14 ≥ 0.40 at Δχ² 6.8 (2.6σ) and ≥ 0.30 at Δχ² 12.7.
  - A dark occulter therefore needs k = Rp/R★ ≤ 0.43, i.e. a depth ≤ ~0.2.
- **Required depths** (A = 14.3 e⁻/s; stellar parameters from Gaia and the Mamajek table):

  | Star | Flux (e⁻/s) | Required depth | Companion radius | Possible? |
  |---|---|---|---|---|
  | M3V (R 0.38, M 0.37) | 837 | 1.7% | 5.4 R⊕ | Yes. Duration fits at b ≈ 0.6–0.7 |
  | pair A (M_G 10.7, R ≈ 0.31) | 611 | 2.3% | 5.2 R⊕ | Yes |
  | pair B (M_G 11.8, R ≈ 0.23) | 246 | 5.8% | 6.1 R⊕ | Yes, but only near b ≈ 0 |
  | G17.5 star at 3.1″ (K dwarf, ~1.5 kpc) | 28 | 51% | k = 0.72 | No: needs a total eclipse of a luminous compact companion |
  | G17.7 star at 5.7″ | 20 | 71% | | No |
  | G17.0 star 7″ from the Mira | 40 | 36% | | No |
  | Mira | 6,000–14,000 (measured) | 0.1–0.24% | | No: a 5.69-d orbit (a ≈ 15 R☉) is inside a star of hundreds of R☉ |

  - Of the 93 stars, 4 would need depth ≥ 100%.
  - Only the M3V, pair A and pair B pass the dark-occultor test (flat-bottomed, T14 1.32–1.64 h, circular orbit).
  - For the bright field stars (23–95″ away; depth 0.3–30%) the test fails on duration. For K/G dwarfs and evolved
    stars, a non-grazing transit at 5.69 d lasts ≥ 2.3 h unless the orbit is eccentric (e ≳ 0.4). So for them the
    physics argument is model-dependent, but the pixels exclude them by Δχ² ≥ 267.
  - Pixels rule out the pair; physics, and less strongly pixels, rule out the faint stars near the M3V.

### 6. Robustness and model mismatch (`v01_variants.py`, `v03_faintcolour.py`, `p04_mira_emp.py`, `r05_mismatch.py`)
Combined Δχ² relative to the M3V across variants (DIFF / AMP):

| Variant | pair A | Mira | G17.5 star |
|---|---|---|---|
| Baseline | 133 / 65 | 275 / 165 | 35 / 19 |
| Fit radius 3.5 / 4.5 / 7 px | 133–142 / 64–66 | 278–283 / 153–167 | 36–39 / 19 |
| Colour curve ×0.7 | 136 / 63 | 274 / 165 | 38 / 21 |
| Colour curve ×1.3 | 129 / 67 | 272 / 162 | 29 / 14 |
| WCS ±0.03 px (4 directions) | 124–141 / 60–71 | 266–284 / 158–172 | 28–42 / 14.5–23 |
| Mira at Gaia position / Gaia + colour curve | | 284, 277 / 149, 171 | |
| Mira with its empirical image as PSF | | 291 / 172 | |
| OOT window 0.15 d; 0.40 d (DIFF only) | 133; 99 | 284; 205 | 40; 29 |
| Spline knots 0.3 d; 1.0 d (AMP only) | 86; 57 | 188; 127 | 16; 11 |
| Drop any one sector | 88–132 / 33–63 | 171–264 / 88–152 | 17–38 / 5–21 |
| G17.5 star assumed BP−RP 2.0 | | | 20 / 12 |
| G17.5 star assumed BP−RP 2.76 (as red as the M3V) | | | 8.3 / 7.1 |
| No background subtraction (DIFF) | 21 | 55 | −0.5 |

- The pair and the Mira are excluded in every variant.
- The G17.5-star exclusion weakens in three cases:
  - without S39: Δχ² drops to 17 (DIFF) and 5 (AMP);
  - if the star were as red as the M3V: Δχ² ~ 7–8;
  - without background subtraction, which the fake-window noise shows to be invalid.

**PSF-mismatch injections** (`r05_mismatch.py` → `out/mismatch_{diff,amp}.json`). Each was 1,000× per source and
injected with a PSF different from the fitted PRF:
- 'alt': the PRF of a CCD location 450 px away;
- 'blur': the PRF smeared by 0.15 px rms;
- 'emp': the PRF plus the residual pattern of the red field star Gaia DR3 5925373341570270208 (BP−RP 2.41, T 10.5);
- 'miraemp': the Mira's own empirical image (Mira only);
- 'ownamp': each hypothesis's own best-fit amplitude instead of the common 14.3 e⁻/s.

Results:
- **Pair and Mira.** For true pair A, pair B or Mira sources, the 99.85% quantile of χ²(true) − χ²(M3V) across all
  variants is ≤ −21 (pair, DIFF), ≤ +16 (pair, AMP), ≤ −59 (Mira, DIFF) and ≤ −3.8 (Mira, AMP). The real values
  (+133/+150/+275 DIFF, +65/+73/+165 AMP) are reached in 0 of 1,000 trials in every variant.
- **The empirical Mira PSF.** Injected on the Mira, the Mira is still recovered as the best star in 100% (DIFF) and
  89% (AMP) of trials.
- **The M3V.** For a true M3V source the M3V stays the best catalogued star in 99% (DIFF) and 91–93% (AMP) of trials.
  The combined free-fit bias changes by ≤ 0.4″.
- **The G17.5 star.** The real Δχ² is reached in 0–0.1% (DIFF) and 1.0–1.4% (AMP) of injections on that star.

### 7. Uncatalogued source (item 5) (`r04_ring.py` → `out/ring_{diff,amp}.json`)
- **No sign of one.**
  - Free-fit χ² minus that of the best catalogued star is 0.28 (DIFF) and 1.17 (AMP).
  - For a source truly on the M3V, P(≥) = 0.78 and 0.63.
- **How close it would have to be.** Injections of the observed signal at r = 1–8″ from the M3V (8 position angles,
  M-dwarf colour, 600 each). Fraction of trials in which the M3V is excluded at 95.4% / 99.73%:

  | r | DIFF | AMP |
  |---|---|---|
  | 1″ | 8% / 0.4% | 5% / 0.2% |
  | 2″ | 29% / 3% | 11% / 0.5% |
  | 3″ | 66% / 13% | 25% / 2% |
  | 4″ | 90% / 42% | 46% / 5% |
  | 6″ | 99.8% / 92% | 84% / 28% |
  | 8″ | 100% / 99.6% | 97% / 67% |

  - 50% exclusion at 99.73% happens at 4.3″ (DIFF) and 7.1″ (AMP); 90% at 5.9″ (DIFF).
  - **A source within ~2–3″ of the M3V cannot be told apart from it with TESS.**
- **Gaia limits.**
  - Such a source would need T ≲ 15.8 (ΔT ≲ 2.8 from the M3V) to give a flat-bottomed transit with depth ≤ 0.2.
    For an exotic total eclipse it would need T ≤ 17.5.
  - Within this field Gaia DR3 resolves a ΔG = 3.2 star at 3.1″ (the G17.5 star).
  - The M3V shows no sign of a companion: RUWE 1.08, ipd_frac_multi_peak 0, ipd_gof_harmonic_amplitude 0.013.
  - So the remaining room is an unresolved companion at ≲ 1″, or a fainter star (ΔG ≳ 4) at 1–2″.
  - The ≲ 1″ / 1–2″ limits are approximate. They come from literature-level Gaia completeness, not a calibration
    done here.

## Caveats
- **Narrow separations.** The M3V–pair separation is 11.7″ = 0.55 px and pair–Mira is 6″. Everything rests on
  sub-pixel PRF modelling next to a star 7–17× brighter whose PSF the model matches only to ~5%. The injections
  include the real systematics (fake windows) and position/colour uncertainties, and Section 6 tests PSF mismatch.
- **The G17.5 star at 3.1″ is excluded only at ~3–4σ by pixels.**
  - It is disfavoured 0/1000 (DIFF) and 10/1000 (AMP).
  - The exclusion relies on Sector 39 and on its catalogue colour, which may be biased by the nearby M3V (BP/RP
    excess factor 1.49).
  - Its stronger exclusion is the 51% depth required.
- **Calibration limits.**
  - Noise draws reuse 259–350 overlapping fake windows per sector, so tail probabilities beyond ~10⁻³ are
    extrapolations.
  - The κ-calibrated posteriors are over-confident at the 10⁻³ level, so read "1 − 10⁻⁵" as "≥ 0.998".
- **Colour inputs.** Pair B's colour is assumed (BP−RP 3.3). The lateral-colour curve beyond BP−RP 2.5 rests on few
  calibrators (±0.5″ used).
- **Which methods are independent.** DIFF and AMP use the same pixels and are not independent of each other. They
  are independent of the other agents' deblended photometry and of the earlier `loc5.py` / referee code; only the
  PRF package and the colour curve are shared.
- **Reliance on background subtraction.** All results assume the background-plane subtraction. Without it the
  discrimination collapses (pair +21, Mira +55, G17.5 star −0.5). The fake-window tests show that subtracting the
  background is the correct treatment.

## Files

### Scripts (`scripts/`)
| Script | What it does | Output |
|---|---|---|
| `a00_gaia.py` | Gaia DR3 (+astrophysical parameters) within 330″; TIC rows | `data/gaia_330.csv`, `data/tic_trio.csv` |
| `l0lib.py` | Library: sector loading, barycentric fix, Gaia → apparent pixel positions (PM, WCS offset, lateral-colour curve with field-angle scaling, boresight from tess-point), PRF, templates | |
| `l1img.py` | Background plane; DIFF and AMP per-window products; real and fake window sets | |
| `l2fit.py` | `SectorModel`: products for real and fake windows, noise inflation, 93 hypotheses, sky grid, closed-form fits, fake draws | |
| `e01_explore.py` | Geometry per sector (camera, field angle, boresight PA, Jacobian), star positions, Mira flux range | stdout |
| `c01_wcs.py` | WCS calibration on normal-colour stars; red-calibrator check | `out/wcs_cal.json` |
| `c02_mira.py` | Mira apparent position, 18 variants, 2-d chunks | `out/mira_pos.json`, `figs/mira_direct_s*.png` |
| `c03_prfcheck.py` | PRF residuals of isolated field stars | `out/prf_match_fieldstars.json` |
| `e02_noise.py` | Noise ratios, pixel correlations, crude totals | stdout |
| `e03_pca.py` | PCA of fake difference images WITHOUT background subtraction (shows why it is needed) | `figs/noise_pca.png` |
| `e04_quick.py` | Quick per-sector look at both methods | stdout |
| `e05_perwindow.py` | Single-transit fits | `out/per_window.csv` |
| `r01_sim.py` | Injection–recovery engine (modes prf/alt/blur/emp/miraemp; targets all/main/mira/ring; `AMPMODE=own`) | `out/sim_*.npz` (large) |
| `run_sims.sh`, `run_sims2.sh` | Simulation batches | |
| `r02_post.py` | κ calibration, reliability, posteriors (model-based and empirical-position), p-values, localisation thresholds and regions, uncatalogued statistic | `out/post_*.json`, `figs/loc_map_*.png` |
| `r03_cubecheck.py` | Cube-level injection check of linearity | `out/cubecheck.json` |
| `r04_ring.py` | Distinguishability of an uncatalogued source vs offset | `out/ring_*.json` |
| `r05_mismatch.py` | PSF-mismatch and own-amplitude injection comparison | `out/mismatch_*.json` |
| `p01_joint.py` | Free-amplitude multi-source fits (3/4/7 sources; non-negative fit with all 93); Mira flux scaling | `out/joint_fits.json` |
| `p02_shape.py` | Matched-filter light curve at the M3V; shape profiles | `out/shape_profile.json`, `figs/shape_fold.png` |
| `p03_physics.py` | Required depth and feasibility for all 93 stars | `out/physics.json`, `out/physics_table.csv` |
| `p04_mira_emp.py` | Mira hypothesis with its empirical image as PSF | `out/mira_emp_fit.json` |
| `v01_variants.py` | Analysis-choice variants | `out/variants.json`, `out/log_variants.txt` |
| `v02_colourmap.py` | Localisation maps for source colour 2.76 and 1.3 | `out/colourmap.json` |
| `v03_faintcolour.py` | Faint-neighbour discrimination vs its assumed colour | `out/faint_colour.json` |
| `f01_figs.py` | Per-sector difference images and residuals | `figs/diff_images.png` |
| `f02_locmap.py` | Localisation maps and per-sector free fits vs injection clouds | `figs/localisation.png` |
| `s01_summary.py` | Collects key numbers | `out/summary.json` |

### Main outputs
- `out/summary.json` holds all key numbers.
- `out/post_diff_prf_all.json` and `out/post_amp_prf_all.json` hold the posteriors, p-values and localisation.
- `out/sim_*.npz` are the injection arrays: `sim_{diff,amp}_prf_all.npz` (baseline, all 93 stars, ~210 MB each),
  `sim_{diff,amp}_ring.npz` (uncatalogued-source rings), `sim_{diff,amp}_{alt,blur,emp}_main.npz`,
  `sim_{diff,amp}_miraemp.npz` and `sim_diff_prf_main_ownamp.npz`. The folder is ~0.8 GB in total; the npz files can
  be deleted and regenerated with `run_sims.sh` / `run_sims2.sh` (about 2 h on 2 threads).
  - Per realisation × sector × hypothesis: χ²; per realisation × sector: free-fit position.
  - Keys: `inj_<star>__chiH`, `__free`, `__cmin`, `__gcomb_pos`, `__dchi_m3v`; `real__chiH`, `real__gsum`,
    `hyp_name`, …
- `figs/`: `diff_images.png`, `localisation.png`, `shape_fold.png`, `mira_direct_s*.png`, `noise_pca.png`.
