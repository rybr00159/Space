# Referee check of the 'loc' attribution (TIC 214311908): which star dims?

Folder: `tic214311908/analysis/referee_loc/`
(AI-run referee analysis, 2026-10-17. Nothing was posted, uploaded or submitted.)

I re-derived the key numbers with my own code. I share only the inputs with the investigator ('loc'): the TESScut
cubes, the Gaia DR3 catalogue (queried again), the 'PRF' package and its SPOC model files, and the TIC 71316387
lateral-colour curve. I read the investigator's scripts to look for errors but did not run or import them, except
`rq01_time.py`, which imports `l0lib` only to compare time stamps.

## Method: differences from the investigator

| Step | Investigator | Referee (this folder) |
|---|---|---|
| Background | Per-cadence plane fitted to the faintest 30% of pixels and subtracted from the cube | None in the time domain. A plane (3 terms) is fitted in every image fit; the noise is estimated after removing that plane |
| Transit images | DIFF: straight lines through out-of-transit (OOT) windows. AMP: global B-spline | 'LQ': per pixel, a quadratic baseline fitted at the same time as the trapezoid template, within ±0.30 d (variants ±0.20 and ±0.45 d) |
| Fit region | 5.5 px around the Mira | 5 px around the M3V (variants 3.5 and 7 px) |
| Noise | Formal errors × per-pixel inflation | Formal errors × the median of \|plane-removed fake residual\| / σ (variant: per-pixel inflation) |
| WCS / Mira | `c01_wcs.py`, `c02_mira.py` | Own global fit of 470–520 stars; Mira fitted with a free position, free blur and the trio flux free or fixed; plus a model-light centroid |
| Colour-error propagation | Independent radial errors per sector | The same, plus a scale error of the colour curve common to all sectors, N(1, 0.2) (variant: 0) |

## Scripts (`scripts/`)

| Script | Purpose | Output |
|---|---|---|
| `rl_lib.py` | Library: cube loading; target-direction barycentric time; Gaia → apparent pixel position (PM, header WCS, offset, colour shift along the tess-point boresight direction); PRF wrapper; trapezoid | |
| `rq00_gaia.py` | Gaia DR3 cone, 400″ | `data/gaia_400.csv` |
| `rq01_time.py` | My TESScut time correction compared with the raw TESScut times and with the investigator's | `out/time_check.json` |
| `rq02_wcs_mira.py` | WCS offset from a global PRF fit to the median image (7 variants); the Mira's apparent position (8 fit variants plus a centroid) | `out/wcs_mira.json` (merged from `wcs_mira_12_39_66.json` and `wcs_mira_93_103_104.json`) |
| `rq03_img.py` | LQ per-window depth images for the real windows and fake windows (0.05-d grid, no overlap with any real transit) | `out/img_W0.30.npz`, `img_W0.20.npz`, `img_W0.45.npz` |
| `rq04_fit.py` | Fits of the hypotheses (93 Gaia stars within 80″ with T < 17.6; the Mira at its measured position), a free-position grid co-moving with the M3V, injection–recovery on 7 sources, and ring injections | `out/fit_<tag>.json`, `out/inj_<tag>.npz`, `out/ring_<tag>.npz` |
| `run_variants.sh` | 23 analysis variants of `rq04_fit.py` | `out/fit_v_*.json`, `out/log_fit_v_*.txt` |
| `rq05_tic.py` | TIC-8 rows near the Mira | `data/tic_near.csv` |
| `rq06_shape.py` | PRF-weighted light curve at the M3V; joint trapezoid fit; profiles in T23/T14 and T14; jackknife over transits with V-shape profiles | `out/shape.json`, `out/shape_lc.npz`, `figs/shape.png` |
| `rq07_perwindow.py` | Single-transit fits; the Mira hypothesis with the Mira's empirical image as its PSF, plus 1,000 injections | `out/perwindow.csv`, `out/mira_emp.json` |
| `rq08_summary.py` | Collects all runs into one table | `out/summary.json` |
| `rq09_locregion.py` | Localisation regions from the combined Δχ² map, with injection-calibrated thresholds | `out/locregion_base.json`, `out/locregion_base_cssd0.json` |
| `rq10_persector.py` | For each sector, the real free-fit position and hypothesis Δχ² compared with the injection clouds | `out/persector.json` |
| `rq11_joint.py` | Joint fits with free amplitudes (3 and 4 sources) | `out/joint.json` |
| `rq12_mirascale.py` | Mira flux-scaling test: amplitude per sector under the Mira and M3V hypotheses against the Mira's fitted TESS flux | `out/mirascale.json` |
| `dbg_time.py` | Timing test of one global WCS fit (development only) | |

Logs: `out/log_*.txt`. `out/log_shape_partial_slow.txt` is an interrupted first jackknife run with a finer grid, kept
for reference; its seven subsamples agree with the final run.

## Main results (baseline run `base`: W = 0.30 d, campaign colour curve, referee WCS and Mira positions)

### Hypothesis tests
Combined χ²(h) − χ²(M3V):

| Hypothesis | Δχ² | Investigator (DIFF / AMP) |
|---|---|---|
| pair A | +126.2 | +132.8 / +65.1 |
| pair B | +145.5 | +150.1 / +72.9 |
| Mira | +284.8 | +274.7 / +164.6 |
| G17.5 star | +35.8 | +35.3 / +18.8 |
| G17.7 star | +86.6 | +79.6 / +51.0 |
| G17.0 star | +189.3 | |

- **Injections on the true source** (1,000 each) never reach the real value for pair A, pair B, the Mira or the
  G17.5 star (0/1000).
  - Under a pair A source, χ²(pairA) − χ²(M3V) has 2.5/50/97.5/99.85% quantiles of −201/−114/−58/−35.
  - Under a G17.5 source: −83/−44/−7/+6.
- **Injections on the M3V reproduce the real values.** The fractions of injections at or below the real value are
  0.12 (pair A), 0.14 (Mira) and 0.31 (G17.5).
- **Across 23 variants** (`out/summary.json`):
  - pair A: +75 to +138; Mira: +183 to +301. In every variant, 0 of the injections on the true source reach the
    real value.
  - G17.5 star: +18 to +42, with p_true ≤ 0.007, except the implausible case where it has BP−RP 2.76 (+5.5,
    p = 0.04).
  - Gaussian log10 likelihood ratio L(M3V)/L(G17.5): 3.9 at baseline, ranging over 1.4–4.2 across the variants.
- **No background subtraction is needed.** My pipeline does no time-domain background subtraction and recovers the
  full discrimination. The investigator's "no background subtraction" variant (pair +21) came from leaving the
  uniform offset in the per-pixel noise, not from the data.
- **Mira's empirical image as its PSF** gives +310.7 combined. 1,000 injections with that PSF reach at most −128
  (99.85%) (`out/mira_emp.json`).
- **Joint free-amplitude fit** (`out/joint.json`): M3V 14.8 ± 1.7, pair −3.0 ± 4.8, Mira +4.7 ± 3.6 e⁻/s. With the
  G17.5 star added: G17.5 1.3 ± 2.5, M3V 14.0 ± 2.3.
- **Per sector** (`out/persector.json`):
  - S12 (pair A −1.0) and S66 (+1.1) are individually ambiguous. Each sits at the 5–7% tail of the M3V injections.
  - S39 (+37) and S104 (+51) dominate.
  - Fisher combination of the per-sector p-values under the M3V hypothesis gives p ≈ 0.09, which is consistent.
- **Goodness of fit.** χ²(M3V) totals 465.8 for 445 degrees of freedom. 79% of M3V injections give a larger value.
- **Single transits** (`out/perwindow.csv`):
  - The pair fits worse than the M3V in 15 of 20 transits; the Mira fits worse in 20 of 20.
  - The amplitudes are consistent with a constant (χ² ≈ 20 for 19 degrees of freedom).
  - S103 epoch 56 is low at 8.4 ± 2.8 e⁻/s.

### Position
- **Free fit** (M-dwarf colour): (+1.6, +0.3)″ from the M3V.
  - M3V injections give a mean of (+0.2, 0.0)″ and a scatter of (1.6, 1.15)″, or (1.3, 1.0)″ without the common
    colour-scale error.
  - The free fit is better than the M3V hypothesis by only Δχ² = 3.1, for which P(≥) ≈ 0.45. There is no evidence
    for an uncatalogued source.
- **Localisation thresholds** (68/95/99.7%): 4.6 / 13.8 / 56.6 ± 26 (bootstrap) at baseline, and 4.5 / 13.5 /
  32.1 without the common colour-scale error.
- **Region sizes** (largest distance from the M3V): 3.2 / 4.6 / 8.1″ at baseline, or 3.2 / 4.6 / 6.3″ without the
  common colour-scale error.
  - The investigator's DIFF regions are 2.2 / 3.2 / 5.1″, which is too small.
  - In the M-dwarf-colour frame, the 95% region contains the sky position of the G17.5 star (map value 7.4). That star
    is excluded only by its own-colour hypothesis test (above) and by physics.
- **Ring injections** (uncatalogued source at r from the M3V; fraction of trials in which the M3V is excluded at
  95.45% / 99.73%):

  | r | Baseline | No common colour error |
  |---|---|---|
  | 2″ | 0.25 / 0.00 | 0.26 / 0.03 |
  | 3″ | 0.60 / 0.01 | 0.61 / 0.11 |
  | 4″ | 0.86 / 0.05 | 0.86 / 0.37 |
  | 6″ | 0.99 / 0.52 | 0.99 / 0.90 |

  - So 50% exclusion at 99.73% needs about 4.2″ (no common colour error) to about 6″ (with it).
- **Colour of the source.**
  - With the source assumed to have BP−RP 1.3, the best free fit is χ² +6.9 worse than the M3V hypothesis.
  - It is about 10 worse than the free fit for an M-dwarf colour.
  - So the sector-to-sector pattern of positions mildly prefers a red source, given the colour curve.

### Calibration inputs
- **WCS offsets** (px, from my global fit): S12 (−0.087, +0.068), S39 (−0.073, +0.108), S66 (+0.025, +0.133),
  S93 (+0.065, +0.185), S103 (+0.027, +0.104), S104 (−0.121, +0.140). The investigator's agree to ≤ 0.013 px.
- **Mira apparent offset from Gaia**:
  - PRF fit: radial −5.0 to −7.5″, tangential ≤ 0.6″. The 8 fit variants scatter by ≤ 1.2″.
  - Centroid lower bound: −3.4 to −6.2″.
  - The investigator's "5–7″, not 16–17″" is confirmed.
- **Time**: my target-direction correction agrees with the investigator's to < 0.04 s (−21 to +33 s relative to the
  raw TESScut times).

### Shape and depth
`out/shape.json`, `out/shape_lc.npz`, `figs/shape.png`. 19 transits.
- Best fit: depth 14.7 e⁻/s, T14 1.51 h, T23/T14 0.72, time offset +1.2 min from the ephemeris.
- T14 and the flat fraction are strongly degenerate:
  - The profile is flat between T14 = 1.3 h (flat fraction ≈ 0.93) and 1.5 h (Δχ² < 1).
  - T14 = 1.6 h gives Δχ² +1.6, and 1.7 h gives +6.5.
  - The jackknife scatter of the best flat fraction is 0.35.
- Profile in T23/T14 (Δχ² against the best fit):

  | T23/T14 | Δχ² |
  |---|---|
  | 0 (V shape) | 37.4 |
  | 0.15 | 29.3 |
  | 0.30 | 17.3 |
  | 0.45 | 8.6 |
  | 0.6 | 1.9 |

- In all 19 leave-one-transit-out subsamples, Δχ² ≥ 22.6 at 0.15 and ≥ 11.9 at 0.30.
- So a dark body with k ≥ 0.6 (maximum T23/T14 ≤ 0.25) is excluded at about 4σ or more.

### Mira flux scaling
`out/mirascale.json`. Amplitudes under the Mira hypothesis, with errors from the formal variance × k:
- Constant: χ² 6.2. Proportional to the Mira's TESS flux (6,100–14,000 e⁻/s): χ² 16.4, for 5 degrees of freedom
  (p ≈ 0.006).
- Proportional scaling is mildly disfavoured but not excluded.
- Under the M3V hypothesis the amplitudes are constant: χ² 5.0 for 5 degrees of freedom.

### Required depth
For A = 15.2 ± 0.6 e⁻/s (combined M3V amplitude):

| Star | Required depth |
|---|---|
| M3V | 1.62% (TIC T = 13.01, 937 e⁻/s) to 1.82% (Gaia relation T = 13.13, 837 e⁻/s) |
| pair A | 2.3–2.5% |
| pair B | 6% (assumed BP−RP 3.3) to 13% (TIC T = 15.27 ± 0.6) |
| G17.5 star | 54% |
| G17.7 star | 75% |
| G17.0 star | 38% |

- The G17.5 star has parallax 0.65 ± 0.11 mas, so at G = 17.5 it is a K-type dwarf. If it were as red as the M3V,
  its parallax would have to be ≥ 4 mas. Contamination by the M3V would make its catalogue colour too red, not too
  blue.

## Files that matter downstream
- `out/summary.json`: all runs, with p-values and likelihood ratios.
- `out/fit_base.json` and `out/inj_base.npz`: the baseline run.
- `out/wcs_mira.json`: WCS offsets and Mira positions.
- `out/locregion_base*.json`: localisation regions.
- `out/shape.json`: transit shape.
- `out/joint.json`, `out/mira_emp.json`, `out/perwindow.csv`, `out/persector.json`.
- The whole folder is about 45 MB.
