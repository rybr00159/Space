# Review log: TIC 214311908 write-up

*2026-10-17. All reviewers were AI agents. Nobody outside this analysis has seen the write-up; nothing was submitted.
Raw output is in [`review_findings.json`](review_findings.json); scratch paths there are shortened to `<scratch>` and the
name of the AI tool is redacted.*

## What was run

| Step | Status |
|---|---|
| Systematics panel: own raw-pixel image fits, shape fits and blind box searches from the TESScut cubes ([`analysis/panel_systematics/`](analysis/panel_systematics/)) | Done; 13 findings (2 major) |
| FPP panel: sign-fixed re-run of the prior variants and 1,500 joint prior draws, follow-up ladder with the adopted inputs, parameter re-derivation ([`analysis/panel_fpp/`](analysis/panel_fpp/)) | Done; 12 findings (3 major) |
| Audit panel: numbers, links, prediction times, manifest, figures ([`analysis/panel_audit/`](analysis/panel_audit/)) | Done; 29 findings (4 major) |
| Applying the review | Done by hand. Each finding was checked against the source files before it was applied. For the FPP panel, 3 named variants and 3 joint prior draws were re-run with its code and reproduced to machine precision, and the draw statistics were recomputed from `p05_draws.csv`. The follow-up ladder was recomputed independently ([`scripts/w05_fpp_ladder.py`](scripts/w05_fpp_ladder.py)); with the referee's own inputs it reproduces the published ladder exactly |

The three panels agree that the central conclusion stands: a strong planet candidate, not statistically validated, which
only a ground-based transit can test against the Mira and the photometric systematics. No finding changed a headline
number of the planet (depth, radius, period, ephemeris) or the FPP at the adopted priors.

## Findings and what was done

S = systematics panel, FPP/PAR/BOOK = FPP panel, F = audit panel. "Fixed" means the change is in this version of the
README, captions, figures, `key_numbers.json`, `corrections_to_repo.md` and the manifest.

| ID | Severity | Issue | Action |
|---|---|---|---|
| SYS-1 | major | Deblending pipelines and the raw 3×3 box were credited with putting the dip on the M3 dwarf; the box also holds the Mira's core | Fixed: they now "confine the dip to the pixels holding the M3 dwarf and the Mira"; only the image fits separate the M3 dwarf from the Mira (summary, localisation caption, captions.md, corrections lines 41 and 61). The panel's own template fit is quoted as support |
| SYS-2 | major | Late-minus-early contrast quoted as +3.3 to +4.1 ppt, omitting REF with E58 (+6.1 ± 2.1, 3.0σ) and mixing analyst and referee values | Fixed: caption gives the referee's INV +4.1 ± 2.6, REF +6.1 ± 2.1 (χ² 15.0/5) and +3.4 ± 2.1 without E58; key_numbers from `q14_contrast_all.json` |
| SYS-3 | minor | "The trend comes from modelling the scene" overstated | Fixed: "most likely", paired difference ~2σ, omnibus p 0.06–0.45, apertures cannot exclude 3–4 ppt. The aperture errors are quoted as ±2.4–3.2 ppt (the referee's apertures), not ±2.4–3.8 |
| SYS-4 | minor | "No field-wide regressors" is false for the raw tests | Fixed: field-median background and, in the aperture fits, a field-star pointing term are named |
| SYS-5 | minor | T14 ± 0.07 h and T23/T14 depend on the density prior; light-curve-only fits differ | Fixed: duration row gives the free-density values (1.65 ± 0.10, 1.75 ± 0.14 h) and the raw-aperture preference (1.3–1.5 h, 1.5–2σ); same qualifier in the corrections file |
| SYS-6 | minor | 2P exclusion "≥ 7.7σ" is a phase-2 number | Fixed with F11: "≥ 6σ in every cleaned photometry (6.0σ, 6.7σ, 7.2σ)". The panel's own differently detrended raw aperture gives 3.8σ; this is quoted too |
| SYS-7 | minor | Raw apertures confirm the dip but are not a blind detection | Fixed: new caveat bullet and a pointer in the first table row |
| SYS-8 | minor | 3.1″ star "~3–4σ" too strong | Fixed with F14: ~2.5–3.5σ, with p-values, and "lies inside the 95% region" |
| SYS-9 | minor | Corrections file overstated the Mira flux-scaling test | Fixed: "disfavours a Mira host at 2.7–3.4σ there (p ≈ 0.006 in the referee's Mira-hypothesis fits)"; the Mira is excluded mainly by the pixel image and the orbit |
| SYS-10 | nit | TTV χ² from the analyst, not the referee | Fixed with F06 |
| SYS-11 | nit | 13.9 e⁻/s is not corrected for the 3×3 box's encircled energy | Partly: the README now says the box holds 76–81% of the light and that light-corrected raw apertures give 16.7–17.5 ppt. The suggested "≈ 16.7 ppt after this correction" was not used: 13.9 e⁻/s ÷ (979 × 0.78) is ≈ 18 ppt, and 16.7 ppt comes from a different (fixed-shape) fit |
| SYS-12 | nit | Source of the 979 e⁻/s; Mira ratio 7–17× | Fixed: Gaia XP calibration named (TIC T in brackets); 6–14× |
| SYS-13 | nit | "Transit unchanged" after the RR Lyrae removal | Fixed: global depth, shape and timing unchanged; per-transit depths moved by 2.7 ppt rms; 0.8-ppt residual in S93 |
| FPP-1 | major | Follow-up ladder used f_BD = 10⁻³ and the uncalibrated shape grid after "now" | Fixed: [`scripts/w05_fpp_ladder.py`](scripts/w05_fpp_ladder.py) keeps the referee's rules and factors and swaps in the adopted composition, over three brown-dwarf calibrations. Speckle + two-band: 3.1–6.1 × 10⁻⁵ (REF) / 3.6–11.8 × 10⁻⁵ (INV); + brown-dwarf channel closed: 3.0–4.9 × 10⁻⁵ / 3.0–8.3 × 10⁻⁵; closing the channel gains at most ~2× (the panel's 1.7× is for one calibration only). Figure panel b redrawn with the published ladder as open markers; captions, key_numbers and the corrections file updated |
| FPP-2 | major | The robustness range ignored prior dimensions that were not re-run after the sign fix | Fixed: joint draws (median 0.0022, 5–95% 0.0002–0.012, P(< 0.015) 0.97) and pessimistic variants (0.015–0.045) added to the FPP row, total row, caveats and "Is it validated?"; companion-limit sensitivity stated. Verified by re-running variants and draws (see above). The panel folders are added to the bundle as the source |
| FPP-3 | major | Comparison with the Giacalone et al. 2021 threshold not like for like | Fixed: "Is it validated?" now gives the strict FPP (0.005 / 0.009; draws median 0.013, below 0.015 in 56%) and adds prior dependence and a reconnaissance spectrum. The paper's sections on cadence were not re-read for this step, so that part of the panel's argument is not repeated |
| FPP-4 | minor | NFPP definitions mixed; "most lenient" not true; variants reach ~1.1 × 10⁻³ | Fixed: NFPP = pair-A term only (3.7 / 4.2 × 10⁻⁴), "most lenient Student-t likelihood", 0.9–1.1 × 10⁻³ cases, typical 10⁻⁹–10⁻⁷; background counted in the FPP; corrections file "~4e-4" |
| FPP-5 | minor | "0.9995 on some star of the bound system" leaves out pair A | Fixed: "on the M3 dwarf or an unresolved companion"; 0.9999 with pair A (README, captions, figure) |
| FPP-6 | minor | The 15–30× circular/eccentric contrast holds for REF only | Fixed with F10: ~30× (REF), ~2× (INV) |
| FPP-7 | minor | The referee recommends the Sagear & Ballard prior | Partly: the table keeps Kipping as its reference row (the ladder factors and the figure are computed for it) but states the recommendation, gives the Sagear & Ballard planet-on-companion values and quotes the total as 0.0013–0.0016 / 0.0031–0.0033. The headline 0.0013–0.0033 already spans both |
| FPP-8 | minor | White-dwarf posterior leaves out the Gaia RV term | Fixed: posterior with and without the RV term (~10⁻⁷ with it) and a note that the shown value is conservative |
| FPP-9 | minor | History of the superseded FPP values inaccurate | Fixed with F13 |
| PAR-1 | nit | Mass, K, TSM, ESM computed for Rp = 4.90 R⊕ | Fixed by a note with scaled values (M ≈ 21.7 M⊕, K ≈ 15.1 m/s, TSM ≈ 110, ESM ≈ 9.3). `r06` was not re-run (it needs a network query); ESM 9.3 is r06's 9.1 scaled by (4.95/4.90)², not the panel's 9.4 |
| PAR-2 | nit | "Neptune-sized" for 1.3 Neptune radii | Fixed: title "super-Neptune", summary "a little larger than Neptune (1.3 Neptune radii)"; the corrections row for the repository README now proposes the same wording |
| BOOK-1 | nit | Rounding (0.0037, 1.15 × 10⁻³), Gaia RV factor in the BD row, twin definition, "PTP" acronym | Fixed (a)–(e); the figure labels the companion channel "STP" with an explanation in the caption |
| F01 | major | X-ray activity still from eRASS1; eRASS:3 makes the attribution a coin flip and halves the flux | Fixed: activity row, rotation advice (12–17 d; 8.5–12.8 d from eRASS1) and caveat |
| F02 | major | Novelty row overclaimed ("not catalogued", "no paper"); novelty/ not linked or bundled | Fixed: row "Not previously reported as a candidate" with the full list of sources, catalogue mentions and the ADS limitation; caveat rewritten; novelty/ linked and in the manifest; corrections line 117 |
| F03 | major | ~420 bundled files hard-code the session scratch path, which names the AI tool | Fixed for the copy step: [`scripts/w06_copy_bundle.py`](scripts/w06_copy_bundle.py) copies the manifest, rewrites the working paths (analysis root → `tic214311908/analysis`) and fails if any absolute path is left; the manifest header says to use it. Tested on a scratch copy (no path left). The analysis folders themselves were not edited. README says paths may need adjusting |
| F04 | major | Corrections file both kept and "corrected" the old prediction file | Fixed: the vetting replacement links the new file and says the old one is kept for the record, uncorrected; "do not overwrite". When the corrections were applied, the old file was instead overwritten with the corrected one (same format and epochs; the uncorrected version is in the git history), so there is only one version to use |
| F05 | minor | Per-transit caption omits REF +6.1 ppt | Fixed with SYS-2. The look-elsewhere p is quoted as 0.31–0.81 (it includes REF without E58, 0.52–0.81), not 0.31–0.64 |
| F06 | minor | Timing χ² not the referee's | Fixed: 8.9/17 and 12.6/19 (analyst's in brackets); key_numbers source `q05_timing.json` |
| F07 | minor | ±10-min margins < 3σ in late 2027 and 2028 | Fixed: ±3σ margins, ±10 min (2027), ±11–14 min (2028), ~1 min pipeline systematic |
| F08 | minor | Summary's post-follow-up FPP used f_BD = 10⁻³ | Fixed with FPP-1 |
| F09 | minor | NFPP values used different definitions | Fixed with FPP-4 |
| F10 | minor | Circular/eccentric contrast REF only | Fixed |
| F11 | minor | "≥ 7.7σ" is a phase-2 number | Fixed (README and FINDINGS entry) |
| F12 | minor | "In two separate codes" overstated independence | Fixed: one framework, re-run independently; `fpp` not re-run |
| F13 | minor | Reason for the superseded 0.011 wrong | Fixed: 0.011 was 93% planet on a companion from the over-luminosity model; the sign error inflated the companion channel ~2× and the brown-dwarf channel 6–7× |
| F14 | minor | "Every other catalogued neighbour" overclaims for the 3.1″ star | Fixed (summary and table) |
| F15 | minor | Two phase-0.5 numbers unlabelled | Fixed: "transit-shaped" vs "occultation-shaped" in the README and the corrections file |
| F16 | minor | followup_2027.png used the analyst's ephemeris; caption's date range wrong | Fixed: figure rebuilt from the write-up's own observability tables (times now match the README table; window 5.0–5.2 h); caption corrected |
| F17 | minor | Mira 7–17× brighter | Fixed: 6–14× (6,100–13,900 against 979 e⁻/s) |
| F18 | minor | Too much jargon | Fixed: "Terms" list after the summary; "hotter than"; position-angle wording |
| F19 | minor | Figures, captions, figscripts and novelty/ not in the manifest; figures not rebuildable | Fixed: all added (including `figscripts/data/rawap_box3.npz`); `common.py` resolves paths in either layout; README says what rebuilding needs |
| F20 | nit | M8 limit marginal; 2,540 K is from the phase-2 bound | Fixed |
| F21 | nit | "Δχ² ≥ 12,600" | Fixed: "≈ 12,600" |
| F22 | nit | White-dwarf K is 35–73 km/s, not 73 | Fixed |
| F23 | nit | fold.png panel c quotes REF but does not plot it; INV 18 vs 19 transits | Fixed: "REF (not shown)"; captions explain 18 (clean) vs 19 (referee) |
| F24 | nit | "Independent" photometries in captions | Fixed: "two separate deblending photometries of the same TESS pixels" |
| F25 | nit | Comment line in the observability CSVs | Fixed: CSVs are header-first; descriptions moved to `observability/README.md`. `w02` re-run: all CSVs identical apart from the removed line |
| F26 | nit | Some replacement numbers not in key_numbers.json | Fixed: sources given inline in the corrections file, and line 7 reworded |
| F27 | nit | novelty/README.md: 64.5 pc, separations without epochs, −3.3 vs −3.4 | Fixed |
| F28 | nit | Pessimistic red-noise case and non-calibrated Δχ² ≥ 108 missing | Fixed |
| F29 | nit | WD Δχ² order; "chance-alignment probability" | Fixed: "+1.5 REF / +1.1 INV"; "chance-alignment parameter R" |

## Still open

* **The repository copy has not been made.** It must be made with `scripts/w06_copy_bundle.py` (or an equivalent path
  rewrite) so that no absolute scratch path is published; search the copy for `/tmp/` once more before committing.
* **The background term (3 × 10⁻⁵) is carried over, not re-derived.** It is the floor of every follow-up step; speckle
  imaging would also constrain it, but that was not modelled.
* **The mass, K, TSM and ESM** are scaled to the adopted radius, not recomputed.
* **The transit duration** from raw-pixel apertures is shorter (1.3–1.5 h at 1.5–2σ) than the density-prior value;
  only a ground-based light curve will settle it.
* **ADS full text** has still not been searched for the star's identifiers.
