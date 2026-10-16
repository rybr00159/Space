# Review log: TIC 71316387 observing note and supporting report

*2026-10-16. All reviewers were AI agents. Raw output is in [`review_findings.json`](review_findings.json); scratch paths
there are shortened to `<scratch>`.*

## What was run

| Step | Status |
|---|---|
| Corpus study of 33,156 ExoFOP observing notes (community, TFOP and closest-analogue styles) and ExoFOP docs | Done (4 agents) |
| Note: fact sheet + independent audit, 3 drafts, 3 judges, synthesis, fact-check 1 | Done |
| Note: fact-check 2 (TFOP / RAVEN sceptic) and automated fix-up | **Not run** (account usage limit). Done by hand instead; see below |
| Report: number audit, adversarial referee, presentation review | Done (3 agents, about 75 findings) |
| Report: reference check and the verification pass on each finding | **Not run** (usage limit). Done by hand instead: each applied finding was re-checked against the source files, and the references were checked (the RAVEN paper by web search) |

## Report findings and what was done

N = number audit, F = referee, P = presentation. "Fixed" means the change is in v1.0 of the PDF.

| ID | Issue | Action |
|---|---|---|
| N1 | "16 variants" range 1.06–1.54″ included the shallow-only fit (2.46″) | Fixed: 15 deep-eclipse variants; shallow-only quoted separately |
| N2, F2 | V depth prediction should be 0.28–0.33 of TESS | Fixed: V 1,600–1,900 ppm; Fig. 3b band spans 0.28–0.33 at V |
| N3 | R★ = 1.07 is TIC, not Gaia | Fixed: 1.08 R☉ (FLAME) |
| N4, F1, P4 | Fig. 3a plotted the narrow calibration but BP−RP 2.31 came from the wide one | Fixed: Fig. 3a now uses the wide sample (5,560 fits); text explains the non-linear red end. K7 numbers kept as in the source (1.35″, 5.7σ) |
| N5, F27, P1 | ±20% infrared band not drawn | Fixed in the figure |
| N6, F8 | Colour systematic is ±0.1–0.2 mag | Fixed |
| N7 | SuperWASP offset is against a TESS-only ephemeris | Reworded ("relative to a TESS-only ephemeris"). The reviewer's −2.6 min against Table 2 was not independently verified, so not quoted |
| N8, F18 | 0.9″ is a 95% bound, not "most likely" | Fixed |
| N9, F30 | Shape-only odds range | Fixed: "about 10³ (roughly 1–4 × 10³)" |
| N10, F4, P3 | S1/S28/S68 are TESS-SPOC FFI products, not SPOC 2-min | Fixed throughout |
| N11, F10 | Calibration "singles" are imaging non-detections; define the four diagnostics; give upper limit | Fixed |
| N12, F15, P5 | Sign of k; describe the joint fit correctly | Fixed |
| N13, F26, P18 | Table 5 caption | Fixed |
| N14 | UTC times truncated, not rounded | Fixed (eclipse_windows.py now rounds) |
| N15, F24, P16 | S1/S28 excess is in the primaries only | Fixed |
| N16 | G/F V-depth bound is ≈4,650 ppm | Fixed: "4,600 ppm or more" |
| N17, F31 | P_out/P_in lower bound | Reworded "of about 3×10⁴ or more" |
| N18 | Parallax zero-point | Fixed: 293.9 pc after the −0.041 mas correction (recomputed with the `zero_point` package) |
| N19 | CTOI date may be last-modified only | Kept: ExoFOP's change log shows the CTOI created on 2026-08-27 |
| N20, F28, P13 | Author/username placeholders | **Left for the uploader** (see README) |
| N21, F29, P8 | BJD_TDB subscript clipped in Table 5 header | Fixed ("BJD (TDB)") |
| F3 | W1 G/F range should reach 1.35; 8 in-eclipse points | Fixed in text and Fig. 3b |
| F5 | "Follows the final referee's corrections" was not yet true | All listed corrections applied |
| F6 | Title asserted M dwarfs as fact | Fixed: "most likely an M-dwarf pair" |
| F7 | Gaia signature does not argue for binding | Fixed |
| F9 | Light-curve shape gives no upper limit on the light fraction | Fixed (the reviewer's specific flat-range numbers were not used; they were internally inconsistent) |
| F11 | "Confirm" too strong for SuperWASP/W1 | Fixed: "independently support" |
| F12 | Gaia DR4 detection marginal | Fixed |
| F13 | "BEB" ambiguous | Fixed: "false positive (blended eclipsing binary) … not a background star" |
| F14 | Why a ≥0.1 M☉ partner? | Fixed: secondary-eclipse argument added (0.7 M☉ partner → K ≈ 95 km/s; recomputed) |
| F16 | "The referee" undefined | Fixed: identified as an AI re-analysis |
| F17 | Inside-0.1″ posterior fractions | Fixed (Table 4, summary) |
| F19 | Early-K exclusion needs the TESS shift | Fixed |
| F20 | VHS residual could carry the Gaia signature | Fixed (Sections 6.3, 8) |
| F21 | Spectroscopy feasibility | Fixed |
| F22 | Metallicity systematic on M1 | Fixed |
| F23 | Total occultation only in the main-sequence solution; flat bottom is the model | Fixed (Table 3, Fig. 1 label and caption) |
| F25 | Odds wording consistent with the note | Fixed: both use the final referee's baseline range and the ~200 worst case |
| F32–F39 | Wording nits (160 combinations, cross-references, uncited references, statistical-only error, p range across variants, ±5% depth-scale systematic, spot modulation, "raw") | Fixed; every reference is now cited where it is used, and 2MASS, APASS, SkyMapper and Legacy Survey references were added |
| F40 | Campaign `partial/lcmodel_result.json` and `astro_result.json` held each other's contents | Fixed in the repo (commit e5d7c8d) |
| P2 | Gaia DR4 units | Fixed: uses the synthesis's explicit DR4 prediction (3.6 / 1.2 mmag). The synthesis's multi-colour table gives 3,600 / 1,240 ppm for G, which is ~8% different; both are within the stated ±15% |
| P6, P7, P10 | Cross-reference and page-break problems | Fixed: Table 4 reference; the PDF now uses floating figures and tables (8 pages, no half-empty pages) |
| P9, P22 | Sub/superscript collisions; awkward wraps in Table 4 | Fixed (table leading 10.8 pt; one band group per line) |
| P11 | Pair-level values under the "Primary" column | Fixed (cells span two columns) |
| P12 | References not cited | Fixed |
| P14 | Kink from the broad Gaia G point | Fixed (G removed from the prediction curve) |
| P15, P17, P19, P20, P21, P24, P25, P27–P30 | Caption scope, a_out/P_out symbols, "transits" wording, units, J in AO row, legend wording, upright T, 2.1–3.7″, T0 precision, TESSCut spelling, circle radius | Fixed |
| P23 | Arrows in Fig. 2a unlabelled | Partly: caption now says arrow length is arbitrary; sector labels sit on the points |
| P26 | Bullet font, heavy ≲ glyph, justified spacing | Bullets now serif and ≲ removed; justification kept |

## Note: verifier issues and what was done

| Issue (fact-check 1) | Action |
|---|---|
| `[file/tag]` placeholder | Replaced with the report's file name (initials `rb`) |
| Colours were credited with the ~M2V + ~M4V types and 2.6% | Fixed: colours give "red (early-M, possibly late-K)"; types and 2.6% attributed to the eclipse shape with MS relations |
| Odds range 1000–10000 hid the worst case | Fixed: "~5000–14000x for baseline priors, ~200x worst case" |
| Basis of the 18 km/s | Fixed: "for any stellar EB partner" (the report gives the full argument) |
| 1.24″ variant range | Left to the report (length) |
| Depth labels | Kept short; the report gives products and the ±8% systematic |
| AI assistance not disclosed | Added a short clause to the note. The report's AI-assistance section was later removed at the account holder's request, and its wording no longer calls the internal re-checks "referee" passes, so it does not suggest human review |
| Length | 1,222 characters including the 50-character file name: inside the range of comparable multi-part EB notes (IQR 804–1,375) |

Hand check in place of fact-check 2: every number in the final note was traced to the campaign results or recomputed
(the CTOI epoch is at phase 0.5005 of the EB ephemeris, i.e. on a secondary eclipse, 1 min from phase 0.5; Gaia RV statistics from
the DR3 row: 15 transits, 12 visibility periods, 945.6 d, p = 0.2751; K = 17.97 km/s for 0.1 M☉).
