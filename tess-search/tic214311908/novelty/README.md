# TIC 214311908: novelty and literature re-check (2026-10-17)

This folder re-checks whether the 5.69-d transit candidate on the M3 dwarf TIC 214311908 has already been reported
anywhere. The star is Gaia DR3 5925373272850784768 = 2MASS J17235292-5119271, at 64.3 pc (zero-point-corrected parallax). The check also covers the Mira
TIC 214311906, which lies 17.8" to the east, and the pair TIC 214311910, 11.7" away.

* All queries were read-only GET requests, run on **2026-10-17 between 09:53 and 10:47 UTC**. The exact time of each
  query is the `queried_utc` field of the matching `out/*.json` file.
* Nothing was submitted, uploaded or posted.
* An AI agent ran this check. No professional astronomer has reviewed it.

An earlier check of the same sources is in `../../ctx/` (sections 4 and 5 of its README), and `../../referee_ctx/`
confirmed it. This folder repeats those queries and extends them.

## Verdict

**The candidate is still unreported.** As of 2026-10-17, none of the stars in the table below has any of the following:

* a TOI or CTOI;
* any ExoFOP community data (observing notes, files, imaging, spectroscopy, time series, stellar companions or planet
  parameters);
* an entry in the NASA Exoplanet Archive (`pscomppars`, `ps`, `stellarhosts` or `toi`);
* a planet-type object in SIMBAD;
* a hit in arXiv or a web search;
* a known transit or eclipsing-binary entry in VSX, the Gaia DR3 vari/NSS tables or the TESS EB catalogue.

| TIC | Star |
|---|---|
| 214311908 | the M3V |
| 214311906 | the Mira |
| 214311910 | the pair |
| 1515709322, 1515709367 | the pair's split TIC entries |
| 1515709327 | the 3.1" neighbour |

**Prior mentions of the M3V.** Only catalogues mention it. No paper discusses it.

* SIMBAD lists it only under its Gaia DR3 name, with a single reference: the Gaia Catalogue of Nearby Stars.
* Several catalogue entries are worth citing; they are listed in section 6.

**What changed since the ctx/ check.** One new catalogue matters, and one entry needs context:

1. **eROSITA DR2 / eRASS:3 (Ramos-Ceja et al. 2026, A&A 712, A171).** This is a new catalogue.
   * It moves the X-ray position to between the M3V and pair A.
   * Its best Gaia counterpart is still the M3V, but only at p_i = 0.54. Pair A gets 0.38, and p_any is 0.52.
   * The flux is about half the eRASS1 value.
   * So the statement that the M3V is X-ray active, close to saturation, is now less secure.
   * Details are in section 7.
2. **Maíz Apellániz et al. (2023) flag the M3V as variable in G, BP and RP** (S_G = 8.0 ± 0.7 mmag).
   * This flag is not informative. 94% of comparable nearby M dwarfs within 15° are also flagged variable in G, and the
     M3V sits at the 42nd percentile of S_G.
   * It must not be cited as evidence of activity or of a blend.

## What was checked

| # | Source | Query | Result | Files |
|---|---|---|---|---|
| 1 | ExoFOP-TESS TOI table (`download_toi.php?sort=toi&output=csv`) | All 6 TIC IDs; any row within 60′; free-text search for the TIC, Gaia and 2MASS IDs | 8151 rows, highest TOI 7927.01 (latest alert 2026-10-08, latest update 2026-10-14, modified 2026-10-15). **0 matches.** Nearest is TOI 4398.01 (TIC 214243287), 19.1′ away (see note A) | `out/exofop_status.json`, `raw/toi_rows_within_60arcmin.csv` |
| 2 | ExoFOP-TESS CTOI table (`download_ctoi.php?sort=ctoi&output=csv`) | Same | 5138 rows (latest `lastmod` 2026-10-08). **0 matches.** Nearest is 55.3′ away | `out/exofop_status.json`, `raw/ctoi_rows_within_60arcmin.csv` |
| 3 | ExoFOP target pages | `target.php?id=…&json`, `download_uploads.php`, `download_filelist.php`, `download_imaging/spect/tseries/stellarcomp.php`, `view_obsnotes.php` for all 6 TIC IDs | **All empty.** Every user-upload count is 0, and no observing notes are matched | `raw/exofop_<TIC>_*.{json,txt}` |
| 4 | NASA Exoplanet Archive TAP | `pscomppars`, `ps`, `stellarhosts` in a ±10′ box; `pscomppars` by TIC/Gaia ID; `toi` by `tid` and in a ±30′ box (note B) | **0 planets, 0 TOIs for these stars.** Only TOI 4398.01 is within 30′. `pscomppars` has 6451 rows; `toi` has 8151 rows, last updated 2026-10-09 | `out/nea_status.json`, `raw/nea_*.csv` |
| 5 | SIMBAD TAP | Identifiers, object types, bibliography and measurement tables for the M3V, pair A, the Mira and the 3.1" star; everything within 60″; Pl/EB types within 30′ | See section 3 | `out/simbad_status.json`, `raw/simbad_*.csv` |
| 6 | arXiv API (`all:"…"`, i.e. title, abstract, comments, journal ref) | 28 strings: TIC, Gaia DR2/DR3, 2MASS, WISE, UCAC4, APASS, 1eRASS and PIC names of the three stars; coordinate strings | **0 hits for every string.** Positive controls ("TIC 172900988", "TOI-3884") return 5–8 | `out/literature_search.json`, `out/literature_controls.json` |
| 7 | Web search (general engine) | "TIC 214311908"; "5925373272850784768"; "J17235292-5119271"; "TIC 214311906"/"TIC 214311910"/"ASASSN-V J172354.64-511925.3"; a descriptive query (M3 dwarf, Neptune, 5.69 d) | **No relevant result.** All hits are other objects | summarised here only |
| 8 | Semantic Scholar snippet search (full text of open-access papers) | Main identifiers plus controls | **Not usable:** the first query (the TOI-3884 control) got HTTP 429 (rate limit) on 5 tries over about 8 min, and a separate single request also got 429, so the run was stopped. See note C | `out/literature_controls.json` |
| 9 | NASA ADS (full text) | — | **Not queried:** no ADS API token is available in this environment | — |
| 10 | Gaia DR3 (ESA TAP) | `gaia_source` flags; all 13 `vari_*` tables; all 4 `nss_*` tables for 7 sources; flagged sources within 60″ | **The M3V is in no vari or NSS table.** Details in section 5 | `out/gaia_vari_nss.json`, `raw/gaia_*.csv` |
| 11 | VSX | VizieR copy `B/vsx/vsx`, 5′ cone. The AAVSO VSX API returned a Cloudflare challenge (HTTP 403) | **Nothing at the M3V or the pair**; the only VSX star within 16″ is the Mira. 15 entries within 5′ | `out/vsx_status.json`, `raw/vsx_vizier_5arcmin.csv` |
| 12 | VizieR, all catalogues | 4″ cones at the J2000 and J2016 positions | 63 catalogues, the same set as the ctx/ list (apart from the GUMS simulation). References and rows of the non-survey catalogues fetched | `out/vizier_catalogues_m3v.csv`, `out/vizier_status.json`, `raw/vizier_rows_*.csv` |
| 13 | VizieR, X-ray/UV/radio within 30″ and named catalogues | eRASS1 main, Freund+24, Salvato+25, eRASS:3 (DR2); M-dwarf, variable and EB catalogues by name | See sections 6 and 7 | `out/xray_mdwarf.json`, `raw/vizier_xray_*.csv`, `raw/readme_*.txt` |
| 14 | MAST CAOM | Observations within 5″ | No SPOC 2-min/20-s data, so no SPOC DV report or TCE. Holdings: QLP S12/39/66/93/103/104; TARS S12/39/66/93; GSFC-ELEANOR-LITE S12; ELEANOR postcards and TICA FFIs (note D) | `out/mast_inventory.csv`, `out/mast_status.json` |
| 15 | Maíz Apellániz+23 comparison sample | VizieR cone 15°, 13.8 < G < 14.8, 2.6 < BP−RP < 2.95, plx > 8 mas | 287 stars; 94% flagged variable in G; median S_G 8.6 mmag | `out/gaia_dispersion_context.json`, `raw/maiz23_comparison_mdwarfs.csv` |

**Notes**

* **A. TOI 4398.01** (TIC 214243287, TFOPWG APC) has P = 5.7525 d, close to our 5.6915 d, and lies 19.1′ away (about 55
  TESS pixels). It is a different signal on a different star: ΔP = 0.061 d is thousands of σ, and the separation rules
  out any shared pixels.
* **B. Archive box queries.** The archive's web firewall rejected ADQL `CONTAINS(POINT…)` cone queries sent through this
  environment's proxy (HTTP 403), so RA/Dec boxes were used instead.
* **C. Semantic Scholar.** The keyword `paper/search` endpoint ignores quotes and returned only unrelated papers
  ("TiC" composites, tic disorders), so it was dropped. The snippet endpoint was rate-limited.
* **D. TGLC.** TGLC light curves for S12 and S39 exist (found by ctx/ with a name query). This 5″ position query did not
  return them.
* **E. Planet Hunters TESS.** This project showed volunteers 2-min cadence light curves. This star never had 2-min data,
  so it was almost certainly never shown there. Any PHT candidate would in any case appear in the CTOI table, which has
  none.

## 1. ExoFOP

* All six TIC pages exist and hold only the TESS-provided stellar parameters and magnitudes.
* ExoFOP's name list for the M3V is: TIC 214311908, 2MASS J17235292-5119271, APASS 24986318, Gaia DR2/DR3
  5925373272850784768, UCAC4 194-162330 and WISE J172352.86-511926.1.
* Its TIC contamination ratio is 2.23. This comes from the Mira.
* The obsnotes pages contain only the "Matched Notes for TIC ID …" header, with no notes.

## 2. NASA Exoplanet Archive

* No confirmed planet, no Planetary Systems row and no stellar host lie within the ±10′ box.
* No TOI-table row matches any of the six TIC IDs.
* The only TOI within ±30′ is 4398.01.

## 3. SIMBAD (identifiers and bibliography)

| Object | SIMBAD main id | Type | Identifiers in SIMBAD | References | Measurements |
|---|---|---|---|---|---|
| M3V | Gaia DR3 5925373272850784768 | PM* | **Only the Gaia DR3 ID.** No TIC, 2MASS or other cross-identification | 1: 2021A&A...649A...6G (GCNS) | Gaia DR3 RV 16.28 ± 2.22 km/s; d = 64.53 pc |
| Mira | 2MASS J17235464-5119254 | LP* | Gaia DR3 5925373268528684288, 2MASS | 1: 2023A&A...674A..15L (Gaia DR3 LPVs) | mesVar: LPV, P = 326.17 d (Gaia DR3) |
| pair A | Gaia DR3 5925373272850782208 | PM* | Only the Gaia DR3 ID | 1: GCNS | d only |

The 3.1″ star (Gaia DR3 5925373272827260544) is not in SIMBAD.

Within 60″ of the M3V, SIMBAD has six objects: the three above, two Gaia DR3 eclipsing-binary candidates (G = 19.5 and
20.0, 28″ and 60″ away, P = 0.316 d and 0.261 d) and TYC 8353-2676-1 (53″).

* Those two EBs were already in the search/ analysis of Gaia variables.
* They are too faint to matter. They give at most about 4 e⁻/s in TESS, against the 13.9 e⁻/s flux loss.
* Their periods are unrelated to 5.69 d.

There is no planet-type (Pl, Pl?) object within 30′.

## 4. Literature

* **arXiv:** 0 hits for all 28 identifier and coordinate strings, while the controls work.
* **Web search:** no relevant result.
* **ADS full text:** not checked (no token).
* **Semantic Scholar full text:** blocked by rate limiting.
* So "no paper mentions the star" rests on the following:
  * arXiv metadata;
  * a general web search;
  * the SIMBAD bibliography (1 reference);
  * the VizieR holdings, which list every published table containing the star.
* A paper that names the star only in its full text, without a VizieR table or a SIMBAD link, cannot be excluded.
  **An ADS full-text search for "214311908" and "5925373272850784768" before submission is recommended.**

## 5. Variability and binarity flags

**Gaia DR3, M3V:**

* phot_variable_flag = NOT_AVAILABLE and has_epoch_photometry = False.
* No entry in any `vari_*` table, including `vari_planetary_transit`, `vari_eclipsing_binary` and
  `vari_rotation_modulation`.
* non_single_star = 0, with no entry in `nss_two_body_orbit`, `nss_acceleration_astro`, `nss_non_linear_spectro` or
  `nss_vim_fl`.
* RUWE = 1.076, ipd_frac_multi_peak = 0, ipd_gof_harmonic_amplitude = 0.013.
* RV = 16.28 ± 2.22 km/s from 21 transits. `rv_chisq_pvalue` and `rv_renormalised_gof` are not published at this
  G_RVS.
* Blended or contaminated BP/RP transits: 2 and 1.

**Gaia DR3, pair A:**

* RUWE = 13.7 and ipd_frac_multi_peak = 20, but no NSS solution and no variability entry.
* Pair B (5925373272827253632) has a 2-parameter solution only.

**Gaia DR3, Mira:**

* `vari_long_period_variable`: frequency 0.003066 ± 0.000429 d⁻¹, i.e. P = 326 (+53/−40) d, with amplitude 0.75 mag.
* Classifier: LPV (score 0.98).
* This agrees with our 313–320 d within its error and, like it, disagrees with the VSX 294 d.

**Within 60″:** only the Mira and the 28″ EB candidate are flagged in Gaia DR3.

**VSX** (VizieR `B/vsx/vsx`, 5′ cone):

* Nothing at the M3V or the pair. The nearest entry is the Mira, ASASSN-V J172354.64-511925.3 (type "M:", P = 294 d,
  V = 12.91–17.01).
* The RRab that leaked into one pipeline is **ASASSN-V J172352.39-512245.7** (RRAB, P = 0.648487 d), 3.31′ from the M3V.
  It is the same star as the Gaia DR3 RRab 5925372615693467904 used in referee2_search/ (P = 0.648528 d): the positions
  agree to under 1″.
* No VSX period within 5′ is close to 5.69 d, 11.38 d or 2.85 d.

**Other variable and EB catalogues** (10″ cones), with no entry for the M3V:

* TESS EBs (Prša et al. 2022, J/ApJS/258/16);
* ASAS-SN variables (J/MNRAS/486/1907, J/MNRAS/491/13, II/366);
* Gaia DR3 LPV and classifier tables in VizieR (I/358).

## 6. Catalogue entries for the M3V worth citing

| Catalogue | What it says about the M3V | Use |
|---|---|---|
| Gaia DR3 (Gaia Collaboration, Vallenari et al. 2023) | Astrometry, photometry, RV 16.3 ± 2.2 km/s, RUWE 1.08, no variability or NSS flags | Primary |
| GCNS (Gaia Collaboration, Smart et al. 2021, A&A 649, A6) | In the 100-pc catalogue. This is the only SIMBAD reference for the star | Cite as the only prior bibliographic mention |
| El-Badry, Rix & Heintz 2021 (MNRAS 506, 2269) | Wide binary M3V + pair A: 11.7″ (≈757 au projected), MS+MS, chance-alignment parameter R = 2.0 × 10⁻⁵ | Cite for the bound wide companion |
| Kervella et al. 2022 (A&A 657, A7) | Common-proper-motion companion of pair A, bound flag set (Δv_tan 0.54 ± 0.16 km/s vs v_esc 1.62 km/s) | Cite with El-Badry |
| Merloni et al. 2024 (A&A 682, A34); Freund et al. 2024 (A&A 684, A121); Salvato et al. 2025 (A&A 704, A344) | eRASS1 source 1eRASS J172352.9−511928, coronal (p = 0.995). Counterparts: M3V p_ij = 0.73, pair A 0.27; Salvato+25 Gaia match is the M3V (p_any 0.73, p_i 0.94) | Cite for X-ray, together with eRASS:3 (section 7) |
| Ramos-Ceja et al. 2026 (A&A 712, A171), eRASS:3 / DR2 | 3eRASS J172353.3−511924, 5.8″ from the M3V at the Gaia epoch J2016 (6.3″ with both moved to the survey epoch, section 7). Best Gaia counterpart the M3V (p_i 0.54) vs pair A (0.38); p_any 0.52 | **New since ctx/.** Cite and weaken the activity claim |
| Medan, Lépine & Hartman 2021 (AJ 161, 234) | In their 1.7-million K/M-dwarf Gaia–2MASS–WISE cross-match (2MASS J17235292−5119271, AllWISE J172352.86−511926.1) | M-dwarf catalogue membership; secures the 2MASS/WISE IDs |
| Hardegree-Ullman et al. 2023 (AJ 165, 267) | In their < 120 pc main-sequence table: Teff 3293 ± 110 K, R = 0.383 ± 0.012 R☉, M = 0.374 ± 0.012 M☉, L = 0.0152 L☉. Flagged "Gaia binary / any binary" (the wide pair), not an AO or DR3 binary | Independent published R and M, which agree with ours (0.381–0.385 R☉, 0.370 M☉) |
| Suissa et al. 2020 (ApJ 891, 58) | In their TIC cool-dwarf exposure-time list (Teff 3301 K, R 0.381 R☉, M 0.369 M☉ from TIC 8) | Optional |
| TIC 8 (Stassun et al. 2019) | TIC 214311908, contamination ratio 2.23 | Cite for the IDs |
| Huson et al. 2025 (ApJ 984, 58) | XP [Fe/H] = −0.03 | Do not cite as a metallicity. host/ showed it carries no information for M dwarfs |
| Maíz Apellániz et al. 2023 (A&A 677, A137) | Flag VVV; S_G = 8.0 ± 0.7, S_BP = 18.6 ± 1.3, S_RP = 6.8 ± 0.4 mmag | Not informative (section 5 and the verdict); at most a footnote |
| Kordopatis et al. 2023 (A&A 669, A104); Boettner et al. 2024 (A&A 692, A150) | Galactic orbit (Z_max ≈ 0.95 kpc, e ≈ 0.10); "thick-disk candidate" (TD/D 0.35) | Consistent with host/'s UVW and TD/D. Optional |
| Lebzelter et al. 2023 (A&A 674, A15); ASAS-SN / VSX (Jayasinghe et al.) | The Mira: Gaia DR3 LPV (P = 326 d), ASASSN-V J172354.64-511925.3 (294 d) | Cite for the Mira |

**Catalogues checked by name that do not contain the M3V:**

* TESS Cool Dwarf Catalog (Muirhead et al. 2018);
* Lépine & Gaidos 2011;
* Frith et al. 2013;
* the 10-pc sample (Reylé et al. 2021);
* Andrae et al. 2023 XP metallicities (CatWISE is missing because of the Mira);
* the TESS EB catalogue.

## 7. X-ray: eRASS1 vs eRASS:3

| | eRASS1 (Merloni+24) | eRASS:3 / DR2 (Ramos-Ceja+26) |
|---|---|---|
| Name | 1eRASS J172352.9−511928 | 3eRASS J172353.3−511924 (linked to the eRASS1 source) |
| Counts; DET_LIKE | 7.1; 6.2 | 10.3; 7.6 |
| Position error (`posErr`) | 5.4″ | 3.6″ |
| Separation from the M3V / pair A, both moved to the survey epoch | 4.2″ / 8.9″ | 6.3″ / 5.5″ |
| Best Gaia counterpart | M3V: p_ij 0.73 (Freund); p_i 0.94, p_any 0.73 (Salvato) | M3V p_i 0.54; pair A p_i 0.38; p_any 0.52 |
| F_X (0.2–2.3 keV, erg s⁻¹ cm⁻²) | 5.4 (±2.5) × 10⁻¹⁴ | 2.5 (±1.0) × 10⁻¹⁴ |
| L_X if from the M3V (64.3 pc) | 2.7 × 10²⁸ erg/s, log L_X/L_bol ≈ −3.4 | 1.2 × 10²⁸ erg/s, log L_X/L_bol ≈ −3.7 |

**eRASS:3 supersedes eRASS1.** It stacks the first three surveys, so it contains the eRASS1 photons.

* Its position lies between the M3V and pair A, which is itself a pair of M dwarfs.
* Taking σ ≈ 3.1″ per axis and both stars moved to the epoch of about 2020.7, the position alone favours pair A by about
  1.6:1. NWAY, which also uses priors, gives the M3V 0.54 against 0.38.
* So the X-ray source is real (DET_LIKE 7.6), but **its attribution to the M3V is uncertain, roughly a coin flip with a
  mild preference for the M3V.**
* If it is the M3V, the flux is half what the eRASS1-based statements assumed: log L_X/L_bol ≈ −3.7 instead of −3.4.
  The star would then be moderately active rather than near saturation.

The following earlier statements should be softened or carry the alternative attribution:

* ctx/ §4: "the M3V is active, close to saturation";
* search/: "X-ray active vs kinematically old" tension;
* host/: "L_X/L_bol = 4e-4 implies Ro ≈ 0.17, P_rot ≈ 11–13 d".

This does not affect the transit analysis.

## 8. Limits of this check

* ADS full text was not searched (no token), and the Semantic Scholar full-text search was rate-limited.
* The AAVSO VSX API was blocked; its VizieR copy was used instead.
* TESS working-group and TFOP-internal material is not public. A candidate under internal review that has not reached
  the TOI table cannot be seen.
* QLP and SPOC TCE lists for FFI stars are not public for this star; no SPOC products exist.
* ExoFOP and NEA tables change daily. This check is valid for 2026-10-17 only and should be repeated just before any
  submission.

## Scripts and outputs

All scripts take the folder path as their first argument and are run with `python -I`.

| Script | What it does | Outputs |
|---|---|---|
| `scripts/n01_exofop.py` | TOI and CTOI tables (search by TIC and position, freshness); target JSON, upload counts, file list, imaging/spectroscopy/time-series/stellar-companion summaries and obsnotes page for 6 TIC IDs | `out/exofop_status.json`, `out/n01.log`, `raw/toi_rows_within_60arcmin.csv`, `raw/ctoi_rows_within_60arcmin.csv`, `raw/exofop_<TIC>_*` |
| `scripts/n02_nea_simbad.py` | NASA Exoplanet Archive TAP (pscomppars, ps, stellarhosts, toi); SIMBAD TAP (identifiers, types, bibliography, measurements, 60″ cone, Pl/EB within 30′) | `out/nea_status.json`, `out/simbad_status.json`, `out/n02.log`, `raw/nea_*.csv`, `raw/simbad_*.csv` |
| `scripts/n03_gaia_vsx.py` | Gaia DR3 gaia_source flags, 13 vari tables, 4 NSS tables, flagged sources within 60″; VSX (VizieR + AAVSO API attempt) | `out/gaia_vari_nss.json`, `out/vsx_status.json`, `raw/gaia_*.csv`, `raw/vsx_*` |
| `scripts/n04_literature.py` | arXiv API, 28 identifier and coordinate strings | `out/literature_search.json`, `out/n04.log` |
| `scripts/n04b_controls.py` | arXiv positive controls; Semantic Scholar snippet search (stopped after HTTP 429; the JSON was rebuilt from the log) | `out/literature_controls.json`, `out/n04b.log` |
| `scripts/n05_vizier.py` | All VizieR catalogues within 4″ (two epochs) vs the ctx/ list; ReadMe references; rows of the citable catalogues; X-ray/UV/radio within 30″ | `out/vizier_catalogues_m3v.csv`, `out/vizier_status.json`, `out/n05.log`, `raw/vizier_rows_*.csv` |
| `scripts/n06_mast.py` | MAST product inventory within 5″ | `out/mast_inventory.csv`, `out/mast_status.json`, `out/n06.log` |
| `scripts/n07_xray_mdwarf.py` | eRASS1 / Freund / Salvato / eRASS:3 rows within 30″; named M-dwarf, variable and EB catalogues (10″) | `out/xray_mdwarf.json`, `out/n07.log`, `raw/vizier_xray_*.csv`, `raw/readme_*.txt` |
| `scripts/n08_gaia_dispersion.py` | M3V Gaia photometry counters; Maíz Apellániz+23 comparison sample | `out/gaia_dispersion_context.json`, `out/n08.log`, `raw/maiz23_comparison_mdwarfs.csv` |
