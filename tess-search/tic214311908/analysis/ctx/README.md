# TIC 214311908: independent data, future data and priors (`ctx`)

Part of the adversarial test of the 5.69-d transit-like signal attributed to the M3 dwarf TIC 214311908
(Gaia DR3 5925373272850784768, 64.5 pc), 17.8" W of the Mira TIC 214311906. Ephemeris used throughout:
`/home/user/Space/tess-search/vetting_2026-10-16/results/214311906/ephemeris_v2.json`
(P = 5.6914968 ± 0.0000145 d, T0 = BTJD 3840.9153 ± 0.0019; covariance already scaled by χ²/dof = 53.4/11).
Template: trapezoid with T14 = 1.63 h and a flat part 70% of T14, integrated over each exposure.

Note: the work was done in two sessions. The first wrote `c01`–`c11` and their outputs; the second checked them and
added `c05e`, `c08b`, `c11b`, `c12`–`c18`, the figure and this README. One claim in the first session was slightly
wrong: the five GSFC-ELEANOR-LITE files are not byte-identical. RAW_FLUX, PCA_FLUX and FLUX_BKG are identical, and
CORR_FLUX differs only at a few corrupted cadences. All five files use the same aperture (L_180 at column 1588,
row 268), so they are one blended light curve.

## Main results

### 1. Independent light-curve products (`c01`–`c06`, `c05b`–`c05e`, `c17`)
- **What exists on MAST:**
  - QLP: Sectors 12/39/66/93/103/104.
  - TGLC: Sectors 12 and 39 only, because the HLSP covers Sectors 1–55.
  - TARS: Sectors 12/39/66/93, for the M3V only.
  - GSFC-ELEANOR-LITE: Sector 12, one blended aperture.
  - TEQUILA: Sector 12, the Mira's pixel.
  - None of the following exist: SPOC 2-min/20-s, TESS-SPOC, CDIPS, PATHOS, DIAmante, TASOC.
- **Depths** in `out/fold_products.csv`. Empirical errors come from 400 random-phase depths; S/N uses the empirical error.

| Product (sectors) | Depth (ppm) | S/N | Odd / even (ppm) | Phase 0.5 (ppm) |
|---|---|---|---|---|
| QLP M3V DET (39–104) | 13351 ± 754 | 17.7 | 14516 ± 842 / 11877 ± 947 (2.1σ) | −485 ± 591 |
| QLP pair 214311910 DET | 17587 ± 999 | 17.6 | 18923 ± 1120 / 15966 ± 1233 (1.8σ) | −1496 ± 762 |
| QLP Mira DET | 2268 ± 136 | 16.7 | 2133 ± 164 / 2435 ± 183 (opposite sign, 1.2σ) | 98 ± 113 |
| TGLC M3V psf_flux, S39 | 17170 ± 3156 | 5.4 | 17424 ± 2396 / 16916 ± 2394 | 2693 ± 1693 |
| TGLC M3V aperture, S39 | 19492 ± 4908 | 4.0 | 21484 / 17503 (± 2850) | 4058 ± 2014 |
| TGLC M3V psf, S12 (30-min) | 3548 ± 12106 | 0.3 | – | – |
| TARS M3V (12–93) | 2950 ± 187 | 15.8 | 3092 ± 226 / 2768 ± 256 (0.9σ) | 93 ± 167 |
| ELEANOR-LITE S12 PCA (blend) | 2260 ± 732 | 3.1 | 1701 / 2657 (± 950) | 484 ± 662 |
| TEQUILA S12 (Mira pixel) | 289 ± 1328 | 0.2 | – | – |

- **The odd/even pattern is not consistent across products.** Odd is deeper than even by 20 ± 10% in the M3V and pair
  QLP light curves. It is the other way round in the Mira's QLP light curve, and equal in TGLC. The differences
  look like systematics, not a property of the eclipses.

**TGLC, the Gaia-based deblending, read carefully** (`c05`, `c05c`, `c05d`, `c05e`):
- **TGLC does not deblend this signal.** In Sector 39 the dip appears in the psf_flux of every star near the Mira,
  with similar flux lost:

  | Star | Flux lost (e⁻/s) |
  |---|---|
  | M3V | 15.0 ± 2.9 |
  | Pair A | 13.9 ± 3.2 |
  | Pair B | 13.5 ± 3.2 |
  | Mira | 12.2 ± 5.2 |

  - The reason: TGLC subtracts every neighbour as a static Gaia-flux model and fits each star on its own 5×5-pixel
    cut. A dimming in one star therefore leaks into its neighbours' light curves.
  - So the 1.7% TGLC depth of the M3V is not evidence that the M3V is the host.
- **A linear model of TGLC's fit** (`c05`; checked against the tglc 0.8.1 source, `effective_psf.fit_lc`) gives the
  expected leak coefficients. Fitting the six-star pattern of flux lost, with the covariance taken from random phases:

  | Source placed on | χ² (5 dof) |
  |---|---|
  | M3V | 21.4 |
  | Pair A | 55.1 |
  | Pair B | 62.9 |
  | Mira | 88.1 |
  | G = 17.5 star 3.1″ NE of the M3V | 13.8 |
  | G = 17.7 star 5.5″ N of the M3V | 9.7 |
  | Free position: best fit (−15, +7)″ from the Mira, 3.4″ E and 6.9″ N of the M3V | 8.1 (3 dof) |

  - **Pair and Mira:** TGLC rules out the M-dwarf pair (Δχ² +34) and the Mira (Δχ² +67).
  - **The M3V is disfavoured relative to a position about 7″ further north**, by Δχ² ≈ 13.
  - **Noise simulations** inject the signal on the M3V into the model, using real random-phase noise vectors:
    - the median offset of the recovered position is 2.5″, and 95% fall within 5.5″;
    - only 1 in 300 (p = 0.003) gives a Δχ² as large as the observed one.
  - **The preference survives every model change tried:**

    | Variant | Δχ² (M3V − best) |
    |---|---|
    | PSF broadened, or wings heavier or lighter (`c05e`) | 12–14 |
    | ±0.15-pixel jitter in the world coordinate system (WCS) | 9–12 |
    | Emulation of TGLC's edge mask | 12–15 |
    | M3V image moved 2–4″ by TESS lateral colour | worse |

  - **Sector 12 TGLC adds nothing.** At 30-min cadence, systematics from the Mira and the background dominate.
- **How to read this:**
  - TGLC points the same way as the referee's difference-image fit, whose position is 4.1″ N of the M3V.
  - The `loc` agent's pixel analyses exclude the G = 17.7 star (Δχ² +80 and +51), and exclude the G = 17.5 star only
    at about 3–4σ.
  - TGLC's fit model is not calibrated against this field, and these are 3 transits in one sector. So this is an
    **open tension, not a refutation**.
  - A seeing-limited transit would settle it, because the M3V and the stars 3–6″ to its north are separately
    resolvable.

**TARS** (`c06`):
- TARS is a normalised light curve. **We assume** a 1-pixel aperture; the TARS aperture was not verified against the
  TARS documentation.
- Its depths are converted to total flux lost using the pixel's TESScut flux and the TESS PRF fraction of each
  candidate host that falls in that pixel.
- Sectors 39/66/93 imply a host flux loss of:

  | Host | Flux lost (e⁻/s) | χ² for a constant (2 dof) |
  |---|---|---|
  | M3V | 16.8, 21.3, 19.0 | 2.0 |
  | Pair | 20.7, 34.6, 31.5 | 10.4 |
  | Mira | 29.8, 59.5, 53.8 | 20.0 |

- So TARS weakly favours a constant-brightness host near the M3V. Including Sector 12 (34 ± 6 e⁻/s) raises the M3V's
  χ² to 7.0 for 3 dof.

### 2. Other archival time series (`c07`, `c09`, `c10`, `c10b`)
- **NEOWISE** 2014–2024: 317 epochs, but not usable.
  - Scatter is 8% in W1 and 5% in W2. All detections carry the halo flag `hh00`, from the Mira.
  - Only 3 points fall in transit, giving W1 −5.3 ± 4.6% and W2 +1.7 ± 2.9%.
  - AllWISE multi-epoch photometry has 28 epochs from 2010, none in transit.
- **Gaia**:
  - DR3 epoch photometry exists only for the Mira.
  - **Scanning law:** the 2.5 GB commanded scanning law was streamed without saving it to disk.
    - Validation: all 39 DR3 G-band epochs of the Mira match predicted crossings within 14 s. 19 of the 58 predicted
      crossings in the DR3 window have no published epoch.
  - **DR4** (December 2026; data up to 2020-01-20): 96 field-of-view crossings, **none within 6.7 h of a transit**.
  - **DR5** (~2030): 110 crossings.
    - **One falls fully in transit:** 2024-10-29 18:20 UTC, 13.9 ± 2.5 min after mid-transit.
    - One falls at ingress: 2020-02-14, −42 ± 6.5 min, about 45% of the depth expected.
    - One falls at the egress edge: 2023-07-26.
  - At G = 14.3 a single in-transit epoch should show a 1.5–1.8% dip at ≳10σ. The faint stars 3–6″ away are separate
    Gaia sources, so this would also test which star dims.
- **ASAS-SN**:
  - The Sky Patrol v2 API timed out from here (port 9006), and the service was reported down on 13 October 2026.
  - It would be unusable anyway. Its 8″ pixels blend the Mira (V 12.9–17.0); per-epoch errors are 0.05–0.1 mag at
    V ≈ 15.6–15.9; and the transit covers only 1% of the orbit.
- **ATLAS:** forced photometry needs an account, so it was not used. Its 1.86″ pixels would resolve the Mira.
- **SuperWASP and KELT:** their apertures of tens of arcseconds blend all three stars, and the target is too faint for
  them.
- **MEarth-South:** the star is not a plausible MEarth target; its target list was not checked.
- **ZTF and Pan-STARRS:** the field is outside their footprints.

### 3. Future data (`c11`, `c11b`)
- **TESS:** tess-point 0.9.5.post1 has pointings through Sector 134, which ends 13 September 2028.
  - Those pointings match the HEASARC sector table to 0.005°.
  - They give **no future sector**: the last data are Sectors 103 and 104.
  - The nearest camera centre in Sectors 105–134 is 20.6° away (Sector 114).
  - So there will be no new TESS data before at least September 2028.
- **Gaia DR4** has no in-transit epoch. **Gaia DR5** has 1 full in-transit epoch plus 1 partial one.

### 4. ExoFOP and literature (`c08`, `c08b`, `c16`)
- **ExoFOP:** no TOI, CTOI, planet parameters, imaging, spectroscopy, time series, files or observing notes for
  TIC 214311906, 214311908, 214311910, 1515709322 or 1515709367.
  - The TOI table (8151 rows) has nothing within 19′; the CTOI table (5138 rows) has nothing within 55′.
- **SIMBAD:** only Gaia catalogue papers cite these stars.
- **VSX:** the Mira only (ASASSN-V J172354.64-511925.3, P = 294 d, V 12.91–17.01).
- **EB catalogues:** no entry for any of the stars.
- **arXiv:** no hits on any identifier.
  - ADS full-text search was not possible, because no ADS token is available.
- **New context:**
  - **eRASS1 X-ray source** 1eRASS J172352.9-511928, 3.6″ from the M3V, classified coronal (p = 0.995) with the M3V
    as the most likely counterpart (p_ij = 0.73).
    - L_X ≈ 2.5×10²⁸ erg/s and log L_X/L_bol ≈ −3.4, so the M3V is active, close to saturation.
    - The pair, 11.7″ away, may contribute.
  - The El-Badry et al. (2021) wide-binary catalogue lists the M3V and the pair as a binary.

### 5. Priors (`c14`, `c14b`, `c15`; abstracts in `ext/lit/arxiv_abstracts*.json`, full-text extracts `ext/lit/*.txt`)
**Planet occurrence**, verified from abstracts or full text:
- **Dressing & Charbonneau 2015** (Kepler, early M):
  - Their radius bins stop at 4 R⊕.
  - At 3.5–4 R⊕: 0.36 (+0.72/−0.18)% for P = 1.7–5.5 d, and 0.87 (+1.65/−0.43)% for P = 5.5–18.2 d.
  - "Hot mini-Neptunes and Neptunes are rare around low-mass stars."
- **Gillis, Cloutier, Pass et al. 2026** (TESS, M < 0.4 M☉): no Neptunes found; 95% upper limit of 0.034 per star for
  3–6.5 R⊕ at P < 30 d.
- **Ment & Charbonneau 2023:** fewer than 0.07 planets larger than 1.5 R⊕ per mid-to-late M dwarf (1σ).
- **Hardegree-Ullman et al. 2019:** 0.86 (+1.32/−0.68) planets of 0.5–2.5 R⊕ at P = 0.5–10 d per M3V star.
- **Hardegree-Ullman et al. 2025:** short-period sub-Neptunes peak at Teff 3750 K.
- **Gan et al. 2023:** hot Jupiters around early M dwarfs, 0.27 ± 0.09%.
- **Bryant et al. 2023:** giant planets at 0.26–0.42 M☉, 0.108 ± 0.083%.
- **Sabotta et al. 2021:** planets of 10–100 M⊕, 0.18 per star.
- **Powers et al. 2023:** about 10 confirmed Neptunes of 4–7 R⊕ around M dwarfs.
- **Our synthesis**, not a published number: f(4–7 R⊕, P = 3–10 d) is about 0.1–1% per M3V star, best guess about 0.3%.

**Known planets** (NASA Exoplanet Archive, `out/known_mdwarf_neptunes.csv`):
- 18 transiting planets of 4–7 R⊕ with P < 10 d orbit stars below 4000 K.
- **Only TOI-3884 b** (3180 K, 0.30 M☉, 6.4 R⊕, P = 4.54 d) orbits a star cooler than 3500 K. Only it and TOI-674 b
  (0.42 M☉) orbit stars below 0.45 M☉.
- **If confirmed, the candidate would be the second such planet.** Its prior is low.

**Binaries:**
- **Winters et al. 2019:**
  - Multiplicity is 26.8 ± 1.4% overall, and 28.2 ± 2.1% for 0.30–0.60 M☉.
  - Separations follow a log-normal peaking at 20 au with σ_log a = 1.16. That puts about 67% of companions inside
    65 au, i.e. 1″ at this distance (our calculation).
- **Shan et al. 2015:** 0.11 (+0.02/−0.04) close companions with P = 1–90 d per M dwarf.
- **Moe & Kratter 2021:** binaries closer than 1 au suppress planets; those wider than 200 au do not.

### 6. Follow-up (`c12`, `c13`)
- **The 2027 transit table checks out.**
  - `predicted_transits_2027.csv` agrees with the ephemeris to 0.07 min in time and 0.05 min in σ.
  - Its σ of 4.2–4.9 min already includes the χ²/dof inflation.
- **Which transits can be observed:**
  - "Observable" means altitude above 30° and the Sun below −12° throughout mid-transit ± (T14/2 + 1 h).
  - On that test, 17 transits between October 2026 and December 2027 are observable from at least one of CTIO,
    La Silla, SAAO and Siding Spring, 5–6 per site.
  - All fall in April–August 2027. The first is 2027-04-11 07:36 UTC from CTIO/La Silla, with the Moon 22% lit and
    155° away.
  - **No transit is usable in late 2026.** The best is 2026-11-02, with only 16% of the window observable.
- **Precision** on the M3V, using SkyMapper magnitudes g 15.82, r 14.91, i 13.39, z 12.77:

  | Telescope, band | Noise per exposure | Exposure | Per-event S/N for 1.5% (nominal) | (poor conditions) |
  |---|---|---|---|---|
  | 0.4 m, i | 2.3 ppt | 120 s | 20 | 12 |
  | 1 m, i | 1.4 ppt | 60 s | 25 | 13 |
  | 0.4 m, g | – | – | 11 | 7 |

  - The noise model includes photon noise, sky, read noise and scintillation (Osborn et al. 2015).
  - **Nominal** means a red-noise floor of 1 ppt per 30 min; **poor** means 2″→3″ seeing, airmass 1.8, bright Moon
    and a 2 ppt floor.
  - **One event detects a 1.5–1.8% transit at S/N > 5** in any band with either telescope size.
- **Contamination from the Mira and the pair:**
  - In 2″ seeing with a 3″ aperture, the Mira (at maximum) and the pair together add less than 0.3% to the M3V's flux in i and 0.5% in z.
  - In 3″ seeing it reaches 2–4% at Mira maximum in i/z (Moffat β = 2.5).
  - The Mira's core saturates in i/z in 60 s on a 1 m telescope. Keep the M3V off the Mira's bleed columns, or use a
    bluer band or shorter exposures.

## Scripts (`scripts/`) and outputs
| Script | What it does | Outputs |
|---|---|---|
| `common.py` | ephemeris, transit template, masked biweight detrending, depth fit, random-phase noise | – |
| `c01_mast_products.py` | MAST time-series products within 60″ | `out/mast_timeseries_60arcsec.ecsv` |
| `c02_mast_hlsp_byname.py` | MAST queries by name and by provenance (TGLC, TESS-SPOC, CDIPS, PATHOS…) | `out/c02_mast_byname.txt` |
| `c03_download.py` | download/copy QLP, ELEANOR-LITE, TARS, TEQUILA files | `lc/*.fits` |
| (curl lines) | TGLC files from the MAST bulk-download scripts (cam1-ccd2) | `out/tglc_s0012_lines.txt`, `out/tglc_s0039_lines.txt` → `lc/tglc/*.fits` |
| `c04_fold_products.py` | fold every product: depth, empirical error, odd/even, phase 0.5, flux lost | `out/fold_products.csv`, `out/fold_events.csv` (per transit), `out/fold_arrays.npz` |
| `c05_tglc_leak.py` | TGLC leak model: source hypotheses, free-position grid | `out/tglc_leak.json`, `out/tglc_leak_grid_{psf,aperture}_flux_s{12,39}.npz` |
| `c05b_tglc_plot.py` | raw TGLC psf_flux light curves (diagnostic) | `fig/tglc_raw_psf_flux.png` |
| `c05c_tglc_leak_stars.py` | leak model at every Gaia star within 30″, variants, noise simulation | `out/tglc_leak_stars.json`, `fig/tglc_leak_map_s39.png` |
| `c05d_tglc_leak_robust.py` | robustness: lateral colour, edge mask, WCS jitter | `out/tglc_leak_robust.json` |
| `c05e_tglc_psf_shape.py` | robustness to PSF shape | `out/tglc_psf_shape.json` |
| `c06_tars_pixel.py` | TARS depth → flux lost under each host hypothesis | `out/tars_pixel.json` |
| `c07_gaia_scanlaw_stream.py` | stream the Gaia commanded scanning law; keep samples within 1.5° | `out/gaia_scanlaw_near.csv`, `out/c07_stream.log` |
| `c08_catalogs.py` | ExoFOP, TOI/CTOI tables, SIMBAD, VizieR within 3–5″ | `out/catalog_status.json`, `out/vizier_hits.csv` |
| `c08b_ids_bib.py` | ExoFOP JSON, TIC within 30″, Gaia DR3 flags, SIMBAD bibliography, arXiv identifier search | `ext/exofop/*.json`, `out/tic_30arcsec.ecsv`, `out/gaia_dr3_flags.json`, `out/simbad_bib.json`, `out/arxiv_id_search.txt` |
| `c09_neowise.py` | NEOWISE and AllWISE single exposures folded | `out/neowise_fold.json`, `lc/neowise_m3v.csv`, `lc/allwise_mep_m3v.csv` |
| `c10_gaia_epochs.py` | Gaia field-of-view crossings, in-transit probabilities, DR4/DR5 | `out/gaia_fov_transits.csv`, `out/gaia_in_transit.json` |
| `c10b_validate_scanlaw.py` | scanning-law validation on the Mira's DR3 epochs | `out/gaia_scanlaw_validation.json`, `lc/gaia_epoch/mira_dr3_epoch_G.csv` |
| `c11_future_tess.py` | past and future TESS sectors (tess-point) | `out/future_tess.json` |
| `c11b_future_tess_crosscheck.py` | HEASARC pointings vs tess-point; camera-centre distances | `out/future_tess_crosscheck.json` |
| `c12_photometry_mags.py` | SkyMapper, 2MASS, VHS, APASS, ASAS-SN, VSX magnitudes; Gaia→VRI | `out/mags.json` |
| `c13_followup.py` | 2026–27 transits, observability at 4 sites, precision and S/N model, contamination | `out/followup_transits.csv`, `out/followup.json` |
| `c14_lit_arxiv.py`, `c14b_lit_arxiv_more.py` | arXiv API: abstracts for the priors | `ext/lit/arxiv_abstracts.json`, `ext/lit/arxiv_abstracts_2.json` |
| (pdftotext) | full-text extracts of 8 papers (PDFs deleted) | `ext/lit/*.txt` |
| `c15_known_mdwarf_neptunes.py` | NASA Exoplanet Archive: transiting 3.5–8 R⊕, P < 15 d, Teff < 4000 K | `out/known_mdwarf_neptunes.csv`, `.json` |
| `c16_vizier_extra.py` | details of the VizieR entries (wide binary, eRASS1, etc.) | `out/vizier_extra.json` |
| `c17_fold_figure.py` | folded light curves in the independent products | `fig/fold_products.png` |
| `c18_summary.py` | collects everything | `out/summary_ctx.json` |

### Other files
- **`ext/`**, external pages fetched on 2026-10-16/17:
  - `heasarc_sector.html`: TESS sector pointings and dates.
  - `TESS_orbit_times.csv` and `tess_mit_observations.html`: from MIT.
  - `cosmos_fmsl.html`: listing of the Gaia scanning-law files.
  - `gost_*`: a failed attempt to use the Gaia Observation Scheduling Tool (GOST) web form, replaced by `c07`/`c10`.
  - `asassn/`: Sky Patrol landing page.
  - `exofop/`: ExoFOP target JSON.
  - `lit/`: literature.
- **`ext_src/tglc_whl/`:** tglc 0.8.1 source, used only to read how TGLC fits its light curves.
- **`lc/`:** all downloaded light curves.
  - TESS: QLP, TGLC, TARS, ELEANOR-LITE, TEQUILA.
  - NEOWISE and AllWISE CSVs.
  - Gaia DR3 epoch photometry of the Mira.

## Caveats
- **TGLC:** the leak-model conclusions assume TGLC's fit is linear in a static neighbour model, as its source code
  shows. They use the lab TESS PRF rather than TGLC's own ePSF, which was varied in `c05e`. They rest on one sector
  of 3 transits. The noise simulation does not include systematics that correlate with the signal.
- **TARS:** the 1-pixel aperture is an assumption.
- **Precision model:** throughputs, sky brightness and the red-noise floor are generic assumptions, uncertain by about
  30%.
- **Gaia scanning law:** it is the commanded law. Some crossings yield no epoch (about 33% in the DR3 window).
- **Occurrence-rate prior:** it is a synthesis; no paper measures 4–7 R⊕ at P = 3–10 d specifically for 0.3–0.4 M☉
  hosts.
- **Not accessed:** ATLAS, ASAS-SN g band (Sky Patrol down), ADS full text.
