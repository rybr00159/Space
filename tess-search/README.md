# TESS search for undiscovered planets and eclipsing binaries (work in progress)

An independent transit / eclipse search of public TESS data, in the spirit of the
[TIC 4206066 Reddit post](../tic4206066/). The search code, intermediate findings and candidate plots are here.
See [`FINDINGS.md`](FINDINGS.md) for the running log of results.

**Status (2026-10-11, in progress):** searching QLP light curves of Sectors 105 and 104 (about 160,000 dwarf and
subgiant stars brighter than TESS magnitude 12, plus cool dwarfs to 13.5). A final report will replace this note.

## Pipeline in brief

1. **Targets:** QLP target lists for S104/S105, with TIC v8.2 parameters from MAST; keep T <= 12 stars with R <= 2.5 Rsun
   (or unknown radius), plus Teff < 4000 K, R < 0.7 Rsun dwarfs to T <= 13.5.
2. **Light curves:** QLP `DET_FLUX`, QUALITY == 0, binned to 10 min; per-camera time windows hit by sector-wide
   systematics are masked (found from thousands of stars at once).
3. **Search:** biweight detrending (0.75 d), box least squares 0.5 d to half the baseline, isolated single events masked,
   best of the top 10 peaks by red-noise S/N, 2P/3P alias check. Validated on known TOIs and with injection-recovery
   (96-99 % recovery for expected S/N > 25 on quiet stars).
4. **Selection:** red-noise S/N >= 9, >= 3 transits, no dominant single event; split into planet-like, eclipsing-binary
   and variable-star tracks (odd/even, secondary, phase-curve uniqueness).
5. **Novelty:** TOI, CTOI, confirmed planets, 341,351 SPOC TCEs, VSX, ASAS-SN, Gaia DR3 EB and variability, TESS EB
   catalogues (Prsa+22, Kostov+25, IJspeert+24), ATLAS, SuperWASP, WISE, SIMBAD and every VizieR catalogue within 15".
6. **Multi-sector test:** every other TESS sector of the star (SPOC 2-min, TESS-SPOC or QLP) is searched near the
   candidate period and compared with the same search at random periods.
7. **Vetting:** transit-model fit and stellar density, held-out-sector predictions, TESScut difference-image
   localisation, Gaia neighbours, TRICERATOPS, Gaia DR3 radial-velocity diagnostics, and adversarial review by an
   independent model (GPT-5.5 via OpenRouter).

## Layout

| Path | What it is |
|---|---|
| `src/` | All pipeline code (Python) |
| `FINDINGS.md` | Running log of candidates and verdicts |
| `candidates/` | Plots and numbers for each vetted candidate |
| `reviews/` | Dossiers sent for independent review, and the reviews |
