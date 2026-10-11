# Review request: candidate discoveries from a TESS search (round 1)

## What was done
- Searched QLP light curves (MIT Quick-Look Pipeline, 200-s FFI cadence, detrended `DET_FLUX`, QUALITY==0, binned to 10 min)
  of TESS Sector 105 (2026) for ~77,000 stars (T<=12 dwarfs/subgiants with R<=2.5 Rsun, plus cool dwarfs T<=13.5).
- Box least squares (0.5-12.6 d), biweight detrend (0.75 d window), sector-wide systematic time windows masked per camera,
  single-event masking, red-noise S/N, odd/even, secondary, phase-curve uniqueness tests.
- Candidates were checked against: ExoFOP TOI and CTOI lists (downloaded 2026-10-11), NASA Exoplanet Archive pscomppars,
  341,351 SPOC TCEs (all public 2-min `dvr-tcestats` tables, S1-S107), VSX, ASAS-SN variables (II/366), Gaia DR3
  eclipsing binaries (I/358/veb) and variability classes (I/358/vclassre), TESS EB catalogue (Prsa+2022, J/ApJS/258/16),
  a VizieR query over *all* catalogues within 15" filtered by variability/eclipse/binary keywords, and SIMBAD object types.
- All other TESS sectors of each star were fetched from MAST (SPOC 2-min > TESS-SPOC FFI > QLP), detrended the same way.
- Difference-image localisation from TESScut FFI cutouts (in-transit vs adjacent out-of-transit frames, flux-weighted
  centroid within 3 px of the target; offset vs catalogue position with 0.15 px systematic floor).

## Candidate A: TIC 143042674 = TYC 7034-792-1 (Gaia DR3 4860406341885274240) — claimed new G + M-dwarf eclipsing binary
- RA 55.15363, Dec -35.17917; T = 10.69, V = 11.47, G = 11.15; TIC Teff 5620 K, logg 4.11, R = 1.46 Rsun, M = 1.0 Msun; d = 269 pc.
- TESS sectors with data: 4 (TESS-SPOC 30-min FFI), 31 (TESS-SPOC 10-min FFI), 105 (QLP, discovery).
- Transit-model fit (batman, quadratic LD fixed u=[0.45,0.2], circular): P = 5.044137 +- 0.000003 d,
  T0 = BTJD 1413.82666 +- 0.00035 (BJD_TDB 2458413.82666), Rc/R* = 0.1384 +- 0.0009, a/R* = 9.57 +- 0.22, b = 0.27 +- 0.08.
  Implied rho* = 0.46 solar (TIC: 0.32). Rc = 0.203 Rsun = 1.97 RJ using TIC R*.
- Primary eclipse 2.2-2.4 %, flat-bottomed, ~4.3 h; odd vs even depths 2.17 % vs 2.24 %.
- Secondary eclipse at phase 0.502, depth 870 ppm (10 sigma).
- Gaia DR3: RUWE 0.97, non_single_star 0, radial_velocity 40.3 +- 3.1 km/s, rv_nb_transits 22,
  rv_amplitude_robust 34.9 km/s, rv_chisq_pvalue 0.0.
- Difference image: S31 offset 2.4" (0.8 sigma), S4 offset 5.0" (1.6 sigma).
- Catalogue checks: no VSX/ASAS-SN/Gaia-EB/TESS-EB entry within 150" with matching period; none of 167 VizieR catalogues within 15"
  is a variability/binary catalogue; SIMBAD type '*'; not a TOI, CTOI, or SPOC TCE.
- Claimed interpretation: an M dwarf (~0.2 Rsun, ~0.18 Msun, Teff ~2800-2900 K from secondary depth) orbiting a G star.

## Candidate B: TIC 123267483 = TYC 7006-2246-1 (Gaia DR3 5017925504414632832) — claimed new detached EB
- RA 28.55235, Dec -33.07382; T = 10.46, G = 10.90; Teff 5878 K; parallax 3.904 mas; RUWE 1.22; no Gaia RV.
- Sectors 3, 30 (TESS-SPOC), 97, 104 (QLP), discovered in S105 (QLP).
- Eclipse timing of 14 primary eclipses over 387 cycles: P = 7.2732209 +- 0.0000020 d, T0 = BTJD 1388.62225 +- 0.00057, rms 1.5 min.
- Primary ~44 % deep in TESS-SPOC, ~37 % in QLP (pipelines dilute differently); secondary ~32 % / ~27 % at phase 0.485 (e cos w ~ -0.024).
- Eclipse duration ~5.5 h total. Difference image: S97 1.8" (0.6 sigma), S104 3.3" (1.1 sigma).
- Only Gaia neighbour within 20": G = 17.8 at 10.6".
- Catalogue checks: none (117 VizieR catalogues within 15", none variability-related; SIMBAD '*').

## Candidate C: TIC 80064194 = TYC 8394-1748-1 (Gaia DR3 6671496606390802176) — claimed new detached EB
- T = 10.55, G = 10.96, Teff 6030 K, parallax 3.626 mas, RUWE 1.23, no Gaia RV, non_single_star 0.
- Sectors 13, 27, 67 (TESS-SPOC/SPOC), 104, discovered in S105.
- Timing of 15 primaries over 438 cycles: P = 5.8106944 +- 0.0000015 d, T0 = BTJD 1658.62910 +- 0.00038, rms 1.4 min.
- Primary ~21 %, secondary ~8 % at phase 0.500; duration ~4.5 h. Difference image: S67 0.4" (0.1 sigma), S104 6.9" (2.2 sigma).
- Catalogue checks: nearest VSX entries are other stars at 85" (RS CVn, P=0.12 d) and 120" (delta Sct); none within 15".

## Candidate D: TIC 72651793 — new periodic transit-like signal (claimed: not a validated planet; most likely an unresolved EB)
- T = 10.49, Teff 6153 K, R = 1.29 Rsun, M = 1.17 Msun, d = 250 pc, contratio 0.059.
- P = 0.6345291 d, depth ~410 ppm, duration 0.83 h, V-ish; independent recoveries in S3, S30, S97 (each FAP < 0.007 vs
  150 random-period nulls); S105 discovery not used in those tests.
- Odd/even 3.5 sigma (odd deeper); fitted rho* 5.7 solar vs TIC 0.54 solar; no secondary.
- TRICERATOPS 1.1 (our folded light curve): FPP = 0.86, NFPP ~ 0; SEBx2P 0.45, SEB 0.26, STP 0.13, TP 0.11.
- Difference image on target to ~1-3" in S97/S30. Gaia star G=19.5 at 6.0".

## Questions
1. For each candidate, is the claim (new EB / new low-mass companion / not-a-planet) supported? What would you dispute?
2. Which catalogues or surveys might already list these eclipsing binaries that the checks above would miss?
3. Is the M-dwarf interpretation of Candidate A sound (radius, secondary depth, Gaia RV amplitude)? Any arithmetic errors?
4. What is the single most valuable additional check for each?
