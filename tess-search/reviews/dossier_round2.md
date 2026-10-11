# Review request (round 2): TESS search results — new eclipsing binaries and two transit-like signals

## Context
An AI agent searched TESS QLP light curves of Sector 105 (~77,000 dwarf/subgiant stars, T <= 12, plus cool dwarfs to
T <= 13.5) with box least squares, masked sector-wide systematics, and required recurrence in every other TESS sector of
each star (SPOC 2-min > TESS-SPOC > QLP). Novelty checks: ExoFOP TOI and CTOI lists, NASA Exoplanet Archive, all 341,351
public SPOC TCEs, VSX, ASAS-SN variables and EBs (Rowan+22), Gaia DR3 EB and variability classes, Prsa+22 TESS EBs,
Kostov+25 "TESS Ten Thousand" EBs (vetted tables 3/4 and the ~873,000-star unvetted neural-network list, table 2),
IJspeert+24, ATLAS (Heinze+18), SuperWASP (Thiemann+21), Gao+25 TESS periodic variables, WISE variables, SIMBAD object types,
and every VizieR table within 15-20". Localisation by TESScut difference imaging (flux-weighted centroid vs catalogue
position, 0.15 px systematic floor). A previous reviewer (GPT-5.5) judged three earlier EBs "very likely real; newness
unproven" and suggested the catalogue sweep that found the Kostov+25 unvetted list.

## Part 1 — nine eclipsing binaries claimed as "first vetted confirmation and characterisation"
All nine: on target in difference images (worst offset <= 2.2 sigma); present ONLY in Kostov+25 table 2 (unvetted NN list,
no periods); no vetted EB/variable catalogue entry within 21"; SIMBAD lists each as a plain star (or nothing).
Ephemerides from timing of individual primary eclipses across 2018-2026 sectors.

| TIC | Name | T | P (d) | T0 (BTJD) | primary | secondary @ phase | Gaia RUWE | note |
|---|---|---|---|---|---|---|---|---|
| 143042674 | TYC 7034-792-1 | 10.69 | 5.044137 +- 0.000003 (transit fit) | 1413.82666 | 2.2 % flat-bottomed | 0.087 % @ 0.502 (10 sigma) | 0.97 | Rc/R* = 0.1384, b = 0.27, a/R* = 9.57; Gaia rv_amplitude_robust 34.9 km/s (22 transits, p=0) |
| 123267483 | TYC 7006-2246-1 | 10.46 | 7.273221 +- 0.000002 | 1388.6223 | ~32-44 % | ~28-32 % @ 0.485 | 1.22 | near-twin detached |
| 80064194 | TYC 8394-1748-1 | 10.55 | 5.810694 +- 0.000001 | 1658.6291 | ~19-21 % | ~7-8 % @ 0.500 | 1.23 | |
| 122543151 | CD-31 1060 | 9.95 | 7.483459 +- 0.000001 | 1386.5798 | 17.6 % | 10.6 % @ 0.452 | 4.43 | eccentric (e cos w ~ -0.075) |
| 179677500 | TYC 7000-227-1 | 10.82 | 7.911417 +- 0.000003 | 1359.2084 | 11.1 % | 3.3 % @ 0.568 | 1.25 | eccentric (e cos w ~ +0.106) |
| 143345010 | TYC 7940-1151-1 | 10.63 | 2.669561 +- 0.000001 | 1660.4571 | 7.0 % | 1.47 % @ 0.500 | 1.93 | Gaia non_single_star=1, rv_amplitude_robust 140 km/s |
| 89040715 | TYC 8003-1087-1 | 11.02 | 4.621055 +- 0.000001 | 1326.3282 | 5.9 % | 3.3 % @ 0.502 | 8.17 | |
| 200247352 | (no SIMBAD entry) | 11.03 | 3.338011 +- 0.000001 | 1659.3141 | 3.7 % | 1.3 % @ 0.500 | - | Gaia 2-parameter solution |
| 231286114 | TYC 8054-890-1 | 10.55 | 2.020458 +- 0.000001 | 1386.4481 | 1.3 % | 0.77 % @ 0.500 | 0.89 | shallow |

Questions: Is "first vetted confirmation and characterisation" the right claim? Is the G + M-dwarf interpretation of
TIC 143042674 sound (secondary depth 870 ppm with k^2 = 0.0192 implies a TESS-band surface-brightness ratio 0.045,
blackbody brightness temperature ~2,880 K vs a typical ~3,100-3,300 K for a 0.2 Rsun M dwarf)?

## Part 2 — TIC 155991187: robust new transit-like signal
- K dwarf, T = 11.81, TIC Teff 4227 K, R = 0.69 Rsun. Gaia DR3 4975693262795009152: G = 12.58, BP-RP 1.57,
  parallax 6.58 +- 0.26 mas, RUWE 21.1, 54 RVS transits, rv_amplitude_robust 14.2 km/s, rv_chisq_pvalue 0.65,
  non_single_star 0, phot_variable_flag VARIABLE. VSX/Gaia classify the star as a rotational / solar-like variable (no period).
- Co-moving neighbour 8.2" away: Gaia DR3 4975693262795009024, G = 15.88, BP-RP 2.72, parallax 6.74 +- 0.11 mas,
  proper motion (114.3, -14.0) vs target (109.8, -12.7) mas/yr; RUWE 3.46. TIC 155991188: T = 14.58, Teff 3329 K.
- Signal: P = 0.5012110 d, depth ~2.8-3.1 ppt in each of S2, S29, S69, S96, S103, S104 (QLP; 2018-2026), ~0.7 h,
  V-ish; odd/even 0.3 sigma; secondary 43 +- 50 ppm; flat outside transit; Lomb-Scargle shows only the transit's own harmonics.
- Transit fit: Rp/R* = 0.0657 +- 0.003, b = 0.91 +- 0.02, a/R* = 3.12 +- 0.33 -> rho* = 1.6 solar (expected ~2.1 for the K dwarf).
- Difference-image offsets 3.6-7.2" (<= 1.5 sigma); target and companion are only 0.31 TESS px apart.
- Not a TOI/CTOI/TCE; not in Kostov+25 table 2; no EB/variable catalogue match.
- TRICERATOPS 1.1 (M = 0.68, R = 0.69, plx = 6.58 supplied): FPP 0.93, NFPP 0.0; SEBx2P 0.60, STP 0.25, PEBx2P 0.08, TP 0.07.
- Claim: real, previously unreported periodic signal in a likely triple system; most likely an eclipsing binary
  (equal eclipses, P_EB = 1.0024 d) on an unresolved companion; a planet is not excluded.

## Part 3 — TIC 72651793 (from round 1): P = 0.6345 d, ~410 ppm, 0.83 h, recovered in S3/S30/S97; odd/even 3.5 sigma;
rho_fit 5.7 vs TIC 0.54; TRICERATOPS FPP 0.86 (SEBx2P 0.45). Claimed: new signal, most likely an unresolved EB.

## Questions
1. Do you agree with each verdict? Where is the reasoning or arithmetic weakest?
2. For TIC 155991187: what does RUWE 21 at ~150 pc imply about the unresolved companion (separation, mass ranges)? Does
   it make the planet interpretation more or less likely? Is the eclipsing-binary-at-2P scenario physically consistent
   with a 0.7 h, 3 ppt, flat-bottomed-or-V dip and no secondary/odd-even difference?
3. Any additional catalogue or archive that could already contain these objects?
4. What is the single most informative follow-up for each?
