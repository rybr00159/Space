# Final review request (round 3): 30 newly characterised eclipsing binaries + planet-track conclusions

## Search summary
TESS QLP light curves (200-s FFI, detrended, binned to 10 min) of Sector 105 (77,153 stars) and Sector 104 (64,121 stars
not also in S105): dwarfs/subgiants T <= 12 (R <= 2.5 Rsun or unknown) + cool dwarfs T <= 13.5. Box least squares 0.5 d to
half the baseline after masking per-camera systematic windows; candidates required red-noise S/N >= 9, >= 3 transits.
Known-signal rejection: ExoFOP TOI + CTOI, NASA Exoplanet Archive, 341,351 SPOC TCEs. Of the S105 candidates, 327 of 1,420
matched known TOIs/TCEs/planets/EBs-with-TCEs (pipeline sanity check).

## EB selection criteria (all required)
>= 2 TESS sectors; TESScut difference-image centroid within 2.5 sigma of the catalogue position in every measured sector;
primary eclipse >= 0.5 %; coherent ephemeris across 2018-2026 (eclipse timing rms < 15 min, or primary >= 1 % with a
>= 10-sigma secondary); no entry within 21" in VSX, ASAS-SN variables, Rowan+22 ASAS-SN EBs, Gaia DR3 EB, Prsa+22,
Kostov+25 vetted tables, IJspeert+24, ATLAS, SuperWASP, Gao+25, WISE variables, and no period-matched entry within 150";
no variability/binary catalogue among all VizieR tables within 15"; SIMBAD type not variable/binary; visual inspection of
every fold (5 sinusoidal variables removed); period twins within 3' (light leaking into a neighbour) removed.
28 of the 30 appear in Kostov+25 table 2 (~873,000 unvetted neural-network "eclipse-like" flags, no periods); 2 do not.

| TIC | Other name | T mag | Period (d) | Primary | Secondary (phase) | Sectors | Type | Notes |
|---|---|---|---|---|---|---|---|---|
| 146708570 | - | 11.49 | 4.732100 | 34.2 % | 22.60 % (0.500) | 5 | detached | Gaia RV amplitude 121 km/s; RUWE 15 |
| 123267483 | TYC 7006-2246-1 | 10.46 | 7.273221 | 31.8 % | 27.79 % (0.485) | 4 | detached | near-twin |
| 80064194 | TYC 8394-1748-1 | 10.55 | 5.810694 | 19.1 % | 7.26 % (0.500) | 4 | detached |  |
| 122543151 | CD-31  1060 | 9.95 | 7.483459 | 17.6 % | 10.61 % (0.452) | 4 | detached | eccentric; RUWE 4.4 (tertiary?) |
| 114917180 | - | 12.00 | 2.406130 | 13.9 % | 7.59 % (0.500) | 4 | detached |  |
| 142274982 | UCAC2  11682533 | 11.47 | 7.93719 | 12.7 % | 1.39 % (0.507) | 4 | detached |  |
| 197958602 | [R97] GRB 790406 2 | 11.93 | 5.417814 | 11.3 % | 8.60 % (0.500) | 8 | detached | Gaia RV amplitude 180 km/s |
| 179677500 | TYC 7000-227-1 | 10.82 | 7.911417 | 11.1 % | 3.31 % (0.568) | 5 | detached | eccentric |
| 156003824 | TYC 8027-291-1 | 11.64 | 4.969439 | 7.5 % | 0.63 % (0.502) | 6 | detached | cool companion (secondary 8 % of primary) |
| 143345010 | TYC 7940-1151-1 | 10.63 | 2.669561 | 7.0 % | 1.47 % (0.500) | 4 | detached | Gaia NSS solution exists (binarity known; eclipses new) |
| 180444443 | TYC 7001-246-1 | 11.19 | 2.497060 | 6.9 % | 1.05 % (0.502) | 6 | detached |  |
| 89040715 | TYC 8003-1087-1 | 11.02 | 4.621055 | 5.9 % | 3.31 % (0.502) | 4 | detached | RUWE 8.2 (tertiary?) |
| 100666720 | - | 13.39 | 7.341134 | 5.8 % | 0.14 % (0.117) | 4 | detached | M-dwarf pair, 122 pc; period 2x search value |
| 48939280 | - | 13.25 | 2.73651 | 5.5 % | 3.72 % (0.500) | 4 | detached | M-dwarf pair |
| 128677131 | - | 11.53 | 1.251952 | 5.4 % | 2.78 % (0.500) | 5 | detached |  |
| 217648205 | TYC 8341-933-1 | 10.54 | 7.61104 | 4.5 % | 0.31 % (0.382) | 6 | detached | eccentric (secondary at phase 0.38); pulsating primary?; RUWE 18 |
| 396957421 | TYC 8352-1599-1 | 10.94 | 6.329213 | 4.2 % | 1.20 % (0.500) | 6 | detached |  |
| 129032726 | TYC 7955-2777-1 | 11.97 | 3.350564 | 4.0 % | 1.78 % (0.492) | 3 | detached | RUWE 19.5 |
| 320327029 | TYC 8786-1621-1 | 10.55 | 4.440732 | 3.9 % | 0.70 % (0.500) | 6 | detached |  |
| 49617965 | TYC 7016-804-1 | 11.35 | 7.53040 | 3.9 % | 3.90 % (0.750) | 3 | detached | equal eclipses (3.765 d spacing) |
| 200247352 | - | 11.03 | 3.338011 | 3.7 % | 1.31 % (0.500) | 5 | detached |  |
| 197985960 | - | 11.74 | 3.329832 | 3.0 % | 2.45 % (0.500) | 7 | detached |  |
| 278802197 | TYC 7999-451-1 | 11.79 | 1.061827 | 2.8 % | 1.35 % (0.502) | 4 | detached | Gaia NSS solution exists (binarity known; eclipses new) |
| 143042674 | TYC 7034-792-1 | 10.69 | 5.044137 | 2.2 % | 0.09 % (0.502) | 2 | detached | G + M-dwarf companion (EBLM-like); Gaia RV amplitude 35 km/s |
| 261500541 | CD-66  2541 | 9.77 | 3.126963 | 1.6 % | 0.09 % (0.502) | 11 | detached | grazing; cool companion (secondary 5 % of primary) |
| 231286114 | TYC 8054-890-1 | 10.55 | 2.020458 | 1.3 % | 0.77 % (0.500) | 4 | detached |  |
| 115242114 | TYC 7527-529-1 | 11.52 | 3.935908 | 1.1 % | 0.17 % (0.500) | 5 | detached | RUWE 46, plx 0.16 mas (astrometry unreliable) |
| 32398505 | TYC 7003-400-1 | 11.50 | 0.601493 | 2.0 % | 0.66 % (0.950) | 4 | contact/ellipsoidal |  |
| 114787953 | TYC 7951-1401-1 | 11.74 | 0.6057 | 1.0 % | 0.67 % (0.500) | 3 | contact/ellipsoidal | not in Kostov+25 table 2 |
| 71315603 | CD-44 14744 | 11.02 | 0.507699 | 0.6 % | 0.30 % (0.500) | 3 | contact/ellipsoidal |  |

## Planet-track conclusions
- No signal could be validated as a planet. ~200 planet-like single-sector signals were tested for recurrence in every
  other TESS sector; most failed, or were blends (difference-image offsets 15-45"), EBs (secondaries, odd/even),
  pulsators (delta Scuti troughs), rotational variables, or giants.
- TIC 155991187 (K dwarf, 152 pc, RUWE 21, co-moving M dwarf 8.2" away): P = 0.5012110 d, ~2.9 ppt, consistent in 6 sectors
  2018-2026; no odd/even, no secondary, no ellipsoidal signal (< 80 ppm, 3 sigma); directional centroid test disfavours
  the 8.2" companion (chi2 21.8 vs 10.9 for 10 dof); TRICERATOPS FPP 0.93 (SEBx2P 0.60, STP 0.25, TP 0.07).
  Reported as a new transit-like signal, most likely an eclipsing binary within the system.
- TIC 72651793: P = 0.6345 d, 410 ppm, 3 independent sectors; odd/even 3.5 sigma, density mismatch; TRICERATOPS FPP 0.86.
- TIC 245976814: P = 0.5070757 d, ~1.1-1.5 ppt, 4 sectors (2021-2026), F/G star T=11.6, R=1.29; no odd/even, no
  secondary; on target in S104 (1.1 sigma) but 20" (2.8 sigma) offset in S66; faint Gaia stars (T~17.6-19) at 4.7-6.2".
  TRICERATOPS pending. Not in Kostov+25, no SIMBAD entry, not a TOI/TCE.

## Questions
1. Is it fair to call these 30 "first vetted confirmation and characterisation" of eclipsing binaries? Any that you would
   remove or downgrade based on the table (e.g. contact/ellipsoidal ones, high-RUWE/odd parallax ones)?
2. Are the planet-track conclusions appropriately cautious? Anything overstated?
3. One-paragraph summary you would give a non-astronomer about what was and was not discovered.
