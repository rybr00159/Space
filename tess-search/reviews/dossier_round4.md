# Review request (round 4): two new planet candidates, and 23 more eclipsing binaries

Context: TESS QLP FFI search of 158,263 dwarf/subgiant stars in Sectors 104/105 (2026), plus a combined S104+S105
search of the 39,999 stars in both (periods 5-25 d). Since your last review:
- All 94 Sector 104 recurrence survivors were vetted (the first pass had only vetted the small-dip dwarf hosts).
- The combined two-sector search finished.
- The EB list grew from 30 to 53 under the same cuts and visual review.

Please be adversarial. For each planet candidate, say whether "new planet candidate (unvalidated)" is a fair label,
what most likely explains it if it is not a planet, and which single follow-up would most change your view. Check the
arithmetic.

## Candidate A: TIC 87055304 = TYC 7912-71-1 (Gaia DR3 6720986204578421888)
- Star:
  - T = 11.72, G = 12.19, BP-RP = 0.86. Gaia DR3 parallax 2.81 +- 0.36 mas, but RUWE = 18.8; ipd_frac_multi_peak = 0;
    non_single_star = 0; no NSS solution.
  - Stellar parameters disagree: TIC Teff 7276 K, R 1.79; Gaia GSP-Phot Teff 5696 K, log g 3.50, R 3.49;
    FLAME R 1.24, M 0.95 (uses the RUWE-degraded parallax).
  - TIC contamination ratio 0.51. The neighbours are a T = 13.5 star 9.0" away and T ~ 17.7-17.9 stars 4.3-5.3" away.
- Transit:
  - P = 4.136195 +- 0.000007 d, T0 = BTJD 1658.4133 +- 0.0027. Sectors 13, 39, 66, 93 and 104 (2019-2026), about 25 transits.
  - QLP SAP flux with a transit-masked 1-day biweight detrend: depth 3,580 +- 120 ppm, flat-bottomed, 4.5 h.
  - Fit: Rp/R* = 0.0668 +- 0.0022 (not corrected for dilution), b = 0.85 +- 0.04, a/R* = 4.30 +- 0.44, so
    rho* = 0.063 solar (circular orbit).
  - Odd 3,529 +- 172 vs even 3,624 +- 169 ppm.
  - Secondary at phase 0.5: 109 +- 203 ppm. Scanning phase 0.05-0.95, the only feature is one event at the start of
    S104 where the sky background is 1.8x normal; it is absent in the other 4 sectors and the other 5 S104 epochs.
  - Out-of-transit harmonics at P: cos 2phi = -267 +- 32 ppm formal, but the same fit at 40 unrelated periods has an
    rms of 209 ppm, so it is not significant.
  - Note: QLP's detrended flux shows +0.5-1 ppt "shoulders" around each transit and 25-45 % shallower depths. This is a
    detrending artefact, absent in SAP flux.
  - Per-sector SAP depths: 4101, 4238, 3792, 3629, 2939 ppm.
- Location:
  - TESScut difference-image centroid 1.1" (0.3 sigma) from the target in S66 and 1.7" (0.3 sigma) in S104.
  - The 9" neighbour is excluded at > 6 sigma. The 4-5" stars would need 81-91 % eclipses.
- Gaia DR3 RVs: 6 transits, rv_amplitude_robust 5.2 km/s, rv_chisq_pvalue 0.93, mean RV -13.3 +- 0.9 km/s.
- TRICERATOPS (QLP DET light curve, TIC star with M = 1.3): FPP 0.037, NFPP < 1e-4 (TP 0.78, PTP 0.16, STP 0.04).
  A rerun with the SAP light curve and R = 2.7 R_sun is pending (stalled on a Gaia archive query).
- Novelty:
  - No TOI, CTOI or confirmed planet within 3' (nearest 24').
  - Not a SPOC TCE. Not in Kostov+25 table 2. Green+23 ellipsoidal score 0.0.
  - None of the 76 VizieR tables within 15" is a variability, EB or planet-candidate table. SIMBAD: '*'.
- Companion radius:
  - 0.067-0.082 R*, depending on dilution correction.
  - 1.2-1.4 R_Jup for R* = 1.8; 1.8-2.2 R_Jup for R* = 2.7.

## Candidate B: TIC 320686138 (Gaia DR3 6641518696737518208)
- Star:
  - T = 11.78, G = 12.19, BP-RP = 0.78. Parallax 2.18 mas (d ~ 458 pc), RUWE 1.21.
  - TIC Teff 6046 K, R 1.37, M 1.12. Gaia GSP-Phot log g 4.19; FLAME R 1.32, M 1.11.
  - Gaia GSP-Spec (Matsuno+24): Teff 6006, log g 4.02, [Fe/H] = +0.41 +- 0.21.
  - TIC contamination ratio 0.05.
- Transit:
  - P = 11.91415 +- 0.00003 d, T0 = BTJD 1668.1611 +- 0.0028. Sectors 13 and 27 (TESS-SPOC), 67, 94 and 104 (QLP);
    9 transits, 2019-2026. Duration 2.2 h.
  - Per-sector depths, SAP/PDCSAP with transit-masked detrend: 2772 (S13, 1 transit, 30-min cadence), 4385, 3542, 3247,
    4078 ppm. Combined 3,782 +- 144 ppm.
  - QLP DET_FLUX gave 1,698 ppm in S104 versus 4,078 in SAP: the same detrending artefact as candidate A.
  - Per-transit depths range from 2.8 to 5.4 ppt, chi2 = 35 for 8 dof against a constant depth, so the errors are
    underestimated about 2x.
  - Within-sector odd minus even: -1128 +- 646, +124 +- 606, -105 +- 532, -2218 +- 597; weighted -782 +- 295 (2.6
    sigma formal, ~1.3 sigma after the 2x error scaling).
  - Secondary at phase 0.5: 4 +- 160 ppm.
  - A central transit would last ~5.8 h for this star, so the orbit is near-grazing. Fits:
    - free, chi2 845: Rp/R* 0.34, b 1.26;
    - a/R* fixed at 16.6 from the TIC density, chi2 851: Rp/R* 0.073, b 0.93 (~1.0 R_Jup);
    - b <= 0.9 forced, chi2 863.
    - Null chi2 1763, 623 points.
- Location:
  - Difference images 3.1" (1.0 sigma, S104), 0.9" (0.3 sigma, S94) and 4.6" (1.1 sigma, S67) from the target.
  - The only Gaia star able to mimic it (T ~ 15.1, 19.1") is excluded at several sigma.
- Gaia DR3 RVs: 12 transits, rv_amplitude_robust 11.8 km/s, rv_chisq_pvalue 0.77.
- TRICERATOPS (clean light curve, TIC star): FPP 0.98 (SEBx2P 0.86, EBx2P 0.10, BEBx2P 0.013, PTP 0.010, TP 0.007).
  We have therefore downgraded B to "most likely a grazing EB" and do not claim it. Do you agree?
- Novelty:
  - No TOI, CTOI or confirmed planet within 3' (nearest 49'). Not a SPOC TCE. Not in Kostov+25 table 2.
  - Green+23 score 0.0. In Doyle+24 (TESS-SPOC main-sequence target list) only as a target.
  - No variability, EB or planet-candidate table among the VizieR entries. No SIMBAD entry.

## Rejected in this round (summary)
- Sector 104: 40 further recurrence survivors.
  - 17 distant giants (Gaia radius 113-278 R_sun, or 3-6.5 kpc).
  - 23 vetted:
    - blends, with 10-46" offsets at 3-11 sigma;
    - EBs (odd/even 4-15 sigma, secondaries at 5-9 sigma, grazing fits, companions > 2 R_Jup);
    - TIC 158159741: P = 1.807 d, V-shaped, b = 0.95, RUWE 6.6 with a Gaia acceleration solution, a 3.6" T = 16.7
      neighbour, TRICERATOPS FPP 0.89.
- Combined search: 4 planet-like signals not recovered in other sectors; 1 EB (odd/even 11 sigma).

## EB additions (23; the list is now 53)
Same cuts as before. The images show all 53 folded at the adopted periods.
Five more sinusoidal or double-wave variables were rejected visually. Period fixes: TIC 406279367, where odd and even
eclipses differ in all 8 sectors, so the period is 2 x 1.562732 = 3.125464 d.
TIC 6598941 has a flat-bottomed 1.7 % eclipse on a 1.86 R_sun star, a 570 +- 47 ppm secondary at phase 0.500, and a
Gaia RV range of 62 km/s, so it is reported as an EBLM-like binary rather than a planet.

| TIC | Other name | T mag | Period (d) | Primary | Secondary (phase) | Sectors | Type | Notes |
|---|---|---|---|---|---|---|---|---|
| 314677221 | - | 11.56 | 3.275130 | 44.4 % | 11.41 % (0.500) | 4 | detached | Gaia RV amplitude 150 km/s |
| 100378891 | 2MASS J20310324-5307568 | 11.99 | 11.20054 | 38.1 % | 8.37 % (0.357) | 5 | detached | eccentric (secondary at phase 0.357, e cos w ~ -0.22); found by the combined two-sector search |
| 217276354 | - | 11.52 | 6.89185 | 27.2 % | 28.70 % (0.526) | 3 | detached | near-twin; eccentric (secondary at phase 0.526); Gaia RV amplitude 116 km/s |
| 234968334 | UCAC2   4235460 | 11.93 | 7.685760 | 23.8 % | 9.23 % (0.500) | 5 | detached | Gaia RV amplitude 108 km/s; RUWE 37.9; found by the combined two-sector search |
| 139256801 | TYC 8448-663-1 | 10.69 | 15.523717 | 21.6 % | 16.10 % (0.722) | 7 | detached | highly eccentric (secondary at phase 0.722, e cos w ~ +0.35); found by the combined two-sector search |
| 76224662 | TYC 8750-903-1 | 11.13 | 4.93228 | 18.8 % | 5.20 % (0.500) | 7 | detached |  |
| 122203391 | - | 11.74 | 4.71438 | 17.3 % | 10.09 % (0.500) | 5 | detached | Gaia RV amplitude 108 km/s; RUWE 10.0 |
| 122852877 | - | 11.16 | 3.053672 | 12.4 % | 11.55 % (0.500) | 5 | detached |  |
| 327599302 | - | 11.77 | 5.020329 | 10.8 % | 2.83 % (0.500) | 10 | detached | Gaia RV amplitude 142 km/s |
| 404124262 | UCAC2  10760330 | 11.45 | 7.59867 | 10.8 % | 0.83 % (0.505) | 5 | detached | flat-bottomed (total) primary eclipse; cool companion (secondary 8 % of primary); found by the combined two-sector search |
| 80453626 | UCAC4 227-000753 | 11.75 | 3.337312 | 10.2 % | 5.72 % (0.500) | 7 | detached | found by the combined two-sector search |
| 304444218 | - | 11.45 | 3.241675 | 10.1 % | 2.92 % (0.500) | 4 | detached |  |
| 218359654 | TYC 7892-1121-1 | 11.10 | 1.391544 | 9.3 % | 7.27 % (0.502) | 3 | detached | plx 0.85 mas: luminous primary (subgiant or giant?); RUWE 7.6 |
| 406279367 | TYC 9119-1863-1 | 11.27 | 3.125464 | 7.8 % | 7.47 % (0.500) | 8 | detached | near-equal eclipses; odd and even depths differ in all 8 sectors, so the orbit is twice the search period; Gaia RV amplitude 175 km/s |
| 372170356 | - | 11.72 | 2.835650 | 7.2 % | 1.31 % (0.500) | 7 | detached | Gaia RV amplitude 159 km/s; RUWE 10.5 |
| 281740175 | UCAC4 155-000744 | 11.79 | 11.05219 | 6.3 % | 0.53 % (0.515) | 10 | detached | flat-bottomed (total) primary eclipse; cool companion; secondary at phase 0.514; RUWE 9.2; found by the combined two-sector search |
| 412733784 | - | 11.02 | 3.854711 | 6.2 % | 0.51 % (0.517) | 6 | detached | slightly eccentric (secondary at phase 0.517); RUWE 6.9 |
| 90919279 | TYC 7914-380-1 | 11.73 | 7.392323 | 4.2 % | 1.45 % (0.512) | 4 | detached | slightly eccentric (secondary at phase 0.513); a T~13.5 star 5.1" away could produce the eclipse and cannot be separated by TESS; Gaia RV amplitude 126 km/s |
| 63594264 | TYC 7907-2230-1 | 11.38 | 1.813884 | 4.2 % | 0.60 % (0.500) | 3 | detached | eclipse 4.3 % deep in S13/S39 but 2.9 % in S104: dilution differs between sectors |
| 314203738 | - | 11.71 | 3.398441 | 4.1 % | 2.50 % (0.500) | 4 | detached | Gaia RV amplitude 149 km/s; RUWE 7.7 |
| 6598941 | TYC 8382-1377-1 | 11.18 | 2.11961 | 1.7 % | 0.06 % (0.500) | 4 | detached | F star + low-mass companion (EBLM-like): flat-bottomed eclipse, companion ~0.25 solar radii; Gaia RV amplitude 62 km/s, as expected for a ~0.25 solar-mass companion (a planet would give < 1 km/s) |
| 33445288 | TYC 6431-1534-1 | 10.99 | 0.7986 | 4.2 % | 1.85 % (0.500) | 4 | contact/ellipsoidal | eclipses on top of ellipsoidal variation; Gaia NSS solution exists (binarity known; eclipses new); Gaia RV amplitude 176 km/s |
| 76440895 | - | 11.26 | 0.677319 | 3.0 % | 1.95 % (0.502) | 6 | contact/ellipsoidal | well-separated eclipses; grouped with close binaries only because P < 0.8 d; RUWE 11.5 |
