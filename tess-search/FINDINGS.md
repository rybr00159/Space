# Running findings log (TESS S104/S105 search)

## New eclipsing binaries (not in VSX / ASAS-SN / Gaia DR3 EB / TESS EB / SIMBAD-variable / any VizieR variability catalogue)

IMPORTANT NUANCE (found after GPT-5.5 review prompted a wider catalogue sweep): all four EBs below appear in
Kostov et al. 2025 (ApJS 279, 50; VizieR J/ApJS/279/50) **table 2**: ~873,000 TESS targets where a neural network saw
"eclipse-like features" (score > 0.9) — listed with TIC/RA/Dec/Tmag only, explicitly *unvetted and unvalidated*, no
periods or depths. None is in that paper's 10,001 vetted EBs (tables 3/4), nor in Rowan+22 ASAS-SN EBs, ATLAS,
SuperWASP (Thiemann+21), Gao+25 TESS periodic variables, IJspeert+24, WISE variables, VSX, Gaia DR3 EB, Prsa+22.
=> Honest claim: first vetted confirmation + characterisation (period, eclipses, on-target localisation), not first
   ever flagging.

### TIC 143042674 = TYC 7034-792-1 (Gaia DR3 4860406341885274240) — G dwarf + M dwarf (EBLM-like)
- RA 55.15363, Dec -35.17917; T=10.69, V=11.47, G=11.15; Teff 5620 K, R=1.46 Rsun, d=269 pc
- Transit-model fit: P = 5.044137 +- 0.000003 d, T0 = BTJD 1413.82666 +- 0.00035; Rc/R* = 0.1384 +- 0.0009,
  a/R* = 9.57 +- 0.22, b = 0.27 +- 0.08 -> rho* = 0.46 solar (TIC 0.32); Rc = 0.203 Rsun; sectors 4, 31 (+ S105)
- Primary: flat-bottomed, 2.2-2.4 % deep, ~4.3 h; odd/even equal (2.17 vs 2.24 %)
- Secondary at phase 0.502, ~870 ppm (10 sigma) -> companion Teff ~2800-2900 K, R ~0.22 Rsun
- Gaia DR3 RV: 22 transits, rv_amplitude_robust = 34.9 km/s, rv_chisq_pvalue = 0 -> stellar-mass companion
  (expected K ~ 22 km/s for 0.2 Msun at 5.04 d); a planet would give < 0.5 km/s
- Difference image: on target (2.4" / 0.8 sigma S31; 5.0" / 1.6 sigma S4)

### TIC 123267483 = TYC 7006-2246-1 (Gaia DR3 5017925504414632832) — detached near-twin EB
- RA 28.55235, Dec -33.07382; T=10.46; G=10.90; Teff 5878 K; d=256 pc (plx 3.904); RUWE 1.22
- P = 7.2732209 +- 0.0000020 d, T0 = BTJD 1388.62225 +- 0.00057 (14 primary eclipses over 387 cycles, rms 1.5 min)
- Primary ~44 % (SPOC) / ~37 % (QLP, different dilution), secondary ~32 % at phase 0.485 (e cos w ~ -0.024)
- Sectors 3, 30, 97, 104; difference image on target (1.8"/0.6 sigma S97, 3.3"/1.1 sigma S104)
- Only neighbour within 20": G=17.8 at 10.6" (cannot produce 40 % eclipses)

### TIC 80064194 = TYC 8394-1748-1 (Gaia DR3 6671496606390802176) — detached EB
- T=10.55, G=10.96, Teff 6030 K, plx 3.626 (276 pc), RUWE 1.23
- P = 5.8106944 +- 0.0000015 d, T0 = BTJD 1658.62910 +- 0.00038 (15 primaries, 438 cycles, rms 1.4 min)
- primary ~19-21 %, secondary ~7.3-8 % at phase 0.500; sectors 13, 27, 67, 104
- Difference image on target in S67 (0.4"), S104 6.9" (2.2 sigma)
- Nearest VSX entries are other stars (85", 120") with unrelated periods

### TIC 146284379 = CD-33 815 — shallow EB (nature unclear: grazing or blended)
- P = 3.85379 d, primary 0.54 %, secondary 0.13 % at phase 0.500; sectors 3, 30, 97, 104

## Rejected / known (examples)
- TIC 281639741 (CD-55 151): Gaia DR3 SB2 orbit -> known binary
- TIC 389527961 (TYC 7964-1312-1): Gaia DR3 SB1 orbit -> known binary
- TIC 410309657 (HD 188319): IJspeert+2021 OBA EB catalogue -> known
- TIC 201251813: VSX KELT KS29C08150 EA, P=1.4091 -> known
- TIC 79271759, 142200409: delta Scuti pulsators (BLS locked onto pulsation troughs)

## Planet-like candidates
### TIC 72651793 (T=10.49, F8, R=1.29) — new periodic signal, most likely an unresolved EB
- P = 0.634529 d, depth ~410 ppm, 0.83 h, V-ish; independent recoveries S3, S30, S97 (each FAP < 0.007)
- odd/even 3.5 sigma (odd deeper); fitted stellar density 5.7 vs TIC 0.54 solar
- TRICERATOPS FPP 0.86 (SEBx2P 0.45, SEB 0.26, STP 0.13, TP 0.11); NFPP ~0
- Difference image on target to ~1-3"
- NOT in Kostov+25 table 2 (the unvetted NN list) nor any EB/variable catalogue
- GPT-5.5 review: (c) likely false positive as a planet; plausible unresolved EB

## Newly characterised eclipsing binaries — summary table (eb_final.py)

All: on target in TESScut difference images; no entry in any vetted EB/variable catalogue within 21"; flagged only in
Kostov+25 table 2 (unvetted NN list); SIMBAD plain star. Depths are from QLP/TESS-SPOC/SPOC data averaged (pipelines dilute differently).

| TIC | Name | T | P (d) | T0 (BTJD) | primary | secondary @ phase | e cos w | Gaia RUWE | note |
|---|---|---|---|---|---|---|---|---|---|
| 143042674 | TYC 7034-792-1 | 10.69 | 5.044136 ± 0.000112 | 1413.8299 | 2.2 % | 0.09 % @ 0.502 | +0.004 | 0.97 | G + M dwarf; Gaia RV amplitude 35 km/s |
| 123267483 | TYC 7006-2246-1 | 10.46 | 7.273221 ± 0.000002 | 1388.6223 | 31.8 % | 27.79 % @ 0.485 | -0.024 | 1.22 |  |
| 80064194 | TYC 8394-1748-1 | 10.55 | 5.810694 ± 0.000001 | 1658.6291 | 19.1 % | 7.27 % @ 0.500 | -0.000 | 1.23 |  |
| 122543151 | CD-31  1060 | 9.95 | 7.483459 ± 0.000001 | 1386.5798 | 17.6 % | 10.61 % @ 0.452 | -0.075 | 4.43 | eccentric; RUWE 4.4 suggests tertiary |
| 179677500 | TYC 7000-227-1 | 10.82 | 7.911417 ± 0.000003 | 1359.2084 | 11.1 % | 3.31 % @ 0.568 | +0.106 | 1.25 | eccentric |
| 143345010 | TYC 7940-1151-1 | 10.63 | 2.669561 ± 0.000001 | 1660.4571 | 7.0 % | 1.47 % @ 0.500 | -0.000 | 1.93 | Gaia NSS=1, RV amplitude 140 km/s (binarity known, eclipses new) |
| 89040715 | TYC 8003-1087-1 | 11.02 | 4.621055 ± 0.000001 | 1326.3282 | 5.9 % | 3.31 % @ 0.502 | +0.004 | 8.17 | RUWE 8.2 suggests tertiary |
| 200247352 |  | 11.03 | 3.338011 ± 0.000001 | 1659.3141 | 3.7 % | 1.31 % @ 0.500 | -0.000 | - | no Gaia parallax (2-parameter solution) |
| 231286114 | TYC 8054-890-1 | 10.55 | 2.020458 ± 0.000000 | 1386.4481 | 1.3 % | 0.77 % @ 0.500 | -0.000 | 0.89 | shallow; grazing or diluted |

### TIC 281705082 (T=11.43, G/K dwarf) — rejected: blend from a neighbour ~40-47" away
- P = 1.457797 d, recovered in 8/8 other sectors (2018-2026), but depth varies with pipeline/aperture (523-1898 ppm)
- Difference image offset 39" (11 sigma) in S104, 29" (2.9 sigma) in S96 -> eclipse is on a neighbour
  (Gaia stars at 46" T~15.6 and 48" T~14.6); fit needs b~1 and rho* 0.11 vs TIC 1.41 -> not this star
- Likely an uncatalogued faint EB on the neighbour, not a planet

### Other planet-track rejections
- TIC 117545527: signal 37" off target; neighbour is in IJspeert+24 TESS OBAF EB catalogue
- TIC 382526005: inconsistent depths (-22..316 ppm), 5.7 sigma secondary, variable star
- TIC 156015928: odd/even 4.8 sigma, secondary 4.3 sigma, variability
- TIC 166810584: inconsistent depths across sectors; EB catalogue entry within 21"

### TIC 155991187 (K dwarf, T=11.81) — new robust transit-like signal in a likely triple system; most likely an EB
- Gaia DR3 4975693262795009152: G=12.58, BP-RP 1.57, plx 6.58 +- 0.26 mas (~152 pc), **RUWE 21.1** (unresolved companion
  very likely), 54 RVS transits, rv_amplitude_robust 14 km/s, rv_chisq_pvalue 0.65 (RV not significantly variable)
- Co-moving M-dwarf companion at 8.2" (Gaia DR3 4975693262795009024, G=15.88, plx 6.74, same proper motion) -> ~1200 au
- P = 0.5012110 d (12.03 h), depth ~2.8-3.1 ppt consistent in S2, S29, S69, S96, S103, S104 (2018-2026); odd/even 0.3 sigma,
  no secondary (43 +- 50 ppm), flat out of transit; transit fit Rp/R* = 0.066, b = 0.91 (grazing), a/R* = 3.1 +- 0.3
  -> rho* 1.6 solar (K dwarf ~2.1; an M-dwarf host would be ~8)
- Difference image offset 3.6-7.2" (<= 1.5 sigma); the 8.2" companion is only 0.31 px away, so TESS pixels cannot separate them
- Not in TOI/CTOI/TCE lists, not in Kostov+25 table 2, no EB/variable catalogue match (VSX/Gaia list the star only as a
  rotational/solar-like variable)
- TRICERATOPS (M=0.68, R=0.69, plx 6.58 supplied): FPP 0.93, NFPP 0.0; SEBx2P 0.60, STP 0.25, PEBx2P 0.08, TP 0.07
- Verdict: real astrophysical signal, previously unreported; most likely an eclipsing binary within the system
  (P_EB = 1.0024 d with equal eclipses) rather than a planet; needs high-resolution imaging + RVs

### More planet-track rejections
- TIC 231663104: distant M giant (d~5.6 kpc), VSX irregular variable; inconsistent depths, odd/even 11.5 sigma
- TIC 139285753: blend (38" offset at 10.5 sigma in S68), inconsistent depths
- TIC 61898320: on target but transit far too short for the R=1.83 Rsun star (rho_fit 8.8 vs 0.20) -> blended EB
- Follow-up tests after Gemini 3.1 Pro review (which argued the 8.2" companion is the source):
  * Directional centroid test with Gaia positions propagated to each sector's epoch, offsets converted to sky coordinates
    (offset_test.py): companion lies at (+7.9", -2.2"); a loss on it should shift the difference-image centroid by
    ~+5..+8" in RA every sector. Measured RA offsets: -2.8, 0.0, +3.4, -3.8, -1.0" (S29, S69, S96, S103, S104).
    Total chi2 over 5 sectors (10 dof): target 10.9, companion 21.8 -> 8.2" companion disfavoured (delta chi2 ~ 11, ~3 sigma).
  * No ellipsoidal/phase-curve modulation: cos(2 pi phase) amplitude -29 +- 25 ppm (3-sigma limit ~80 ppm), cos(4 pi) -12 +- 23 ppm.
    Rules out a bright near-equal-mass EB companion (would give hundreds of ppm); a faint late-M + late-M EB contributing
    < ~3 % of the light (eclipses >~10 %, (R1+R2)/a ~ 0.08 -> ~0.6 h eclipses) remains allowed.
  * Pair is a known co-moving wide binary (SUPERWIDE, Hartman & Lepine 2020, J/ApJS/247/66); not in WDS; no Gaia NSS solution.
- Updated verdict: real, previously unreported transit-like signal on the K dwarf or an unresolved companion of it;
  planet (TP 7 %, STP 25 % in TRICERATOPS) not excluded; SEB with a faint late-M binary companion most likely.
- TIC 101820625: 0.85 % deep, secondary 20.7 sigma, Rp 24 Re -> EB (possibly on a neighbour: 17-22" offsets at ~2.5 sigma)
- TIC 100566758: Rp 11.7 Re, secondary 5.6 sigma, VSX RS CVn variable -> EB / active binary, not a planet

### Sector 104 planet-track survivors (all rejected so far)
- TIC 371800580, 218560701: M giants at 2.6-3.1 kpc (Teff ~3250 K, no TIC radius) -> not planet hosts of interest
- TIC 65096543: offsets 16" at 3.1-3.6 sigma in two sectors -> blend (crowded field with several VSX EBs)
- TIC 200133198: depth varies 306-581 ppm, 15.7" offset (2 sigma), Gaia/VSX EB 31" away -> likely contaminated
- TIC 123507075: duration too long for the star (rho_fit 0.06 vs 1.41), anti-secondary, variable -> not a transit on target
- TIC 218690437: Rp 28 Re, b 1.14, secondary 5.5 sigma -> grazing EB

## EB curation (eb_select.py + visual review of fold plots)
- 123 EB analyses -> 35 pass cuts (>= 2 sectors, on target <= 2.5 sigma, primary >= 0.5 %, coherent ephemeris, no catalogue)
- Visual review removed 5 (pure sinusoids = rotation/pulsation: TIC 237329004, 166225041, 60829240, 422003560; inconsistent fold: 387291911)
- Period fixes: TIC 100666720 (M dwarf, alternating 6.4 % / 5.0 % eclipses in every sector) P = 7.3411 d, not 3.67;
  TIC 49617965 four equal ~4 % eclipses per 15.06 d -> P = 7.5304 d (equal eclipses; 3.7652 d spacing)
- => 30 newly characterised EBs (to be finalised after Sector 104 completes)
