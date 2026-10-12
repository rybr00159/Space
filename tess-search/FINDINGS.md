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
- Not in the final table: eclipse timing is incoherent (rms 92 min) because the 0.5 % eclipses are too shallow to time well

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
  TIC 49617965: P = 3.765198 d (primary 3.7 %, secondary 0.22 % at phase 0.5; BLS first found 4x the period)
- => 30 newly characterised EBs (first batch)
- Second batch, after the Sector 104 EB triage finished (193 analyses -> 52 pass cuts, 22 new):
  - Visual review removed 5 more: TIC 351349985, 206360252, 119045231 (sinusoids; BLS locked onto a harmonic of a
    rotational or pulsating variable) and TIC 314437829, 219931386 (low-amplitude double waves with no distinct eclipse:
    W UMa or ellipsoidal or two-spot rotation cannot be told apart, so not claimed)
  - TIC 406279367: odd eclipses shallower than even ones in all 8 sectors (by 0.25-0.33 %, ~7 sigma combined) and
    nothing at phase 0.5 of the search period (< 0.01 %) -> orbit is 2 x 1.562732 = 3.125464 d with near-equal eclipses
  - TIC 6598941 (TYC 8382-1377-1): flat-bottomed 1.7 % eclipse, P = 2.1196 d, on a 1.86 R_sun F star (TIC) -> companion
    ~0.25 R_sun (~2.4 R_Jup), too large for a planet; secondary 0.08 % at phase ~0.5 in all 4 sectors, far deeper than a
    hot Jupiter's dayside would give (~0.02 %) -> low-mass star (EBLM-like)
  - TIC 63594264: eclipse depth 4.3 % in S13/S39 but 2.9 % in S104 (secondary scales with it) -> dilution differs between
    sectors; centroid within 1.9 sigma. Kept with a caveat
  - => 17 added, 47 in total

### TIC 245976814 (F/G star, T=11.6, crowded bulge-direction field) — new transit-like signal; most likely a blend
- P = 0.5070757 d, 1.1-1.5 ppt consistent in S39, S66, S93, S104 (2021-2026); odd/even 0.7 sigma, secondary 0.2 sigma;
  rho_fit 0.66 vs ~0.52 expected; not a TOI/CTOI/TCE, not in Kostov+25, no SIMBAD entry
- TRICERATOPS (M=1.12, R=1.29, plx=2.87 supplied): FPP 0.19, NFPP 0.0; TP 0.54, BEBx2P 0.17, DTP 0.15, PTP 0.12
  (formally "likely planet" by Giacalone+21 criteria, FPP < 0.5 and NFPP < 1e-3)
- BUT directional difference-image test (epoch-corrected positions, on-sky offsets): S39 and S66 both put the flux loss
  ~21" to the south-east (+13.9/-16.3" and +13.5/-15.1"), matching a G=15.3 Gaia star 16.6" away (needs only ~2.3 %
  eclipses): chi2 on-target 13.7 and 4.9 vs 2.2 and 1.6 for that star; S104 prefers the target (0.9 vs 8.8).
  Total (6 dof): target 19.6, neighbour 12.6. ~20 faint Gaia stars (G 19-20.4) also lie within 8".
- Verdict: real, previously unreported periodic signal, but most likely an eclipsing binary on the G=15.3 neighbour;
  not claimed as a planet. Seeing-limited ground photometry separating the 16.6" star would settle it.

### Sector 104 recurrence survivors not covered above (40 with p_comb < 1e-3; found in the final bookkeeping pass)
The first S104 pass vetted only the small-dip dwarf hosts. The remaining survivors:
- Distant giants (Gaia radius_gspphot 113-278 R_sun, or 3-6.5 kpc away, where T = 10-12 means absolute magnitude
  -1 to -3, i.e. a giant), so a periodic dip cannot be a planet transit on them: TIC 382930116, 60304199, 23594093, 326418285, 121581606, 217175116,
  332194501, 57250866, 348041216, 348168019, 213895533, 56614513, 332189784, 469404629, 61797685
- TIC 201598744 (G star, 0.99 %): secondary 489 ppm at 9 sigma, b = 0.91, Rp >= 1.75 R_Jup -> grazing EB
- TIC 339451799 (G dwarf, 3.0 %): odd/even 5.2 sigma, fitted density 35x the star's -> EB
- TIC 234511777 (M0 dwarf at 129 pc, Gaia RUWE 8.5; 0.65 %): odd/even 14.8 sigma, strongly variable -> EB at 2 x 0.96462 d
- TIC 367745648: 27" offset at 8.5 sigma in S103, depth varies 7x between sectors -> blend
- TIC 263447699 (TIC giant, Gaia RUWE 22.5): odd/even 3.8 sigma, fitted density 0.05 solar -> EB or giant
- TIC 97012825: 32" offset at 9.1 sigma (S103) and 24" at 5.2 sigma (S93); TIC contamination ratio 1.4 -> blend
- TIC 1542186739 (Gaia G = 11.06 at ~830 pc, so M_G ~ 1.5: evolved): grazing fit (b = 1.17), depth differs 2.4x between
  sectors, and a 0.75 d orbit would graze the surface of a ~3 R_sun star -> not a planet
- TIC 122017198 (F star, 0.38 %): secondary 637 ppm at 5.3 sigma, 12" offset at 3.1 sigma in S104, TIC contamination
  ratio 0.72 -> EB, probably blended
- TIC 30775686: 30-35" offsets at 7.7-10.7 sigma in three sectors, odd/even 4.3 sigma -> blend of a neighbouring EB
- TIC 119943912 (2.3 R_sun F star): on target, but the fit needs a 2.6 R_Jup companion on a grazing orbit (b = 0.93),
  with a 230 ppm secondary at 3.1 sigma -> stellar companion (grazing EB), not a planet
- TIC 112079520: offsets 16" (4.6 sigma, S104) and 22" (5.1 sigma, S66), depth varies 2.5x, secondary 6.2 sigma -> blended EB
- TIC 213812460: odd/even 8.9 sigma, 29" offset at 8.8 sigma (S103), depth varies 2.7x -> blended EB at twice the period
- TIC 123912679: secondary 483 ppm at 4.8 sigma, b = 0.95, fitted density 7x too low, offsets 10-16" at ~3 sigma,
  depth varies 2.6x -> blended EB
- TIC 232873972 (0.538 d): neighbours contribute ~2x the target's flux (TIC contamination ratio 1.94), 30" offset at
  3.9 sigma (S104), strong phase variation (brightening at phase 0.5, 5 sigma) -> blended variable
- TIC 200959803: odd/even 15.5 sigma, 42-46" offsets at ~4-5 sigma -> blended EB at twice the period
- TIC 64352125: 25-41" offsets at 3.5-6.8 sigma, depth varies 6x -> blend
- TIC 49442138: grazing fit (b = 1.08), odd/even 4.3 sigma, 27" offsets at ~4 sigma -> blended EB
- TIC 469039599 (0.603 d): ~3 sigma offsets of 12-20" in both measured sectors, neighbours outshine the target in the
  aperture (TIC contamination ratio 1.24), depth varies 2x -> blend
- TIC 247736519 (11,100 K star): only a grazing fit works (b = 1.2) -> grazing EB
- TIC 216975471: odd/even 5.0 sigma, 28" offset at 5.7 sigma (S104), other sectors only a quarter as deep -> blend
- TIC 217581961 (K dwarf): fitted density 50x too low for the star, 34" offset at 10.5 sigma (S93) -> blend

### Novelty check added late: Green+23 ellipsoidal candidates (J/MNRAS/522/29)
- Table 2 lists every TESS target they processed, with a candidate score; the all-VizieR keyword sweep cannot flag it
  because its table description ("Results for all TESS input targets") names no variability.
- 33 of the 47 EBs are in it, all with Score <= 0.19 (none is one of their ~15,000 candidates). TIC 87055304: Score 0.0.
- TIC 158159741 (G star, T = 11.1; P = 1.80689 d, 1.7 ppt, 5 sectors): on target in 3 difference images and no
  odd/even or secondary, but the fit is near-grazing (b = 0.95, V-shaped); Gaia RUWE 6.6 with an astrometric
  acceleration solution and a 34 km/s RV range over 19 epochs (hierarchical system); a T ~ 16.7 star 3.6" away could
  produce the dip with a ~32 % eclipse. TRICERATOPS FPP 0.89 (SEBx2P 0.57) -> most likely an EB, not a planet

## Combined Sector 104+105 search (39,999 stars in both sectors; periods 5-25 d)
- 410 candidates after rescoring (>= 3 transits, red-noise S/N >= 10): 340 EB-like, 42 variable, 28 planet-like.
  After known-signal and VizieR checks, 6 new planet-like and 118 new EB-like signals.
- Planet-like: TIC 128972240, 183588924, 220015693, 231063191 not recovered in any other sector (systematics);
  TIC 387292667 odd/even 11.2 sigma and grazing (b = 1.2) -> EB at 2 x 9.645 d; TIC 320686138 -> candidate (below).
- EB-like: 53 passed broad novelty; 6 passed every cut and visual review and were added to the EB table (53 total),
  all detached with periods 3.3-15.5 d, including two strongly eccentric systems (TIC 100378891, 139256801).

## QLP DET_FLUX depth bias (found while vetting)
- For transits lasting several hours, QLP's detrended flux (DET_FLUX / KSPSAP_FLUX) overshoots around each transit
  (+0.5-1 ppt "shoulders") and underestimates the depth, by up to 60 % in some sectors (TIC 320686138, S104: DET 1,698
  ppm vs SAP 4,078 ppm). Depth-consistency checks across sectors must therefore use SAP flux with a transit-masked
  detrend; the candidate dossiers below do.

## New planet candidates (not validated)
### TIC 87055304 = TYC 7912-71-1 (T = 11.7, F/G star; Gaia DR3 6720986204578421888)
- P = 4.136195 +- 0.000007 d, T0 = BTJD 1658.4133 +- 0.0027; 5 sectors (S13, S39, S66, S93, S104; 2019-2026)
- SAP-flux transit: 3,580 +- 120 ppm (TIC contamination ratio 0.51, so ~5,400 ppm undiluted), flat-bottomed,
  b = 0.85 +- 0.04, a/R* = 4.3 +- 0.4 -> host density 0.063 solar (circular): an evolved F/G subgiant (Gaia GSP-Phot:
  R = 3.5, log g = 3.5; TIC R = 1.79; FLAME R = 1.24 uses a parallax degraded by RUWE 18.8)
- Odd/even 3,529 vs 3,624 ppm (0.4 sigma); secondary at phase 0.5: 109 +- 203 ppm; no secondary anywhere in phase
  (a phase-0.39 feature is one scattered-light event at the start of S104); no ellipsoidal or beaming signal above the
  ~200 ppm noise floor set by fitting 40 unrelated periods
- Difference images on target (S66 1.1", S104 1.7", both 0.3 sigma); the T = 13.5 star 9" away is excluded at > 6
  sigma by S66; the T ~ 17.8 stars at 4-5" would need 81-91 % eclipses
- Gaia DR3 RVs: 6 epochs, robust amplitude 5.2 km/s, constant (p = 0.93) -> a companion orbiting this star with
  K >~ 5 km/s (>~ 50 M_Jup) is disfavoured
- Companion radius 0.067 R* (diluted) to 0.082 R* (undiluted): 1.2-1.4 R_Jup if R* = 1.8, 1.8-2.2 R_Jup if R* = 2.7
- TRICERATOPS (QLP DET light curve, TIC star): FPP 0.037, NFPP < 1e-4 (TP 0.78, PTP 0.16, STP 0.04) -> "likely planet"
  by its criteria (validation needs FPP < 0.015)
- TRICERATOPS rerun with the clean SAP light curve and R* = 2.7 R_sun (the transit density): FPP 0.52, NFPP < 1e-5
  (TP 0.38, EB 0.31, PEB 0.18, PTP 0.08, DEB 0.03). The larger star makes the transiting body star-sized, so on-target
  EB scenarios gain. TRICERATOPS ignores RVs; Gaia's flat RVs (6 epochs, p = 0.93) argue against an EB on this star
  (a stellar companion at 4.1 d gives K >~ 10-20 km/s). FPP is therefore 0.04-0.52 depending on R*, and the candidate
  is not validated
- Not a TOI/CTOI (none within 3'), not a SPOC TCE, not in Kostov+25, Green+23 score 0, no variability catalogue
- Risks: RUWE 18.8 (unresolved companion that could host an EB); host size uncertain, so the companion could be > 2 R_Jup
- Follow-up: RVs (a 1 M_Jup planet gives K ~ 106 m/s); high-resolution imaging

### TIC 320686138 (G star, T = 11.8, [Fe/H] +0.41) — rejected: eclipsing binary (secondary eclipse found)
- P = 11.91415 +- 0.00003 d, ~3,700 ppm (SAP/PDCSAP), 2.2 h, 9 transits in 5 sectors 2019-2026; on target (0.3-1.1
  sigma); no secondary (4 +- 160 ppm); odd/even ~1.3 sigma after allowing for the excess per-transit scatter
- The duration is 0.37x that of a central transit, so the orbit is near-grazing. The free fit prefers Rp/R* 0.34 with
  b 1.26 (chi2 845) over a b = 0.93, 1.0 R_Jup planet at the stellar density (chi2 851)
- Gaia RV constant over 12 epochs (p = 0.77): no massive companion around the G star itself
- TRICERATOPS (clean light curve): FPP 0.98 (SEBx2P 0.86, EBx2P 0.10, TP 0.007) -> most likely two near-equal grazing
  eclipses at 23.8 d on an unresolved companion; not a planet candidate worth claiming
- Settled by the round-4 review's prompt to inspect the phase-0.5 panel: a secondary eclipse of 1,531 +- 155 ppm at
  phase 0.511 (3.1 h after half phase, so e cos w ~ 0.017), present in all 10 epochs from 2019 to 2026. The original
  secondary check sat exactly at phase 0.5 with a +-0.65 h window and missed it. -> EB (primary 0.37 % is below the
  0.5 % cut for the EB table)

### Final Sector 104 driver pass (18 more recurrence survivors with p_comb < 1e-3)
- Distant giants: TIC 194792481, 214494772, 213237492 (Teff ~3,200 K at 2.5-9.5 kpc), 111366822 (Gaia R 148 R_sun),
  234754204, 196375549, 125767306 (BP-RP 5-7 with parallax ~0: reddened giants)
- TIC 218663210 (M dwarf at 86 pc, P = 0.789 d, 0.8 %): 41-46" offsets at 10 sigma (S12, S39), fitted density 20x too
  low -> blend
- TIC 57645079 (M dwarf at 136 pc, P = 4.929 d, 0.8 %): transit 2.2x too long for the star (density 11x too low), 9"
  offset at 2.7 sigma (S66) -> most likely a blend

### EB batch 5 (final Sector 104 triage pass: 456 candidates, 131 pass broad novelty, 64 pass cuts)
- Visual rejects: TIC 218655462 (inconsistent fold), 268886954 (noisy fold, eclipse timing rms 234 min)
- Period doubled from odd/even depths: TIC 212554303 (6/6 sectors, 13 sigma) -> 1.198154 d; TIC 316096209 (3/3, 8 sigma)
  -> 2.296088 d (both M-dwarf pairs, each with an equally bright star 2-3" away)
- Deeper eclipse moved to primary: TIC 363482012, 351474229
- TIC 379508029: no secondary (< 0.04 %), odd/even inconclusive -> period ambiguous (P or 2P), kept with a note
- 9 added -> 62 EBs (56 detached, 6 contact/ellipsoidal); 59 in Kostov+25 table 2; TIC 351474229 is not, but has a
  Gaia NSS solution
- TIC 218815113 (F star, 0.5 %): secondary 798 ppm at 4.3 sigma, 17" offsets at 2.6-4.0 sigma (S66, S104), neighbours
  outshine the target in the aperture (TIC contamination ratio 1.36) -> blended EB
- TIC 351474229: the window-averaged depth of the broad ellipsoidal minimum (2.8 %) exceeds that of the sharp eclipse
  (2.3 %), so the automatic swap in eb_merge.py put phase 0 on the wrong feature; T0 was reset by hand to the sharp eclipse
- TIC 365836059: only a grazing fit works (b = 1.2), odd/even 3.8 sigma -> grazing EB
- TIC 365086729 (0.650 d): odd/even 9.5 sigma, grazing -> EB at twice the period
- TIC 246721728 (9,200 K star): 41" offset at 7.5 sigma (S104), fitted density 13x too low, depth varies 2.6x -> blend
- TIC 32899487: 33" offset at 7.4 sigma (S66), near-grazing, depth varies 2.5x -> blend
- TIC 63610172 (0.556 d): 48" offset at 9.4 sigma, neighbours give 2.5x the target's flux, depth varies 8x -> blend
- TIC 143766264 (0.511 d, 430 ppm): secondary 270 ppm at 8.6 sigma, fitted density 23x too high -> EB or variable
- TIC 21880444: offsets 21-26" at 3.5-5.2 sigma (S13, S39), fitted density 36x too high -> blend
- TIC 172722554: secondary 1,683 ppm at 12.6 sigma, 44" offset at 14 sigma (S103) -> blended EB
- => all 18 rejected; every Sector 104 recurrence survivor (p_comb < 1e-3) has now been vetted

## CORRECTION (2026-10-12): TIC 87055304 withdrawn; the signal is on TIC 87055315 (EBLM-type eclipsing binary)
Full account: tess-search/tic87055304/README.md and tic87055304/analysis/.

**Localisation.**
- First-round two-source fits, pixel-depth maps and injections (tic87/transit) all put the signal on the 9.33"
  neighbour Gaia DR3 6720986208895208960 = TIC 87055315 (G 14.02; 1.29 kpc).
- Independent PRF localisation (loc_independent): joint position 0.66" from the neighbour, 8.8" from the target;
  injection-calibrated z_T = 8.6-10.1 (count floor 3.9 sigma); all 26 transits point to the neighbour.
- Simultaneous PSF photometry (psf_photometry): f_target = 0.16 [0.07, 0.25] +- 0.075 sys, so f = 1 is excluded at
  7.5 sigma; target PSF depth 0.64 +- 0.36 ppt.
- Referee (referee_audit):
  - no shared WCS, registration, jitter or timing systematic;
  - model-free pixel test 6.4 sigma (4.9 conservative; 3.1 most conservative envelope).
- Adopted: target excluded at >= 5 sigma.

**The neighbour-hosted eclipse (neighbour/).**
- Depth 21.4 +- 0.6 ppt on a 1.93 +- 0.09 R_sun star; companion 0.269 +- 0.013 R_sun = 2.61 R_Jup.
- Twin EB at 2P excluded (4.7-5.8 sigma); occultation T_b < 3510 K.
- TRICERATOPS on TIC 87055315: FPP = 1.000.
- Not in any catalogue searched.

**Errors found in our own pipeline.**
- src/centroid.py is biased toward the target, and S66/S93 were dominated by scattered light.
- The ">6 sigma exclusion of the 9" neighbour" divided by the floor-free err_px; the code's own sigma gives ~2.4. It was
  invalid.
- src/run_tri.py passed tdepth ~0.33 (np.convolve 'same' edge effect), so every nearby-star scenario was dropped and
  NFPP was 0 in every run. Fixed; corrected runs are in tri/fixed/:
  - 87055304: FPP 0.977, NFPP 0.977
  - 72651793: 0.918 / 0.399
  - 158159741: 0.873 / 0.162
  - 320686138: 0.981 / 0.250
  - 155991187 and 245976814: pending
