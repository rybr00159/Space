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
