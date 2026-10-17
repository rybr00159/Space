# Corrections to earlier repository statements about TIC 214311906 / TIC 214311908

*Prepared 2026-10-17 from the follow-up campaign in `tic214311908/`, and applied the same day to the three files
below (with small wording changes; the old prediction file was overwritten rather than kept). This page is the record
of what changed. Each item quotes the text as it stood, with its line number then, and gives the replacement. AI-run analysis; not
reviewed by a professional astronomer.*

Sources for the replacement numbers are in `tic214311908/key_numbers.json`, or cited inline where a number appears only
here. Line numbers refer to the files as they stood on 2026-10-17. Three files were checked: `vetting_2026-10-16/README.md`, `README.md` and `FINDINGS.md`. Other
mentions of these stars (`reviews/`, `vetting_2026-10-16/referee/`, `reaudit/`) are records of past reviews and runs;
they need no change.

## Short paragraph pointing to the campaign folder

Suggested for the top of `vetting_2026-10-16/README.md` (below the status line), and in shortened form as a dated
note at the top of `README.md`:

> **Update (2026-10-17).** The planet candidate on TIC 214311908 has been tested further in a separate campaign:
> eleven investigations in three phases, each checked by an adversarial referee. See [`../tic214311908/`](../tic214311908/).
> The signal is on the sky and on the M3 dwarf. The transiting body is 4.95 ± 0.24 (stat) +0.25/−0.18 (sys) R⊕, and
> the false-positive probability is 0.0013–0.0033. It is a strong candidate but not statistically validated: a
> ground-based multi-band transit and high-resolution imaging are still needed. Several numbers on this page are
> superseded (listed in `../tic214311908/corrections_to_repo.md`). The 2027 predictions now use an updated ephemeris
> and correct barycentric UTC times.

## `vetting_2026-10-16/README.md`

| Line | Current text | Problem | Replacement |
|---|---|---|---|
| 21–22 | "It is unconfirmed. One ground-based transit would confirm the host; predictions for 2027 are below." | Outdated. One transit would confirm the host, but validation also needs imaging | "It is a strong candidate but not yet validated (false-positive probability 0.0013–0.0033). A ground-based multi-band transit and high-resolution imaging are needed; see [`../tic214311908/`](../tic214311908/)." |
| 26 | "an M3 dwarf 64 pc away: a ~5–6 R⊕ body every 5.6915 d" | Radius superseded | "an M3 dwarf 64 pc away: a ~5 R⊕ body (4.95 ± 0.24 (stat) +0.25/−0.18 (sys) R⊕) every 5.6915 d" |
| 45 | "folded on 5.6915 d; flat-bottomed, 1.6 h." | Duration superseded | "folded on 5.6915 d; flat-bottomed, 1.7 h." |
| 49 | "Gaia LPV, VSX and ASAS-SN P ≈ 294 d" | Wrong. Gaia DR3 gives 326 ± 46 d (`tic214311908/analysis/followup/README.md`; `tic214311908/novelty/out/gaia_vari_nss.json`). NEOWISE and Gaia epochs give 318.0 ± 2.6 d and exclude 294 d (Δχ² 82) | "VSX and ASAS-SN list P ≈ 294 d, but NEOWISE and Gaia epochs give 313–320 d" |
| 52–53 | "Gaia DR3 5925373272850784768, T = 13.01, 64.5 pc, Teff 3,300 K, R 0.38 R☉, M 0.37 M☉, RUWE 1.08." | Imprecise: the parallax had no zero-point correction, and the star's over-luminosity was not noted | "Gaia DR3 5925373272850784768, T = 13.01, 64.31 ± 0.15 pc, Teff 3,350 ± 55 K, R 0.381 ± 0.012 R☉, M 0.370 ± 0.011 M☉, RUWE 1.08. It is 0.75 mag brighter than the main sequence at its colour, most likely because it is metal-rich ([Fe/H] ≈ +0.2)." |
| 54–55 | "shares its proper motion: an M3–M4 dwarf at 61 pc (G 14.65)." | Imprecise. 61 pc comes from A's own parallax, which is unreliable (RUWE 13.7). The pair is a catalogued bound companion of the M3 dwarf | "is a bound companion of it (El-Badry et al. 2021; chance-alignment probability 2 × 10⁻⁵), 754 au away: an M3–M4 dwarf (G 14.65). Its own parallax (61 pc) is unreliable because of its RUWE." |
| 57 | "The M dwarfs form one bound system about 62 pc away." | Distance superseded | "The M dwarfs form one bound system 64.3 pc away." |
| 60–61 | "Dips every **P = 5.691497 ± 0.000014 d**, T0 = BTJD 3840.9153 ± 0.0019 (BJD_TDB 2460840.9153). This comes from 13 transits in TIC 214311908's QLP light curve" | Ephemeris superseded | "Dips every **P = 5.691506 ± 0.000011 d**, T0 = BJD_TDB 2460840.9150 ± 0.0013. This comes from 21 transit times measured in two deblended photometries of the TESS pixels ([`../tic214311908/ephemeris_adopted.json`](../tic214311908/ephemeris_adopted.json))" |
| 62 | "Flat-bottomed: T14 = 1.61–1.65 h, the flat part about 60–80% of the transit; a V shape is worse by Δχ² = 8." | Superseded by physical transit models and injection calibration | "Flat-bottomed: T14 = 1.70 ± 0.07 h, the flat part 71–74% of the transit (fits with the star's density as a prior; without it, 1.65–1.75 h). Simulated V-shaped (grazing) eclipses reproduce the observed flat-bottom preference in only 0 of 25 and 1 of 25 trials." |
| 63 | "Odd and even depths agree within 1.9σ." | Imprecise | "Odd and even depths agree in every photometry (raw 3×3 pixel box: −1.1 ± 2.2 ppt)." |
| 64–65 | "No dip at phase 0.5 (−700 ± 600 ppm), so half the period is ruled out. No secondary eclipse at any phase: the one feature in the scan, at phase 0.825, is a single artefact in Sector 103." | Superseded by tighter limits | "A transit-shaped dip at phase 0.5 measures −0.98 ± 0.78 ppt, against 14.66 ± 0.78 ppt at phase 0 with the same fixed shape, so half the period is excluded at 14σ. No secondary eclipse: an occultation-shaped dip at phase 0.5 measures +0.3 ± 0.8 ppt (3σ limit 2.8–3.1 ppt)." |
| 66–67 | "The flux lost per transit is 14–17 e⁻/s and stays constant while the Mira's flux changes by a factor of ~2." | Imprecise | "The flux lost per transit is ≈ 15 e⁻/s (13.9 ± 1.4 e⁻/s in a raw 3×3 pixel box that holds ~78% of the M3 dwarf's light; χ² 4.2/5 across sectors). The PRF photometries show no dependence on the Mira's flux, which disfavours a Mira host at 2.7–3.4σ there (p ≈ 0.006 in the referee's Mira-hypothesis fits; `tic214311908/analysis/clean/out/depth_epoch.json`, `tic214311908/analysis/referee_loc/out/mirascale.json`). The Mira is excluded mainly by the pixel image (Δχ² +285) and because it cannot host a 5.69-d orbit." |
| 73–83 | "**Which star dims.** Three independent pixel tests agree on TIC 214311908: …" (the χ² values 66 / 151 / 204, the referee position (−17.2 ± 1.7, +3.7 ± 1.0)″ and the injection χ² 1.6 / 28.9 / 56) | Superseded by injection-calibrated fits for all 93 nearby Gaia stars in two separate codes | "**Which star dims.** Injection-calibrated PRF fits of the transit images (six sectors, two separate codes) put the dimming 0.7–1.6″ from TIC 214311908, with a 95% region of 4.6″. Relative to it, the M-dwarf pair is worse by Δχ² +126 and +146 and the Mira by +285. None of 1,000 injections on these stars reaches those values. A model-free test finds the dip only in the raw pixels around the M3 dwarf and the Mira (0.9 px apart); no pixel more than 2.5 px away dims." |
| 87–89 | "Its extreme colour (BP−RP 4.7) shifts its TESS image, and the calibration was off by 16–17″, about the distance to TIC 214311908." | Wrong explanation. The Mira's colour displaces its TESS image by 5–7.5″. The 16–17″ error came from a degenerate fit of the Mira's direct image: the TESS PRF does not describe the Mira, and the pair is 6″ away | "The calibration was off by 16–17″, about the distance to TIC 214311908. The Mira's colour (BP−RP 4.7) displaces its TESS image by only 5–7.5″; the rest came from a degenerate fit of its direct image, which the TESS PRF does not describe, with the pair 6″ away." |
| 95 | "On TIC 214311908 (about 940 e⁻/s in TESS), the transit is 1.5–1.8% deep." | Superseded normalisation and depth | "On TIC 214311908 (979 e⁻/s in TESS, measured on the same scale as the fits), the transit is 15.8 ± 0.9 (stat) +1.5/−1.1 (sys) ppt deep." |
| 96–97 | "The transiting body is about **5–6 R⊕**, between Neptune and Saturn in size." | Radius superseded | "The transiting body is about **5 R⊕** (4.95 ± 0.24 (stat) +0.25/−0.18 (sys) R⊕), a little larger than Neptune." |
| 98 | "The duration fits the star's density with an impact parameter of about 0.6." | Superseded | "The duration fits the star's density with an impact parameter of about 0.5." |
| 99 | "receiving ~8 times Earth's sunlight; its equilibrium temperature is near 460 K." | Superseded | "receiving 8.2 times Earth's sunlight; its equilibrium temperature is about 470 K." |
| 101–106 | "**What is not yet excluded.** An eclipsing binary unresolved from TIC 214311908 … The PRF forward model fits only roughly (χ² 66 for 19)." | Outdated: the false-positive scenarios have since been modelled one by one | "**What is not yet excluded.** Every eclipsing-binary scenario that was modelled is strongly disfavoured, and the false-positive probability is 0.0013–0.0033. The largest remaining term is a planet on an unresolved ~0.2 M☉ companion, which would still be a planet. Validation still needs a ground-based detection and high-resolution imaging ([`../tic214311908/`](../tic214311908/))." |
| 110 | "On TIC 214311908, the transit should be 1.5–1.8% deep and 1.6 h long." | Superseded | "On TIC 214311908, the transit should be ~1.6% deep and 1.7 h long." |
| 111–112 | "[`results/214311906/predicted_transits_2027.csv`] lists 42 predicted transits for February–October 2027, each with ±4–5 min (1σ) timing." | Wrong. Its `utc_mid` column has no barycentric correction (−3.4 to +7.3 min). The file also uses the superseded ephemeris (σ 4.2–4.9 min) | "[`../tic214311908/predicted_transits_2027.csv`](../tic214311908/predicted_transits_2027.csv) lists 42 predicted transits for February–October 2027, each with ±2.9–3.5 min (1σ) timing, in barycentre-corrected geocentric UTC (updated 2026-10-17). Full observing windows exist only in April–August 2027." Applied by overwriting the old file with the corrected one (same format and epochs); the uncorrected version is in the git history. |
| 113–114 | "If the host is confirmed, the next steps are multi-colour photometry (a planet's transit depth is the same in all colours) and high-resolution imaging." | Outdated plan | "Observe the transit in two bands at once (g and z, or g and i), which tests the host and the colour in one event. Add 8-m speckle imaging and low-precision radial velocities." |
| 117 | "(checked 2026-10-12 and 2026-10-16)" | Outdated | "(checked 2026-10-12, 2026-10-16 and 2026-10-17; ADS full text not searched)" |
| 304–308 | "`results/214311906/`: `ephemeris.json` (uncorrected QLP times) is replaced by `ephemeris_v2.json`; …" | Incomplete | Add: "`ephemeris_v2.json` is superseded by `../tic214311908/ephemeris_adopted.json`, and `predicted_transits_2027.csv` was recomputed from it on 2026-10-17 (the uncorrected version is in the git history)." |
| 314–316 | "The host rests on TESS pixel tests of modest S/N next to a bright, very red, variable star. An unresolved eclipsing binary at the same position is not excluded. The planet's size depends on the star's TESS brightness and on dilution, both estimated." | Outdated | "The host and the false-positive probability (0.0013–0.0033) rest on TESS pixels next to a bright, very red, variable star; no ground-based detection or high-resolution imaging exists yet. The planet's radius carries a +0.25/−0.18 R⊕ systematic from the depth normalisation and the photometry method." |

## `README.md`

| Line | Current text | Problem | Replacement |
|---|---|---|---|
| 19–21 | "However, the signal found on the Mira TIC 214311906 comes from the neighbouring M dwarf TIC 214311908, where it is a **Neptune-sized planet candidate**." | Right, but "Neptune-sized" is loose (4.95 R⊕ is 1.3 Neptune radii), and the campaign result is missing | Change "**Neptune-sized planet candidate**" to "**planet candidate a little larger than Neptune**", and add a new dated note at the top: "**Update (2026-10-17).** A follow-up campaign on TIC 214311908 finds a 4.95 R⊕ planet candidate with a false-positive probability of 0.0013–0.0033. It is strong but not yet validated; see [`tic214311908/`](tic214311908/)." |
| 34 | "A ~5–6 R⊕ body transits the M3 dwarf TIC 214311908 (64 pc) every 5.691497 d, for 1.6 h each time." | Radius, period and duration superseded | "A ~5 R⊕ body (4.95 ± 0.24 (stat) +0.25/−0.18 (sys) R⊕) transits the M3 dwarf TIC 214311908 (64 pc) every 5.691506 d, for 1.7 h each time." |
| 35–36 | "Three pixel-level tests place the source on TIC 214311908." | Superseded | "Injection-calibrated pixel fits place the source on TIC 214311908; two deblending pipelines and a model-free raw-pixel test confirm the dip in every sector at the M3 dwarf/Mira position." |
| 38–39 | "It is unconfirmed. One ground-based transit would confirm the host. See [`vetting_2026-10-16/`]." | Outdated | "It is a strong candidate (false-positive probability 0.0013–0.0033) but not validated. A ground-based multi-band transit and high-resolution imaging are needed. See [`tic214311908/`](tic214311908/)." |
| 454 | "Not on the Mira. A ~5–6 R⊕ planet candidate transiting the M3 dwarf TIC 214311908, 17.8″ away (host still to be confirmed from the ground)" | Radius superseded | "Not on the Mira. A ~5 R⊕ planet candidate transiting the M3 dwarf TIC 214311908, 17.8″ away (false-positive probability 0.0013–0.0033; not yet validated; [`tic214311908/`](tic214311908/))" |
| 472–475 | "**Observe one transit of the TIC 214311908 planet candidate from the ground.** … It should dim by 1.5–1.8% for 1.6 h. The 42 predicted transits for February–October 2027 are in [`vetting_2026-10-16/results/214311906/predicted_transits_2027.csv`]." | Depth and duration superseded; the prediction file had uncorrected UTC times | "**Observe a transit of the TIC 214311908 planet candidate from the ground, in two bands at once.** Seeing-limited images separate it from the Mira (17.8″) and the M-dwarf pair TIC 214311910 (11.7″). It should dim by ~1.6% for 1.7 h. Full windows exist only in April–August 2027, for example 2027-06-01 12:56 UTC (Siding Spring), 2027-06-07 05:32 and 2027-08-03 03:32 UTC (Chile). All predictions are in [`tic214311908/predicted_transits_2026_2028.csv`](tic214311908/predicted_transits_2026_2028.csv). High-resolution imaging is the next step." |
| after 500 (end of the Files table) | (no entry) | Missing folder | Add the row: "\| `tic214311908/` \| Planet candidate TIC 214311908 (found on the Mira TIC 214311906): full adversarial campaign, adopted parameters, false-positive analysis, 2026–2028 transit predictions and observing plan \|" |

## `FINDINGS.md`

`FINDINGS.md` is a dated working log, so its old entries should stay as they are. Each item below says what is now
superseded; the replacements are collected into one new dated entry at the end of the list.

| Line | Current text | Problem |
|---|---|---|
| 425–426 | "the sources are within ~10" and ~5-8" respectively, and unidentified; … 214311906 candidate: an M-dwarf binary 6.5" away (Gaia DR3 5925373272850782208, 61 pc, RUWE 13.7)." | Superseded: the 214311906 source is the M3 dwarf TIC 214311908, 17.8″ from the Mira |
| 431–432 | "TIC 214311910 (Gaia DR3 5925373272850782208; 61 pc; RUWE 13.7; ipd_multi 20%) … TIC 214311908 (M3V, 64.5 pc) is a CPM companion 17.8" W." | Distances imprecise: the system is at 64.31 ± 0.15 pc, and the pair's own parallax is unreliable (RUWE 13.7) |
| 435–437 | "Real: Mira/pair/M3V/bright41 = 0/15/57/122. … So the source is at the Mira/pair position, not the M3V" | Wrong (already corrected at line 452); the Mira-based calibration was degenerate |
| 438–439 | "Shape: flat-bottomed (T23/T14 0.78 +- 0.09), T14 1.61 h, so k <~ 0.2 … Implied host density ~12-15 rho_sun (M dwarf)." | Superseded: T14 1.70 ± 0.07 h, T23/T14 0.71–0.74 (with the stellar-density prior; 1.65–1.75 h without); the host density is 6.7 ρ☉ (9.4 g/cm³; `tic214311908/analysis/star/star_params.json`), consistent with the transit |
| 440–441 | "Ephemeris P = 5.691471 +- 0.000018, T0 = BTJD 3840.9101 +- 0.0022 (14 events). 2027 predictions are in results/214311906/predicted_dips_2027.csv." | Superseded ephemeris (uncorrected QLP S103/S104 times). The file `predicted_dips_2027.csv` is not in the repository |
| 452 | "the M3V TIC 214311908 (64.5 pc, R 0.381)" | 64.31 ± 0.15 pc (zero-point-corrected parallax) |
| 453–454 | "Our earlier map used a WCS correction fitted on the Mira (BP-RP 4.7). It was off by 16-17" in S93/S103/S104" | The offset is right, but its cause is not colour: the Mira's colour shift is 5–7.5″; the rest came from a degenerate direct-image fit |
| 456–460 | "Forward aperture model chi2 = 66 (M3V) / 151 (pair) / 204 (Mira) … Delta F = 16.5 +- 1.1 e/s" | Superseded by the campaign's injection-calibrated localisation and photometry |
| 461–462 | "Corrected ephemeris: P 5.691497 +- 0.000014, T0 3840.9153 +- 0.0019 BTJD. 42 predicted transits Feb-Oct 2027 (results/214311906/predicted_transits_2027.csv)." | Superseded ephemeris; the prediction file's UTC column lacked the barycentric correction (−3.4 to +7.3 min) |
| 463 | "No phase-0.5 dip (-712 +- 591 ppm) … Odd/even within 1.9 sigma." | Superseded by tighter limits |
| 464 | "Rp ~5-6 Re (depth 1.5-1.8% on the M3V; about 940 e/s for T 13.0); b ~0.6; Teq ~460 K." | Superseded |

Suggested new entry, appended after line 471:

```
**Follow-up campaign on TIC 214311908 (2026-10-17).** Full write-up: tess-search/tic214311908/.
- Eleven investigations in three phases (photometry, localisation, star, context; transit model, planet search,
  host, FPP; cleaned photometry, brown-dwarf channel, follow-up design), each checked by an adversarial referee.
- Signal on sky: raw-pixel dip 13.9 +- 1.4 e/s in the 3x3 box holding the M3V and the Mira's core, all six sectors
  (chi2 4.2/5); no pixel beyond 2.5 px dims. On the M3V: injection-calibrated PRF fits give pair A +126, pair B +146,
  Mira +285, 5.5" star +87, 3.1" star +36 (Delta chi2 vs M3V); 95% region within 4.6" (contains the 3.1" star, which
  colours and shape exclude). An RRab 200" away leaked into one pipeline; removed at source, global transit unchanged.
- Raw pixels confirm the dip at the known ephemeris but are not an independent blind detection.
- Depth 15.8 +- 0.9 (stat) +1.5/-1.1 (sys) ppt of 979 e/s; Rp 4.95 +- 0.24 (stat) +0.25/-0.18 (sys) Re; b ~0.5;
  T14 1.70 h, T23/T14 0.71-0.74 (density prior; 1.65-1.75 h without); Teq ~470 K; a 0.0448 au.
- Star: 64.31 +- 0.15 pc; R 0.381-0.385, M 0.370, Teff 3350 K; 0.75 mag over-luminous, explained by [Fe/H] ~+0.2
  (colours, XP); P(twin) ~0.3-1%.
- Ephemeris (conservative): T0 = BJD_TDB 2460840.9150 +- 0.0013, P = 5.691506 +- 0.000011 d; 1 sigma ~3 min in 2027.
  The old predicted_transits_2027.csv had UTC without barycentric correction (-3.4 to +7.3 min); superseded by
  tic214311908/predicted_transits_2027.csv, which now also replaces the old file (uncorrected version in git history).
- P/2 excluded at 14 sigma; 2P main-sequence EB excluded at >= 6 sigma in every cleaned photometry; no secondary
  (3 sigma < 2.8-3.1 ppt); no TTVs.
- Both FPP codes had a photo-eccentric sign error (rho_c = rho g^-3 instead of g^+3); fixed. FPP 0.0013 (REF) to
  0.0033 (INV); 64 variants 0.0007-0.0037; 1,500 joint prior draws: median 0.0022, 95th percentile 0.012,
  pessimistic combinations up to 0.045. NFPP ~4e-4 (most lenient pixel likelihood). Earlier 0.011 and 0.003-0.007
  superseded. With planets on unresolved companions counted as false positives (as TRICERATOPS does): 0.005-0.009.
- After speckle + a two-band ground transit the FPP would be ~3-12e-5 (adopted brown-dwarf inputs); floor 3e-5 from
  the carried-over background term.
- Not statistically validated: no ground-based detection, no high-resolution imaging, TRICERATOPS unstable
  (0.002-0.04). First full windows April-August 2027.
- Side results: Mira period 313-320 d (not 294 d); Mira TESS image displaced 5-7.5" by colour (not 16-17");
  the eRASS:3 X-ray source may be the M3V or pair A (p 0.54 vs 0.38).
```
