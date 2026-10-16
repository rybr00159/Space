# Vetting the five recurring signals (2026-10-16)

*AI-run analysis; no human has reviewed it yet. Nothing has been submitted anywhere.*

On 2026-10-16 we found that our recurrence test had been too coarse for stars with few TESS sectors (see the
[search README](../README.md#caveats)). With a finer test, five more Sector 104 signals recur. This page vets them as
we vetted the earlier candidates:
- the stars themselves;
- the light curves in every sector and product;
- which star in the image is dimming;
- nearby stars;
- earlier catalogue entries.

It also reflects two external reviews (GPT-5.5 and Gemini 3.1 Pro; [`../reviews/`](../reviews/), round 6).

**None of the five is a planet orbiting the star we searched. One signal, however, found on the Mira-type giant
TIC 214311906, most likely comes from a neighbouring M dwarf 6″ away. There it would be a Neptune-sized transiting
planet candidate.** It is unconfirmed, and identifying its host needs one ground-based observation.

| TIC | Signal as first found | Verdict | Decisive evidence | New? |
|---|---|---|---|---|
| 214311906 | 5.69 d, 1,900 ppm | **Possible planet candidate on a neighbour**: a ~5–7 R⊕ body transiting an M dwarf in the 61 pc pair TIC 214311910, 6″ from the Mira | The Mira cannot host it. Calibrated pixel analysis puts the source at the Mira/pair position, not the other bright stars. The flat-bottomed 1.6 h transit needs a small occulter and a dense host, which only the M dwarfs provide | Yes: no TOI, CTOI, alert or ExoFOP record for any of the three stars |
| 90086266 | 2.36 d, 0.7% box | **Eclipsing binary** coincident with an F star; true period 2.33947 d (or 4.679 d) | The dimming is centred on the F star. A 3% eclipse of it would need a 3 R_Jup companion, i.e. a star, which Gaia RVs and the missing secondary eclipse rule out. So a fainter, unresolved binary is eclipsing | Yes, as an eclipsing binary |
| 61332759 | 0.740 d, 550 ppm | **Not a transit**: 1.4806 d variability of a neighbouring star | Every sector traces the variability to the catalogued RS CVn variable Gaia DR3 6711181657530577280, 57″ away | No period found for the neighbour in VSX, ASAS-SN or Gaia DR3 |
| 61141412 | 6.525 d, 1,100 ppm | **Not a planet**: dips every 6.5285 d in or near a 5.45 d spectroscopic binary, with timing irregularities | The star is a Gaia SB1, so a 6.53 d planet orbit is impossible. Dip times do not follow a constant period | No: a TESS-SPOC alert at 13.054 d (2×) in 2023 |
| 55829034 | 0.616 d, 1,000 ppm | **Blended eclipsing (or transiting) source** in the bulge | The target is a red giant; a 0.6 d orbit would lie inside it. The source is unidentified in a crowded field | Signal not catalogued |

## TIC 214311906: a possible Neptune-sized planet on a neighbouring M dwarf

![The two giant-star signals](figures/giants.png)

*Left: the TIC 214311906 dips, folded on 5.6915 d; flat-bottomed, 1.6 h. Right: TIC 55829034 (next section).*

**The stars.**
- **The target, TIC 214311906**, is a Mira-type long-period variable: Gaia LPV, VSX and ASAS-SN P ≈ 294 d, BP−RP 4.7,
  parallax 0.15 ± 0.09 mas. It is a distant giant hundreds of R☉ across, and cannot host 1.6 h dips.
- **TIC 214311910**, 6.2″ west of it, is a nearby M dwarf: Gaia DR3 5925373272850782208, G 14.65, 61 pc. Gaia sees a
  double image in 20% of scans and RUWE 13.7, so it is a close pair. The TIC has no stellar parameters for it.
- **TIC 214311908**, 17.8″ west, is an M3V dwarf: G 14.28, 64.5 pc, 0.37 M☉, 0.38 R☉. It shares the pair's proper
  motion, so the three M dwarfs form one bound system about 62 pc away.

**The signal.**
- Dips every **P = 5.691471 ± 0.000018 d** (T0 = BTJD 3840.9101 ± 0.0022), from 14 individual events in Sectors 39,
  66, 93, 103 and 104. It is also seen weakly in Sector 12.
- Flat-bottomed, T14 = 1.61 ± 0.09 h, with the flat part 78 ± 9% of the transit. Odd and even depths agree, and there
  is no secondary eclipse.
- About 2,000–2,400 ppm deep in the Mira's light curve, varying between sectors as the Mira's own brightness changes
  the dilution.
- Our search also detected it on TIC 214311908 in Sector 104. That detection was dropped because, after masking
  instrumental time windows, one of its three transits carried 58% of the S/N (the limit was 50%).

**Where the dimming is.**
- Free fits of the source position fail in this field: fake signals injected at known stars come back scattered by
  1–2′.
- A joint χ² map over Sectors 39, 66, 93 and 104 does work, once calibrated with injections:

  | Signal | Map minimum | Δχ² at the Mira / pair / M3V / bright star 41″ away |
  |---|---|---|
  | Real | (+3, −3)″ | 0 / 15 / 57 / 122 |
  | Injected on the pair (6.2″) | (0, 0)″ | 0 / 17 / 65 / 138 |
  | Injected on the M3V (17.8″) | (−21, +3)″ | 52 / 38 / 0 / 154 |
  | Injected on the bright star (41″, T 11.1, RUWE 15) | (+36, −24)″ | 284 / 297 / 300 / 0 |

- The real signal matches a source at the Mira/pair position, which TESS cannot separate. It does not match the M3V
  or the bright star.

**Why the host is probably the M-dwarf pair.**
1. The Mira cannot host the dips: a 5.7 d orbit would lie deep inside it.
2. The flat bottom limits the occulter to at most ~0.2 of the host's radius. So the host must supply at least ~7% of
   the light in the aperture (T ≲ 14). The only such stars at the source position are the Mira and the pair.
3. If the transit is central, the shape implies a host density of ~12–15 ρ☉: an M dwarf. A Sun-like background star
   would need a grazing, V-shaped transit, which is not what is seen.
4. A background M dwarf bright enough (T ≲ 14) would have to be nearby and within ~0.5″ of the pair. That is very
   unlikely.

**What it would be.**
- On the pair, the dips are 1.5–3% deep in the pair's light, depending on the Mira's brightness.
- If the pair's light is dominated by one star (R ≈ 0.35 R☉), the transiting body is 4.5–6.5 R⊕. If the pair is two
  similar stars (R ≈ 0.27 R☉ each), it is 5–7 R⊕.
- Either way it is smaller than any star or old brown dwarf: **a Neptune- or sub-Saturn-sized planet**.
- It would be at about 0.04 au from its star, with an equilibrium temperature near 450 K.

**What is not yet excluded.**
- An eclipsing binary hidden inside the pair. This needs a stellar occulter much smaller than the host, which the size
  limit rules out, or an unlikely chance alignment.
- A brown dwarf with a grazing path. The flat bottom disfavours this at about 2σ.
- A host other than the pair if the TESS calibration is wrong. The map uses signals of only S/N 6–7 per sector and
  one injection per hypothesis.

**The decisive test.** One transit observed from the ground in seeing-limited images separates the Mira, the pair
(6.2″) and the M3V (17.8″). On the right star the dip should be 1.5–3% deep and last 1.6 h.
[`results/214311906/predicted_dips_2027.csv`](results/214311906/predicted_dips_2027.csv) lists 38 predicted transits
for March–September 2027, each with ±5–6 min (1σ) timing. The field is at Dec −51°.

**Novelty.**
- None of the three stars has a TOI, CTOI, SPOC or TESS-SPOC alert, QLP-alert (Tey+23) entry, VSX period or any ExoFOP
  record (checked 2026-10-12 and 2026-10-16).
- In the TOI process, this would be a planet candidate on TIC 214311910 with an uncertain host.

## TIC 90086266: a 3% eclipsing binary coincident with an F star

![TIC 90086266](figures/tic90086266.png)

*Left: the eclipse in 2019 and 2026. Middle: depths of successive eclipses; alternate eclipses are equal in 2019 but
differ by 7.5% in 2026. Right: signal S/N per pixel in Sector 104, centred on the target (blue cross); circles are Gaia
stars.*

**The star.** Gaia DR3 6722830841491505536; T = 11.10; F-type subgiant, Teff 6,160–6,370 K, R 1.63–1.71 R☉,
M 1.2–1.3 M☉, 413–423 pc. Gaia RUWE 0.96 and no double-peaked scans. Gaia RV −7.1 ± 1.5 km/s from 5 transits, with a
peak-to-peak range of only 5.0 km/s (constancy p = 0.22; RV S/N 23 per transit).

**The signal.**
- Eclipses 2.7–3.0% deep in raw (SAP) flux and 3.2% in crowding-corrected PDCSAP flux. Duration 3.1 h, nearly
  V-shaped (flat bottom 12 ± 1% of the eclipse).
- **The period is 2.339469 d.** Sector 13 gives 2.33940 ± 0.00006 d and Sector 104 2.33936 ± 0.00005 d. Only one
  cycle count (1,078) between them agrees with both. Our search had reported 2.3606 d and 2.35477 d, both artefacts of
  a too-wide box on detrended QLP flux.
- No secondary eclipse at phase 0.5: −60 ± 45 ppm.

**Where the dimming is.** Matched-filter TESS PRF fits put the source (+0.2, −0.1)″ ± 0.2″ from the target in Sector
13 and (−1.2, −0.2)″ ± 0.2″ in Sector 104. Fake eclipses injected on the target give (+0.2, −0.2)″ and (−1.3, −0.2)″.
So the source is coincident with the target to within TESS's systematic accuracy (of order an arcsecond). The four
Gaia stars within 6″ (T = 15.3–18.7) each supply at most 2% of the light, too little for a 3% dip.

**Why it is not a planet, and not an eclipse of the F star itself.**
1. A 3% eclipse of a 1.7 R☉ star needs a companion of about 0.3 R☉ (3 R_Jup). No planet or brown dwarf is that
   large, so the companion would be a star.
2. A ~0.3 M☉ star in a 2.34 d orbit would swing the F star's velocity by about ±30–35 km/s. Gaia's 5 RV transits vary
   by only 5 km/s.
3. Such a star would also show a secondary eclipse. For an M dwarf of 3,300–3,500 K, the expected depth is roughly
   1,000–3,000 ppm, depending on its temperature. None is seen (−60 ± 45 ppm).
4. A twin F-star binary at 4.68 d fits the brightness, but would give double lines moving by ±86 km/s. Gaia's RVs rule
   that out too, unless all its transits fell near conjunction (probability of order 10⁻³ or less).

So the eclipses most likely come from a **fainter eclipsing binary unresolved from the F star**, either bound to it
(a hierarchical triple) or behind it. It must supply at least 3.2% of the light (ΔT ≲ 3.7 mag), or at least ~6% if
its eclipses are 50% deep. "Within about 1″" is an inference: such a star farther out would be bright enough to
appear in Gaia as a separate source. The TESS pixels alone only show coincidence to about an arcsecond.

**Open issue: alternating depths.**
- In Sector 104, odd and even eclipses differ by 7.5% (5.6–8.9σ). The difference holds in QLP flux and in our own
  TESScut photometry with three aperture sizes.
- In Sector 13 they agree to within 0.5 ± 1.5%.
- It is not an aperture or dilution effect: those would change all eclipses in a sector equally, not alternate ones.
- A binary with two unequal eclipses at 2 × 2.339 = 4.679 d would show the difference in both sectors. So either the
  eclipses changed between 2019 and 2026, or the 2026 difference has a cause we have not found.

**What would settle it.**
- High-resolution imaging or spectroscopy, to find the companion pair or its moving lines.
- Two-colour eclipse photometry, to measure the eclipsed star's colour.

**Novelty.**
- Not a TOI, CTOI or SPOC/TESS-SPOC alert. No variability entry in VSX, ASAS-SN or Gaia DR3.
- Green+23 gives this star a Score of 0 (not flagged). Its listed best period, 0.2355 d, is close to one-tenth of the
  eclipse period, so their periodogram probably picked up a harmonic of the eclipses. They did not classify the star
  as variable.
- TARS has no period. This is a new identification as an eclipsing binary.

## TIC 61332759: variability of a neighbour 57″ away

![TIC 61332759](figures/tic61332759.png)

*Left: signal S/N per pixel in Sector 67; the source (orange circle) is the RS CVn variable, not the target (blue
cross). Right: the variability in the QLP aperture folded on 1.4806 d.*

- **Not a transit.** The "dip" is a smooth variation covering about half of each 0.74 d cycle (trapezoid fits give an
  8 h duration). The true period is **1.4806 d**, with two unequal minima and unequal maxima per cycle. That is typical
  of a spotted (RS CVn-type) or ellipsoidal binary.
- **It comes from a neighbour.** The amplitude depends on the aperture. In Sectors 13 and 27, QLP's raw flux shows 440
  and 1,340 ppm, but TESS-SPOC's smaller apertures show only 50–80 ppm. In Sector 67, QLP shows 2,700 ppm.
- **Every sector traces it to the same place.** Pixel-level fits in Sectors 13, 27, 67, 94 and 104 put the source at
  (−36.5, −44.3)″ ± 0.3″ from the target. That is the catalogued variable Gaia DR3 6711181657530577280 (G 14.4, Gaia
  class "RS"). The target is excluded by Δχ² ≈ 195,000 and the next-best star by Δχ² ≈ 45,000.
- The target itself is an F star with constant Gaia RV (7 transits).
- VSX and Gaia DR3 list the neighbour as an RS CVn variable but give no period, and ASAS-SN has no entry. We found no
  catalogued 1.48 d period.

## TIC 61141412: dips with irregular timing near a spectroscopic binary

![TIC 61141412](figures/tic61141412.png)

*Left: dip times minus a constant-period prediction, 2019–2026, with a quadratic fit. Right: the dips in Sectors 27
and 67, folded on that fit.*

**The star.** A subgiant: T = 10.0, Teff ≈ 5,900–6,000 K, R 2.4–2.7 R☉, 354–384 pc. It is a **Gaia DR3 spectroscopic
binary** (SB1): P = 5.452 d, e = 0.01, significance 51; its Gaia RVs span 35 km/s. TESS shows rotation at 5.4–5.8 d
(TARS, Gao+25, Green+23), which matches tidally locked rotation in that orbit.

**The signal.**
- About 1,000–2,000 ppm, roughly 1 h long, every **6.5285 d**. It is seen in Sectors 13, 27, 67, 94 and 104, with 11
  single events at S/N > 4.
- Our Sector 104 detection (6.5252 d) was the same signal.
- Depths agree between products (QLP, SPOC, TESS-SPOC) within each sector.

**Timing.** The 11 event times do not fit a constant period with their formal 3–5 min errors (χ² = 243 for 9 degrees
of freedom; residuals from −25 to +29 min). A steadily changing period fits better (χ² = 22 for 8 degrees of
freedom). Even with all errors inflated by 1.66 so that the quadratic fit has χ² = 8, the constant period still gives
χ² ≈ 88. A third body pulling the eclipsing pair is one explanation. But the times come from three pipelines, and the
dips are shallow, short and on a spotted star, so the errors may be underestimated. **We treat the timing variations
as suggestive, not established.**

**Where the dimming is.** Our joint pixel map puts it on the target. Every Gaia star 17″ or more away fits worse than
the target by Δχ² ≥ 88. The TESS-SPOC pipeline's own Sector 67 analysis gives an offset of 3.3 ± 2.9″.

**Why it is not a planet around the target.**
- A 6.53 d orbit around one star of a 5.45 d binary is dynamically unstable.
- A planet orbiting both stars would need a period of at least ~25 d.
- The ~1 h duration would also need an almost grazing path across a 2.5 R☉ star.

The likeliest picture is an eclipsing pair within a few arcseconds: either in the same multiple system, or an
unresolved background binary. If Gaia's 5.45 d period were an alias of 6.53 d, the dips would instead be grazing
eclipses within the SB1 itself.

**Not new.** The TESS-SPOC pipeline flagged this signal in Sector 67 (2023), as a threshold-crossing event at
13.054 d (twice the true period). Its report is public. The 2019 SPOC alerts at 16.6 d and 20.6 d each combined one
of these dips (BTJD 1661.01) with unrelated dips at 1677.42 and 1677.64. SPOC's own centroids place those unrelated
dips 18–62″ away.

## TIC 55829034: a blended source in the bulge

- **The target** is a red giant towards the Galactic bulge: semi-regular variable, P = 42 d; parallax 0.07 ± 0.03 mas;
  TIC Teff 3,260 K. A 0.6 d orbit would lie deep inside it.
- **The signal** is a V-shaped, ~2 h, ~1,000–2,300 ppm dip every 0.6161 d, with equal odd and even depths and no
  secondary eclipse.
- **The source** is within several arcseconds of the target, in a field with about 6,000 Gaia stars within 2.5′. TESS
  cannot identify it.
- The V shape and short period make an eclipsing binary likely. A transiting planet around a faint neighbour is not
  formally excluded.

## How it was done

| Step | What | Script (in `scripts/`) |
|---|---|---|
| Stars | TIC and Gaia DR3 parameters, RV statistics, NSS orbits, variability classes | `stars.py` |
| Earlier catalogue entries | TOI/CTOI lists, SPOC and TESS-SPOC pipeline reports, every VizieR table within 30″, VSX, ASAS-SN, Gaia variables and EBs within 4–5′, TARS, ExoFOP | `products.py`, `dv.py`, `allviz5.py`, `nbr.py`, `nbr2.py` |
| Light curves | Every product and aperture (QLP SAP/DET, SPOC and TESS-SPOC SAP/PDCSAP), transit-masked detrending, ephemeris, depths, odd/even, secondary scan, phase curve, shape | `lcan.py` |
| TIC 214311906 field | Depths in the three stars' light curves; per-event times and 2027 predictions; host hypotheses with injections; joint-map calibration | `trio214.py`, `eph214.py`, `hyp214.py`, `jointinj214.py` |
| TIC 90086266 | Per-sector timing and cycle count; alternate-eclipse depths in pipeline flux and in our own TESScut photometry | `eph90.py`, `oddeven90.py`, `ap90.py` |
| TIC 61332759 | Period and shape of the variability | `var613.py` |
| TIC 61141412 | Single-event search; per-event times | `events614.py`, `etv614.py` |
| Localisation | TESScut cutouts (21×21 px); per-pixel regression on a signal template (eclipse shape or free-phase harmonics) with spline baselines and jitter regressors; TESS PRF fits to the amplitude images, with the WCS calibrated on the target's direct image; injected fake signals; joint χ² maps across sectors evaluated at every Gaia star | `dl_cutouts.py`, `loc5.py`, `locmulti.py` |
| Figures | | `figs.py` |

The scripts import the search pipeline from `/home/user/work/src`, which is identical to [`../src/`](../src/), and
read cutouts and light curves from local caches. Results are in [`results/`](results/), one folder per star.
- `results/61141412/out_61141412.txt` used the wrong period (6.5458 d) and is superseded by
  `out_61141412_fixed.txt`.
- `mira214.py` tested whether the Mira's brightness changes dilute the eclipse. It gave depths inconsistent with QLP
  and is not used.
- The free-fit positions for TIC 214311906 (`loc/s*.json`, `loc_214.txt`) are unreliable, as the injections show.

## Caveats

- **TIC 214311906/214311910 is a candidate, not a planet.**
  - The host identification rests on one injection-calibrated TESS map plus physical arguments.
  - The planet's size depends on the Mira's brightness at each epoch and on how light is shared within the M-dwarf
    pair.
  - Statistical validation would need the host confirmed first.
- **TIC 90086266's companion binary is inferred, not seen.** Its separation, brightness and period (2.34 or 4.68 d)
  are unknown.
- **TIC 55829034:** the conclusion that the eclipsed star is not the giant rests on physics; the pixel data cannot
  resolve the field.
- **TIC 61141412:** the timing variations need a homogeneous re-measurement before they can be called real.
- An independent referee re-check is in progress and will be added here.
