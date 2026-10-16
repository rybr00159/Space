# Review request (round 6): vetting verdicts for five recurring TESS signals

Context: an AI-run search of TESS Sectors 104/105 for new planets. A coarse recurrence test (150 random periods per
sector, so p could never go below 1/151 per sector) had dropped five signals; with a finer test they recur, and they
were vetted on 2026-10-16. The agent concluded that **none is a planet** and that one (TIC 90086266) is a **newly
identified eclipsing binary** blended within ~1" of an F star. Please try to find errors in each verdict: wrong
physics, wrong arithmetic, a plausible planet scenario that was dismissed too quickly, or an overclaim (in
particular the "new eclipsing binary within ~1 arcsec" claim for TIC 90086266 and the eclipse-timing-variation claim
for TIC 61141412). Four figures are attached (TIC 90086266; TIC 61332759; TIC 61141412; the two giant-star signals).

The full write-up follows.

# Vetting the five recurring signals (2026-10-16)

*AI-run analysis; no human has reviewed it yet. Nothing has been submitted anywhere.*

On 2026-10-16 we found that our recurrence test had been too coarse for stars with few TESS sectors (see the
[search README](../README.md#caveats)). With a finer test, five more Sector 104 signals recur. This page vets them in
the same way as the earlier candidates: the stars themselves, the light curves in every sector and product, which star
in the image is dimming, nearby stars, and earlier catalogue entries.

**None of the five is a planet.**

| TIC | Signal as first found | Verdict | Decisive evidence | New? |
|---|---|---|---|---|
| 90086266 | 2.36 d, 0.7% box | **Eclipsing binary** blended within ~1″ of an F star; true period 2.33947 d (or 4.679 d) | 3% eclipses would need a 3 R_Jup companion, i.e. a star. Gaia RVs and the missing secondary eclipse rule that out for the F star itself, yet the dimming is centred on it | Yes: no earlier catalogue entry or pipeline alert |
| 61332759 | 0.740 d, 550 ppm | **Not a transit**: 1.4806 d variability of a neighbouring star | Every sector traces the variability to the catalogued RS CVn variable Gaia DR3 6711181657530577280, 57″ away | The neighbour's 1.48 d period is not in VSX or Gaia |
| 61141412 | 6.525 d, 1,100 ppm | **Eclipsing pair with eclipse-timing variations** in or near a spectroscopic binary | The star is a Gaia SB1 with a 5.45 d orbit, so a 6.53 d planet orbit is impossible; dip times drift by ±25 min | No: a TESS-SPOC alert at 13.054 d (2×) in 2023 |
| 214311906 | 5.69 d, 1,900 ppm | **Blended eclipsing binary** | The target is a Mira-type giant, far too big for 1.6 h eclipses; the source is within ~10″ but unidentified | Signal not catalogued |
| 55829034 | 0.616 d, 1,000 ppm | **Blended eclipsing binary** | The target is a bulge red giant; a 0.6 d orbit would lie inside it; the source is within ~5–8″ in a very crowded field | Signal not catalogued |

## TIC 90086266: a 3% eclipsing binary hidden by an F star

![TIC 90086266](figures/tic90086266.png)

*Left: the eclipse in 2019 and 2026. Middle: depths of successive eclipses; alternate eclipses are equal in 2019 but
differ by 7.5% in 2026. Right: signal S/N per pixel in Sector 104, centred on the target (blue cross); circles are Gaia
stars.*

**The star.** Gaia DR3 6722830841491505536; T = 11.10; F-type subgiant, Teff 6,160–6,370 K, R 1.63–1.71 R☉,
M 1.2–1.3 M☉, 413–423 pc. Gaia RUWE 0.96 and no double-peaked scans. Gaia RV −7.1 ± 1.5 km/s from 5 transits, with a
peak-to-peak range of only 5.0 km/s (constancy p = 0.22).

**The signal.**
- Eclipses 2.7–3.0% deep in raw (SAP) flux and 3.2% in crowding-corrected PDCSAP flux. Duration 3.1 h, nearly
  V-shaped (flat bottom 12 ± 1% of the eclipse).
- **The period is 2.339469 d.** Sector 13 gives 2.33940 ± 0.00006 d and Sector 104 2.33936 ± 0.00005 d. Only one
  cycle count (1,078) between them agrees with both. Our search had reported 2.3606 d and 2.35477 d, both artefacts of
  a too-wide box on detrended QLP flux.
- No secondary eclipse at phase 0.5: −60 ± 45 ppm.

**Where the dimming is.** Matched-filter TESS PRF fits put the source (+0.2, −0.1)″ ± 0.2″ from the target in Sector
13 and (−1.2, −0.2)″ ± 0.2″ in Sector 104. A fake eclipse injected on the target gives (+0.2, −0.2)″ and (−1.3,
−0.2)″, so both are fully consistent with the target. The four Gaia stars within 6″ (T = 15.3–18.7) each supply at
most 2% of the light, too little to cause a 3% dip.

**Why it is not a planet, and not an eclipse of the F star itself.**
1. A 3% eclipse of a 1.7 R☉ star needs a companion of about 0.3 R☉ (3 R_Jup). No planet or brown dwarf is that
   large, so the companion would be a star.
2. A ~0.3 M☉ star in a 2.34 d orbit would swing the F star's velocity by about ±30 km/s. Gaia's 5 RV transits vary by
   only 5 km/s.
3. Such a star would also show a secondary eclipse of about 2,000 ppm. None is seen (−60 ± 45 ppm).
4. A twin F-star binary at 4.68 d fits the brightness, but would give double lines moving by ±86 km/s. Gaia's RVs rule
   that out too, unless all its transits fell near conjunction, which has a probability of about 10⁻³ or less.

So the eclipses most likely come from a **fainter eclipsing binary within about 1″ of the F star**, either bound to it
(a hierarchical triple) or behind it. It must supply at least ~6% of the light, i.e. be no more than ~3 mag fainter
in the TESS band.

**Open issue: alternating depths.** In Sector 104, odd and even eclipses differ by 7.5% (5.6–8.9σ). The difference
holds in QLP flux and in our own TESScut photometry with three aperture sizes. In Sector 13 they agree to within
0.5 ± 1.5%. A binary with two unequal eclipses at 2 × 2.339 = 4.679 d would show it in both sectors. So either the
eclipses changed between 2019 and 2026 (for example through spots, or orbital changes driven by the F star), or the
2026 difference has a cause we have not found.

**What would settle it.** High-resolution imaging or spectroscopy would reveal the companion pair or its moving
spectral lines. Two-colour eclipse photometry would show the colour of the eclipsed star.

**Novelty.** Not a TOI, CTOI or SPOC/TESS-SPOC alert. No variability entry in VSX, ASAS-SN or Gaia DR3. Green+23 lists
0.2355 d (no match), and TARS has no period. This is a new identification as an eclipsing binary.

## TIC 61332759: variability of a neighbour 57″ away

![TIC 61332759](figures/tic61332759.png)

*Left: signal S/N per pixel in Sector 67; the source (orange circle) is the RS CVn variable, not the target (blue
cross). Right: the variability in the QLP aperture folded on 1.4806 d.*

- **Not a transit.** The "dip" is a smooth variation covering about half of each 0.74 d cycle (trapezoid fits give an
  8 h duration). The true period is **1.4806 d**, with two unequal minima and unequal maxima per cycle. That is typical
  of a spotted (RS CVn-type) or ellipsoidal binary.
- **It comes from a neighbour.** The amplitude depends on the aperture. In Sectors 13 and 27, QLP's raw flux shows
  440 and 1,340 ppm, but TESS-SPOC's smaller apertures show only 50–80 ppm. In Sector 67, QLP shows 2,700 ppm.
- **Every sector traces it to the same place.** Pixel-level fits in Sectors 13, 27, 67, 94 and 104 put the source at
  (−36.5, −44.3)″ ± 0.3″ from the target. That is the catalogued variable Gaia DR3 6711181657530577280 (G 14.4,
  Gaia class "RS"). The target is excluded by Δχ² ≈ 195,000 and the next-best star by Δχ² ≈ 45,000.
- The target itself is an F star with constant Gaia RV (7 transits).
- VSX and Gaia DR3 list the neighbour as an RS CVn variable but give no period, so the 1.4806 d period appears to be
  new.

## TIC 61141412: timed eclipses in a spectroscopic binary

![TIC 61141412](figures/tic61141412.png)

*Left: dip times minus a constant-period prediction, 2019–2026, with a quadratic fit. Right: the dips in Sectors 27
and 67, folded on that fit.*

**The star.** A subgiant: T = 10.0, Teff ≈ 5,900–6,000 K, R 2.4–2.7 R☉, 354–384 pc. It is a **Gaia DR3 spectroscopic
binary** (SB1): P = 5.452 d, e = 0.01, significance 51; its Gaia RVs span 35 km/s. TESS shows rotation at 5.4–5.8 d,
which matches tidally locked rotation in that orbit.

**The signal.**
- About 1,000–2,000 ppm, roughly 1 h long, every **6.5285 d**. It is seen in Sectors 13, 27, 67, 94 and 104, with 11
  single events at S/N > 4.
- Our Sector 104 detection (6.5252 d) was the same signal.
- Depths agree between products (QLP, SPOC, TESS-SPOC) within each sector.

**Eclipse-timing variations.**
- A constant period does not fit the 11 event times: χ² = 243 for 9 degrees of freedom. The residuals form an arch
  from −25 to +29 min.
- A period that changes steadily fits much better: χ² = 22 for 8 degrees of freedom.
- The dipping pair is therefore being pulled or moved by a third body.

**Where the dimming is.** Our joint pixel map puts it on the target. Every Gaia star 17″ or more away fits
worse than the target by Δχ² ≥ 88. The TESS-SPOC pipeline's own Sector 67 analysis gives an offset of 3.3 ± 2.9″.

**Why it is not a planet.**
- A 6.53 d orbit around one star of a 5.45 d binary is dynamically unstable.
- A planet orbiting both stars would need a period of at least ~25 d.
- The ~1 h duration would also need an almost grazing path across a 2.5 R☉ star.

The likeliest picture is an eclipsing pair in the same multiple system (for example a 2+2 quadruple, with the timing
drift caused by the two pairs' mutual orbit), or an unresolved background binary with its own companion. If Gaia's
5.45 d period were an alias of 6.53 d, the dips would instead be grazing eclipses within the SB1 itself. That would
still need a third body to explain the timing drift.

**Not new.** The TESS-SPOC pipeline flagged this signal in Sector 67 (2023), as a threshold-crossing event at
13.054 d (twice the true period). Its report is public. The 2019 SPOC alerts at 16.6 d and 20.6 d each combined one
of these dips (BTJD 1661.01) with unrelated dips at 1677.42 and 1677.64. SPOC's own centroids place those unrelated
dips 18–62″ away.

## TIC 214311906 and TIC 55829034: eclipses blended with giant stars

![The two giant-star signals](figures/giants.png)

**TIC 214311906.**
- The target is a Mira-type long-period variable: Gaia LPV; VSX and ASAS-SN P ≈ 294 d; BP−RP 4.7; parallax
  0.15 ± 0.09 mas. It is a distant giant hundreds of R☉ across.
- The signal is a flat-bottomed, 1.6 h, ~2,000–2,800 ppm dip every 5.6915 d, with equal odd and even depths and no
  secondary eclipse.
- A star in a 5.7 d orbit around this giant would be inside it, so the eclipsed star is something else.
- The pixel data put the source within ~10″ of the Mira but cannot single out a star.
- One candidate is a nearby M-dwarf pair 6.4–6.6″ away (Gaia DR3 5925373272850782208; 61 pc; RUWE 13.7, which marks
  it as a binary; T ≈ 13.5). It would need ~2.3% eclipses.
- One eclipse observed with a small telescope could check this, since 6.6″ is easily resolved from the ground.

**TIC 55829034.**
- The target is a red giant towards the Galactic bulge: semi-regular variable, P = 42 d; parallax 0.07 ± 0.03 mas;
  TIC Teff 3,260 K.
- The signal is a V-shaped, ~2 h, ~1,000–2,300 ppm dip every 0.6161 d.
- A 0.6 d orbit would lie deep inside the giant.
- The source is within ~5–8″ of the target, in a field with about 6,000 Gaia stars within 2.5′. It cannot be
  identified with TESS.

## How it was done

| Step | What | Script (in `scripts/`) |
|---|---|---|
| Stars | TIC and Gaia DR3 parameters, RV statistics, NSS orbits, variability classes | `stars.py` |
| Earlier catalogue entries | TOI/CTOI lists, SPOC and TESS-SPOC pipeline reports, every VizieR table within 30″, VSX, ASAS-SN, Gaia variables and EBs within 4–5′, TARS | `products.py`, `dv.py`, `allviz5.py`, `nbr.py`, `nbr2.py` |
| Light curves | Every product and aperture (QLP SAP/DET, SPOC and TESS-SPOC SAP/PDCSAP), transit-masked detrending, ephemeris, depths, odd/even, secondary scan, phase curve, shape | `lcan.py` |
| TIC 90086266 | Per-sector timing and cycle count; alternate-eclipse depths in pipeline flux and in our own TESScut photometry | `eph90.py`, `oddeven90.py`, `ap90.py` |
| TIC 61332759 | Period and shape of the variability | `var613.py` |
| TIC 61141412 | Single-event search; per-event times and timing variations | `events614.py`, `etv614.py` |
| Localisation | TESScut cutouts (21×21 px); per-pixel regression on a signal template (eclipse shape or free-phase harmonics) with spline baselines and jitter regressors; TESS PRF fits to the amplitude images, with the WCS calibrated on the target's direct image; injected fake signals; joint χ² maps across sectors evaluated at every Gaia star | `dl_cutouts.py`, `loc5.py`, `locmulti.py` |
| Figures | | `figs.py` |

The scripts import the search pipeline from `/home/user/work/src`, which is identical to [`../src/`](../src/), and
read cutouts and light curves from local caches. Results are in [`results/`](results/), one folder per star.
`results/61141412/out_61141412.txt` used the wrong period (6.5458 d) and is superseded by `out_61141412_fixed.txt`.
`mira214.py` (a test of whether the Mira's brightness changes dilute the eclipse) gave depths inconsistent with QLP
and is not used.

## Caveats

- **TIC 90086266's companion binary is inferred, not seen.** Its existence rests on the arguments above. Its
  separation, brightness and period (2.34 or 4.68 d) are unknown.
- **TIC 214311906 and TIC 55829034:** the conclusion that the eclipsed star is not the giant rests on physics, not on
  the pixel data, which cannot resolve the field within ~10″.
- **TIC 61141412:** the timing variations rest on 11 events with 3–5 min uncertainties from three different
  pipelines. A quadratic fit is only a description; the true cause could be periodic.
- An independent referee re-check is in progress and will be added here.
