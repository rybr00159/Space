# Would our TESS search have found TIC 4206066?

*Checked 2026-10-16. The unchanged search code in [`../../tess-search/src/`](../../tess-search/src/) was run on this
star's TESS light curves. AI-run; no human has reviewed it yet.*

## Answer

- **Our search never saw this star.**
  - TESS has observed TIC 4206066 (RA 05h23m26s, Dec −01°19′38″) only in Sectors 6 (Dec 2018–Jan 2019), 32
    (Nov–Dec 2020) and 98 (Nov 2025–Jan 2026). Sector 110 starts on 31 October 2026.
  - Our search covered Sectors 104 and 105 (Apr–Jun 2026), which pointed at a different part of the sky. The star is
    on neither sector's target list.
- **Given the same data, our search would have found the 3.18-day signal, then dropped it.**
  - On the Sector 98 light curve, the search detects the signal, finds it in no catalogue and classes it as
    planet-like.
  - The recurrence test in Sectors 6 and 32 then gives p ≈ 0.004. Our cut for vetting is p < 0.001, so the signal
    would never have been vetted.
- **The drop is a flaw in our test, not a weakness of the signal.** The test compares each sector with only 150
  random periods, so no sector can score better than 1 in 151. With 3,000 random periods, the same data give
  p = 0.0003, which passes.
- **The search would not have reported the 11.13-day signal.** It keeps one signal per star and sector, and needs at
  least three transits in one sector.

## Step by step

| Step | Our cut | TIC 4206066 | Outcome |
|---|---|---|---|
| Search of the Sector 98 QLP light curve (the product type we used) | Red-noise S/N ≥ 9, ≥ 3 transits, no dominant single event | P = 3.1811 d, 477 ppm, 9 transits, S/N 11.8 | Detected |
| Known-signal check | Not in our TOI, CTOI, confirmed-planet or 341,351-entry SPOC TCE lists | Not listed | New |
| Classifier | Odd/even depths, secondary eclipse, phase curve, duration | Odd/even 0.3σ; phase-0.5 depth 0.8σ; 1.5 R⊕; duration 0.74 of the central-transit maximum | Planet-like |
| Recurrence in Sectors 6 and 32 (discovery sector excluded) | Combined p < 0.001 | S6 false-alarm probability 0.007, which is the 1/151 floor; S32 0.06. Combined p = 0.0035–0.0042 (two random seeds) | **Dropped** |
| The same test with 3,000 random periods (added here) | Combined p < 0.001 | S6 0.0003 (again the floor); S32 0.065. Combined p = 0.00025 | Would pass |

The joint fit across all three sectors gives P = 3.182789 d, against the published 3.182784 ± 0.000006 d. Depths at
that ephemeris are 444 ± 53 ppm (S6), 419 ± 49 ppm (S32) and 408 ± 40 ppm (S98), consistent with each other and with
the paper.

**Each sector searched on its own:**

| Sector and product | Strongest signal | Passes our cuts? |
|---|---|---|
| S98 QLP, 200 s | 3.1811 d, 477 ppm, S/N 11.8 | Yes |
| S6 TESS-SPOC, 30 min | 3.1825 d, 427 ppm, S/N 10.7 | Yes |
| S6 QLP, 30 min | 11.12 d, two dips, S/N 9.8 | No (only 2 transits) |
| S32 SPOC 2 min / TESS-SPOC 10 min / QLP 10 min | 3.18 d, S/N 6.9 / 8.2 / 7.0 | No (S/N < 9) |

Sector 32 alone falls short. That fits the author's note that it sits at SPOC's 7.1σ alert threshold.

## Notes

- **Thresholds compared.**
  - The author's shortlist needed a multiple-event statistic of 15. The signal scored 10.6–11.0 in Sector 98, so he
    looked at it only after the fact.
  - Our first cut, S/N ≥ 9, is lower, so our search would have flagged it directly.
  - That cut does not pick out real signals on its own. In Sector 104, 92 other new planet-like signals were at least
    this strong, and most turned out to be noise, systematics or eclipsing binaries. The recurrence test has to do
    that job, which is why its flaw matters.
- **The 11.13-day signal.**
  - In Sector 98 it is the strongest raw peak. Our search, however, ranks peaks by red-noise S/N and requires that no
    single event dominate. One of the signal's three events carries 57% of its S/N (the limit is 50%), and the
    3.18-day signal scores higher anyway.
  - A blind search of the Sector 6 QLP light curve does make two dips, at BTJD 1473.74 and 1484.85, its strongest
    signal. They fall within 4 and 17 minutes of the published 11.13-day ephemeris. The manuscript says an unanchored
    search failed to find them.
  - This is not independent confirmation. The published 11.13-day ephemeris was fitted partly to these Sector 6 dips,
    and the detection is weak (two events, SDE 4.0).
- **Limits.**
  - Our search also masks time windows hit by sector-wide systematics. That step needs every star in a sector, so it
    was not reproduced for Sector 98. Masking one of the nine transits would lower the S/N by about 6%.
  - This check only asks whether our pipeline would detect the signal. It does not vet the candidate.

## What this means for our own search

The same flaw affected our own search of Sectors 104 and 105. Five Sector 104 signals missed the p < 0.001 cut with at
least one sector at the 1/151 floor. Re-tested with 3,000 random periods, all five pass (p = 0.00002–0.0003), so they
should have been vetted. They have since been vetted, and none is a planet: see
[`../../tess-search/vetting_2026-10-16/`](../../tess-search/vetting_2026-10-16/). The re-test is
[`../../tess-search/reaudit/recurrence_retest.py`](../../tess-search/reaudit/recurrence_retest.py).

## Files

| File | What it does |
|---|---|
| `run_all.sh` | Runs the four scripts below (`PY=<python> ./run_all.sh`) |
| `would_detect.py` | Search of the Sector 98 QLP light curve, with the ranking of the raw peaks |
| `would_detect2.py` | Search of every sector and product, the classifier and the recurrence test; writes `would_detect2.json` |
| `null_hires.py` | Recurrence test with 150 and with 3,000 random periods; writes `null_hires.json` |
| `s6_b_check.py` | Compares the Sector 6 QLP detection with the 11.13-day ephemeris |
| `output.txt` | Output of `run_all.sh`, reproduced from this folder |
