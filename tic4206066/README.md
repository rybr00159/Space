# TIC 4206066: a TESS planet candidate found with Claude Code

In September 2026, Pavel Rabtsevich, a product manager and amateur astronomer in Spain, reported two transit-like signals
in public TESS data of the nearby K dwarf TIC 4206066. He did most of the analysis with Claude Code and used Codex and
fresh read-only agent sessions as reviewers. TESS was approved to observe the star again at 2-minute cadence in
**Sector 110 (2026-10-31 to 2026-11-26)**. Before that data existed, he published the exact transit times and pass/fail
rules, so the November data is a direct test.

**Status (2026-10-11):** planet candidate, not confirmed or validated. Not yet peer reviewed.

## Key numbers

From the preprint ([doi:10.5281/zenodo.22967456](https://doi.org/10.5281/zenodo.22967456)) and the pre-registration
([doi:10.5281/zenodo.23175179](https://doi.org/10.5281/zenodo.23175179)).

**Star:** TIC 4206066 = StKM 1-561 = Gaia DR3 3220388198192519424. Late-K dwarf, TESS mag 9.68, 35.6 pc (116 ly).
R★ = 0.623 ± 0.019 R☉, M★ = 0.606 ± 0.015 M☉.

| | Signal A | Signal B |
|---|---|---|
| Period | 3.182785 ± 0.000006 d | 11.13274 ± 0.00003 d |
| Depth | ≈500 ppm (0.05%) | ≈1000 ppm |
| Duration | 1.87 ± 0.05 h | ≈1.7 h |
| Radius if it orbits the target | 1.4 ± 0.1 R⊕ | ≈2.2 R⊕ |
| Archival transits | 23, in Sectors 6 (2018), 32 (2020), 98 (2025/26) | 2 clean in S98, dips in S6, marginal in S32 |
| Strength | False-alarm probability ≲0.1% across sectors | ≈3σ, in a test defined after looking (post hoc) |
| Expected S110 depth (frozen) | 442 ± 52 ppm | 586 ± 137 ppm |
| Predicted S110 transits | 9 (one in a scheduled data gap) | 3 |

**Vetting summary (preprint):**
- Passed all LEO-Vetter flux-level tests.
- Odd and even transits have the same depth, and there is no secondary eclipse (phase-0.5 depth −20 ± 32 ppm).
- A V-shaped dip is disfavoured (Δχ² = 32), and there is no sign of transit timing variations.
- The difference-image centroid is (+5 ± 3, 0 ± 3)″ from the target, about 2σ. Every Gaia star bright enough to fake
  the dip is at least 43″ away.
- TRICERATOPS false-positive probability is 0.03–0.04, above the 0.015 validation bar, so he claims no validation.
  The main remaining alternative is a planet orbiting an unresolved bound companion star. Ruling that out needs
  high-resolution imaging and radial velocities.

## What to watch for in November

[`predicted_transits_S110.csv`](predicted_transits_S110.csv) lists the 12 predicted transits: 9 of signal A and 3 of
signal B. Each has a centre time known to about ±4 min (1σ). The table comes from the frozen
[`sources/ephemeris_final.json`](sources/ephemeris_final.json). The times are calendar labels of barycentric (TDB)
times; light from the star reaches Earth 5.6–7.2 min earlier in UTC.

Under the pre-registration, each signal gets exactly one outcome from the first released SPOC 2-minute light curve:

| Outcome | Meaning |
|---|---|
| O1 recovered | Depth above the detection threshold, inside the allowed range, and on time (within 3σ) |
| O2 inconsistent | Detected, but the depth or timing disagrees with the prediction |
| O3 excluded | Not detected even though sensitivity was at least 95%. Counts against the signal |
| O4 insufficient | The data cannot decide, for example too few usable transits |

Even an O1 only shows the signal is real and periodic. It does not prove a planet, because an unresolved eclipsing
companion would also pass. The author will post the frozen scripts' outputs as new versions of the
[Zenodo record](https://doi.org/10.5281/zenodo.23175179).

## Feedback from the Reddit thread

See [`reddit-post/comments.md`](reddit-post/comments.md) for the saved comments.

- **u/ErrorlessQuaak (works on exoplanets):** ran it through LEOVetter and TRICERATOPS and called it a "totally
  reasonable candidate". Said the S/N is borderline by SPOC's noise estimate. OP added that Sector 32 alone sits right
  at SPOC's 7.1 detection cutoff, so no SPOC alert was raised. Suggested submitting it to ExoFOP as a Community TOI. OP replied that ExoFOP requires a
  peer-reviewed paper first, and that the full paper is in progress.
- **u/twotootstotucson (runs a similar Claude Code TESS survey):** worries most about a blend at ≈500 ppm on 21″
  pixels, because a background eclipsing binary would also pass the held-out-year test. Asked about high-resolution
  imaging and the 30-minute cadence of Sector 6. Their own best candidate fell apart once they derived the stellar
  radius themselves.
- **u/cyril_zeta:** if both signals are real planets, transit timing variations of 10–30 min could help validate the
  system.
- **u/Natural-Lifeguard-38:** could it be a variable star? The preprint's tests argue against this: odd and even
  depths match, there is no secondary eclipse, and the dip has a flat bottom rather than a V shape.

**Audit lessons for agent-driven science,** from OP and other commenters:
1. Ask the agent to disprove the result, not confirm it.
2. Make every check fail once on purpose. A guard that has never failed might just be `return True`.
3. Work out the false-alarm rate from your own null distribution, and recheck it as that distribution grows.
4. Re-derive headline numbers from the raw output, never from a summary the agent wrote.
5. Watch which way a check fails. An empty list wired to "pass" looks identical to a clean result.

## Would our own search have found it?

Checked 2026-10-16; details in [`pipeline_check/`](pipeline_check/).

- **No, because the star was never in our data.** TESS observed it only in Sectors 6, 32 and 98, and our
  [search](../tess-search/) covered Sectors 104 and 105.
- **Given the Sector 98 data, our search would have detected the 3.18-day signal** (S/N 11.8, 9 transits), found it
  in no catalogue and classed it as planet-like.
- **Our recurrence test would then have dropped it,** with p ≈ 0.004 against a cut of 0.001. The test was too coarse:
  with finer random trials it gives p = 0.0003. The same flaw left five of our own recurring signals unvetted.
- **The 11.13-day signal would not have been reported.** Our search keeps one signal per star and sector and needs
  three transits in a sector.

## Files

| Path | What it is |
|---|---|
| [`pipeline_check/`](pipeline_check/) | Our search pipeline run on this star's light curves: scripts, outputs, write-up |
| [`reddit-post/post.md`](reddit-post/post.md) | The r/ClaudeAI post (2026-10-07), verbatim, with its 3 images |
| [`reddit-post/comments.md`](reddit-post/comments.md) | 100 of the thread's 732 comments, threaded, as saved on 2026-10-11 |
| [`predicted_transits_S110.csv`](predicted_transits_S110.csv) | Predicted Sector 110 transit windows, derived from the frozen ephemeris |
| [`sources/preprint_zenodo-22967456.pdf`](sources/preprint_zenodo-22967456.pdf) | Preprint: *Two Transit-like Signals in TESS Photometry of the Nearby K Dwarf TIC 4206066* |
| [`sources/PREREG_S110_SHORT.pdf`](sources/PREREG_S110_SHORT.pdf) | Pre-registration summary (start here) |
| [`sources/PREREG_S110_FINAL.pdf`](sources/PREREG_S110_FINAL.pdf) | Full pre-registration (governs) |
| [`sources/ephemeris_final.json`](sources/ephemeris_final.json) | Frozen ephemeris with per-transit predicted times |
| [`sources/expected_depths.json`](sources/expected_depths.json) | Frozen expected depths |
| [`sources/prereg_zenodo-23175179_README.md`](sources/prereg_zenodo-23175179_README.md) | The deposit's README: how to reproduce, file guide |

The files in `sources/` are unmodified copies from Zenodo. Their sha256 values match the deposit's `SHA256SUMS`. The
full 13.8 MB deposit, with the frozen code and calibration runs, is at
[doi:10.5281/zenodo.23175179](https://doi.org/10.5281/zenodo.23175179). The supporting data and code for the
preprint are at [doi:10.5281/zenodo.22967099](https://doi.org/10.5281/zenodo.22967099).

## Links

- Reddit post: <https://www.reddit.com/r/ClaudeAI/comments/1wzw8zd/i_think_i_found_a_planet_nobody_knew_existed_i/>
- Author's X thread: <https://x.com/p_rabtsevich/status/2107803696058822742>
- Interactive 3D model of the system: <https://tic4206066.pages.dev>
- TESS Director's Discretionary Targets (DDT 100): <https://tess.mit.edu/science/ddt/>
- Gizmodo coverage: <https://gizmodo.com/did-claude-help-an-amateur-astronomer-discover-a-new-planet-2000824393>

## Credit and licences

Research, preprint, pre-registration, figures and post images © Pavel Rabtsevich. The Zenodo materials are CC BY 4.0
(cite as in their `CITATION.cff`). Reddit comments belong to their authors and are kept here for reference only.
