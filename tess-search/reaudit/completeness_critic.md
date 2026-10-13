The bundle's main answer still holds: there is no planet candidate, and the only ExoFOP-worthy item is TIC 87055315 as an eclipsing-binary/false-positive note. But the bundle overlooked one source of prior detections, got the novelty verdict on TIC 158159741 wrong, and missed a whole category of findings about objects already listed on ExoFOP.

**1. Checks I ran on the five priority objects (2026-10-13)**

- **MAST cone search, 4″, all collections** (`mast_cone.txt`, `dvcheck.txt`):
  - TIC 87055315 has only QLP S13, GSFC-eleanor-lite S13, eleanor postcards and raw FFIs.
  - None of the five objects has TGLC, CDIPS (except 100666720), PATHOS, DIAmante or TESS-SPOC DV/TCE products. TESS-SPOC light curves exist for 143042674 (S4, S31) and 404124262 (S13, S27), with no DV.
- **TEQUILA** (Ogunwale+2025, A&A in press, arXiv:2511.19313; variable sources found by difference imaging in S1–S26, no periods). I read the S2, S4 and S13 catalogues:
  - TIC 87055315: no source within 90″ in S13.
  - TIC 143042674: source 16.6″ away (one TESS pixel is 21″), magnitude 10.74, classed astrophysical by TEQUILA's neural net (is_systematic=0, confidence 0.93). This is a prior "variable" flag with no period.
  - TIC 281740175 (1.9″) and 404124262 (2.5″): sources exist but are classed as systematics (is_systematic=1, confidence 0.60 and 0.53), so they are not credible prior detections.
- **CDIPS header for TIC 100666720** (`cdips_out.txt`): CDIPSREF = HATSandHATNcandidates20210505, CDREFBIB = JoelHartmanPrivComm. The star is in CDIPS because it was on a HATSouth/HATNet candidate list by May 2021. The HAT team probably knew of the signal; its period is not public. This weakens the bundle's "strongest survivor" verdict.
  - Side finding: TIC 212554303 is in CDIPS as a Kounkel+2019/2020 group member (log age 8.3), so it is possibly a young M-dwarf EB.
- **Zooniverse Talk**, Planet Hunters TESS (project 7929) and Eclipsing Binary Patrol (project 22939), searched for all 62 EBs plus the key objects:
  - 0 hits for 87055315, 87055304, 100666720, 143042674, 281740175 and 404124262.
  - Hits on other EBs:
    - TIC 63594264, EB Patrol, 2025-02-15: "#eclipsingbinary … 3.6278/2 = 1.8139 days" (that is our period).
    - TIC 114917180 (2024-09-23) and 200247352 (2024-12-11): "#possible-new-eb".
    - TIC 80453626, PHT, Aug–Sep 2026: "Maybe EB".
  - Limitation: Talk search only finds comments that contain the TIC number. EB Patrol subject metadata could not be searched.
- **Live VSX** (vsx.aavso.org API, 1′, all 71 objects): nothing new compared with the VizieR copy. No entry within 21″ of any object except the already-known ROT entry for 155991187.
- **TAPVizieR sweep of every table with a TIC column** (414 tables; all succeeded after retries; the TIC catalogue itself was skipped). This covers tables that have no coordinates, which the cone sweeps could not reach. The only new variability source is **Tey et al. 2023, AJ 165, 95 (J/AJ/165/95/table2)**, a set of 24,926 QLP threshold-crossing events with human labels.
  - In Tey's scheme, "E" means a periodic eclipsing signal, planet or EB.
  - It lists 7 of our EBs with matching periods: 80064194, 63594264, 304444218, 76224662, 396957421, 200247352, and 231286114 at half our period.
  - It also lists **TIC 158159741 at P = 1.807655 d, depth 1720 ppm, epoch BTJD 1658.3672** (ours: 1.8068919 d, 1658.3739), labels E,E,N,E,J,E.
  - None of 87055315, 87055304, 100666720, 143042674, 281740175, 404124262, 320686138, 155991187, 72651793 or 245976814 is in it.
- **ExoFOP re-check today:** the three target pages are still empty, and the TOI list is unchanged (8151 rows, max TOI 7927.01).

**2. Corrections to the bundle**

- **TIC 158159741 is not "novel_as_claimed".** Its 1.807-d dip was already a published QLP threshold-crossing event (Tey+23) that passed triage but never became a TOI. Why it was dropped is not public.
- **More EBs already have public periods.** Counting DV reports, Green+23, TARS, KELT and Tey+23 together, 34 of the 62 have a public period matching P, 2P or P/2. The bundle said 30; Tey+23 adds 76224662, 200247352, 231286114 and 396957421. 28 have none found (`prior_union.txt`).
- **The four "Gaia NSS solutions" are not binary detections** (143345010, 278802197, 33445288, 351474229). All four are astrometric accelerations (Acceleration7 or 9) with no orbit and no period (`nss4.txt`). They point to a long-period third body, not the eclipsing pair, so the README's "binarity known" is wrong for all four, not just 351474229.
- **Internal inconsistency in the bundle.** The novelty caveats say the TOI list stops at sector 99, but the policy section shows TOIs from S100–S102 exist. The policy section is right.
- **Wording for the four priority EBs:**
  - 100666720: "on a HATS/HATNet candidate list, period not public".
  - 143042674: "TEQUILA variable flag, no period".
  - 281740175 and 404124262: unchanged.

**3. ExoFOP-relevant items the bundle missed**

- **Existing CTOIs that our pipeline shows as eclipsing binaries.** 37 stars carry planet-candidate CTOIs where our search recovered the signal with a secondary eclipse or odd/even difference above 5σ (`ctoi_matches_detail.csv`).
  - 25 of these are lipponen's 2019 bulk uploads, plainly EBs (depths 1–99%). Low value.
  - 12 are published candidates:
    - RAVEN (Lafarga+2026, uploaded 2026-08-27): 71316387.01 (odd/even 13.9σ) and 147285336.01.
    - DTARPS (Melton+2023): 88571421.01, 92352621.01, 102064145.01, 115693377.01, 206599822.01.
    - DIAmante (Montalto+2020): 65450588.01, 278536652.01.
    - Bryant & Van Eylen 2025: 218856662.01.
    - Feinstein+2019: 234503282.01.
    - chiarello: 382360740.01.
  - Also flagged: devadas 321044007.01, and TOIs 6787.01, 3288.01, 1124.01 and 138.01.
  - A crude re-check at the catalogue ephemerides (`oddeven_out.txt`) supports three:
    - 71316387.01: odd 3282 ± 200 vs even 1816 ± 113 ppm, 6.4σ.
    - 88571421.01: phase-0.5 dip of 1025 ppm, 6.0σ.
    - TOI 3288.01: 2031 ppm, 3.6σ.
  - The ephemerides no longer line up for 147285336.01, 206599822.01, 65450588.01 and 218856662.01 (no transit seen at the catalogue times), so those are untested.
  - **Caveat:** the same flags fire on 5 of 42 known or confirmed planets (TOI 107.01, 117.01, 185.01, 190.01, 403.01). These are leads that need proper vetting, not findings.
  - If confirmed, they concern objects already on ExoFOP, where an observing note would help. No paper is needed to add a parameter set or disposition to someone else's CTOI only if ExoFOP allows that, and its rules don't say (not verified).
- **TIC 320686138** is a fully characterised, on-target EB (secondary 1531 ± 155 ppm). It was left out of the 62-EB table only because its 0.37% primary is under the 0.5% cut. It could go to VSX.
- **TIC 146284379** (CD-33 815; 3.85379 d, 0.54%/0.13%) was dropped for incoherent eclipse timing and never novelty-checked.
- **Eclipses on neighbouring stars:** about 20 signals fall 10–46″ off target (for example the neighbours of TIC 281705082). None is pinned to a specific star, so none can be submitted. Cross-matching all 60 against TOIs/CTOIs found only TIC 139285753, which is the known false positive TOI 377.01 (`blend_xmatch.txt`). No new TOI contamination finding.
- **Single transits:** `bls.py` masks single events and never saves them, so there are no single-transit candidates to submit. That is a limit of the search, not a missed result.

**4. Still under-checked**

- **Localisation of the 62 EBs.** They were placed with `src/centroid.py`, the same estimator shown to be biased in the TIC 87055304 case. For 22 of the 62, a neighbour is bright enough to cause the eclipse, so "on target" rests on that biased test. VSX moderators will ask which star varies; a pixel-level PSF check is needed before submitting those 22.
- **Timing risk from Gaia DR4, due 2 December 2026** (https://www.cosmos.esa.int/web/gaia/data-release-4). It will release epoch photometry and could publish these eclipses first. Re-check novelty after that date.
- **Not verified:**
  - ADS full-text search.
  - HATSouth, WASP, NGTS and KELT internal candidate lists.
  - QLP and TOI-team internal vetting of TIC 87055304 and 158159741.
  - EB Patrol subject metadata.
  - TEQUILA for the other 58 EBs.
  - VizieR tables keyed only by Gaia, TYC or 2MASS ID, with no coordinates.
  - Whether the RAVEN paper itself flags 71316387.01 or 147285336.01.
- **Open policy questions** (ask exofop-support@ipac.caltech.edu):
  - Can a third party add a false-positive parameter set or observing note to another user's CTOI without a paper?
  - Does an RNAAS that cites published tools (BLS, TRICERATOPS, PRF fitting) meet the "methodology previously published" rule? Our pipeline as a whole is unpublished.
  - Is the user eligible for an account? They need an affiliation or "Citizen Scientist" approval; unknown.
- **Venues the bundle did not consider:**
  - The HATSouth team (Hartman/Bakos) about TIC 100666720.
  - The TESS team's contact page (https://tess.mit.edu/contact/; contents not verified), to flag TIC 87055304 ahead of any QLP alert.
  - A Zenodo/GitHub data release with a DOI: citable, but it would not count for ExoFOP.
  - arXiv, which needs an endorsement.

Everything is in /tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/exofop/agents/critic2/; `notes.txt` indexes the outputs and `scripts/` holds the code. Nothing under /home/user/Space or /home/user/work was modified, and nothing was submitted anywhere.