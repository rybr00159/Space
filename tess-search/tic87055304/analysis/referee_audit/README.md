# Referee audit: is the "neighbour-hosted" localisation of the TIC 87055304 signal a shared systematic? And what did src/centroid.py really measure?

Folder: `/home/user/work/tic87/referee_audit/`. All code here was written for this audit. It reads the shared data and
(read-only) the other investigators' JSON outputs. Python: `/home/user/work/venv/bin/python`, `OMP_NUM_THREADS=1`, niced. No
TRICERATOPS runs and no paid APIs were used. All numbers are collected in `results.json`.

Ephemeris for Task A: P = 4.136195 d, T0 = 1658.41325 BTJD, T14 = 0.1875 d. Task B uses the vet run's own ephemeris and
duration (see Task B).

## Verdicts

**Task A: adversarial referee of the localisation.**
- I found **no shared systematic** that could make a signal hosted by the target look as if it were hosted by the neighbour.
  The conclusion of the three analyses survives: the host is TIC 87055315 (Gaia DR3 6720986208895208960), or a source
  within about 1" of it.
- **WCS and parity.** The WCS orientation and parity are correct in all five sectors, including S104 on CCD 3. This was
  checked three independent ways: against tess-point, against PSF-free centroids of field stars, and against mirrored or
  transposed alternatives.
- **Registration, jitter, timing.** Registration, pointing jitter and timing/quality masking cannot produce the offset.
- **Own PRF-free test** (the single-pixel depth ratio). It rejects a target host. The size of the rejection depends on how
  much is assumed:
  - **6.4σ** (Gaussian-equivalent; 5 sectors combined, primary a-priori pixel pair). This assumes the most target-friendly
    of 180 contamination models.
  - **≥ 4.2σ** by pure counting: none of 2×10^5 fake-transit combinations reaches the observed value.
  - **4.9σ** if the noise is 30 % larger than measured.
  - **4.0–4.3σ** from a zero-PSF-model flux inequality in S93 alone.
- **Verdict scale: supports false positive** (for a planet on TIC 87055304).
- **Caveat.** Several of the quoted significances are overstated or not defensible as written; see section 6.

**Task B: audit of src/centroid.py.**
- The earlier claim "on target; the 9" neighbour excluded at >6σ" was **invalid**. src/centroid.py is biased toward the
  target by construction.
- Its S66 result, the basis of the ">6σ" claim, is dominated by two scattered-light windows. Its quoted σ contains a
  0.15-px floor.
- The ">6σ" came from dividing the target–neighbour distance by the floor-free per-axis bootstrap error `err_px` (0.04–0.07
  px). The same code's own σ is about 0.16 px (3.4") per axis, which is what produced "0.3σ".
- **Injection test.** For signals injected on the neighbour, centroid.py reports offsets of median 6.3–7.3" at a median
  "offset_sigma" of 1.6–1.9. Only 6 % of these reach 3σ, so it would have called a neighbour-hosted signal "on target"
  about 94 % of the time.

---

## TASK A

### 1. Gaia positions vs the TESScut WCS; could a rotation or flip mirror the source? (`r1_wcs_check.py` → `r1_wcs_check.json`)

**Method:**
- **Gaia positions.** Taken from the shared 90" CSV, plus my own VizieR copy of Gaia DR3 for the full stamp (`fetch_gaia.py`,
  G < 18 within 330"; it matches the CSV to 2e-8").
- **Propagation.** Positions propagated with proper motion to each sector's median epoch.
- **WCS.** From ext 2 (APERTURE), and cross-checked against the FLUX-column WCS in ext 1.

| Sector | cam/ccd | det(CD) | PA of +x / +y | Neighbour vector (px) | Sep., PA | tess-point − WCS vector (px) |
|---|---|---|---|---|---|---|
| 13 | 1/1 | + | 92° / 1° | (+0.317, +0.332) | 9.37", 44.9° | (0.0000, −0.0001) |
| 39 | 1/2 | + | 84° / 355° | (+0.367, +0.290) | 9.39", 44.9° | (+0.0003, −0.0001) |
| 66 | 1/2 | + | 85° / 356° | (+0.343, +0.292) | 9.41", 44.8° | (0.0000, +0.0001) |
| 93 | 1/2 | + | 86° / 356° | (+0.342, +0.294) | 9.44", 44.8° | (−0.0001, 0.0000) |
| 104 | 1/3 | + | 224° / 134° | (−0.449, +0.006) | 9.45", 44.8° | (+0.0001, 0.0000) |

- **Parity.** It is the same (det > 0) in all sectors.
- **S104.** The pixel grid is rotated by about 135°, so the neighbour lies along −x there and along the +x+y diagonal in
  S13–S93.
- **tess-point check.** tess-point uses the spacecraft focal-plane model, not the TESScut WCS. It gives the identical
  target→neighbour pixel vector, to ≤ 3e-4 px.
- **The two WCS headers agree** exactly (ext 2 and the ext-1 column WCS).
- **Field-star test (no PSF model).** I used 7–10 isolated stars with T < 15 per sector and measured flux-weighted
  centroids in 3×3 boxes.
  - Under the actual WCS they lie a median 0.08–0.24 px from their predicted positions.
  - Under the alternatives they lie 1.25–2.15 px away. The alternatives tested were: mirror in x or in y about the target,
    transpose about the target, 180° rotation, and swapped array axes.
- **Proper motions.** Changing the proper motions moves the target→neighbour vector by < 0.1". The target's RUWE of 18.8
  implies mas-level errors, which are irrelevant here.

**Why a flip cannot be the explanation:**
- A mirror about the target leaves a target-hosted signal at the target.
- Any mirror or rotation that does not pass through the target would also misplace every field star by more than 1 px.
- A pixel-frame artefact rotates with the CCD, but the source stays at sky PA about 45° in S104 as in the other sectors.

**Verdict:** no WCS, parity or rotation error. Neutral for the planet, but it removes this route to a target host.

### 2. Is "registration on field stars" independent of the pair? (`r2_registration.py` → `r2_registration.json`)

**Method:** Fisher information for a global shift along u (the target→neighbour axis). I used a SPOC-PRF scene of all Gaia
stars with T < 18 and computed the fraction of that information carried by the target+neighbour pair.

| Pixel set (whose registration) | Pair share of the shift information along u |
|---|---|
| r ≤ 4 px circle (loc_independent "local" fit) | **85–97 %**, i.e. *not* independent: this is essentially registration on the pair |
| annulus 2.5–8 px (loc_independent "external") | 0.1–0.4 % |
| r > 3 px (psf_photometry) | ≤ 0.04 % |

- The "field-star" registrations of psf_photometry and of loc_independent's external variant **are effectively independent
  of the pair**.
- The local-fit registration is not independent. It is also the one that can be biased by a mis-modelled flux ratio, which
  is what happened in `transit/prf_localize.json` with its negative fluxes. Local and external agree to ≤ 1.3".
- **It does not matter in any case.** Turning a target host into the observed answer needs a registration error of
  0.38–0.39 px (8") along u in *every* sector. The WCS alone, with no registration, puts field stars within −0.075 to
  +0.091 px along u.
- My pixel test (section 5) uses no registration at all: only which pixel is which. Its model predictions include ±0.1-px
  shifts.

### 3. Is the pointing-jitter leverage argument right? (`r5_jitter.py` → `r5_jitter.json`)

**Leverage:**
- If the image shifts by s during transit, D picks up F_pair·s·∇PSF. A free-position fit then moves the source by
  −(F_pair/A)·s.
- F_pair/A = 220–258, where F_pair = 3100 e/s (Gaia T mags, zero point calibrated on the stamp) and A = 12–14 e/s
  (measured).
- **Checked numerically** with the PRF: injecting a 1 mpx shift and fitting the position gives a leverage of −201 to −239.
  The argument is right in size and sign: the apparent source moves opposite to the image shift.
- **Shift needed to fake the neighbour:** s_u = −1.8 to −2.1 mpx along u in every transit, i.e. the image moving away from
  the neighbour.

**Measured pointing:**
- **Method.** The in-transit shift was measured from the difference images themselves, without a PSF model. I regressed
  D in field pixels (r > 3 px) on the gradient of the mean image, calibrated by a known 0.01-px spline shift (the factor is
  1.24–1.43).
- **Per-transit shifts:** |s_u| ≤ 0.64 mpx.
- **Sector stacks:** −0.44, −0.07, +0.06, −0.17, +0.20 mpx, giving spurious s_par = +0.21, +0.04, −0.03, +0.09, −0.11 (mixed
  sign).
- **Fake-window scatter:** 0.18–0.31 mpx (MAD).
- **Conclusion.** Jitter explains at most about 0.2 of the observed s_par ≈ 0.8–1.0, and only in S13. Jitter is not locked
  to a 4.136-d ephemeris over 7 years, so it is not a viable shared systematic.

**Caveat:** the per-transit sign test (26/26 transits on the neighbour side) protects against jitter and noise only. A
shared offset (WCS, registration, PRF) would move every transit together. Those are excluded by sections 1, 2 and 5, not by
the sign test.

### 4. In-transit selection, baselines and quality masking in BTJD (`r4_timing.py` → `r4_timing.json`)

- **Time system.** TESScut has TIMEREF = SOLARSYSTEM, TIMESYS = TDB and BJDREFI = 2457000 in all sectors, so TIME is BTJD.
  QLP is the same.
  - The nearest-cadence TESScut−QLP time difference is −5 to −20 s (S13–S93) and +54 s (S104).
  - This is negligible against the ±70-min in-transit windows.
- **Transit time in the pixels.** I used a pair aperture (r ≤ 2 px, background-subtracted), a per-transit linear baseline
  and a fixed-shape trapezoid.
  - dt0 = −4.3 ± 4.3, +2.9 ± 3.2, +6.5 ± 2.9, −4.3 ± 3.2 and −4.3 ± 4.0 min. The ephemeris σ is 3.9–7.3 min.
  - Scanning an artificial ephemeris shift from −0.4 to +0.4 d, the depth peaks at 0.00 d in every sector.
  - The depth in this aperture is 2.3–2.9 ppt, consistent with the other analyses' background-subtracted values.
- **Quality flags.** Flagged fractions are 0–17 % in transit and 6–20 % in the nearby out-of-transit cadences, so in-transit
  cadences are not preferentially removed.
  - The bits are 2048 (stray light), 128, 16384, 8 and 32.
  - S104 has no flagged in-transit cadences. Its scattered-light episode (BTJD 4178.8–4179.2) falls at phase 0.39, not in
    transit.
- **Windows and baselines of the three analyses** are correctly placed:
  - loc_independent: in-transit |dt| < 0.35 T14; out-of-transit 1–2.5 T14, linear;
  - psf_photometry: ±0.45 d with a quadratic baseline and a batman template;
  - this audit: in-transit |dt| < 0.375 T14; out-of-transit 0.75–2.0 T14.
- **centroid.py** used the 3.0-h BLS duration, so its out-of-transit window starts at 0.094 d, right at egress.
  - Only 0–2 % of its out-of-transit cadences fall inside the true transit. This is harmless.
- **Verdict:** no timing or masking error.

### 5. Decisive PSF-model-free test: single-pixel depths (`r3_pixel_ratio.py`, `r3c_flux_inequality.py`; `fig_pixel_test.png`)

**Statistic.**
- Per transit, D = baseline − in-transit (per-pixel linear baseline from both sides), with a far-field plane removed.
  Transits are stacked with one weight each. O is the out-of-transit image minus its 10th percentile.
- The relative depth is δ(p) = D(p)/O(p).
- **If the target hosts the transit, δ(p) = δ_T·f_T(p)** in every pixel, where f_T is the target's share of that pixel's
  light (0 ≤ f_T ≤ 1). So for any two pixels:

  Q = δ(x)/δ(y) = f_T(x)/f_T(y) ≤ 1/f_T(y)

  No PSF enters this relation. With x on the neighbour side, f_T(x) < f_T(y) and Q < 1 is expected.
- If the neighbour hosts, Q = f_N(x)/f_N(y), which is ≫ 1.

**Pixels**, fixed a priori:
- a = the pixel nearest x_t + 1 px·u (neighbour side);
- b = the pixel nearest x_t − 1 px·u (far side);
- t = the target's pixel;
- n = the neighbour's pixel (distinct from t only in S93 and S104).

**Noise.** I used 4000 stacks of fake transits at phases 0.20–0.80, with the same estimator and weights.

**Predictions.** A PSF is needed only to say how small f_T(y) could be. I used the most target-friendly value over 180
scene models: SPOC PRF or Gaussian σ = 0.5/0.65/0.8 px; ±0.1-px shifts; contaminant fluxes ×0.7/1/1.3; background
±20 e/s.

**Primary pixel pair a/b (the requested test).** Values are listed in sector order S13, S39, S66, S93, S104.

| | S13 | S39 | S66 | S93 | S104 |
|---|---|---|---|---|---|
| δ(a) / δ(b) (ppt) | 3.3 ± 1.7 / 0.3 ± 5.8 | 3.6 ± 0.6 / 1.0 ± 0.7 | 3.1 ± 0.8 / 2.9 ± 1.3 | 4.4 ± 0.4 / 1.0 ± 0.4 | 4.4 ± 0.5 / 2.0 ± 1.1 |
| Q observed | (13; δ(b) ≈ 0) | 3.8 | 1.05 | 4.3 | 2.2 |
| Q if target hosts (max of 180 models) | ≤ 0.65 | ≤ 0.69 | ≤ 1.32 | ≤ 0.69 | ≤ 1.14 |
| Q if neighbour hosts (range) | 1.9–24 | 2.0–35 | 2.6–63 | 1.7–10 | 2.1–17 |
| z against a target host | 0.6 | 4.8 | −0.4 | 7.9 | 1.5 |

- **Combined (Stouffer): 6.4σ** (6.3σ with symmetric scattered-light cleaning).
- With a 30 % larger noise σ: 4.9σ.
- **Empirical.** None of 2×10^5 combined fake-null draws reaches the observed value (p < 1.5e-5, > 4.2σ). The uncleaned
  null has heavy tails (its maximum is 6.40), so the Gaussian translation is somewhat optimistic.
- **A neighbour host is consistent:** the observed Q lies inside or at the edge of the neighbour-host range in every sector.
- S13 (δ(b) ≈ 0, noisy) and S66 (pixel b holds a T = 14.7 star, so f_T(b) is unconstrained) carry no information. I report
  them as they are.

**Other a-priori pairs** (same data, not independent of a/b, so not combined with it):

| Pair | Per-sector z against a target host (S13, S39, S66, S93, S104) | Combined |
|---|---|---|
| a/t | 1.1, 2.8, 1.4, 6.3, 3.9 | 6.9σ |
| half-planes Ga/Gb (s ≥ +0.5 vs ≤ −0.5 px, r ≤ 2 px) | 0.7, 6.0, 2.1, 6.3, 1.6 | 7.5σ |
| n/t (S93, S104 only) | 6.0, 3.9 | 7.0σ |

Example: in S93, Q(n/t) = 2.08 against ≤ 0.91 for a target host. In S104 it is 2.03 against ≤ 0.80.

**Increasingly assumption-free versions:**
- **Conservative envelope instead of a model** for the contamination bound (every non-target star puts the maximum
  pixel fraction of any of the four PSFs into the denominator pixel, ×1.25, plus 20 e/s).
  - n/t: z = 2.8 (S93) and 1.6 (S104), combined **3.1σ**.
- **Zero PSF model** (`r3c_flux_inequality.py`). For any symmetric, unimodal PSF, the pixel containing a star gets at least
  as much of its light as the edge-adjacent pixel. So a target host requires ΔF(n) ≤ ΔF(t).
  - **S93:** ΔF(n) = 3.56 > ΔF(t) = 1.96 e/s, even though the target pixel is brighter (O = 972 vs 850 e/s). This gives
    **z = 4.0–4.3 from S93 alone**, and every one of its 5 transits shows ΔF(n) − ΔF(t) = 1.2–2.2 e/s.
    - It requires the target to be in pixel t. It sits 0.146 px from the boundary, against a WCS/registration accuracy of
      about 0.04 px.
  - **S104:** 0.32 ± 0.59 e/s, i.e. 0.5σ.

**The single most decisive number** that I recomputed myself:
- **Q(a/b) = 4.3 in S93 and 3.8 in S39.** A target host allows at most about 0.7 under any of 180 contamination models,
  while a neighbour host predicts 1.7–35.
- Over 5 sectors this gives **6.4σ** (Gaussian-equivalent) against a target host, with no PSF fit and no registration.

### 6. Statistical independence and the quoted significances (`r6_stats.py` → `r6_stats.json`)

**Per-sector independence:**
- The noise is independent between sectors (separate data, different pointings, and S104 rotated).
- All analyses nevertheless share the TESScut cutouts, the Gaia positions/WCS, the ephemeris and (two of three) the same
  SPOC start_s0004 PRF files.
- Combining sectors is valid only once common-mode errors are excluded. Sections 1, 2 and 5 do that: the WCS is verified
  independently, the registration error needed is 0.38 px, and the pixel test is PRF-free.

**Injection calibrations cannot detect a WCS or registration error, by construction.** The signal is injected at the
WCS-predicted position and recovered there. The ">8σ" therefore depends on section 1 being right, and it is.

**Fake and control windows overlap, so the effective numbers of independent windows are small:**

| Analysis | Effective independent windows |
|---|---|
| loc_independent | about 3.2 per transit epoch (11 nominal), i.e. about 16–19 per sector |
| this audit | about 4.3 per epoch (16 nominal) |
| psf_photometry control phases (each a joint fit over all epochs) | about 4 (65 nominal) |

Consequences:
- The null σ values are uncertain by about 20–35 %.
- Tail probabilities beyond about 3σ are parametric extrapolations, not measurements.

**Assessment of the quoted significances:**
- **loc_independent ">8σ (injection-calibrated)".**
  - z_T = 9.1 assumes a Gaussian null with σ = 0.103. With a 30 % larger σ it is 7.0σ.
  - Empirically, 0 of 20,000 simulations reach the real value. That gives only p < 5e-5 (3.9σ) by counting, and the
    simulations are combinations of about 16–19 independent windows per sector.
  - The Fisher "≥ 5.3σ conservative" is floor-limited by 600 resamples that are strongly correlated. It is not an
    assumption-free bound.
  - The "10^-31.6" per-transit sign test protects only against noise and jitter.
  - **Defensible wording:** "target host rejected at ≳5–7σ, with the exact value depending on tail assumptions".
- **psf_photometry Δχ² = 129.** This is a stat-only likelihood ratio, so it is not a significance. Its errors (max of four
  estimators) partly rest on about 4 independent control windows.
  - With the adopted 0.07 systematic, f_target = 0.163, stat error 0.087 → f = 1 is excluded at **7.5σ**, and that is the
    number to quote.
  - A 30–40 % error underestimate would still leave Δχ² ≳ 75.
- **Overall.** The qualitative conclusion is robust. Several partly independent methods on the same pixels agree, and two of them (the
  Gaussian PSF and this audit's pixel test) do not use the shared PRF. The quoted σ values should be read as "≳5σ", not
  as precise numbers.

---

## TASK B: src/centroid.py and the ">6σ" claim (`b1_centroid_audit.py` → `b1_centroid_audit.json`; `fig_centroidpy_injections.png`)

**Replication:**
- `centroid_copy.py` is a byte-identical copy of `src/centroid.py` (md5 1053541e…).
- I used the central 15×15 pixels of the 21×21 cutouts; the target pixel agrees with report.json to < 0.001 px.
- Inputs as in the vet run: P = 4.136190053729093 d, t0 = 1658.411589268281 BTJD, **dur = 0.1249 d** (the BLS duration, not
  T14 = 0.1875 d), and the TIC J2000 position.
- The results are identical to report.json: S66 1.10" with err_px (0.039, 0.069) at 0.33σ; S104 1.69" with err_px (0.265,
  0.089) at 0.27σ; S93 NaN.

**What centroid.py measures:**
1. Per transit, it builds a per-pixel linear baseline from out-of-transit cadences at 0.75–1.75 dur. D = baseline −
   in-transit mean over |dt| < 0.4 dur.
2. It takes an **equal-weight mean** of the per-transit D.
3. It computes a **flux-weighted centroid of max(D, 0) inside a circle of r = 3 px centred on the target's catalogue
   pixel**.
4. The offset is measured from that catalogue pixel, not from the out-of-transit centroid. It is then multiplied by
   21"/px.
5. `err_px` is the per-axis bootstrap standard deviation over 5–6 transits, in **pixels**, with no floor.
6. `offset_sigma` = hypot(dx/hypot(ex, 0.15), dy/hypot(ey, 0.15)), i.e. it includes a **0.15-px floor in quadrature on each
   axis**.

**Why the offsets were small:**
- **Aperture-centre bias.** Positive-clipped noise and any positive pedestal add weight that is centred on the aperture,
  i.e. on the target.
  - A noiseless neighbour-hosted D gives s_par = 0.70–1.11 (0 = target, 1 = neighbour).
  - Adding a 0.3 e/s/px pedestal gives 0.43–0.76; adding 1 e/s/px gives 0.25–0.49.
- **S66 was dominated by scattered light.** Two of the six windows carry far-field pedestals of +5.0 and +29.8 e/s/px. Their
  clipped aperture sums are 147 and 816 e/s, against 13–16 e/s for the clean windows.
  - The four clean windows point toward the neighbour individually (s_par 0.30–0.88).
  - The equal-weight mean gives s_par 0.12, i.e. 1.1".
  - Subtracting each window's far-field median raises it to s_par 0.42.
- **S104** is pulled by one window. Five of six windows point to the neighbour (s_par 0.28–0.91). One gives s_par −1.41,
  which makes the stack −0.19 and makes err_px(x) = 0.27 px.
- **S93** is NaN because two windows have pedestals of −150 and +87 e/s/px. Without clipping, or with the pedestal
  subtracted, it is computable (s_par 1.6, noisy).

**Why the σ values were large.** The 0.15-px floor is 3.2" per axis. In S66, σ = hypot(0.039–0.069, 0.15) = 0.155–0.165 px
= 3.3–3.5" per axis, so 1.1" corresponds to 0.33σ. In S104 the bootstrap term dominates: σ_x = 0.30 px = 6.4". **The
"0.3σ" therefore implies σ ≈ 3.4" (S66) and ≈ 6" (S104), not 0.04–0.07 px.**

**Why the ">6σ" inference was invalid.** Yes, it was a misread, of the error definition more than of the unit:
- **Wrong error.** It divided the 0.37-px gap between the S66 centroid and the neighbour by the floor-free per-axis
  bootstrap error (err_px = 0.04/0.07 px = 0.8"/1.5"). My replication of that arithmetic gives 8.2. The code's own σ,
  including the floor, gives **2.4σ**: even an unbiased estimator could not have excluded a 0.45-px (9.4") neighbour at > 3σ
  with that σ.
- **Wrong expectation.** It assumed that a neighbour-hosted signal would put the centroid at the neighbour's catalogue
  position. Calibrated on injections, centroid.py returns s_par = 0.63 [0.24, 1.10] (PRF) for a neighbour host, not 1.0.
- **One bad sector.** It relied on one sector (S66) whose stack is dominated by scattered light. A 6-transit bootstrap also
  understates the error.

**Quantified with injections.**
- **Set-up.** Transits of the real amplitude (12–14 e/s) were injected into the real cubes at 11 fake phases per sector
  (55 per host and PSF; 53 valid). Two PSFs were used: the SPOC PRF and a Gaussian with σ = 0.6. The code was then run
  exactly as in the vet.

| Injected on | Reported offset (") median [16–84 %] | "offset_sigma" median (max) | Fraction with offset_sigma ≥ 3 | s_par median [16–84 %] | Lead's statistic > 6 |
|---|---|---|---|---|---|
| Target (PRF) | 2.2 [1.0, 2.8] | 0.63 (1.8) | 0 % | 0.04 [−0.17, 0.25] | 72 % |
| Target (Gauss) | 2.0 [0.6, 2.6] | 0.58 (1.8) | 0 % | 0.05 [−0.12, 0.22] | 83 % |
| Neighbour (PRF) | 6.3 [2.6, 10.3] | 1.60 (3.4) | 5.7 % | 0.63 [0.24, 1.10] | 7.5 % |
| Neighbour (Gauss) | 7.3 [2.7, 9.9] | 1.85 (3.5) | 5.7 % | 0.75 [0.29, 1.07] | 0 % |

The "lead's statistic" is |centroid − neighbour| / err_px with no floor; the earlier claim used it at > 6.

**What this means:**
- By its own σ, centroid.py would have reported a **neighbour-hosted signal as "consistent with the target"** (< 3σ) **94 %
  of the time**, at 3–10".
- The floor-free statistic used for ">6σ" is wrongly calibrated in the other direction. It "excludes" the true neighbour
  host in up to 7.5 % of clean cases, and the real S66 value of 8.2 is a scattered-light artefact.
- **The real centroid.py numbers carry no host information.**

**Comparison with `../audit/`** (another investigator's centroid audit, read only after my numbers were final). Their
replication is byte-identical to mine. Their neighbour-injection σ (median 1.5–1.9) and their pedestal diagnosis agree.

---

## Corrections for CONTEXT.md and the record

1. The "S66 1.1" (0.3σ), S104 1.7" (0.3σ)" values are correct as outputs of centroid.py but carry no localisation power. The
   σ behind them is about 3.4" (S66) and about 6" (S104) per axis. **"9" neighbour excluded at > 6σ" is wrong**, both in its
   arithmetic and in its premise.
2. Quote the localisation significance as "target host rejected at ≳5σ by several partly independent pixel-level methods".
   Avoid precise ">8σ", "Δχ² = 129" or "10^-31.6" statements.
3. Use PSF photometry's f_target exclusion of 7.5σ, which includes systematics, rather than Δχ² = 129.
4. The pixel-level evidence now includes a PSF-free relative-depth test (this folder) and a zero-model flux inequality
   (S93, 4σ).

## Open issues

- **Unresolved sources.** A source within about 1" of the neighbour cannot be separated from it with TESS pixels (unchanged).
  Ground-based seeing-limited photometry at a predicted transit settles the host directly.
- **S13 and S66** are weak in the single-pixel test (noise in S13; a T = 14.7 star in the far-side pixel in S66).
- **Heavy tails.** The fake-transit null has heavy tails from scattered-light windows. Significances quoted in σ assume
  near-Gaussian stacks. Empirical tail probabilities are limited by about 16–26 independent fake windows per sector.
- **Model dependence of the bound.** The contamination bound for Q relies on Gaia fluxes and a broad PSF grid. The fully
  conservative envelope gives only 3.1σ (n/t).
- **TESScut − QLP time offset** of 5–54 s is unexplained, but irrelevant at this level.

## Files

**Code:**
- `rcommon.py`: I/O, Gaia, WCS and windows.
- `fetch_gaia.py`: VizieR Gaia DR3, saved to `gaia/`.
- `r1_wcs_check.py`, `r2_registration.py`, `r3_pixel_ratio.py`, `r3c_flux_inequality.py`, `r4_timing.py`, `r5_jitter.py`,
  `r6_stats.py`.
  - `r3_pixel_ratio.py` takes the environment variables `R3_CLEAN` = none/mad10 and `R3_FLOOR` = plane/median.
- `b1_centroid_audit.py`, with `centroid_copy.py` (a verbatim copy of `src/centroid.py`).
- `figures.py` and `make_results.py`.
- Diagnostics: `diag_null.py`, `diag_null2.py`.

**Results:**
- `results.json` (master summary).
- `r1_wcs_check.json`, `r2_registration.json`, `r4_timing.json`, `r5_jitter.json`, `r6_stats.json`.
- `r3_pixel_ratio_{none_plane,mad10_plane,mad10_median}.json` and `r3_nulls_*.npz`.
- `r3c_flux_inequality_{none,mad10}.json`.
- `b1_centroid_audit.json`.
- `log_*.txt`.

**Figures:** `fig_pixel_test.png` and `fig_centroidpy_injections.png`.

**PRF files:** `prf_files/` (SPOC start_s0004, cam 1, CCDs 1–3; copied from loc_independent and used only for the model
context).

**Reproduce** (about 2 min in total):
1. `fetch_gaia.py`
2. `r1`
3. `r3` (with `R3_CLEAN=none` and with `mad10`), then `r3c`
4. `r2`, `r4`, `r5`
5. `b1`
6. `r6`
7. `figures.py`
8. `make_results.py`
