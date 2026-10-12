# Review by openai/gpt-5.5

## Bottom line

- **Candidate A (TIC 87055304):** A real, repeatable transit-like signal is plausible, and **“new planet candidate (unvalidated)” is a defensible label only with strong caveats**. It is **not close to validation**. The main red flags are the **very poor Gaia astrometry/RUWE = 18.8**, highly inconsistent stellar parameters, **transit-inferred low stellar density**, likely **large companion radius if the low density is real**, and sector-to-sector depth changes. Most likely non-planet explanation: **a bound/unresolved eclipsing binary or hierarchical triple**, possibly involving the source responsible for the RUWE excess.  
  **Verdict for A: (b) plausible but unproven.**

- **Candidate B (TIC 320686138):** I agree with your downgrade. The signal is real-looking, but the **near-grazing geometry, depth variability, TRICERATOPS FPP = 0.98, and acceptable grazing-EB fits** make “new planet candidate” too generous.  
  **Verdict for B: (c) likely false positive / grazing EB.**

- **EB additions:** The 23 additions are plausibly EBs/contact binaries, and the count 30 → 53 is arithmetically correct. But I would not call all “new” without an explicit cross-match to VSX/ASAS-SN/ZTF/TESS EB catalogs. Also, your text has an internal inconsistency: the combined search is said to have found **1 EB**, but the EB table flags **six** as “found by the combined two-sector search,” including one with **P = 3.337 d**, outside the stated 5–25 d combined-search period range.

---

## Candidate A: TIC 87055304

### What is solid

The signal is likely astrophysical or at least repeatable:

- Period over 2019–2026: **P = 4.136195 d**, about **25 transits** in 5 sectors.
- Odd/even depths:  
  - odd = 3529 ± 172 ppm  
  - even = 3624 ± 169 ppm  
  - difference = 95 ± 241 ppm → **0.4σ**, so no meaningful odd/even evidence.
- Secondary at phase 0.5: **109 ± 203 ppm**, insignificant.
- Difference images in S66/S104 are consistent with the target and exclude the 9″ T≈13.5 neighbor at high significance.
- The formal ellipsoidal/harmonic term, cos 2φ = −267 ± 32 ppm, would be alarming if taken literally, but your empirical comparison against 40 unrelated periods gives rms = 209 ppm, so the real significance is only **~1.3σ**. Good correction.

### Arithmetic checks

The circular transit density is consistent with your number:

\[
\rho_\star/\rho_\odot \simeq \left(\frac{a/R_\star}{4.21}\right)^3 P_{\rm d}^{-2}
\]

For A:

\[
(a/R_\star)=4.30,\quad P=4.1362\,{\rm d}
\]

\[
\rho_\star/\rho_\odot \simeq (4.30/4.21)^3 / 4.136^2 \simeq 0.062
\]

So your quoted **0.063 ρ⊙** is right.

But that density is a major issue. It implies, approximately,

\[
R_\star \simeq (M_\star/\rho_\star)^{1/3}
\]

For plausible masses:

- if \(M_\star = 0.95\,M_\odot\): \(R_\star \approx 2.47\,R_\odot\)
- if \(M_\star = 1.3\,M_\odot\): \(R_\star \approx 2.74\,R_\odot\)

That is much larger than the FLAME radius and larger than the TIC radius. It is more consistent with an evolved/subgiant star, but the Gaia astrometry is bad enough that none of the stellar parameters are secure.

The companion radius arithmetic is also important:

- Undiluted \(R_p/R_\star \approx 0.067\)
- Dilution-corrected \(R_p/R_\star \approx 0.082\), if the TIC contamination ratio corresponds to contaminating flux/target flux ≈ 0.51.

Then:

- For \(R_\star = 1.8\,R_\odot\):  
  \(R_p \approx 1.2\)–\(1.4\,R_{\rm Jup}\), as you state.
- For \(R_\star = 2.7\,R_\odot\):  
  \(R_p \approx 1.75\)–\(2.15\,R_{\rm Jup}\).

So if the transit-derived density is correct and the star is ~2.5–2.8 R⊙, the object is **very large for a planet**. A 1.8–2.2 RJ object at 4.1 d is not impossible, but it pushes strongly toward a brown dwarf / low-mass star / blend interpretation unless RVs say otherwise.

### Main red flags

#### 1. RUWE = 18.8 is severe

RUWE = 18.8 is not a small astrometric warning; it is extreme. It strongly suggests unresolved multiplicity, bad astrometric modeling, or a problematic Gaia solution. That directly attacks the clean single-star planet interpretation.

The combination of:

- RUWE = 18.8,
- inconsistent TIC/Gaia/FLAME radii,
- transit-inferred \(\rho_\star = 0.063\,\rho_\odot\),
- TIC contamination ratio = 0.51,

means the host-star identity and radius are not securely known.

#### 2. The transit-inferred density and catalog parameters do not agree

Using TIC \(R_\star = 1.79\,R_\odot\) and \(M_\star \sim 1.3\,M_\odot\), the stellar density would be roughly:

\[
\rho_\star \sim 1.3/1.79^3 \simeq 0.23\,\rho_\odot
\]

That is about **3.6× larger** than the transit density. Using FLAME \(R_\star = 1.24\,R_\odot\), \(M_\star = 0.95\,M_\odot\):

\[
\rho_\star \sim 0.95/1.24^3 \simeq 0.50\,\rho_\odot
\]

That is about **8× larger** than the transit density.

So either:

- the star is actually evolved and larger, around 2.5–2.8 R⊙;
- the orbit is eccentric and the transit occurs near apastron;
- the signal is not on the assumed star;
- the transit model is biased by blending/systematics/grazing geometry;
- or the system is a multiple.

For a 4.1 d hot Jupiter, a large eccentricity solution is possible but not the default, because tidal circularization should be non-negligible.

#### 3. Sector-to-sector depth changes are non-trivial

Per-sector SAP depths:

- S13: 4101 ppm
- S39: 4238 ppm
- S66: 3792 ppm
- S93: 3629 ppm
- S104: 2939 ppm

The mean is about **3740 ppm**, and the sample rms is about **510 ppm**. If the combined depth uncertainty is 120 ppm, then the sector scatter is much larger than white-noise expectations. Even if per-sector uncertainties are roughly \(120\sqrt{5} \approx 270\) ppm, the depth scatter gives approximately:

\[
\chi^2 \sim 14 \text{ for 4 dof}
\]

which is mildly significant.

This could be aperture/dilution differences between sectors, but it is also what one worries about for blends, changing contamination, or imperfect detrending. It does not kill the candidate, but it should be explicitly treated as a systematic uncertainty. The quoted combined depth uncertainty of ±120 ppm is probably too optimistic if the sector scatter is astrophysical or aperture-dependent.

#### 4. TRICERATOPS FPP = 0.037 is not validation-level

An FPP of 3.7% is not negligible. Also, the run used:

- QLP DET light curve affected by a known detrending artifact,
- TIC stellar assumptions,
- \(M_\star = 1.3\,M_\odot\),
- not the pending large-radius SAP rerun.

If the star is actually ~2.5–2.8 R⊙ and the companion is ~2 RJ, the planet prior should drop and the false-positive/brown-dwarf/stellar-companion prior should rise. I would not lean heavily on the current TRICERATOPS result.

### Contamination / neighbor arithmetic

The 9″ T≈13.5 neighbor is only ΔT ≈ 1.8 mag fainter, so its flux ratio is roughly:

\[
10^{-0.4 \times 1.8} \approx 0.19
\]

If it were the source, the required eclipse depth on that neighbor would be only a few percent, i.e. astrophysically easy. Your difference imaging exclusion is therefore important. If the centroid exclusion is really >6σ, that largely removes this object.

The T≈17.7–17.9 stars at 4–5″ are much fainter. ΔT ≈ 6 mag gives flux ratios of order 0.004, so to produce a 3580 ppm aperture signal they need enormous eclipses, roughly **80–90%**, consistent with your statement. That is unlikely but not impossible for an EB. Also, because they are only 4–5″ away, the TESS difference images probably do **not** exclude them strongly by centroid position alone.

### Most likely non-planet explanation

If not a planet, my leading explanation is:

> **A bound or unresolved eclipsing binary / hierarchical triple associated with the RUWE-excess system.**

A background EB on one of the very faint 4–5″ sources is possible but less likely because of the required near-total eclipse. The 9″ neighbor would be an easy contaminant astrophysically, but your centroid result argues against it.

An on-target M-dwarf companion would likely produce RV semi-amplitudes of order km/s to >10 km/s. For example, a 0.1 M⊙ companion at P = 4.136 d around a ~1.3 M⊙ star gives \(K\sim 10\) km/s. Gaia’s non-variable RVs weaken that scenario, but do not eliminate hierarchical/blended configurations.

### Is “new planet candidate, unvalidated” fair?

**Yes, but barely, and only if stated as an unvalidated giant-companion candidate.** I would avoid language implying a high-confidence planet. The cleanest label would be something like:

> “TIC 87055304: unvalidated transiting giant-companion candidate; planet interpretation plausible but limited by RUWE = 18.8, uncertain stellar radius, and possible unresolved multiplicity.”

### Single follow-up that would most change my view

For A, the highest-impact follow-up is:

> **Multi-epoch high-resolution spectroscopy/RVs at the 4.136 d quadratures, with line-bisector/CCF checks.**

A few spectra should distinguish:

- hot Jupiter: tens to few hundreds m/s;
- brown dwarf / very low-mass star: km/s;
- stellar EB: many km/s;
- blended/hierarchical EB: line asymmetries, secondary CCF, or no RV on the bright star.

Speckle/AO imaging would also be valuable because of RUWE = 18.8, but RV+bisectors would more directly attack the planet-vs-EB question.

---

## Candidate B: TIC 320686138

I agree with your downgrade. This should not be claimed as a planet candidate except perhaps in a “rejected/likely EB” appendix.

### What is solid

The transit-like signal is real enough:

- P = 11.91415 d over 2019–2026.
- 9 observed transits.
- Combined depth = 3782 ± 144 ppm.
- Null \(\chi^2 = 1763\), transit model \(\chi^2 \approx 845\), so Δχ² ≈ 918 — a strong detection.
- Difference imaging is consistent with the target and excludes the obvious 19.1″ T≈15.1 contaminant.
- No significant phase-0.5 secondary by your quoted number: 4 ± 160 ppm.

But the astrophysical interpretation is poor for a planet.

### Arithmetic checks

The stellar density/a/R arithmetic is consistent.

Using TIC/FLAME-like values:

\[
M_\star = 1.12\,M_\odot,\quad R_\star = 1.37\,R_\odot
\]

\[
\rho_\star/\rho_\odot = 1.12/1.37^3 \simeq 0.44
\]

For P = 11.914 d:

\[
a/R_\star \simeq 4.21\,P^{2/3}\rho^{1/3}
\]

\[
P^{2/3} \simeq 5.22,\quad \rho^{1/3}\simeq 0.76
\]

\[
a/R_\star \simeq 4.21 \times 5.22 \times 0.76 \simeq 16.7
\]

So your fixed \(a/R_\star = 16.6\) is right.

Central duration:

\[
T_{\rm dur,central} \sim \frac{P}{\pi}\frac{1+k}{a/R_\star}
\]

For \(k \sim 0.073\):

\[
T \sim \frac{11.914}{\pi}\frac{1.073}{16.6}
\approx 0.245\,{\rm d} \approx 5.9\,{\rm h}
\]

So your statement that a central transit would last ~5.8 h is correct. The observed duration is only **2.2 h**, so the event must be high-impact / near-grazing if on the target.

### Main red flags

#### 1. The geometry is grazing or nearly grazing

A non-grazing planet solution requires \(b\) very close to 1. Your own fits show the problem:

- free fit: \(R_p/R_\star = 0.34\), \(b = 1.26\), \(\chi^2 = 845\)
- fixed \(a/R_\star = 16.6\): \(R_p/R_\star = 0.073\), \(b = 0.93\), \(\chi^2 = 851\)
- forced \(b \le 0.9\): \(\chi^2 = 863\)

The free fit preferring \(R_p/R_\star = 0.34\) is exactly the kind of solution expected for a **grazing stellar eclipsing binary**: a large companion partially clips the primary and mimics a shallow transit.

If \(R_\star = 1.37\,R_\odot\), then \(R_p/R_\star = 0.34\) corresponds to:

\[
R_2 \approx 0.47\,R_\odot
\]

That is an M-dwarf-sized object, not a planet.

#### 2. Depth variability is significant

Per-transit depths range from **2.8 to 5.4 ppt**.

You quote:

\[
\chi^2 = 35 \text{ for 8 dof}
\]

That gives an error inflation factor:

\[
\sqrt{35/8} \simeq 2.1
\]

So your “errors underestimated about 2×” statement is correct.

Depth variability can come from systematics/aperture changes, but for a near-grazing event it is also a classic sign that the inferred radius and impact parameter are unstable and the planet interpretation is fragile.

#### 3. Odd/even evidence is weak after rescaling, but not comforting

Your formal weighted odd-even difference:

\[
-782 \pm 295\,{\rm ppm}
\]

is 2.65σ formally. After multiplying errors by ~2.1, it becomes about **1.3σ**, so not significant.

But the fact that TRICERATOPS strongly favors **SEBx2P** means the actual binary period may be twice the detected period. The absence of a strong odd/even difference does not rescue the planet interpretation.

#### 4. TRICERATOPS is strongly against the planet interpretation

FPP = 0.98, with:

- SEBx2P = 0.86
- EBx2P = 0.10
- BEBx2P = 0.013
- PTP = 0.010
- TP = 0.007

That is a decisive model-prior statement: under TRICERATOPS assumptions, a planet is a ~1% possibility, while EB-like explanations dominate.

### Plot concern

In the supplied Candidate B plot, the “phase 0.5 secondary” panel appears to show a negative excursion of order ~1–2 ppt over several points away from exact center. If that is a single-sector artifact, fine. But if it recurs at fixed phase, it would be additional EB evidence. I would explicitly document whether that feature is one-sector/systematic or coherent.

### Most likely non-planet explanation

Most likely:

> **A grazing eclipsing binary, probably at twice the detected period, either on the target or in an unresolved/bound component.**

A background EB is less likely because the obvious Gaia neighbor is centroid-excluded, but unresolved companions remain possible.

### Is “new planet candidate, unvalidated” fair?

**No.** Given FPP = 0.98 and the grazing-EB morphology, I would not label B a planet candidate. Your downgrade to “most likely grazing EB” is appropriate.

### Single follow-up that would most change my view

For B:

> **Reconnaissance high-resolution RVs at the predicted quadratures for P and 2P.**

A stellar EB should show km/s to tens-of-km/s RV motion or CCF/bisector anomalies. A genuine Jovian planet should be at the tens-to-hundreds m/s level. Even modest-precision spectroscopy should be decisive.

---

## EB additions

### Count check

The table contains **23** entries. If the previous list was 30, then:

\[
30 + 23 = 53
\]

So the “list is now 53” arithmetic is correct.

### Internal inconsistency: combined-search yield

You state:

> Combined search: 4 planet-like signals not recovered in other sectors; 1 EB.

But the EB table marks the following as “found by the combined two-sector search”:

1. TIC 100378891  
2. TIC 234968334  
3. TIC 139256801  
4. TIC 404124262  
5. TIC 80453626  
6. TIC 281740175  

That is **six**, not one.

Also, TIC 80453626 has **P = 3.337312 d**, while the combined search was described as covering periods **5–25 d**. That needs explanation. Either the “found by combined search” notes are inherited from another run, the combined search period range was broader than stated, or the summary “1 EB” is wrong.

### EB classification generally looks plausible

Most of the listed objects are obviously not planet candidates:

- primary depths of 4–44%;
- secondaries at percent level;
- Gaia RV amplitudes of 100–175 km/s in several cases;
- high RUWE in many systems;
- contact/ellipsoidal modulation for the short-period systems.

The period doubling for TIC 406279367 is also sensible if odd/even eclipses differ in all 8 sectors.

### TIC 6598941 EBLM-like classification

Your arithmetic is broadly right.

Primary eclipse depth = 1.7%. If the companion is dark and the eclipse is central:

\[
R_2/R_1 \approx \sqrt{0.017} \approx 0.130
\]

For \(R_1 = 1.86\,R_\odot\):

\[
R_2 \approx 0.24\,R_\odot
\]

So “companion ~0.25 R⊙” is reasonable.

A 0.25 M⊙ companion at P = 2.12 d should produce a primary RV semi-amplitude of order ~30 km/s, i.e. a peak-to-peak range around 60 km/s. The Gaia RV range/amplitude of 62 km/s is therefore exactly in the stellar-companion regime. Calling this EBLM-like rather than a planet is correct.

### Novelty of EBs

The EB identifications look credible, but “new” needs care. For EBs, the relevant missed catalogs are not only TOI/CTOI/SIMBAD/VizieR near the planet candidates, but also:

- VSX,
- ASAS-SN variable star catalog,
- ZTF variable catalogs where available,
- Gaia variability tables,
- TESS EB catalogs,
- TIC contamination/nearby-source notes,
- possibly ATLAS/OGLE depending on sky position.

I would phrase them as **“53 EB/eclipsing-source detections in our search, 23 newly added to our list”** unless you have done a full external novelty cross-match.

---

## Final verdict

- **Candidate A:** Real-looking, repeatable transit-like signal; possible giant planet, but the host-star ambiguity and RUWE = 18.8 prevent any strong planet claim. **Plausible but unproven.**
- **Candidate B:** Agree with downgrade. **Likely grazing EB / false positive.**
- **EB list:** Mostly credible as eclipsing binaries, but the claimed combined-search accounting needs correction, and “new” should be used cautiously.

Overall: **not a validated new-planet discovery. At most, one plausible unvalidated giant-planet/giant-companion candidate plus one likely EB.**
