# Review by openai/gpt-5.5

## Bottom line

The EB light curves mostly look like real astrophysical eclipsing/close-binary signals, especially the deep detached systems. But I would **not** phrase the result as “30 newly discovered eclipsing binaries” without qualification. A safer claim is:

> “We provide periods, ephemerides, morphology classifications, and vetting for 30 TESS eclipsing/close-binary candidates, most of which were previously only unvetted eclipse-like machine-learning flags.”

That distinction matters because **28/30 = 93%** are already in Kostov et al. table 2 as eclipse-like objects, and at least two have Gaia NSS binarity solutions. The planet-track conclusions are appropriately conservative overall: no planet is validated; TIC 155991187 and TIC 72651793 are likely false positives; TIC 245976814 remains an unvalidated transit-like candidate with significant blend risk.

My overall verdict: **(b) plausible but unproven** for “first vetted EB characterisation”; **(c) likely false positive / unvalidated** for the planet-like signals.

---

## 1. Are the 30 fair to call “first vetted confirmation and characterisation” of EBs?

### Mostly fair, but wording should be tightened

I would avoid “newly characterised eclipsing binaries” if the implication is discovery from scratch. Since **28 of 30 are already in Kostov+25 table 2**, they were previously flagged as eclipse-like. It is fairer to say:

- “first periods and ephemerides,” if no period was published before;
- “first manual vetting/characterisation,” if the Kostov entries are indeed unvetted and periodless;
- “newly confirmed as TESS eclipsing/close binaries,” with the caveat that most were previously machine-flagged.

I would not call all 30 “confirmed” in the strong astrophysical sense unless you have either resolved photometry, spectroscopy, or very robust centroid/difference-image evidence for each. TESS photometry plus phase coherence is enough for “very likely EB,” especially for the deep systems, but less secure for shallow/contact/ellipsoidal cases.

### Arithmetic/checks

- Total searched unique stars:  
  **77,153 + 64,121 = 141,274**.
- Known-signal sanity check:  
  **327 / 1,420 = 23.0%** matched known TOIs/TCEs/planets/EBs-with-TCEs.
- Prior Kostov overlap:  
  **28 / 30 = 93.3%** already eclipse-like flagged. This weakens “new discovery” but not necessarily “first vetted characterisation.”

### Objects I would downgrade or flag

#### Contact/ellipsoidal systems: keep separate from detached EBs

These three should be reported as “likely close/contact binaries or ellipsoidal variables,” not as clean detached EBs:

- **TIC 32398505**, P = 0.601493 d  
- **TIC 114787953**, P ≈ 0.6057 d  
- **TIC 71315603**, P = 0.507699 d  

They are likely real close-binary variables, but the distinction matters: ellipsoidal variables need not be eclipsing. For the shallowest one, TIC 71315603, the primary depth is only **0.6%**, and the morphology is closer to low-amplitude contact/ellipsoidal modulation than a detached eclipse. I would not mix these unqualified with detached EBs.

#### High-RUWE / bad-parallax systems: do not remove, but avoid physical interpretation

High RUWE is not a reason to reject an EB; it often supports multiplicity. But it means the Gaia single-star solution is poor, so distances, radii, luminosities, and companion inferences should be treated carefully.

Notable cases:

- **TIC 146708570**: RUWE 15, Gaia RV amplitude 121 km/s. Strongly supports binarity. Keep.
- **TIC 122543151**: RUWE 4.4, eccentric; possible tertiary. Keep, but note triple possibility.
- **TIC 89040715**: RUWE 8.2. Keep, possible tertiary/blend.
- **TIC 217648205**: RUWE 18, eccentric, possible pulsating primary. Keep only if eclipse centroids and secondary are robust.
- **TIC 129032726**: RUWE 19.5. Keep, but physical parameters suspect.
- **TIC 115242114**: RUWE 46 and parallax 0.16 mas. This one needs a strong caveat: Gaia astrometry is essentially unreliable. The EB signal may be real, but distance/radius/absolute-luminosity conclusions should not be used.

#### Gaia NSS systems: not new binaries

- **TIC 143345010**
- **TIC 278802197**

These already have Gaia NSS evidence for binarity. It is fair to say “eclipses newly characterised,” but not “binary newly discovered.”

#### Period/alias concerns

A few systems need period-convention clarification.

- **TIC 100666720**: note says “period 2× search value,” secondary at phase 0.117. That is unusual enough that I would explicitly show the adopted ephemeris and explain why the orbital period is 7.341134 d rather than the BLS period.
- **TIC 49617965**: table gives P = 7.53040 d, secondary at phase 0.750, “equal eclipses, 3.765 d spacing.” But the figure appears to show **P ≈ 15.0608 d** for this TIC. This is a real inconsistency. If equal eclipses occur every 3.765 d, the photometric period could be 3.765 d, 7.530 d, or 15.061 d depending on which events are primary/secondary/tertiary-like. This needs to be fixed before publication.
- Contact/ellipsoidal systems near 0.5–0.6 d are especially vulnerable to half-period/double-period ambiguity. For contact binaries the listed period should be explicitly defined as orbital period or photometric period.

#### Figure/table mismatch

The supplied figures include objects not in the 30-object table, for example:

- TIC 166225041
- TIC 237329004
- TIC 387291911
- TIC 60829240
- TIC 422003560

Some look sinusoidal/ellipsoidal rather than detached. If these were rejected objects, they should not appear in the final figure set without explanation. If they are part of the final sample, the table is incomplete. This mismatch weakens confidence in the final bookkeeping.

---

## 2. Are the planet-track conclusions appropriately cautious?

Mostly yes. The key conclusion — **no planet can be validated** — is appropriate and should remain the headline.

### TIC 155991187

The cautious conclusion is justified. Important numbers:

- P = 0.5012110 d
- depth ≈ 2.9 ppt
- six sectors over 2018–2026
- RUWE 21
- co-moving M dwarf 8.2″ away
- no secondary/ellipsoidal signal to <80 ppm at 3σ
- TRICERATOPS FPP = 0.93  
  - SEBx2P = 0.60  
  - STP = 0.25  
  - TP = 0.07

This is not a validated planet. The statement “most likely an eclipsing binary within the system” is reasonable, though I would phrase it slightly more broadly as “most likely a stellar false positive or planet/EB scenario involving another component of the system.”

The centroid comparison is useful but not decisive. For the 8.2″ companion:

- χ² = 21.8 for 10 dof gives p ≈ 0.016.
- Target model χ² = 10.9 for 10 dof gives p ≈ 0.37.
- Δχ² = 10.9, which corresponds to a likelihood ratio of exp(−10.9/2) ≈ 0.004 if the assumptions are ideal.

So the companion is disfavoured at roughly the few-sigma level, but not eliminated, especially given TESS systematics, crowding, and unresolved/bound companions implied by RUWE 21.

### TIC 72651793

This is appropriately treated as non-planetary/unvalidated.

- P = 0.6345 d
- depth = 410 ppm
- three independent sectors
- odd/even difference = 3.5σ
- TRICERATOPS FPP = 0.86

A 3.5σ odd/even difference is a strong EB warning, though not absolutely conclusive by itself. The “density mismatch” is supportive but should not be overemphasised; transit-derived density can be misleading for grazing events, eccentric orbits, blends, or imperfect detrending.

### TIC 245976814

This one should be described as a **transit-like candidate requiring follow-up**, not as a probable planet.

Key red flags:

- P = 0.5070757 d, near the common contact-binary/ultra-short-period false-positive regime.
- Depth = 1.1–1.5 ppt.
- F/G star, T = 11.6, R = 1.29 R☉. If planetary and undiluted, radius would be approximately:  
  \[
  R_p/R_\star \approx \sqrt{0.0011\text{–}0.0015} = 0.033\text{–}0.039
  \]  
  \[
  R_p \approx 0.043\text{–}0.050 R_\odot \approx 4.7\text{–}5.5 R_\oplus
  \]  
  That would be a very hot Neptune/sub-Saturn on a 0.507 d orbit — possible, but rare.
- S104 centroid is on target at 1.1σ, but S66 is offset by 20″ at 2.8σ. A 20″ offset is about one TESS pixel, so this is a substantial warning.
- Nearby Gaia stars at 4.7–6.2″ and T ≈ 17.6–19 can plausibly mimic the observed depth if they are deep eclipsing binaries. For a contaminant ΔT = 6.0 mag fainter than the target, the flux ratio is:
  \[
  10^{-0.4 \times 6.0} \approx 0.004
  \]
  A 100% eclipse on that star would produce a diluted signal of about **4 ppt**, enough to explain a 1.1–1.5 ppt event. Even ΔT ≈ 7.4 gives:
  \[
  10^{-0.4 \times 7.4} \approx 0.0011
  \]
  which is still comparable to the observed signal if the eclipse is deep.

So TIC 245976814 is not validated and has a credible nearby-EB explanation. The most decisive next test would be ground-based time-series photometry with seeing good enough to separate the 4.7–6.2″ neighbours.

---

## Additional tests that would most change my assessment

For the EB sample:

1. **Per-sector difference images with actual centroid offsets in arcsec**, not only “within 2.5σ.” A 2.5σ agreement can still allow a large physical offset if the centroid uncertainty is large.
2. **Gaia-neighbour contamination analysis** for each object: list all Gaia sources within at least 1–2 TESS pixels, their ΔT/ΔG, and whether a deep eclipse on them could reproduce the observed depth.
3. **Independent extraction check** using different apertures/pipelines: QLP vs eleanor vs custom FFI photometry, especially for shallow/contact objects.
4. **Resolve the period aliases**, especially TIC 49617965 and TIC 100666720.
5. **Catalogue reproducibility table**: for each TIC, list exact nearest entries in VSX/ASAS-SN/Gaia/ATLAS/etc., angular separations, and periods if any. This would make the “first vetted” claim much more defensible.

For the planet-like objects:

1. **Ground-based on/off-target photometry** for TIC 245976814 and TIC 155991187.
2. **High-resolution imaging** or Gaia DR3/DR4 companion checks for TIC 155991187, given RUWE 21.
3. **Reconnaissance spectroscopy/RVs** for TIC 155991187 if a planet interpretation is still being considered. A 0.5 d EB or hierarchical binary should often reveal large RV or line-profile changes.

---

## Recommended wording changes

Instead of:

> “30 newly characterised eclipsing binaries”

I would use:

> “30 TESS eclipsing/close-binary systems with newly derived periods, ephemerides, and vetting; 28 were previously only unvetted eclipse-like machine-learning flags.”

Instead of:

> “confirmed eclipsing binaries”

I would use:

> “vetted likely eclipsing binaries,” or “confirmed photometric eclipsing/close-binary variables,” with contact/ellipsoidal cases separated.

Instead of:

> “TIC 155991187 is most likely an eclipsing binary within the system”

I would say:

> “TIC 155991187 is a persistent transit-like signal, but the high FPP, RUWE 21, and bound companion make a stellar false-positive or hierarchical-system scenario more likely than a planet.”

---

## Non-astronomer summary

We did not find a new confirmed planet. The search did, however, identify and measure 30 stars whose brightness repeatedly drops in a way consistent with eclipsing or very close binary stars. Most of these were already flagged by a previous machine-learning catalogue as possible eclipse-like objects, but this work appears to provide the first detailed periods, eclipse depths, and vetting for them. A few planet-like signals were investigated, but each has warning signs — such as nearby contaminating stars, odd/even transit differences, or high false-positive probabilities — so none should be announced as a planet without further observations.
