# Review by openai/gpt-5.5

## Overall assessment

The withdrawal of the TIC 87055304 planet candidate is well supported. I do **not** see a credible way for TIC 87055304 to host the full 4.136 d signal given the pixel-level evidence shown. The conclusion that the signal is on/very near TIC 87055315 is much stronger than the original planet claim.

The interpretation of TIC 87055315 as an eclipsing F/G + M-dwarf system is also plausible and probably correct, but it is not yet as directly confirmed as the host-localisation result. Without RVs or resolved ground-based eclipse photometry, I would phrase it as “a likely EBLM-type eclipsing binary on TIC 87055315 or an unresolved blend very close to it,” not yet a fully confirmed EBLM.

**Verdict:** for the planet candidate on TIC 87055304: **(c) likely false positive**.  
For the new EB interpretation: **plausible and probably right, but still awaiting standard TFOP-style confirmation.**

---

## 1. Could TIC 87055304 still host the signal?

Only in a very contrived way. The evidence against the target hosting the signal is strong and internally consistent.

### Flux/separation arithmetic

- Separation: 9.33″.
- TESS pixel scale ≈ 21″ px⁻¹, so separation ≈ 0.44–0.45 px. Correct.
- TESS magnitudes: T = 11.72 vs 13.54, ΔT = 1.82 mag.
- Flux ratio neighbour/target:

\[
10^{-0.4 \times 1.82} = 0.187
\]

The quoted 0.188 is correct.

Both stars are deeply blended in TESS, so ordinary aperture photometry cannot identify the host. Pixel-level localisation is decisive.

### PRF / difference-image localisation

The difference images strongly favour the neighbour. The shown residuals are qualitatively persuasive: fitting the transit source at TIC 87055304 leaves structured residuals; fitting it at TIC 87055315 largely removes them.

The reported joint source location:

- Offset from target: +6.37″ E, +6.08″ N.
- Distance from target:

\[
\sqrt{6.37^2 + 6.08^2} = 8.80″
\]

- Distance from neighbour: quoted 0.66″.

That is exactly the expected geometry for a source on TIC 87055315.

The sector-by-sector position fractions along the target–neighbour axis are 1.37, 0.99, 0.84, 1.03, 0.83. These are all on the neighbour side. The scatter is not negligible, but the sign and direction are consistent.

The audit point is important: to fake this with registration error would require a roughly 0.38 px systematic shift, i.e. about 8″, aligned in the right direction across multiple sectors. The quoted WCS/registration errors are ≤0.09 px along that axis. That is a large gap.

### Simultaneous PSF light curves

The target PSF light curve has a measured dip of only

\[
0.64 \pm 0.36 \ {\rm ppt}
\]

whereas if TIC 87055304 hosted the signal, the expected dip in its own PSF flux would be roughly 4 ppt. The difference is about 3.5 ppt, or nearly 10σ using the quoted statistical uncertainty alone. Systematics in deblending two stars separated by only 0.45 px are real, but not large enough to hide a full 4 ppt target-hosted event if the injection tests are trustworthy.

The reported target fraction,

\[
f = 0.16 \ [0.07, 0.25] \ {\rm statistical}, \ \pm 0.075 \ {\rm systematic},
\]

is also strongly inconsistent with \(f=1\). If I approximate the statistical uncertainty as ~0.09 and combine with 0.075 in quadrature,

\[
\sigma_f \simeq \sqrt{0.09^2 + 0.075^2} \simeq 0.117
\]

then

\[
(1 - 0.16)/0.117 \simeq 7.2\sigma.
\]

So the quoted “7.5σ” is arithmetically reasonable, though probably too formal given shared modelling assumptions.

### Model-free pixel test

The model-free test is useful because it does not rely on the detailed TESS PRF. The observed neighbour-side/far-side depth ratios of 3.8, 4.3, and 2.2 in the informative sectors are hard to reconcile with a target-hosted signal, for which ratios of order 0.7–1.3 are expected.

The conservative version gives only 3.1σ, while the nominal gives 6.4σ. I would not advertise this alone as decisive, but as an independent sanity check it supports the same conclusion.

### Bottom line on target hosting

For TIC 87055304 still to host the full signal, one would need all of the following:

1. PRF difference imaging biased by nearly the full 0.45-pixel target–neighbour separation.
2. The bias to persist across multiple sectors and spacecraft orientations.
3. The simultaneous PSF light-curve extraction to move most of the target’s dip into the neighbour’s light curve.
4. The model-free pixel-depth asymmetry to point the wrong way.
5. Injection tests to fail in the same direction.

That is not impossible in an abstract sense, but it is not a credible remaining explanation. The target-hosted planet should be withdrawn.

---

## 2. Could the TIC 87055315 object still be a planet?

If the eclipse is on TIC 87055315 itself, no. The radius arithmetic rules out a planet.

### Companion radius arithmetic

Quoted host radius:

\[
R_1 = 1.93 \pm 0.09 \ R_\odot
\]

Quoted radius ratio:

\[
k = R_2/R_1 = 0.1394 \pm 0.0031
\]

Then

\[
R_2 = 1.93 \times 0.1394 = 0.269 \ R_\odot.
\]

Correct.

Converting to Jupiter radii using \(1 R_\odot \simeq 9.73 R_{\rm Jup}\):

\[
0.269 R_\odot \times 9.73 = 2.62 R_{\rm Jup}.
\]

So the quoted \(2.61 \pm 0.13 R_{\rm Jup}\) is correct.

The relative uncertainty is approximately

\[
\sqrt{(0.09/1.93)^2 + (0.0031/0.1394)^2}
= \sqrt{0.0466^2 + 0.0222^2}
= 0.052,
\]

so

\[
0.052 \times 0.269 R_\odot = 0.014 R_\odot,
\]

consistent with the quoted ±0.013 \(R_\odot\).

A radius of 2.6 \(R_{\rm Jup}\) is not a plausible planet or brown dwarf radius, especially for a several-Gyr F/G turn-off host. Brown dwarfs are roughly Jupiter-sized. Even the most inflated hot Jupiters do not reach 2.6 \(R_{\rm Jup}\), and this system’s equilibrium temperature is likely only of order 1600 K, not extreme enough to rescue such a radius.

### How wrong would the host radius need to be?

For the companion to be ≤2.0 \(R_{\rm Jup}\), i.e. ≤0.206 \(R_\odot\), with \(k = 0.1394\),

\[
R_1 < 0.206/0.1394 = 1.48 R_\odot.
\]

That is ~0.45 \(R_\odot\) below the quoted 1.93 \(R_\odot\), about 5σ using the quoted 0.09 \(R_\odot\) uncertainty, and inconsistent with the Gaia/SED/isochrone picture. For a normal Jupiter-sized planet, the required host radius would be completely impossible.

### Effect of partial target contribution

The report says if 16% of the dip were on TIC 87055304, the TIC 87055315 companion radius would become 2.40 \(R_{\rm Jup}\). That scaling is reasonable:

\[
2.61 \times \sqrt{0.84} = 2.39 R_{\rm Jup}.
\]

Still too large for a planet or brown dwarf.

### Unresolved blend near TIC 87055315

The one caveat is not a planet on TIC 87055315 itself, but an unresolved star within ~1″ of TIC 87055315 hosting the eclipse. The report acknowledges this. However, that does not easily make the object planetary.

The observed missing light is about 2.14% of the catalogued TIC 87055315 flux. If an unresolved subcomponent contributes a fraction \(f\) of TIC 87055315’s total flux, the intrinsic eclipse depth is

\[
\delta_{\rm intrinsic} \simeq 0.0214/f.
\]

To keep the intrinsic depth planet-like, say below a few percent, the unresolved host must contribute a large fraction of the light. For example:

- If \(f = 1\), \(\delta = 2.14\%\).
- If \(f = 0.6\), \(\delta = 3.57\%\).
- If \(f = 0.3\), \(\delta = 7.13\%\).
- If \(f = 0.1\), \(\delta = 21.4\%\).

A faint hidden star cannot host a shallow planet-like transit; the eclipse would be very deep. A bright hidden star would have to be a substantial contributor to the Gaia source and would still imply a large eclipsing object unless the stellar parameters are very unusual. The report’s claim that an unresolved host would need to supply ≥60% of the light is directionally right.

Still, TFOP would ask for high-resolution imaging or reconnaissance spectroscopy to rule out unresolved blends if the paper wants to claim a specific EBLM architecture rather than simply “not on TIC 87055304.”

---

## 3. Arithmetic and consistency checks

### Correct / credible arithmetic

- Flux ratio from ΔT = 1.82 mag: correct, ~0.188.
- Separation in TESS pixels: 9.33″ / 21″ ≈ 0.44 px, correct.
- Companion radius: \(0.269 R_\odot = 2.61 R_{\rm Jup}\), correct.
- Radius uncertainty: ±0.013 \(R_\odot\), plausible.
- RV semi-amplitude estimate: quoted \(K \approx 27\) km/s for a ~0.25 \(M_\odot\) companion is plausible.

Using

\[
K_1 \simeq 213.9 \ {\rm km/s}
\left(\frac{M_2}{M_\odot}\right)
P_{\rm d}^{-1/3}
(M_1+M_2)^{-2/3},
\]

with \(M_1 = 1.17 M_\odot\), \(M_2 = 0.25 M_\odot\), \(P=4.136\) d:

\[
K_1 \approx 26 \ {\rm km/s}.
\]

So the proposed RV confirmation is well chosen.

### Minor inconsistency: number of transits

The report alternately says:

- 25 transits in the summary,
- 24 transits in the timing table,
- 26 transits in the difference-image section.

This needs clarification. It may simply reflect different quality cuts — e.g. 26 imaged events, 24 with reliable timing, 25 used in one signal count — but as written it looks inconsistent. A TFOP reviewer will ask for a per-sector event table with which transits were used for timing, localisation, and modelling.

### Nominal duration / stellar density tension

The neighbour-host solution quotes:

- \(R_1 = 1.93 R_\odot\),
- \(M_1 = 1.17 M_\odot\),
- \(\rho_\star = 0.172 \rho_\odot\),
- \(b = 0.49\),
- \(k = 0.139\),
- \(T_{14} = 4.35\) h.

For a circular orbit with \(P=4.136\) d, \(M_{\rm tot} \sim 1.4 M_\odot\), and \(R_1=1.93 R_\odot\), the semi-major axis is roughly

\[
a \simeq 12.2 R_\odot,
\]

so

\[
a/R_1 \simeq 6.3.
\]

With \(b=0.49\), \(k=0.139\), the expected duration is closer to ~5.1 h than 4.35 h. Conversely, \(T_{14}=4.35\) h and \(b=0.49\) imply \(a/R_\star \sim 7.5\), corresponding to \(\rho_\star \sim 0.29 \rho_\odot\), not 0.17 \(\rho_\odot\).

This is not fatal. The report says the extraction-dependent \(T_{14}\) range is 4.2–4.9 h and \(\rho_{\rm transit}=0.19–0.39\), and modest eccentricity could also help. But the phrase “consistent with the star” is a little too smooth. The nominal values are mildly tense; the consistency relies on the longer-duration extractions and/or eccentricity/uncertainty.

### Occultation / M-dwarf consistency

The occultation value is

\[
0.79 \pm 0.77 \ {\rm ppt}.
\]

A 2σ upper limit is about

\[
0.79 + 2 \times 0.77 = 2.33 \ {\rm ppt}.
\]

The report says \(T_b < 3510\) K and that an M dwarf should produce 1.7–2.0 ppt. That is consistent, but only weakly. The non-detection does not strongly confirm the M dwarf; it merely does not rule it out. Also, “excludes any luminous companion” should be softened to “excludes a substantially hotter/brighter companion producing a large secondary eclipse.”

### Ellipsoidal amplitude

The quoted measured ellipsoidal amplitude,

\[
0.87 \pm 0.66 \ {\rm ppt},
\]

is only a ~1.3σ detection. It is consistent with the expected ~0.67 ppt, but not evidence by itself. This should be presented as a consistency check, not a confirmation.

### Twin binary at 2P

The report says a twin binary at twice the period is excluded by shape at 4.7–5.8σ. This is plausible but model-dependent. Equal odd/even depths alone do not exclude a 2P binary; in fact they are what a nearly twin EB would produce. The exclusion apparently comes from the flat-bottomed shallow shape and duration constraints.

I would ask to see the exact model family used for the 2P/twin-EB rejection. If the test only compares against a restricted “twin” model, then the conclusion should be limited to that. A different stellar binary configuration at 8.27 d may be disfavoured but not necessarily excluded at the same formal significance. This does not affect the withdrawal of TIC 87055304, but it matters for the precise EBLM interpretation.

---

## 4. Statistical claims: what is overstated or understated?

### “≥5σ exclusion of TIC 87055304”

The conclusion is probably right, but the exact sigma value is not clean. The analyses share data, Gaia positions, PRF assumptions, detrending choices, and sector-level systematics. The report acknowledges this and uses injection calibration, which is good. But the quoted significances are still partly model-dependent.

Concrete issue: the injection calibration says Gaussian-tail approximations give 7–10σ, while direct simulation counting gives only

\[
p < 5 \times 10^{-5},
\]

i.e. about 3.9σ one-sided as a finite-simulation floor. Combining with other tests may justify ≥5σ, but only if the correlations are treated carefully.

I would phrase this as:

> “The target-host scenario is excluded at very high confidence; conservative tests give at least ~3–5σ depending on assumptions, and the combined evidence strongly favours TIC 87055315.”

Rather than making the exact “≥5σ” the headline.

### Absolute χ² values look over/under-dispersed

Some quoted χ² values are suspiciously low or high without the degrees of freedom clearly stated.

Examples:

- Apertures centred on each star: χ² = 63 versus 3.5 for 15 dof. A χ² of 3.5 for 15 dof is unusually low, implying overestimated errors or correlated data.
- Per-pixel depth map: χ² = 188 versus 53 for 96 dof. Again, the preferred model χ² = 53 for 96 dof is low.
- Difference-image residuals: neighbour fits have χ² values 14–51, target fits 82–156, but the relevant dof are not given in the figure/report.

The model comparisons are still persuasive, but formal sigma conversions from these χ² values should be treated cautiously.

### “TRICERATOPS gives FPP = 1.000”

This is not very meaningful as a physical proof. Once the event is known to be a 2% eclipse on a 1.9 \(R_\odot\) star, planet scenarios are naturally killed by radius/depth priors. TRICERATOPS is useful as a consistency check, but the decisive evidence is the pixel localisation plus the companion radius.

Also, given the earlier depth-unit/runner bug, any TRICERATOPS result should be reported with explicit input depth, dilution assumptions, target star, neighbour list, and version/settings.

### Archival-sector false-alarm probability

The note says archival sectors were recovered at FAP ≤0.0066 each, with combined \(p = 3 \times 10^{-6}\). This is not central to the withdrawal, but the combination is not obviously rigorous. If 0.0066 is the floor from 150 random-period trials, multiplying such floors across sectors can overstate the significance. The signal is real enough from the multi-sector phase coherence; no need to lean on that combined FAP.

---

## 5. Claims I would soften

### “The signal is real and astrophysical”

Likely yes, but “astrophysical” is inferred from coherent repetition, localisation on a star, transit shape, and lack of sector-orientation dependence. That is strong. Still, without ground-based resolved confirmation or RVs, I would say:

> “The signal is almost certainly astrophysical and associated with TIC 87055315.”

### “No eclipse is seen at any other phase”

The secondary limit is only a few hundred ppm in the diluted target aperture, and ~0.8 ± 0.8 ppt in the neighbour frame. That rules out large secondaries but not all secondary eclipses. An M-dwarf secondary of ~1.7–2.0 ppt is apparently still allowed. So “none at any other phase” should be qualified as “no significant secondary is detected.”

### “Companion is most likely an M dwarf”

This is likely, but mass is not measured. Radius strongly excludes planet/brown dwarf under normal assumptions, but the definitive EBLM label requires RVs or a spectroscopic orbit. Better:

> “The radius is stellar; if the eclipse is on TIC 87055315, the companion is most naturally a ~0.25 \(M_\odot\) M dwarf.”

### “Novelty”

The catalogue search is good, but absence from TOI/CTOI/ExoFOP/SIMBAD/VSX/etc. is not a proof of novelty. I would state:

> “No prior report was found in the searched catalogues.”

Also consider checking or explicitly mentioning ATLAS, ASAS-SN light curves, SuperWASP/KELT if applicable, ZTF if coverage exists, AAVSO VSX aliases, and updated TESS EB lists. The object is at Dec −44°, so survey coverage may be uneven.

---

## 6. What would a sceptical TFOP SG1 vetter still ask for?

### Highest-priority request: resolved ground-based time-series photometry

Because the stars are separated by 9.3″, this is easy by TFOP standards. A 1–2 m telescope, or even smaller if conditions are good, with FWHM <4″ should separate them.

Predictions:

- TIC 87055315 should drop by ~2.1%.
- TIC 87055304 should be constant to <0.5 ppt, or at least show no ~4 ppt event.

A single well-covered ingress/egress/full eclipse would essentially settle the host question observationally.

### Reconnaissance RVs of TIC 87055315

Two spectra near quadrature should see tens of km/s motion if the companion is ~0.25 \(M_\odot\). The predicted \(K \sim 27\) km/s is large enough for modest-resolution reconnaissance spectroscopy at \(V \approx 14.2\).

This would confirm:

- the eclipsing object is stellar mass,
- the orbital period/phase,
- whether TIC 87055315 is single-lined or blended,
- whether the systemic velocity is distinct from TIC 87055304.

### High-resolution imaging or Gaia/blend check

TESS cannot exclude an unresolved source within ~1″ of TIC 87055315. A TFOP vetter would ask for:

- speckle/AO/lucky imaging if feasible,
- Gaia DR3 neighbour and duplicated-source checks,
- RUWE/BP-RP excess/astrometric excess noise,
- archival imaging/proper-motion check,
- reconnaissance spectra for secondary lines.

This is not needed to withdraw TIC 87055304 as a planet host, but it is needed before claiming a clean TIC 87055315 + M-dwarf binary architecture.

### Reproducibility package

Before accepting the withdrawal formally, I would want:

- per-sector camera/CCD, cadence type, aperture masks, quality flags;
- list of all individual events and which ones enter each analysis;
- exact TESScut download parameters;
- Gaia positions propagated to each sector epoch;
- PRF version and interpolation details;
- injection-recovery tables, not only plots;
- exact TRICERATOPS inputs, especially depth units and dilution.

The report is already unusually transparent, but the inconsistent transit counts need cleaning.

---

## 7. Final judgement

The adversarial campaign appears to have found the correct failure mode of the original planet claim: a close, fainter neighbour inside the same TESS pixels. The withdrawal is justified.

The strongest evidence is not any single formal FPP number, but the agreement between:

1. PRF difference imaging,
2. simultaneous two-star PSF photometry,
3. per-pixel depth maps,
4. a model-free pixel asymmetry test,
5. injection experiments,
6. the physically too-large radius if the event is on TIC 87055315.

The main remaining caveats are:

- the exact “≥5σ” localisation significance is not perfectly clean because of shared systematics and correlated tests;
- the number of transits used is reported inconsistently;
- the nominal transit duration/density values show mild tension that should be explained;
- the EBLM interpretation lacks RV confirmation;
- an unresolved blend within ~1″ of TIC 87055315 is not excluded by TESS pixels.

None of these resurrects TIC 87055304 as a viable planet host. They only affect how confidently one can label the neighbour system as a confirmed EBLM rather than a likely stellar eclipsing system.
