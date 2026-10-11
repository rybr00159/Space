# Review by openai/gpt-5.5

## Bottom line verdicts

| Candidate | My verdict | Main reason |
|---|---|---|
| **A TIC 143042674** | **Plausible but unproven as a new G+M eclipsing binary** | The eclipses, secondary, and Gaia RV variability strongly support a real stellar companion. The radius arithmetic is mostly OK. But “new” is not proven, and the mass/temperature interpretation needs real RVs, not Gaia summary statistics. |
| **B TIC 123267483** | **Very likely a real detached EB; newness unproven** | The primary/secondary eclipses are enormous and repeat over years. A background contaminant would need to be nearly as bright as the target, so the EB is almost certainly the target or an unresolved bound companion. But such a deep bright EB could easily already exist in ground-based variable-star surveys not fully covered by the stated checks. |
| **C TIC 80064194** | **Very likely a real detached EB; newness unproven** | Same as B, though the eclipses are less extreme. Difference imaging is mostly consistent with the target, but S104 has a 2.2σ offset and neighbour constraints were not described. |
| **D TIC 72651793** | **Likely not a planet; astrophysical/instrumental false-positive still not settled** | The odd/even difference, too-short duration / high fitted density, short period, and high TRICERATOPS FPP all argue against a planet. The detection itself is plausible but the statistical treatment is weak for a 410 ppm signal found in a large search. |

---

# 1. Candidate-by-candidate review

## Candidate A — TIC 143042674

### What is well supported

This looks like a real eclipsing binary or transiting low-mass stellar companion.

Evidence in favour:

- Primary eclipse repeats in multiple sectors: S4, S31, S105.
- Period is coherent: **P = 5.044137 d**.
- Primary depth is large: **2.2–2.4%**, too deep for a planet around a 1.46 R☉ star unless the object is near/beyond the deuterium-burning regime.
- Secondary eclipse at phase **0.502**, depth **870 ppm**, apparently significant.
- Gaia DR3 RV diagnostics are strongly suspicious:
  - `rv_amplitude_robust = 34.9 km/s`
  - `rv_chisq_pvalue = 0.0`
  - RV uncertainty on the mean, **3.1 km/s**, is high for a G ≈ 11 star.
- Difference imaging offsets of **2.4″** and **5.0″** are consistent with the target at only **0.8σ** and **1.6σ**.

So the claim “this is not a planet; it is probably a stellar/substellar companion” is well supported.

### Arithmetic checks

**Radius ratio and companion radius**

Given

\[
R_c/R_\star = 0.1384
\]

and

\[
R_\star = 1.46 R_\odot,
\]

then

\[
R_c = 0.1384 \times 1.46 = 0.202 R_\odot.
\]

So the quoted **0.203 R☉** is correct.

In Jupiter radii,

\[
1 R_\odot \simeq 9.73 R_J,
\]

so

\[
0.203 R_\odot \simeq 1.98 R_J.
\]

The quoted **1.97 RJ** is also correct.

**Transit density**

Using the fitted

\[
a/R_\star = 9.57,\quad P = 5.044137 \rm \ d,
\]

the usual circular transit-density estimate gives roughly

\[
\rho_{\rm tr} \simeq 0.46 \rho_\odot,
\]

so that number is also correct.

However, for a stellar companion this is not exactly the primary-star density. The transit-derived quantity is proportional to \(M_1+M_2\), so

\[
\rho_{\rm tr} \approx \rho_1(1+q).
\]

If \(q = M_2/M_1 \sim 0.18\), then

\[
\rho_1 \approx 0.46/1.18 \approx 0.39 \rho_\odot.
\]

That is closer to the TIC value of

\[
\rho_{\rm TIC} = M/R^3 = 1.0/1.46^3 \approx 0.32 \rho_\odot.
\]

So the density agreement is not bad.

### M-dwarf interpretation

The radius is very M-dwarf-like: **0.18–0.20 R☉**, depending on whether the primary radius is closer to the TIC value or the transit-density-implied value.

The secondary eclipse depth also roughly supports a cool stellar companion. The observed flux ratio is

\[
F_2/F_1 \approx 870 \ {\rm ppm} = 8.7 \times 10^{-4}.
\]

The area ratio is

\[
(R_2/R_1)^2 = 0.1384^2 = 0.0192.
\]

Therefore the TESS-band surface-brightness ratio is

\[
J_T \approx \frac{8.7\times10^{-4}}{0.0192} \approx 0.045.
\]

A TESS-band surface brightness about **4–5%** of a 5600 K primary is plausible for a late-M dwarf. A temperature of **2800–3200 K** is plausible, though the exact number is model-dependent and should not be overclaimed from a single TESS-band occultation.

I would not claim a precise **Teff ≈ 2800–2900 K** without using a proper bandpass-integrated stellar atmosphere model. Also, a **0.18 M☉** main-sequence star is often somewhat hotter than 2800–2900 K; 3000–3300 K may be more typical depending on metallicity and age. The photometry is consistent with an M dwarf, but not enough for a precise spectral subtype.

### Gaia RV amplitude check

For a circular orbit with

- \(P = 5.044\) d,
- \(M_1 = 1.0 M_\odot\),
- \(M_2 = 0.18 M_\odot\),
- \(i \approx 90^\circ\),

the primary RV semi-amplitude is approximately

\[
K_1 \approx 20 \ {\rm km\,s^{-1}}.
\]

That implies a peak-to-peak RV range of about

\[
2K_1 \approx 40 \ {\rm km\,s^{-1}}.
\]

So Gaia’s `rv_amplitude_robust = 34.9 km/s` is very compatible **if** that quantity is closer to a robust peak-to-peak or percentile amplitude.

But if someone interprets `rv_amplitude_robust = 34.9 km/s` as the semi-amplitude \(K\), then the implied companion mass would be much larger, roughly **0.3–0.4 M☉**. The Gaia summary field alone is not enough to measure the mass.

So: Gaia RVs strongly support binarity, but they do not by themselves prove \(M_2 \approx 0.18 M_\odot\).

### Things I would dispute

1. **“New” is not demonstrated.**  
   The catalogue checking is good, but not exhaustive. A bright 2%-deep 5-d EB could exist in ASAS, SuperWASP, KELT, ATLAS, etc., even if not in VSX or ASAS-SN.

2. **The companion mass is not measured.**  
   A radius of 0.20 R☉ suggests an M dwarf, but the mass estimate is model-based. Gaia RV diagnostics are supportive but not a mass solution.

3. **The quoted Teff is too precise.**  
   The secondary eclipse gives a TESS-band flux ratio. Converting that into 2800–2900 K requires assumptions about radius, primary temperature, bandpass, metallicity, and atmosphere models.

### Verdict for A

Likely a real low-mass stellar eclipsing companion. Plausibly new, but not proven new. The M-dwarf interpretation is broadly sound, but mass and Teff need RVs / SED modelling.

**Overall: (b) plausible but unproven.**

---

## Candidate B — TIC 123267483

### What is well supported

This is almost certainly a real detached EB.

The evidence is very strong:

- Primary eclipse depth: **37–44%**, depending on pipeline.
- Secondary eclipse depth: **27–32%**.
- Period coherent over 14 primary eclipses and 387 cycles:
  \[
  P = 7.2732209 \pm 0.0000020 \ {\rm d}
  \]
- Timing RMS only **1.5 min**, which is excellent.
- Difference imaging offsets:
  - S97: **1.8″**, 0.6σ
  - S104: **3.3″**, 1.1σ
- Only Gaia neighbour within 20″ is **G = 17.8 at 10.6″**.

A background contaminant cannot plausibly explain a 37–44% eclipse unless it is very bright. For an observed depth \(d\), a totally eclipsed contaminant must contribute at least \(d\) of the total light.

For \(d = 0.37\),

\[
\frac{F_{\rm cont}}{F_{\rm total}} \geq 0.37,
\]

so

\[
F_{\rm cont}/F_{\rm target} \geq \frac{0.37}{0.63} = 0.59.
\]

That corresponds to a magnitude difference

\[
\Delta m \leq -2.5\log_{10}(0.59) \approx 0.57 \ {\rm mag}.
\]

For a 44% eclipse,

\[
F_{\rm cont}/F_{\rm target} \geq 0.79,
\]

or

\[
\Delta m \lesssim 0.26 \ {\rm mag}.
\]

So a contaminating EB would need to be within about **0.3–0.6 mag** of the target and unresolved or nearly coincident. The known G = 17.8 neighbour is far too faint.

Thus this is either the target itself or a near-equal-brightness unresolved bound companion.

### Things I would dispute

1. **The depth numbers are internally inconsistent.**  
   The text says primary depth is **44% in TESS-SPOC** and **37% in QLP**, but the figure labels the folded all-sector primary as **31.80%** and the secondary as **27.794%**. The plotted model minimum also appears closer to 0.65–0.70 than 0.56. This needs clarification:
   - Are the plotted depths after dilution correction?
   - Are sector-to-sector normalizations different?
   - Is the “44%” depth from raw SAP, PDCSAP, or DET_FLUX?
   - Are the grey points from a different pipeline than the blue binned curve?

2. **“New” is not established.**  
   A T ≈ 10.5 star with 30–40% eclipses and P ≈ 7.27 d is exactly the kind of object that ground-based surveys often detect. Absence from VSX/ASAS-SN/Gaia EB/TESS EB is useful, but not definitive.

3. **The eccentricity statement should be treated cautiously.**  
   Secondary phase 0.485 implies a small eccentricity. Depending on the approximation used,
   \[
   e\cos\omega
   \]
   is of order −0.02 to −0.05. The quoted −0.024 is plausible for the EB convention, but I would not overinterpret it without a full binary light-curve solution.

### Verdict for B

The EB nature is very strongly supported. The “new” claim is plausible but requires more exhaustive archival checking.

**Overall: (b) plausible but unproven as new; physically very likely a real EB.**

---

## Candidate C — TIC 80064194

### What is well supported

This also looks like a real detached EB.

Evidence:

- Primary eclipse depth: roughly **19–21%**.
- Secondary eclipse depth: roughly **7–8%**.
- Period coherent over 15 primaries and 438 cycles:
  \[
  P = 5.8106944 \pm 0.0000015 \ {\rm d}
  \]
- Timing RMS: **1.4 min**.
- Difference imaging:
  - S67: **0.4″**, 0.1σ, excellent.
  - S104: **6.9″**, 2.2σ, mildly concerning but not decisive.

A background EB causing a 21% eclipse would need to be fairly bright. A totally eclipsed contaminant would require

\[
F_{\rm cont}/F_{\rm target} \geq \frac{0.21}{0.79} = 0.266,
\]

which corresponds to

\[
\Delta m \leq -2.5\log_{10}(0.266) \approx 1.44 \ {\rm mag}.
\]

So any contaminating EB would need to be within about **1.4 mag** of the target if fully eclipsed, or brighter if the eclipse is partial. That should be obvious in Gaia unless it is unresolved or very close.

### Things I would dispute

1. **The S104 centroid offset deserves follow-up.**  
   A 2.2σ offset is not proof of contamination, but it is also not negligible. Since the S67 localization is excellent, the signal is probably on-target, but I would still want to inspect:
   - individual sector difference images,
   - Gaia neighbours out to at least 1 arcmin,
   - aperture masks,
   - whether the S104 offset direction points to a specific star.

2. **Neighbour constraints are incomplete.**  
   For B, the nearest Gaia neighbour was listed. For C, no Gaia neighbour information is given except catalogue-variable matches at 85″ and 120″. For a 21% eclipse, any contaminant must be bright; this should be explicitly quantified.

3. **“New” remains unproven.**  
   Like B, this is a bright, deep EB and could be present in ground-based time-domain catalogues that are not covered by the stated checks.

### Verdict for C

Very likely a real EB. The target association is likely but should be checked against nearby Gaia sources and the S104 offset. Newness is plausible but not proven.

**Overall: (b) plausible but unproven as new; likely real EB.**

---

## Candidate D — TIC 72651793

### What is well supported

The “not a validated planet” statement is strongly supported.

Reasons:

1. **Odd/even difference is significant.**  
   The quoted odd/even difference is **3.5σ**, with odd events deeper. That immediately suggests the adopted period may be half the true EB period:
   \[
   P_{\rm true} \approx 2 \times 0.6345291 = 1.2690582 \ {\rm d}.
   \]

2. **Transit-derived density is inconsistent.**  
   The fit gives
   \[
   \rho_{\rm fit} \approx 5.7 \rho_\odot,
   \]
   while the TIC value is
   \[
   \rho_{\rm TIC} \approx 0.54 \rho_\odot.
   \]
   That is a factor of
   \[
   5.7/0.54 \approx 10.6
   \]
   discrepancy.

3. **Duration is short for the stellar parameters.**  
   For \(P = 0.6345\) d, \(M_\star = 1.17 M_\odot\), \(R_\star = 1.29 R_\odot\), the orbital separation is only about
   \[
   a \approx 3.3 R_\odot,
   \]
   so
   \[
   a/R_\star \approx 2.5.
   \]
   A central transit would last of order **2 h**, whereas the observed duration is **0.83 h**. A grazing geometry can shorten it, but then the V-shape and odd/even behaviour become even more suspicious.

4. **TRICERATOPS strongly disfavors the planet interpretation.**  
   FPP = **0.86**, with large contributions from EB-like scenarios:
   - SEBx2P: 0.45
   - SEB: 0.26
   - STP: 0.13
   - TP: 0.11

That is not a validation-quality result; it is a false-positive warning.

### Things I would dispute

1. **The false-alarm test is too weak.**  
   “FAP < 0.007 vs 150 random-period nulls” is not very strong. With 150 trials, the minimum resolvable p-value is about
   \[
   1/151 \approx 0.0066.
   \]
   So the statement mostly says “none of the 150 random periods was better.” That is not enough for a 410 ppm signal found in a search over **77,000 stars** and many trial periods/durations.

2. **Look-elsewhere effects are not adequately handled.**  
   A per-sector p ≈ 0.007 would be unimpressive in a 77,000-star search. Requiring recovery in S3, S30, and S97 helps substantially, but those tests are conditioned on a period/epoch found elsewhere. A proper global false-alarm estimate should include:
   - number of stars searched,
   - period grid,
   - duration grid,
   - phase grid,
   - detrending choices,
   - odd/even selection effects,
   - post-selection testing.

3. **No secondary at phase 0.5 is not decisive.**  
   If the adopted period is half the true binary period, then the “odd” and “even” events are the primary and secondary eclipses. The relevant folded period may be **1.2690582 d**, not 0.6345291 d.

4. **The Gaia neighbour at 6″ is probably too faint, but unresolved blends remain possible.**  
   A G = 19.5 star near a T = 10.49 target is about 9 mag fainter, corresponding to a flux ratio of roughly
   \[
   10^{-0.4 \times 9} \approx 2.5 \times 10^{-4}.
   \]
   A fully eclipsed such star could produce only about **250 ppm**, and usually less. So that particular neighbour probably cannot explain a 410 ppm signal unless colours/bandpasses change the ratio substantially. But a closer unresolved companion or bound EB is still plausible.

### Verdict for D

The planet interpretation is not supported. The signal may be real, but the current evidence points to an EB/blend/systematic rather than a planet.

**Overall: (c) likely false positive as a planet; plausible unresolved EB or systematic.**

---

# 2. Catalogues/surveys that may have been missed

The stated checks are good but not exhaustive. Specific catalogues/surveys I would still check explicitly:

### Ground-based time-domain surveys

- **ASAS-3** / Pojmański ASAS catalogue, separate from ASAS-SN.
- **SuperWASP** variable catalogues.
- **KELT** variable-star catalogues.
- **ATLAS** variable-star catalogue, e.g. Heinze/Tonry ATLAS periodic variables.
- **CRTS** / Catalina variable catalogues, where sky coverage permits.
- **NSVS**, for northern or equatorial fields; probably less relevant for Dec ≈ −33 to −35 but still worth checking depending on field.
- **ZTF variable catalogues**, for objects near the southern limit. Dec −33 to −35 is mostly outside normal ZTF, but some coverage may exist.
- **SkyMapper Southern Survey** variability products, if available.
- **OGLE**, if near OGLE footprint. These coordinates do not obviously scream OGLE field, but for southern targets it is worth checking.
- **HATSouth** and **HATNet** transit/variable catalogues.
- **MASCARA / bRing** for bright stars, though T ≈ 10–11 is often too faint.

### TESS-specific catalogues

- The online **Villanova/TESS Eclipsing Binary Catalog**, not only the static Prša et al. 2022 VizieR table.
- Other TESS EB catalogues from machine-learning searches, e.g. sector-limited catalogues not mirrored cleanly in VizieR.
- **TESS-Gaia / FFI variable-star catalogues**, e.g. Oelkers/Stassun-style or later FFI variability catalogues.
- QLP/TESS-SPOC threshold-crossing-event lists if public separately from the 2-min SPOC DVR TCE tables.
- MAST HLSPs such as **TGLC**, **PATHOS**, **CDIPS**, **eleanor**-derived catalogues, where relevant.

### Spectroscopic binary catalogues

Especially for Candidate A:

- **SB9** spectroscopic binary catalogue.
- **Gaia DR3 NSS** tables beyond just the `non_single_star` flag, including acceleration / spectroscopic solutions if present.
- **RAVE**, **GALAH**, **LAMOST**, **APOGEE**, **Gaia-ESO**, where sky coverage permits.
- Literature RV catalogues via SIMBAD references, not just SIMBAD object type.

### General caution about the VizieR search

A VizieR query over all catalogues within 15″ filtered by variability/binary keywords can miss things because:

- catalogue titles/columns may not contain the expected keywords,
- some variable catalogues use survey-specific names,
- coordinates may be old epoch or rounded,
- bright-star catalogues sometimes have arcsecond-to-tens-of-arcsecond position errors,
- online catalogues are not always mirrored to VizieR,
- the relevant object may be listed under a cross-ID rather than the TIC/Gaia source.

For B and C, a manual search in ASAS-SN, ATLAS, SuperWASP, KELT, and the online TESS EB catalogue is still warranted.

---

# 3. Candidate A: is the M-dwarf interpretation sound?

Broadly, yes — but with caveats.

### Sound parts

- \(R_c = 0.203 R_\odot\) is correctly calculated.
- \(R_c \approx 2 R_J\) is correctly calculated.
- The radius is stellar, not planetary.
- The 870 ppm secondary eclipse is consistent with a cool M dwarf.
- Gaia RV variability is quantitatively consistent with a stellar companion.

### Caveats

- The companion radius inherits the primary-radius uncertainty. If the primary radius is closer to 1.3 R☉ than 1.46 R☉, then
  \[
  R_c \approx 0.18 R_\odot
  \]
  rather than 0.203 R☉.
- The transit-derived density should include the companion mass:
  \[
  \rho_{\rm tr} \approx \rho_1(1+q).
  \]
- The secondary eclipse constrains the **TESS-band flux ratio**, not directly Teff.
- Gaia’s `rv_amplitude_robust` is not an orbital semi-amplitude. It supports binarity but does not by itself give a precise mass.
- A few proper RVs are needed to distinguish, for example, a 0.16–0.20 M☉ M dwarf from a more massive unresolved/binary configuration.

My assessment: the interpretation “late-M dwarf eclipsing a G/subgiant primary” is plausible and probably correct, but the exact mass and Teff are not established.

---

# 4. Single most valuable additional check for each

## A — TIC 143042674

**Obtain phase-resolved radial velocities.**

Even 3–5 reconnaissance spectra at well-chosen orbital phases would be decisive. Expected RV semi-amplitude for a 0.18 M☉ companion is roughly:

\[
K_1 \sim 20 \ {\rm km\,s^{-1}}.
\]

That is easy to detect. RVs would:

- confirm the companion is bound to the target,
- measure the mass function,
- test the 0.18 M☉ interpretation,
- reveal any double-lined or hierarchical-triple scenario.

This is the single most important follow-up.

## B — TIC 123267483

**Check and/or extract archival ground-based light curves from ASAS-SN, ATLAS, SuperWASP, and KELT.**

The EB is already obvious from TESS. The main uncertainty is whether it is truly new. A 30–40% eclipsing binary at T ≈ 10.5 should often be recoverable from ground-based data.

If the goal is astrophysical characterization rather than newness, then the best next step is phase-resolved RVs. But for the specific discovery claim, archival survey checking is the highest value.

## C — TIC 80064194

**Perform a resolved on/off-eclipse image or difference image from ground-based data.**

The EB is likely real, but the S104 offset of **6.9″ / 2.2σ** and the lack of explicit neighbour constraints leave some localization uncertainty. A modest ground-based observation during eclipse should confirm whether the target itself fades.

Second priority: archival survey search for prior discovery.

## D — TIC 72651793

**Test the doubled-period EB hypothesis with independent data.**

Specifically fold and model at

\[
P = 1.2690582 \ {\rm d}
\]

and compare odd/even events as primary/secondary eclipses. The most valuable practical checks would be:

- independent FFI extractions from another pipeline/aperture,
- a proper multi-sector bootstrap including the full search look-elsewhere effect,
- high-resolution spectroscopy or reconnaissance spectra looking for line-profile/RV variability on the 1.269 d period.

If forced to choose one: I would do **independent light-curve extraction plus a doubled-period EB fit** before spending telescope time. The current evidence already argues strongly against a validated planet.
