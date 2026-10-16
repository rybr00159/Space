# Review by openai/gpt-5.5

I broadly agree that these are **not validated new planets around the named TIC targets**, but several verdicts are stronger than the evidence shown. The biggest issues are:

1. **TIC 90086266:** “new eclipsing binary within ~1″” is plausible but not demonstrated; the TESS data only show the dimming is coincident with the target to TESS/PRF precision. Also the Green+23 0.2355 d entry should not be dismissed too quickly.
2. **TIC 61141412:** the signal is probably real and not a normal 6.5 d planet around the Gaia SB1, but the **ETV/third-body claim is overinterpreted** from sparse, heterogeneous timings whose errors may be underestimated.
3. **TIC 214311906 and TIC 55829034:** the eclipsed source is not the giant, but calling both specifically “blended eclipsing binaries” is not proven. TIC 214311906 in particular could plausibly be a planet transiting a nearby dwarf if the proposed 6.5″ M-dwarf pair is the source.

Below is my object-by-object review.

---

## TIC 90086266 — 3% eclipses near an F subgiant

### What looks solid

The “not a planet transiting the F star” conclusion is strong.

Arithmetic check:

- Depth in corrected flux: roughly 3.2%.
- If on the F star with \(R_\star \simeq 1.7 R_\odot\),

\[
R_{\rm occulter} \simeq R_\star \sqrt{\delta}
= 1.7 R_\odot \sqrt{0.032}
\simeq 0.30 R_\odot
\simeq 3.0 R_{\rm Jup}.
\]

That is too large for a planet or brown dwarf.

The RV argument is also broadly correct. For a \(0.3 M_\odot\) companion in a 2.339 d circular orbit around a \(1.2\)–\(1.3 M_\odot\) F star,

\[
K_1 \approx 35\ {\rm km\,s^{-1}},
\]

so Gaia’s reported 5-transit peak-to-peak span of only 5 km/s is difficult to reconcile with that scenario unless the Gaia sampling was extremely unlucky or the RV solution is compromised.

A twin-F binary at 4.68 d would indeed have very large RV amplitudes, of order 80–90 km/s per component, and should usually be obvious in Gaia spectra/RVs.

### Issues / overclaims

#### 1. “Within about 1″” is not directly shown by the TESS centroiding

The PRF fits being “consistent with the target” do **not** demonstrate a source within 1″. The quoted formal offsets, e.g.

- Sector 13: \((+0.2,-0.1)″ \pm 0.2″\)
- Sector 104: \((-1.2,-0.2)″ \pm 0.2″\)

are much smaller than the likely systematic floor for TESS PRF centroiding in a crowded field, especially when the injection test reproduces similar offsets for a signal placed on the target. That injection result is important: it says the apparent \(-1.2″\) offset in Sector 104 is likely a PRF/WCS/systematic bias, not a measured physical displacement.

A better statement would be:

> The source is unresolved from the target in TESS and not attributable to known Gaia neighbours within several arcseconds. If it is a distinct star bright enough to supply the eclipse, it is likely unresolved or marginally resolved by Gaia, hence probably within roughly an arcsecond or otherwise hidden by contrast/crowding.

But “within ~1″” is an inference from Gaia non-detection/contrast, not a direct TESS localization measurement.

#### 2. The “must supply at least ~6% of the light” bound is not a strict lower limit

Observed depth is about 3.2%. The absolute hard lower limit on the eclipsed source’s flux fraction is approximately the observed depth:

\[
f_{\rm source} \ge \delta_{\rm obs} \simeq 0.032,
\]

if the eclipse removes nearly all light from the eclipsed component.

The 6% number corresponds to assuming a maximum intrinsic eclipse depth of about 50%, as in an equal-brightness eclipsing binary. That is reasonable for many EBs, but not a model-independent lower limit. A hierarchical EB could in principle have intrinsic eclipse depths larger than 50% if the eclipsed component dominates the EB light.

So the brightness constraint should be phrased as approximately:

- hard lower limit: \(f \gtrsim 3.2\%\), i.e. \(\Delta T \lesssim 3.7\) mag;
- for a roughly equal-luminosity EB with 50% intrinsic eclipses: \(f \gtrsim 6.4\%\), i.e. \(\Delta T \lesssim 3.0\) mag.

#### 3. The Green+23 0.2355 d period is suspicious

The write-up says Green+23 lists 0.2355 d “no match.” But the originally detected periods were near 2.35 d, and 0.2355 d is almost exactly a factor of 10 smaller. That could be a decimal-place error, a harmonic/alias, or a catalogue-period failure on the same variability.

If Green+23 is associated with this target or a close neighbour, it should not be dismissed as “no match” without showing:

- the Green+23 source position,
- its aperture/contamination situation,
- the periodogram around 0.2355, 2.355, and 2.339 d,
- whether the Green+23 signal could have come from the same eclipses.

This matters for the “new identification” claim. It may still be newly recognized as an EB, but perhaps not a wholly new variable signal.

#### 4. The missing secondary eclipse is supportive, but the predicted 2000 ppm is model-dependent

For a \(0.3 R_\odot\), \(T_{\rm eff}\sim 3300\)–3500 K M dwarf next to a 6200 K, \(1.7 R_\odot\) F subgiant, a secondary eclipse of order a few \(10^{-3}\) in the red optical is plausible. But the exact TESS-band secondary depth depends strongly on the companion temperature and radius.

The non-detection at \(-60\pm45\) ppm is powerful, but the write-up should not imply a uniquely predicted 2000 ppm without showing the assumed companion model.

### What would most change my assessment

High-resolution imaging or reconnaissance spectroscopy. Specifically:

- speckle/AO imaging to \(\Delta T \sim 3\)–4 mag inside 1″;
- multi-epoch spectra looking for moving secondary/tertiary lines at the 2.339 or 4.679 d period;
- two-colour eclipse photometry to determine whether the eclipsed source is redder than the F star.

### Verdict for TIC 90086266

- **Not a planet around the F star:** strong.
- **Eclipsing binary blended with the F star:** plausible to likely.
- **New EB within ~1″:** plausible but not proven; newness weakened by the Green+23 0.2355 d entry.

Overall: **(b) plausible but unproven** as a newly identified unresolved EB; **false positive as a planet**.

---

## TIC 61332759 — variability of neighbour 57″ away

This one looks the cleanest.

### What looks solid

The offset arithmetic is consistent:

\[
\sqrt{36.5^2 + 44.3^2} \simeq 57.4″.
\]

The folded curve shown is not transit-like. It is broad, smooth, and has two unequal minima/maxima per 1.4806 d cycle. The 0.740 d search period is naturally half the true variability period.

The aperture-dependence argument is also strong:

- QLP large apertures show hundreds to thousands of ppm;
- smaller SPOC/TESS-SPOC apertures show only tens of ppm;
- pixel-level localization identifies the Gaia DR3 RS CVn-class variable, not the target.

If the reported \(\Delta\chi^2\) values are accurate, the target is decisively excluded.

### Minor caveat

The “period appears to be new” claim is lower importance and should be phrased conservatively. The neighbour is already catalogued as an RS-type variable, so the astrophysical source is not new. Only the period may be newly measured, and that should be cross-checked against ASAS-SN, ZTF, ATLAS, TESS EB/variable-star catalogues, Gaia DR3 variability tables, and any published TESS variable catalogues.

### Verdict for TIC 61332759

- **Not a planet:** secure.
- **Source is known nearby RS-type variable:** secure.
- **1.4806 d period new:** possible but not central.

Overall: **(c) likely false positive / already known variable source**.

---

## TIC 61141412 — 6.5285 d dips in/near a Gaia SB1

I agree this is unlikely to be a normal planet transiting the catalogued subgiant, but I would weaken the ETV and “third body” claims.

### What looks solid

The Gaia SB1 orbit is a major problem for a 6.53 d planet around one component.

For a 5.452 d stellar binary with total mass of order \(2 M_\odot\), the binary semimajor axis is roughly:

\[
a_{\rm bin} \sim 0.07\text{--}0.08\ {\rm AU}.
\]

A 6.53 d circumstellar orbit around one component would have a semimajor axis comparable to the binary separation, hence dynamically unstable. A circumbinary planet would generally need to orbit outside roughly 2–4 binary separations, giving periods more like tens of days, not 6.5 d.

The duration argument is also reasonable. For \(R_\star \sim 2.5 R_\odot\) and \(P=6.53\) d, a central transit duration would be many hours. A 1 h dip would require an extremely grazing geometry or a different/smaller eclipsed star.

The TESS-SPOC alert at 13.054 d also means this is not a new, completely unreported TESS event.

### Issues / overclaims

#### 1. The ETV claim is statistically impressive only if the timing errors are trustworthy

The write-up says:

- constant period: \(\chi^2 = 243\) for 9 dof;
- quadratic ephemeris: \(\chi^2 = 22\) for 8 dof.

If the quoted 3–5 min timing errors are correct, the constant-period model is decisively rejected, and the quadratic model is a huge improvement.

But the figure raises concerns:

- Several points within apparently similar observing epochs scatter by \(\sim 10\)–20 min.
- The quadratic still has \(\chi^2=22\) for 8 dof, corresponding to a poor fit, \(p\sim0.005\)–0.01.
- The events are shallow, short, and measured in different sectors/products/cadences, likely on top of rotational modulation from the SB1.
- Starspot variability at 5.4–5.8 d can distort local baselines and shift fitted mid-times of shallow 1 h dips.

So the evidence supports “timings are not described well by a simple linear ephemeris using the quoted errors,” but not yet a secure physical ETV detection.

A more cautious statement would be:

> The event times show apparent O–C structure at the tens-of-minutes level, but the significance depends on the adopted timing uncertainties and treatment of correlated stellar variability. A third body is one possible explanation, not established.

#### 2. “The dipping pair is therefore being pulled or moved by a third body” is too strong

A quadratic O–C over sparse epochs can arise from:

- underestimated timing errors;
- sector-dependent detrending or cadence differences;
- spot-induced timing biases;
- fitting a variable or asymmetric dip with a fixed template;
- incorrect cycle counts if there are gaps and ETVs;
- blending with another unresolved variable;
- a wrong assumed period for part of the data.

The improvement from linear to quadratic is suggestive, but “therefore” is not justified without injection/recovery and uniform pixel-level modelling.

#### 3. The Gaia SB1 period should be treated as strong but not infallible

A Gaia NSS SB1 with significance 51 and RV span 35 km/s is compelling. Still, if the TESS dips are truly on the target, one should check whether:

- the Gaia SB1 period has aliases near 6.53 d or 13.05 d;
- the Gaia RV phases are consistent with the TESS dip ephemeris;
- the TESS photometric rotation period at 5.4–5.8 d is independent or simply broad/ambiguous.

The write-up mentions this possibility, but the planet-impossibility argument depends on the 5.452 d SB1 orbit being the relevant physical binary.

#### 4. An unresolved third star hosting the dips is not excluded

The TESS localization rules out catalogued Gaia stars beyond 17″, but not an unresolved or very close companion. A planet transiting a bound third star is unlikely given the timing behaviour and depth/duration, but not fully excluded by the evidence shown. The statement “likeliest picture is an eclipsing pair” is reasonable; “not a planet” is robust for the named Gaia SB1, less absolute for any unresolved contaminant.

### What would most change my assessment

The most important test is a **single homogeneous remeasurement of all 11 event times** from the same pixel-level photometry, with:

- local stellar-variability model fitted simultaneously;
- cadence integration included;
- per-event template freedom or shape diagnostics;
- bootstrap or injection-recovery timing errors;
- comparison between SAP, PDCSAP, and independent TESScut apertures.

Also very valuable: phase Gaia RVs to both 5.452 d and 6.5285/13.054 d, and obtain a few ground-based spectra near predicted RV extrema.

### Verdict for TIC 61141412

- **Not a normal 6.5 d planet around the Gaia SB1:** strong.
- **Eclipsing pair in/near a multiple system:** plausible.
- **ETVs caused by a third body:** suggestive but unproven.
- **Newness:** not new as a TESS event because of the public TESS-SPOC TCE.

Overall: **(c) likely false positive / already alerted**, with the ETV interpretation **(b) plausible but unproven**.

---

## TIC 214311906 — 5.6915 d dips blended with a Mira

The conclusion that the Mira itself is not the eclipsed star is strong. The conclusion that the source is specifically an EB is not.

### What looks solid

A 294 d Mira with BP–RP = 4.7 is almost certainly a large giant, plausibly hundreds of \(R_\odot\). A 5.69 d orbit has a semimajor axis only of order:

\[
a \sim 0.06\ {\rm AU} \sim 13 R_\odot
\]

for a solar-mass system. That is far inside a Mira envelope. A 1.6 h eclipse/transit across a hundreds-\(R_\odot\) star is physically impossible. So the dip is blended from another star.

### Major issue: the proposed nearby M-dwarf source could host a planet

The write-up says a nearby M-dwarf pair at 6.4–6.6″ would need about 2.3% eclipses. That depth is **not diagnostic of an eclipsing binary**.

For an M dwarf of radius \(R_\star \sim 0.3 R_\odot\),

\[
R_p \simeq R_\star \sqrt{0.023}
\simeq 0.046 R_\odot
\simeq 0.45 R_{\rm Jup}
\simeq 5 R_\oplus.
\]

That is a Neptune/sub-Saturn-sized planet, not a star.

The duration is also plausible for a planet around an M dwarf. For \(P=5.69\) d around an \(M\sim0.3 M_\odot\), \(R\sim0.3 R_\odot\) dwarf, a central transit duration is roughly 1–1.5 h, comparable to the reported 1.6 h.

So if the 6.5″ M-dwarf pair is indeed the source, a planet scenario is not absurd. The high RUWE suggests binarity, and an EB may be more likely, but the data shown do not rule out a planet around the contaminating dwarf.

### Better wording

Instead of “blended eclipsing binary,” I would say:

> The signal is not on the Mira. It is a blended eclipsing/transiting source within about 10″. A nearby M-dwarf pair is a plausible host; the required intrinsic depth is only a few percent, so both an EB and a planet around the contaminant remain possible without resolved follow-up.

### What would most change my assessment

A single ground-based time-series observation resolving the 6.5″ M-dwarf pair from the Mira. This is an easy, decisive test. If the M dwarf dims by 2–3%, the next step would be multicolour photometry and RV/imaging to distinguish planet vs EB.

### Verdict for TIC 214311906

- **Not a planet around the Mira:** secure.
- **Blended EB:** not proven.
- **Possible planet around a resolved contaminant:** still viable.

Overall: **(b) plausible but unproven** as an EB; secure false positive for the Mira target.

---

## TIC 55829034 — 0.6161 d dips blended with a bulge giant

Again, the target giant is almost certainly not the eclipsed star, but the exact nature of the contaminant is not demonstrated.

### What looks solid

A 0.616 d orbit has a semimajor axis of only a few stellar radii for a solar-type star and far inside any red giant. Thus the bulge giant cannot be the eclipsed object.

The field is extremely crowded, and TESS cannot identify the source within 5–8″. So this is definitely a blend or unresolved nearby source.

### Planet scenario?

A 0.616 d signal is very short-period. The V-shaped morphology and crowded bulge field make an EB more likely than a planet. However, “EB” is still an inference.

A hot Jupiter or USP planet around a blended dwarf/subgiant is not impossible in principle. The observed depth of 1000–2300 ppm, after dilution by the giant and neighbours, could correspond to a percent-level transit on a fainter star. The 2 h duration is somewhat long for a compact M dwarf at 0.616 d, but could fit a larger/less dense contaminant or be broadened by cadence/noise/detrending. The V-shape favours grazing EB, but is not decisive from the shown data alone.

### Better wording

> The signal cannot originate on the red giant. It is a blended short-period eclipsing/transiting source within several arcseconds in a very crowded bulge field. The morphology and period make an EB likely, but the source is unidentified.

### What would most change my assessment

High-resolution or seeing-limited ground photometry during an event, preferably with image subtraction, to identify the variable star. Given the 0.616 d period, events are frequent, so this should be feasible, though crowding may be challenging.

### Verdict for TIC 55829034

- **Not a planet around the red giant:** secure.
- **Blended EB:** likely but not proven.
- **New signal:** possible, but source unidentified.

Overall: **(b) plausible but unproven** as an EB; secure false positive for the giant target.

---

## Summary table of my revised assessment

| TIC | Original verdict mostly correct? | Main concern | My verdict |
|---|---:|---|---|
| 90086266 | Mostly | “Within ~1″” and “new” are overclaimed; Green+23 0.2355 d suspicious; flux-fraction lower limit overstated | Not planet on F star; unresolved EB plausible but unproven/newness uncertain |
| 61332759 | Yes | Only novelty of 1.4806 d period needs caution | Clear neighbour variable false positive |
| 61141412 | Partly | ETV/third-body claim too strong; timing errors likely underestimated; already TESS-SPOC TCE | Not normal planet around SB1; EB/multiple plausible; ETV unproven |
| 214311906 | Partly | Source not Mira, but could be planet around nearby M dwarf; EB not demonstrated | Blend; EB vs planet around contaminant unresolved |
| 55829034 | Partly | Source not giant, but unidentified; EB likely not proven | Blend; likely EB but unproven |

Final bottom line: **none of the five is a validated new planet around the searched TIC target**, but the write-up should soften several classifications from “eclipsing binary” to “blended eclipsing/transiting source,” and should especially weaken the TIC 90086266 “new EB within ~1″” and TIC 61141412 “ETV caused by a third body” claims.
