import json
B = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/'
F = {}
F['headline'] = (
"CTOI 71316387.01 is not a planet. Its 0.6518-d signal is half the period of an eclipsing binary (EB) with P = 1.3036723 d and unequal eclipses of 5770 and 2600 ppm. "
"The G star (TIC 71316387, 297 pc) is not one of the eclipsing stars, because its Gaia RV is constant. "
"The eclipsing stars are a faint early-M + mid-M pair: about M2V + M4.5V, 0.45 + 0.20 Msun, masses not dynamical. "
"The pair supplies 2.6 +- 0.4% of the TESS flux and lies within 1.24\" of the G star (95%), most likely 0.4-0.65\" away (120-195 au projected). "
"A bound hierarchical triple (G2V + M-dwarf EB) is strongly favoured over a chance-aligned unrelated M-dwarf EB. "
"Posterior odds are about 8e3 (5e3-1.2e4; at least 1.5e3 under conservative assumptions), i.e. P(bound) about 0.9999. "
"This is a probabilistic inference, not a demonstration: there is still no resolved image, common-proper-motion test or radial velocity of the M dwarfs.")

F['bound_vs_background'] = (
"HYPOTHESES. H_T: the EB is a physically bound companion of the G star (hierarchical triple). "
"H_C: the EB is an unrelated M-dwarf EB, in the foreground or background, that happens to lie in the localisation region.\n\n"
"EXCLUDED BEFORE THE ODDS (these define the hypothesis class and are not likelihood factors):\n"
"- The G star as an eclipsing component. Gaia RVS is constant (15 transits, p = 0.275, robust amplitude 9.9 km/s), but K1 >= 18 km/s would be needed even for a 0.1 Msun companion.\n"
"- Resolved neighbours. Colour-corrected TESS centroids put the source within 1.24\" at 95% (1.48\" at 99.7%); separations >= 1.5\" are excluded at p = 0.002, and the nearest Gaia neighbour is 18.2\" away.\n"
"- G/F/K-coloured eclipsed stars. The TESS colour shift gives BP-RP about 2.3 and excludes K7 at 5.7 sigma (statistical). The SuperWASP/TESS depth ratio is 0.36 +- 0.11 (referee) or 0.42 +- 0.15 (chromatic), against >= 0.82 for G/F.\n"
"- Giants. The eclipse shape gives rho1 of about 1-4 rho_sun.\n\n"
"PRIOR ODDS O_0 (like-for-like: position r <= 1.24\", eclipsed-star colour and light fraction matched):\n"
"- H_T: CSF 0.58 x f_sep(0.1-1.24\") 0.272 x q_tot window 0.55-0.70 (0.167 of flat 0.1-1) = 0.0263, times r_EB.\n"
"- H_C: 4.1e-6 unrelated stars with BP-RP 2.1-2.5 supplying 2.2-6.5% of the TESS light (T 14.33-15.50) per 1.24\" circle, times r_EB. This is my recount, 77 Gaia stars in 7.07 deg2 (synthesis/s04); the referee found 3.3e-6 from 7 stars in 0.785 deg2.\n"
"- r_EB cancels.\n"
"- O_0 = 6.5e3 (referee 7.9e3; Poisson +-11% with 77 hosts). Widening to BP-RP 1.8-2.7 with q_tot 0.48-0.78 gives 5.2e3.\n\n"
"LIKELIHOOD FACTORS:\n"
"(1) Position, F_pos = 0.87-1.14. This is one joint separation likelihood applied identically to both hypotheses: TESS colour-corrected localisation map x Gaia IPD window (soft or hard 0.48-1.2\") x VHS Ks non-detection probability for dKs = 2.62 +- 0.18. It replaces the hard 1.24\" window in O_0 rather than multiplying it. The bound prior's preference for small separations and the background's area weighting toward large ones nearly cancel once IPD (>~0.45\") and VHS (<~0.6\") are imposed.\n"
"(2) Light fraction / photometric distance, F_LC = 1.1-1.6. This uses the referee's independent light-curve + main-sequence chi2 profiles L(f). H_T requires the f at which a main-sequence pair sits at 297.4 pc: f = 2.26% (P&M) or 2.47% (Mann), with 0.15 mag scatter, giving L/L_max = 0.35-0.50. H_C averages L(f) over the actual Gaia host population, where 31% of hosts in the window carry high L.\n"
"(3) Colour, F = 1. It is already used to match the host class in O_0.\n"
"(4) Gaia IPD/C* signature, F = 1 beyond its role in F_pos. Both hypotheses predict the same unresolved companion at the same separation and contrast.\n"
"(5) ETV, PMa/RUWE, G-star RV and SED, F = 1. None is sensitive at 110-250 au, and the SED limits (>= 5.8%) lie above both.\n"
"(6) Not applied: an enhancement of close binaries in triples (E ~ 1-3+, which favours bound). Left out to stay conservative.\n\n"
"POSTERIOR ODDS: 6.1e3-1.2e4 (4.9e3-9.6e3 with the wide colour window), central about 8e3, so P(bound) is about 0.9999.\n"
"- Conservative floor: the lowest value divided by about 3 for prior-model systematics (period distribution, q shape, CSF) gives >= 1.5e3, P >= 0.9993.\n"
"- Matching any host bright enough for the depth, with no colour matching: about 230-265 (referee's ~230 at f_EB = 2.8% x F_pos).\n"
"- Most pessimistic any-host prior choices (astro): 25-36, P >= 0.96. This route ignores the colour evidence, which removes most such hosts.\n\n"
"AVOIDING DOUBLE COUNTING:\n"
"- Colour and light fraction are used once, to define the matched host class.\n"
"- TESS position, Gaia IPD and VHS are combined into one separation likelihood that replaces the window.\n"
"- Distance enters only through L(f).\n"
"- astro's 'distance coincidence' factor of 5-7 is NOT applied. The like-for-like window already selects stars at comparable distances: the window hosts have a median distance of 259 pc, and 30% lie at 265-345 pc, not 15-20%.\n"
"- TRICERATOPS (SEBx2P 0.973) is not used. It equals its own prior ratio of 158 times a light-curve likelihood ratio of 0.22 from a single folded eclipse, so it adds nothing.\n"
"- astro's central 338 is superseded: it was computed at f = 5%.\n"
"- C* is corrected to +0.0181 (correct Riello branch, 2.5 sigma, about 2 sigma against unselected twins) and is used only qualitatively.\n\n"
"WHAT WOULD BE DECISIVE:\n"
"- A common-proper-motion test. Bound gives dmu < 2.5 mas/yr (95%). Among unrelated window M dwarfs, 0 of 77 lie within 5 mas/yr of the G star (5th percentile 9.5, median 35 mas/yr), so a test gives a likelihood ratio >~ 100-200.\n"
"- The M-dwarf systemic RV. Bound predicts -6.1 +- 2.0 km/s, while field stars fall within +-3 km/s with probability 7-9%: a likelihood ratio of about 10-13.")

F['eb_components'] = (
"Adopted values combine the eclipse model (lcmodel v6 + v7rep, plus the referee's independent model), the colours (localise TESS colour shift; chromatic and referee SuperWASP/W1 ratios) and the localisation. "
"The masses come from main-sequence relations, not dynamics.\n\n"
"PRIMARY (deep-eclipsed): M1 = 0.45 +- 0.03 Msun (range 0.40-0.50), drawn from:\n"
"- lcmodel free main-sequence fit 0.440 +- 0.025;\n"
"- lcmodel bound fit 0.470 +- 0.014;\n"
"- referee 0.446-0.474 at f 2.5-2.8%;\n"
"- BP-RP 2.31 (2.22-2.40, +-0.1-0.2 systematic), which maps to M2-M2.5V, 0.40-0.44 Msun.\n"
"R1 = 0.48 +- 0.02 Rsun, Teff1 = 3550 +- 100 K, M2V (M1.5-M2.5V), log g 4.73, L 0.033 Lsun. The primary is about 6% larger than a main-sequence star of its mass (dR1 +6 +- 3.5%).\n\n"
"SECONDARY: M2 = 0.20 +- 0.02 Msun (0.17-0.23), R2 = 0.21 +- 0.01 Rsun, Teff2 = 3160 +- 100 K (T2/T1 = 0.89 +- 0.01), M4-M4.5V, BP-RP about 3.0-3.2 (localise: >~3), log g 5.09, L 0.004 Lsun. "
"It supplies about 10% of the pair's TESS light. The secondary eclipse is a total occultation, so L2/L_total = 2600 ppm directly.\n\n"
"ORBIT:\n"
"- P = 1.30367228 +- 0.00000019 d; T0 = BJD_TDB 2460165.670416 +- 0.00018.\n"
"- a = 4.35 +- 0.09 Rsun (0.0202 au).\n"
"- With the main-sequence prior: i = 89.1 (+0.6/-0.65) deg (95% > 87.9), b0 = 0.10 (+0.07/-0.06) in (R1+R2) units, k = 0.44 +- 0.02 (shape alone 0.18-0.47), r1 = 0.112, r2 = 0.049.\n"
"- Robust to the prior: rsum = 0.1616 +- 0.0021; J = 0.52-0.57 (0.44-0.52 over the whole shape-allowed range).\n"
"- e cos w = 0.0006 +- 0.0005 (light curve, and from the ETV secondary offset +46 +- 33 s minus about 8 s Roemer delay). e sin w = 0.01 +- 0.01 (prior-sensitive). Both are consistent with a circular orbit.\n"
"- Synchronised rotation: spot modulation at P_EB in 5 of 6 sectors (p ~ 9e-9), semi-amplitude 280-480 ppm.\n"
"- Ellipsoidal amplitude < 74 ppm (95%), against 20-21 ppm predicted.\n"
"- Predicted K1 = 52 +- 4 km/s and K2 = 117 +- 3 km/s.\n\n"
"LIGHT AND DEPTHS: f_EB(TESS) = 2.6 +- 0.4%. Contributing estimates:\n"
"- lcmodel main-sequence chain 2.82 +- 0.20%;\n"
"- referee soft main-sequence prior: minimum 2.5%, 2.2-3.4%;\n"
"- main-sequence pair at 297 pc: 2.0-2.8% (Mamajek), 2.3-2.5% (referee), 2.70% (lcmodel bound fit).\n"
"Hard limits:\n"
"- Shape-only lower limit f_EB > 2.2% (dchi2 9), reproduced by the referee's independent model.\n"
"- The upper side is limited by main-sequence consistency (chi2_MS >= 12 at f >= 4%), the ellipsoidal limit and the colour.\n"
"- SED limit <= 5.8-11.6%.\n"
"Undiluted depths are 22% (18-26%) for the primary and 10% (8-12%) for the secondary. Main-sequence photometric distance 245-330 pc, against 297.4 +- 1.5 pc for the G star.\n\n"
"COLOUR CONSISTENCY: the M2V+M4.5V model predicts SuperWASP(V)/TESS = 0.33 for the deep eclipse (measured 0.36 +- 0.11 or 0.42 +- 0.15) and W1/TESS = 3.3 (measured 3.7 +- 1.1). The chromatic Teff of 4040 (+330/-280) K is 1.6-1.9 sigma hot; it is driven by the higher SuperWASP ratio and eased by the referee's re-reduction.\n\n"
"CAVEATS:\n"
"- k, b, f_EB, the masses and the radii are set by the main-sequence prior.\n"
"- The main-sequence solutions need 10-15% pulls in J and k: a secondary brighter than main-sequence values, or a primary inflated by about 6%. This is common in active short-period M-dwarf EBs.\n"
"- The uniform-prior geometric branch (f_EB about 7%) rests on e sin w = +0.03. It needs non-main-sequence stars (chi2_MS 24-39), and its main-sequence versions violate the ellipsoidal limit.")

F['separation_and_position'] = (
"SEPARATION POSTERIOR for a bound companion (synthesis/s02). Inputs: Raghavan period prior with isotropic orbits x TESS colour-corrected localisation map x Gaia IPD window x VHS Ks miss probability for dKs = 2.62 +- 0.18.\n"
"- Soft IPD window: median 0.48\", 68% 0.41-0.57\", 95% 0.37-0.79\".\n"
"- Hard IPD window (0.48-1.2\"): 0.54\", 68% 0.49-0.65\", 95% 0.48-0.84\".\n"
"- Without VHS: 0.58\" (95% 0.39-1.0\").\n"
"ADOPTED: 0.40-0.65\" (68%) and 0.37-0.85\" (95%), i.e. 120-195 au (68%) and 110-250 au (95%) projected at 297.4 pc.\n\n"
"LOWER BOUND: Gaia IPD (ipd_frac_multi_peak 41, BP/RP blended 47/47 and 48/48; 0 of 2650 speckle singles match; matching calibrator pairs at 0.48-1.2\"). This bound is soft because it rests on 4-11 calibrators. PMa/RV/RUWE independently exclude most orbits below ~0.03-0.1\".\n\n"
"UPPER BOUND:\n"
"- VHS miss probability at dKs 2.6: 31% at 0.5\", 6% at 0.6\", <= 7% at 0.7-0.9\", 0 at >= 1.0\". At dKs 2.8 it is 42%, 17%, 6-12% and 0. These injections are optimistic because they use the fit PSF.\n"
"- TESS: 1.24\" (95%), 1.48\" (99.7%).\n"
"The referee's 0.48-0.8\" (140-240 au) is consistent with this. Mine sits slightly closer because the VHS miss probability is applied explicitly.\n\n"
"POSITION ANGLE: only a weak NNW preference (PA ~ 340 +- 40 deg), not a measurement.\n"
"- The TESS best PA is 333-345 deg at every radius and in every localisation variant.\n"
"- On the 0.45-0.6\" ring the baseline gives 68% PA 305-12 deg and 95% 268-38 deg, with 89% posterior probability in 300-20 deg (prior 22%).\n"
"- With S105 dropped or errors x1.5 this falls to 55-77%, and the 95% range covers most PAs. p(delta = 0) = 0.09.\n"
"- Gaia gives no PA (harmonic amplitude 0.029, Rayleigh p 0.83).\n\n"
"VHS FEATURE: there is an unexplained VHS Ks residual at 1.12-1.14\", PA 345-350 deg, dKs 4.0-4.4 after bias correction.\n"
"- It cannot be the EB. The EB must have dKs <= 2.9 (f_EB > 2.2%), and injections at 1.0-1.2\" recover dKs 2.5-3.0 100% of the time and measure them 0.1-0.2 mag too bright.\n"
"- Its PA coincides with the TESS hint. This should be disclosed: it may be a fourth object or a fitting artefact.")

F['outer_orbit'] = (
"MONTE CARLO (synthesis/s02): 2e6 bound orbits with log-normal P (Raghavan; mu 5.03, sigma 2.28 in log d), e ~ U(0, 0.8), isotropic orientation, M_G 1.03 + M_EB 0.66 = 1.69 Msun, weighted by the separation likelihood. Medians, 68% and 95% ranges:\n"
"- a = 166 au (113-290; 84-660); r_3D = 174 au.\n"
"- P_out = 1600 yr (920-3800; 590-13000).\n"
"- e_out is unconstrained.\n\n"
"HIERARCHY: P_out/P_in ~ 5e5, so the system is dynamically stable. The Kozai-Lidov timescale is ~1e8-1e9 yr, while GR apsidal precession of the inner orbit takes ~4e3 yr, so Kozai-Lidov cycles are quenched. A circular, synchronised inner orbit is expected, as observed.\n\n"
"EXPECTED SIGNATURES:\n"
"(a) Relative line-of-sight velocity gamma_EB - gamma_G: 0 +- 1.8 km/s (68%); |Delta| < 3.0 km/s (95%); median |Delta| 1.25 km/s. Gravitational-redshift difference ~0.\n"
"(b) Relative proper motion: 1.4 mas/yr (0.9-2.0; < 2.5 at 95%). Sky-plane relative velocity 2.0 km/s. PA change 0.12 deg/yr (0.04-0.22).\n"
"(c) Eclipse timing of the EB: LTTE semi-amplitude ~11 h over ~1600 yr. Over decades it appears only as curvature, 0.002 min/yr^2 (68% 0.0008-0.0034; 95% < 0.0047), with |dP/dt| <= 0.002 s/yr. Measured limits are |c| < 0.077 min/yr^2 and |dP/dt| < 0.032 s/yr, 15-40x larger. Detection would need several decades of timing.\n"
"(d) G-star RV trend: 1.5 m/s/yr (68% 0.6-2.5; 95% < 3.4). This is detectable only with a decade of precision RVs, after modelling M-dwarf line blending at P_EB.\n"
"(e) G-star astrometric acceleration: 0.0025 mas/yr^2 (95% < 0.006). This is below Gaia DR3 sensitivity (TGAS-DR3 PMa p = 0.46, RUWE 1.026, no NSS solution) and probably below DR4 (~0.01 mas/yr^2).\n\n"
"CONSISTENCY WITH CURRENT DATA: the non-detections so far are what this orbit predicts: no ETV curvature, no LTTE (best p = 0.37), no PMa, constant Gaia RV, low RUWE. They exclude only much tighter orbits: ETV <~ 10-15 au (>= 84-99% of orbits); PMa/RV/RUWE <~ 5-10 au.")

F['system_table'] = """| Parameter | Value | Source / note |
|---|---|---|
| **G star**: TIC 71316387 = TYC 8001-154-1 = Gaia DR3 6571343840265900544 | RA 330.887948, Dec -42.835511 | |
| G, BP-RP, Tmag, V, Ks | 11.775, 0.786, 11.358, 11.857 (APASS), 10.341 | blended; the pair supplies ~1.4% (G), 2.6% (T), 8.5% (Ks) |
| Parallax, distance | 3.3622 +- 0.0170 mas, 297.4 +- 1.5 pc | Gaia DR3 |
| Proper motion | (38.50, -3.13) mas/yr | Gaia DR3 |
| RV | -6.07 +- 0.86 km/s, constant (15 transits over 946 d, p = 0.275) | excludes the G star as an EB component |
| Teff, R, M, age | 5880 K, 1.07 Rsun, 1.03 (0.99-1.07) Msun, 5.1 (4.0-6.2) Gyr | GSP-Phot/FLAME; GSP-Spec MG 6235 K is discrepant |
| RUWE, PMa, NSS | 1.026; TGAS-DR3 p = 0.46; none | no acceleration |
| **EB ephemeris** P | 1.30367228 +- 0.00000019 d | lcmodel (chromatic 1.30367220 +- 0.00000022); CTOI P = P_EB/2 |
| T0 (primary) | BJD_TDB 2460165.670416 +- 0.00018 | BTJD 3165.670416 |
| Diluted depths (TESS) | 5770 +- 100 ppm (deep), 2600 +- 55 ppm (shallow) | +-8% between data products |
| f_EB (TESS flux fraction) | 2.6 +- 0.4% (hard lower limit 2.2%) | MS-prior; SED <= 5.8-11.6% |
| Undiluted depths | 22 (18-26)%, 10 (8-12)% | |
| M1, R1, Teff1, type | 0.45 +- 0.03 Msun, 0.48 +- 0.02 Rsun, 3550 +- 100 K, M2V (M1.5-M2.5V) | MS relations, not dynamical; BP-RP ~2.3 |
| M2, R2, Teff2, type | 0.20 +- 0.02 Msun, 0.21 +- 0.01 Rsun, 3160 +- 100 K, M4-M4.5V | BP-RP ~3.0-3.2 |
| a, i, b0 | 4.35 +- 0.09 Rsun (0.0202 au), 89.1 (+0.6/-0.65) deg, 0.10 (+0.07/-0.06) | i and b0 prior-driven |
| k, r1+r2, J, T2/T1 | 0.44 +- 0.02, 0.1616 +- 0.0021, 0.52-0.57, 0.89 +- 0.01 | rsum and J robust |
| e cos w, e sin w | 0.0006 +- 0.0005, 0.01 +- 0.01 | circular within errors |
| Rotation, ellipsoidal | synchronised (spots at P_EB, 280-480 ppm); ellipsoidal < 74 ppm (pred. 21) | |
| K1, K2 (predicted) | 52 +- 4, 117 +- 3 km/s | |
| vsini (synchronised) | 18.6 +- 1.0, 8.2 +- 0.5 km/s | |
| Pair apparent mags | T 15.3 (15.0-15.5), G 16.4 (16.2-16.6), Ic 15.0-15.2, Ks 13.0 (12.8-13.2) | |
| MS photometric distance of the pair | 245-330 pc | matches 297.4 pc |
| **Outer orbit (if bound)**: separation | 0.48\" (68% 0.40-0.65\"; 95% 0.37-0.85\"); <= 1.24\" (95%, TESS) | |
| Projected separation | 120-195 au (68%), 110-250 au (95%) | |
| PA | weak NNW preference, ~340 +- 40 deg | TESS; not significant |
| a, P_out | 166 (113-290) au, 1600 (920-3800) yr | Monte Carlo |
| Contrasts (pair - G star) | dKs 2.4-2.8, dH 2.6-2.9, dJ 2.9-3.2, d832nm 3.7-4.0, dIc 3.7-4.1, dG 4.4-4.8, d562nm 5.1-5.6, dV 5.2-5.6 | +-0.2 mag systematic |
| gamma_EB - gamma_G, relative PM | 0 +- 1.8 km/s (95% < 3.0); 1.4 (95% < 2.5) mas/yr | |
| Posterior odds bound : chance-aligned | ~8e3 (5e3-1.2e4); conservative >= 1.5e3; pessimistic any-host floor 25-36 | synthesis s02/s04 |"""

F['predictions_for_followup'] = (
"(a) SPECKLE / AO. The companion is the unresolved EB (the inner pair spans 0.07 mas).\n"
"- Expected position: 0.40-0.65\" (95% 0.37-0.85\"), PA weakly ~340 deg.\n"
"- Contrasts: d832nm = 3.7-4.0, dIc (SOAR I) = 3.7-4.1, d562nm = 5.1-5.6, dJ = 2.9-3.2, dH = 2.6-2.9, dKs = 2.4-2.8 (referee Mann 2.57-2.77; lcmodel 2.62-2.80; +-0.2 mag systematic).\n"
"- Expect an 8-m speckle detection at 832 nm but probably not at 562 nm. SOAR I-band is borderline. Near-IR AO should detect it easily.\n"
"- Disconfirming outcomes: a companion with dKs <= 2.5 beyond 0.6\" would contradict VHS. Deep non-detections (dKs ~3.5 at >= 0.15\", d832 ~4.5 at >= 0.2\") would push the EB inside ~0.15\" or refute the pair model.\n"
"- Common proper motion: bound predicts dmu < 2.5 mas/yr (95%), i.e. <= 5 mas of relative motion over 2 yr. Unrelated window M dwarfs: 0 of 77 within 5 mas/yr, 5th percentile 9.5 mas/yr, median 35 mas/yr. A second epoch 1-2 yr later therefore gives a likelihood ratio of ~100-200.\n\n"
"(b) SPECTROSCOPY. Use red or near-IR. The pair supplies ~0.7% (V), 3% (I), 5.5% (J), 7% (H) and 8.5% (Ks) of the light, and the secondary ~10-15% of the pair.\n"
"- Expect an SB2: K1 = 52 +- 4 km/s, K2 = 117 +- 3 km/s, q = K1/K2 = 0.43 +- 0.03, (M1+M2) sin^3 i ~ 0.65 Msun.\n"
"- vsini1 = 18.6 +- 1, vsini2 = 8.2 +- 0.5 km/s (synchronised, i ~ 89 deg).\n"
"- Systemic gamma_EB = -6.1 +- 2.0 km/s if bound (95%: -10 to -2). Field stars fall within +-3 km/s of -6.1 with probability 7-9%, so the likelihood ratio is ~10-13.\n"
"- The G-star lines should stay at -6.1 km/s with no 1.3-d signal. The M-dwarf H-alpha emission is diluted below detectability in the blend.\n\n"
"(c) MULTI-COLOUR ECLIPSE DEPTHS (deep / shallow, ppm, scaled to the TESS 5770/2600; +-15% optical, +-20% NIR):\n"
"- g' 1150/200; V 1900/370; 562 nm 2000/400; r' 2800-3350/580-890; Gaia G 3600/1240; i' 5600-6700/2200-3000; TESS 5770/2600; z' 7200/3800; J 12100/8900; H 15800/11400; Ks 18000/13800; W1 18900/15500.\n"
"- The shallow/deep ratio rises from ~0.2 (g') to 0.45 (TESS) to ~0.8 (Ks).\n"
"- Shapes: the secondary eclipse should be total, with a flat bottom of ~36 min in a ~1.6-h eclipse. The primary eclipse has a curved bottom.\n"
"- Timing: secondary eclipse at phase 0.5 + 46 +- 33 s.\n"
"- A G/F eclipsed star would instead give V/TESS >= 0.82 with near-achromatic depths. Chromatic depths cannot tell a bound M-dwarf EB from a chance-aligned one.\n\n"
"(d) GAIA DR4 (2 Dec 2026).\n"
"- G star: no significant acceleration (predicted <= 0.006 mas/yr^2 at 95%), RUWE ~1, no orbit or acceleration solution, parallax and PM unchanged.\n"
"- The IPD multi-peak and BP/RP blending flags should persist.\n"
"- The companion (G ~16.2-16.6, dG 4.4-4.8 at ~0.5\") is unlikely to be catalogued. If it is, expect parallax ~3.36 mas and PM within 2.5 mas/yr of (38.5, -3.1).\n"
"- Epoch photometry should contain the eclipses at P = 1.3036723 d: ~3.6/1.2 mmag in G, ~6.0/2.6 in RP, ~1.8/0.3 in BP. Detection is marginal, with only a few in-eclipse transits.\n"
"- RVS epoch RVs should be constant, with no 1.3-d variation larger than ~1 km/s.\n"
"- Disconfirming outcome: a DR4 acceleration >~0.03 mas/yr^2 or an SB1 orbit would require an additional, closer body.")

F['claim_wording'] = (
"STRONGEST DEFENSIBLE CLAIM: The TESS signal on TIC 71316387 (CTOI 71316387.01, P = 0.6518 d) is an eclipsing binary at twice that period, P = 1.3036723 d, with unequal eclipses of 5770 and 2600 ppm. "
"The G star is not an eclipsing component. Its Gaia RVS velocity is constant, but at least 18 km/s of variation would be required. "
"After correcting for the colour-dependent PSF shift of the TESS optics, the eclipsing source lies within 1.24\" of the G star (95%). "
"Its eclipses are redder than the G star in the TESS centroid colour shift and in the SuperWASP and NEOWISE depth ratios. "
"Under main-sequence assumptions, the eclipse shape and colours describe an early-M (about M2V, ~0.45 Msun) + mid-M (about M4.5V, ~0.2 Msun) pair. "
"The pair supplies ~2.6% of the TESS flux, and its photometric distance (~245-330 pc) matches the G star's 297 pc. "
"Gaia image parameters also indicate an unresolved companion at ~0.5\". "
"Given these, an unrelated M-dwarf EB this close to the G star is improbable, and a bound hierarchical triple (G2V + M-dwarf EB at ~120-200 au projected) is ~10^3-10^4 times more probable than a chance alignment. "
"The association is not yet confirmed: there is no resolved image, common-proper-motion measurement or radial velocity of the M dwarfs, and the component masses are not dynamical.\n\n"
"EXOFOP NOTE (480 characters): CTOI 71316387.01 is an EB, not a planet: P=1.3036723 d (2x CTOI period), unequal eclipses 5770/2600 ppm. The G star is not the EB (Gaia RV constant; K>=18 km/s needed). Colour-corrected TESS centroids put the source within 1.24\" (95%); eclipse colours + shape imply an ~M2V+M4.5V pair giving ~2.6% of TESS flux (dKs~2.6, d832nm~3.8, likely 0.4-0.65\"). Probably a bound triple (odds ~1e3-1e4 vs chance alignment) but unconfirmed: speckle/AO + 2nd epoch (CPM) or M-dwarf RVs needed.")

F['unresolved_issues'] = (
"1. No direct evidence of association. The bound case is a prior-odds argument plus consistency (position, colour, light-fraction distance). It needs imaging with common proper motion, or the M dwarfs' systemic RV.\n\n"
"2. f_EB and everything that follows from it are prior-driven: masses, k, b and the photometric distance. The light-curve shape gives only f_EB > 2.2%. With e sin w free, uniform priors favour f_EB ~ 7% (geometric L3 0.932). That branch is rejected only because it requires non-main-sequence or strongly inflated M dwarfs (chi2_MS 24-39), because of the ellipsoidal limit, and because of the colour. A sceptic invoking young or inflated M dwarfs would place a main-sequence-luminosity pair at 160-200 pc, i.e. not bound.\n\n"
"3. All main-sequence solutions need 10-15% pulls in J and k: a secondary brighter than main-sequence values, or a primary inflated by ~6%.\n\n"
"4. The chromatic Teff of the eclipsed star, 4040 (+330/-280) K, sits 1.6-1.9 sigma above the adopted 3550 K.\n\n"
"5. The NEOWISE W1 shallow-eclipse ratio (0.06 +- 1.88) is ~3 sigma below the predicted ~6, and W1 and W2 disagree on the shallow eclipse. The SuperWASP shallow eclipse (-622 +- 590 ppm) is ~1.7 sigma below the predicted +370 ppm. All rest on sparse archival data.\n\n"
"6. The separation bounds are soft. The lower bound rests on 4-11 Gaia IPD calibrators. The upper bound rests on optimistic VHS injections that used the fit PSF, so the true separation could reach ~1.0\".\n\n"
"7. The NNW PA hint depends mainly on S105 and on the red end of the TESS colour calibration (only 4 stars with BP-RP > 2.7). The VHS feature at 1.12\", PA 345-350 deg (dKs ~4.2) is unexplained.\n\n"
"8. Prior-model systematics of ~x2-3: the Raghavan period distribution and a flat q applied to a companion that is itself a close binary; an assumed equal EB rate with no tertiary enhancement; Gaia field counts. The referee's 7-host count and my 77-host recount differ by 20%.\n\n"
"9. Depth scales differ by 8% between SPOC S1/S28 and S68 (crowding), adding about +-5% to f_EB.\n\n"
"10. ETV errors needed a 12.5% inflation (possibly spot-induced timing noise). e sin w is marginally non-zero.\n\n"
"11. The Gaia GSP-Spec MatisseGauguin Teff (6235 K) conflicts with ANN (5843 K) and GSP-Phot (5883 K). It is not used.\n\n"
"12. TRICERATOPS FPP values are not meaningful here: it fits a single folded eclipse and ignores RV, colour and localisation. They should not be quoted.\n\n"
"13. The masses are not dynamical. Spectroscopic K1/K2 would give (M1+M2) sin^3 i and q directly.\n\n"
"FILES (all new, under /tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/synthesis/):\n"
"- s01_contrast_depths.py / .json: contrasts and multi-band depths;\n"
"- s02_odds_sep_orbit.py / .json: position factor, L(f) factor, odds, separation posterior, outer-orbit Monte Carlo;\n"
"- s03_field_pm.py and field_Mdwarfs_1p5deg.ecsv: Gaia query of field M dwarfs;\n"
"- s04_field_stats.py / .json: host densities, PM and distance statistics;\n"
"- s05_pa.py / .json: PA on the separation ring;\n"
"- s06_write_final.py: writes this answer.\n"
"Draft: /tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/synthesis_draft.json")

note = F['claim_wording'].split('EXOFOP NOTE (480 characters): ')[1]
assert len(note) < 600, len(note)
F['status'] = 'DRAFT v3: s06 text written by previous run; being verified by current run'
json.dump(F, open(B + 'synthesis_draft.json', 'w'), indent=1)
json.dump(F, open(B + 'synthesis/synthesis_final.json', 'w'), indent=1)
print('exofop note chars', len(note))
for k, v in F.items():
    print(k, len(v))
