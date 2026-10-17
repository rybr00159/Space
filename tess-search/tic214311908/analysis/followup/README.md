# Follow-up observing design for TIC 214311908.01 (folder `followup/`)

Folder: `tic214311908/analysis/followup/`

AI-run analysis, 2026-10-17. Nothing was posted, uploaded or submitted, and no one was contacted.

Network use was read-only:
* the Gaia archive (ADQL);
* SVO filter curves (GET);
* VizieR (Claret & Bloemen 2011 limb-darkening table);
* the NASA Exoplanet Archive TAP (GET);
* arXiv PDFs and ESA, ESO, LCO and Gemini web pages, used to verify instrument and literature numbers.

Python: `/home/user/work/venv/bin/python`, with OMP/MKL/OPENBLAS threads set to 1. The `../referee2_fpp` code was
imported read-only (no bytecode was written). No background jobs are left running. The folder is about 11 MB.

**Inputs used:**
* `../clean/out/ephemeris_final.json` (adopted "combined" block). The scripts fall back to `../phot/out/timing.json`
  and accept any ephemeris file as an argument.
* `../clean/out/final_depth.json`: depth 15.2 ppt, Rp 4.90 R_E.
* `../clean/out/per_transit_correlates.csv`: the Mira's TESS flux at each transit.
* `../referee2_fpp/out/r04_companions.npz` and `scripts/r05_fpp.py`: companion draws and the FPP framework.
* `../bd/` (b10): the photo-eccentric sign fix, and f_BD = 1.3e-4.
* `../star/` (catalogue downloads, BT-Settl models, filters), `../shape/data` (CK04, Claret 2017), `../ctx/` (Gaia
  crossings, Mira epochs).

## Bottom line for an observer

1. **The decisive first observation is a seeing-limited, multi-band ground transit on the M3V.**
   * Instrument: preferably MuSCAT4 on the 2-m FTS at Siding Spring, with simultaneous g′ r′ i′ z_s. MuSCAT4 was
     installed there in October 2023 (LCO pages).
   * What one transit does:
     * confirms that the 15-ppt dip is astrophysical and sits on the M3V. This is the risk the FPP does **not**
       cover: Mira, scene and pipeline systematics, the INV odd-even of 2.1σ, and a 1.6σ depth-epoch contrast;
     * resolves pair A/B (11.7″) and the 3.1″ and 5.5″ stars;
     * reaches about 5% on the r/z or g/z depth ratio.
   * Best 2027 MuSCAT4 event: E126, 2027-06-01 12:55 UTC (Moon 14% lit).
   * In Chile, the five sites see full windows on E127 (2027-06-07 05:30 UTC, Moon 8%) and E137 (2027-08-03
     03:30 UTC, Moon 1%).
   * **No transit is observable before April 2027.**
2. **Speckle imaging with Gemini-S/Zorro (832 + 562 nm)** excludes about 97% of the residual planet-on-companion
   (PTP) posterior. SOAR/HRCam excludes 50–65%; VLT/ERIS LGS AO in Ks about 50–77%. SPHERE is not suited at G = 14.3.
3. **Two to four low-precision RV epochs (σ ≤ 0.5 km/s)** close the brown-dwarf channel (K = 2.8–15.6 km/s) and the
   white-dwarf channel (K ≈ 73 km/s). The FPP is then about 3e-5. That floor is the background term, which r05
   carries at a fixed 3e-5 and does not re-derive.
4. **A precise RV mass:**
   * expected mass 21 (16–27) M_E, K = 14.8 (11.5–19.0) m/s;
   * ESPRESSO: 9–26 epochs for a 5σ K at jitter 5–10 m/s;
   * NIRPS+HARPS: 16–33 epochs, scheduled June–November, because January–May puts the star's lines on the OH lines.

**FPP ladder** (r05 framework, sign-fixed per `../bd`; REF / INV light curve):

| Step | FPP (REF / INV) |
|---|---|
| Now | 1.4e-3 / 3.8e-3 (as published: 3.1e-3 / 6.6e-3) |
| + Zorro speckle | 6.0e-4 / 1.0e-3 |
| + two-band transit at 5% (pair resolved) | 1.5e-4 / 4.8e-4 |
| + both | 1.4e-4 / 4.5e-4 |
| + low-precision RV | 3e-5 |
| + precise RV mass | ~4e-5 (bookkeeping floor) |

## 1. Transit predictions (f01)

**Method.**
* Ephemeris from `clean/out/ephemeris_final.json`:
  * T0(E = −26) = 2460692.935829 ± 0.000912 BJD_TDB;
  * P = 5.6915020 ± 0.0000078 d;
  * covariance −1.4e-11 (ρ = −0.002).
* 1σ and 3σ timing errors come from the full covariance.
* BJD_TDB is converted to UTC at each site by iterating t_bjd = t_obs + ltt_bary(t_obs, site) in astropy, with the
  target position propagated using the Gaia DR3 proper motion.

**Checks.**
* The barycentric light-time ranges from −432.7 to +449.6 s over 2026–2028. The ecliptic latitude is −28.1°, so the
  expected amplitude is 499 s × cos β = 440 s.
* TDB − UTC = 69.185 s. This assumes no new leap second before 2029.
* Proper motion changes the light-time by at most 3.5 ms.
* Topocentric terms are below 0.03 s.

**Predictions.**

| Quantity | Value |
|---|---|
| Transits | 139 (epochs 89–227) between 2026-11-01 and 2028-12-31 |
| σ (1σ), 2027 | 1.9–2.5 min |
| σ (1σ), 2028 | 2.5–3.1 min |
| 3σ | 5.6–9.4 min |
| Full window | T14 1.70 h + 2 × 3σ + 2 × 1.5 h = 4.88–5.01 h |

**Comparison with `ctx/out/followup_transits.csv`.**
* ctx used ephemeris v2: P = 5.6914968 ± 0.0000145, T0 = 3840.9153 ± 0.0019 (χ²-inflated).
* Our BJDs are later by 0.04–0.62 min, which the change of ephemeris explains (expected 0.07–1.1 min over 2026–2028).
* Our σ is 0.47–0.49 times ctx's.
* Converting ctx's own BJDs with this code reproduces its `utc_mid` to within its minute truncation (0.0–1.0 min). So
  ctx's barycentric UTC is confirmed. The naive no-barycentre conversion would be off by −7.1 to +8.3 min.
* We find fewer "full" windows: ctx required Sun < −12° with ±1 h baselines and T14 1.63 h; we require Sun < −18°
  with ±1.5 h + 3σ.
  * ctx epoch 130 drops out at CTIO/La Silla, and epochs 121 and 125 at SAAO.

**Observable full windows** (airmass < 2 and Sun < −18° over the whole window; `out/obs_summary.csv`):

| Site | 2026-11 to 2027-12 | 2028 |
|---|---|---|
| CTIO / Cerro Pachón (SOAR, Gemini-S) / Las Campanas / La Silla | 3 (E124 May 21, Moon 99%; E127 Jun 7, Moon 8%; E137 Aug 3, Moon 1%) | 4 |
| Paranal | 2 (E127, E137) | 3 |
| SAAO Sutherland | 4 (E118 Apr 17; E128 Jun 12; E131 Jun 29, Moon 22%; E138 Aug 8) | 3 |
| Siding Spring (FTS/MuSCAT4, LCOGT coj) | 3 (E126 Jun 1, Moon 14%; E129 Jun 18, full Moon 24° away; E139 Aug 14, Moon 93%) | 2 |
| Dome C (ASTEP+, target circumpolar at 36–66° altitude) | 6 (Sun < −18) / 15 (Sun < −12) | 7 / 10 |

* "Good partial" windows add 0–4 per site per year. These cover the transit ± 3σ with at least 45 min of baseline on
  each side.
* No window is usable in November–December 2026 (fraction 0).
* The first usable event is 2027-04-05 15:02 UTC (partial, Siding Spring and Dome C); the first full one is 2027-04-17
  00:12 UTC (SAAO).

**Gaia epochs.**
* With the final ephemeris, the Gaia field-of-view crossing at 2024-10-29 18:20 UTC falls +14.9 ± 1.3 min after
  mid-transit, inside the flat bottom (ctx: +13.9 ± 2.5 min).
* The 2020-02-14 crossing falls at −39.2 ± 3.8 min, on the ingress slope near the bottom (~85% of the depth).
* Both are in the DR5 window, not DR4.

**Files:**
* `out/transits_all.csv`: every epoch, with BJD, σ, geocentric UTC and barycentric correction;
* `out/transits_obs_all_sites.csv`: every epoch × site, with window, flags, baselines, altitude, Sun and Moon;
* `out/obs_<site>.csv`: per site, epochs with any usable time;
* `out/obs_summary.csv`, `out/f01_transits.json`, `figs/observability.png`.

## 2. Magnitudes, the Mira, contamination, precision (f02–f06, f10)

### Magnitudes (`out/mags_table.csv`, `out/mags_adopted.json`)

| Star (2027.5 separation, PA) | G | g | r | i | z (AB) | V | Ic | J |
|---|---|---|---|---|---|---|---|---|
| M3V | 14.28 | 16.36 | 14.98 | 13.57 | 12.81 | 15.59 | 12.86 | 11.36 |
| Pair A (11.7″, 94°) | 14.65 | 16.44 | 15.31 | 13.83 | 12.88 | 15.95 | 13.09 | ~11.6 |
| Pair B (11.4″, 98°) | 15.74 | | | | | | | |
| 3.1″ star (3.21″, 54°) | 17.48 | 18.16 | 17.35 | 17.10 | 16.90 | 17.83 | 16.64 | ~16.0 |
| 5.5″ star (5.32″, 33°) | 17.68 | 18.09 | 17.59 | 17.41 | 17.30 | 17.87 | 17.02 | ~16.4 |
| Mira (19.08″, 92°) | 12.66 (mean) | 13.4–17.1 | ~12.6–15.6 | 10.3–12.1 | ~9.3–10.6 | 12.91–17.01 | ~9.9–11.6 | ~7.6–8.8 |

**Sources and caveats.**
* M3V and pair A: XP synthetic photometry (gaiadr3.synthetic_photometry_gspc).
* 3.1″ and 5.5″ stars: DECaPS2 griz, and the Riello et al. (2021) G→V/Ic relations (no XP photometry).
* Pair B: Gaia G only. In ground seeing it is blended with A; the measured SkyMapper/DECaPS/VHS rows are A+B.
* M3V J: 2MASS. VHS saturates for the M3V, the Mira and pair A. VHS J/Ks for the 3.1″ and 5.5″ stars are flagged
  unreliable (halo); J is estimated as z − 0.9 ± 0.3.
* The M3V moves 1.6″ between 2016 and 2027. The 2027.5 separations are in `out/geometry_2027.csv`.

**Mira ranges.**
* V: VSX 12.91–17.01.
* G / BP / RP: Gaia DR3 39 epochs, 11.35–13.32 / 13.48–17.19 / 9.98–11.73.
* SkyMapper max − min: g ≥ 2.98, i ≥ 1.14, z ≥ 0.66.
* DENIS 1998 versus 2MASS 1999: ΔJ ≥ 0.48, ΔKs ≥ 0.53.
* The r, z, J and Ks amplitudes are assumptions, interpolated in wavelength.

### Mira phase (f03)
* Sinusoid fitted jointly to the Gaia G epochs (2014–2017) and the TESS Mira flux and in-sector slopes at 21 transits
  (2019–2026): **P = 313.5 ± 0.6 d** (formal), G semi-amplitude 0.87 mag.
  * For comparison: VSX 294 d; Gaia DR3 LPV 326 ± 46 d.
  * Other periods are excluded: 267.75 d at Δχ² = 322, 249.5 d at Δχ² = 458.
* Last maximum: BJD 2461074 (2026 February).
* The phase uncertainty at the 2027 events is ±0.06–0.07. This includes an **assumed** 0.05-per-cycle random-walk
  jitter.
* **2027 April to mid-July:** the Mira is near minimum (phase 0.36–0.69, G 12.5–13.0). It brightens to G ≈ 12.3–12.4
  in late July and 12.0–12.2 by mid-August.
* **2028 June–August:** the transits fall near maximum (phase 0.74–0.99, G 11.4–12.3).
* Files: `out/mira.json`, `out/mira_phase_at_transits.csv`, `figs/mira_fit.png`.

### Contamination (f05; Moffat β 2.5–4, FWHM 1–2.5″, aperture radius 1.5–4″)
**Contaminant budget** (c = contaminating flux / M3V flux in the aperture):
* **The Mira's wings are negligible.**
  * 2e-5 to 3e-5 at FWHM 1″ with r = 1.5″.
  * 1e-4 to 2e-4 at FWHM 1.5″ with r = 2″.
  * At most 4.1e-3 (J) at FWHM 2.5″ with r = 4″, at Mira maximum.
  * A pessimistic 5% r⁻³ aureole raises totals by at most ~1%.
* **Pair A/B:** at most 2.3e-3.
* **The 3.1″ star dominates**, especially in g:

  | FWHM / aperture radius | g | i | z |
  |---|---|---|---|
  | 1.5″ / 2″ | 1.5% | 0.30% | 0.18% |
  | 2″ / 3″ | 7.3% | 1.5% | 0.9% |
  | 2.5″ / 4″ | 13% | 2.7% | 1.6% |

**Differential dilution between bands** (bias on the depth ratio before correction):
* ≤ 1.4% (g/z) for r ≤ 2″ and FWHM ≤ 1.5″;
* up to 14% (g/z) for r = 4″ at 2.5″ seeing.

The contamination is correctable from the known magnitudes. For the chromatic test, use small apertures or PSF
photometry. Files: `out/contamination_grid.csv`, `out/contamination_summary.json`.

**Geometry warning.** The Mira, pair A/B and the M3V lie on a nearly E–W line (PA 92–98°). Any E–W diffraction spike
or row-wise bleed from the saturated Mira runs through the pair towards the M3V. Orient the CCD columns or the spider
so that bleeds and spikes avoid PA 90° ± 10°.

### Saturation (f06, FWHM 1.0″ worst case, 70% of full well; `out/saturation.csv`)
* The Mira at maximum saturates in i/z within:
  * 1.1–1.9 s on the 1–2-m telescopes (SPECULOOS, MuSCAT4, TRAPPIST I+z);
  * 4.9–5.6 s on the LCOGT 1-m;
  * 1.9–9 s on the 0.4-m telescopes.
* So the Mira will be saturated in any transit exposure. This is harmless only if bleeds and spikes avoid the M3V.
* The M3V itself reaches 70% of full well, in focus, in:
  * 25–38 s in i/z/I+z on 1–2-m telescopes;
  * 70–140 s on the LCOGT 1-m.
* Defocus to 2–3″ in i/z, or keep exposures at or below those values.

### Precision and transit S/N per event (f06, calibrated by f10; `out/precision.csv`)
**Noise model.**
* Photon, sky and read noise.
* Scintillation from Osborn et al. (2015) eq. 7, verified: C_Y = 1.42 (Tololo) and 1.56 (Paranal); 1.5 assumed
  elsewhere; 1.0 assumed at Dome C; quartile variant ×1.3.
* Per-event depth error: the white term × 1.42 (f10: a fitted quadratic baseline inflates the depth error by 1.40–1.48×
  over the analytic white term), added in quadrature to an **assumed** red-noise term of 0.5 ppt (1–2 m) or 0.8–1.0
  ppt (≤ 0.6 m).
* f10 shows that an Ornstein–Uhlenbeck process with τ = 30 min and rms s adds 0.9–1.0 s to the per-event depth error.
* Throughputs, sky and some detector parameters are **assumptions, ±30%**.

**Results per event:**

| Telescope / band | σ per 2 min (ppt) | Depth error (ppt) | S/N per event (red ×2) |
|---|---|---|---|
| MuSCAT4 g / r / i / z_s | 2.0 / 1.16 / 0.80 / 0.83 | 0.73 / 0.59 / 0.54 / 0.55 | 21 / 26 / 28 / 28 (13–15) |
| LCOGT 1-m g / r / i / z_s / V / Ic | 4.4 / 2.2 / 1.5 / 1.6 / 3.8 / 1.5 | | 12 / 20 / 24 / 23 / 14 / 24 (10–14) |
| SPECULOOS-S I+z | 0.96 | | 27 (15) |
| TRAPPIST-S I+z / z / V | | | 17 / 15 / 8 |
| ExTrA (8″ / 4″ fibre) | | | 20 / 22. The fibre includes the 3.1″/5.5″ stars and cannot test host identity |
| ASTEP+ blue / red | | | 9 / 14 |
| Amateur 0.4-m Ic / r | | | 11 / 8 |

**One transit gives at least 8σ with any of these.**

### Which band and telescope for which test (f04, f06, f07b, f09)

**(a) Host identity.** Any telescope with FWHM ≲ 2″ and pixels ≤ 0.65″ resolves the 3.1″ and 5.5″ stars and the pair:
MuSCAT4, LCOGT 1-m, SPECULOOS, TRAPPIST.
* If the 3.1″ star were the source, it would show a **52%** eclipse of its own (grey EB assumed).
* A 1.5″ aperture on the M3V would then show at most ~0.1% instead of 1.5%.
* Blended in a large aperture, that source would give depths of 8.4% (g), 5.3% (r), 2.0% (i), 1.2% (z), so it is also
  strongly chromatic.
* ASTEP+ (1.39″/px) and ExTrA (8″ fibre) cannot resolve 3.1″.

**(b) Achromaticity.** Predicted depth ratios relative to TESS:

| Hypothesis | Predicted ratio |
|---|---|
| Planet on the M3V | k achromatic. Raw mid-transit depth/TESS: 1.065 (g), 1.058 (r), 1.024 (i), 1.03 (z_s), 0.96 (J), from limb darkening (Claret & Bloemen 2011 PHOENIX quadratic; Claret 2017 for TESS). Fit k with LD, do not compare raw depths. |
| Unocculted cool spots (real planet, systematic) | g/z_s 1.01–1.09 for filling factor 5–20% and spot Teff 2950–3150 K. Opposite sign to PTP. |
| PTP (planet on an unresolved ~0.2 Msun companion) | Posterior (sign-fixed REF, 16–84%): g/z_s 0.785 (0.72–0.85); r/z_s 0.75 (0.68–0.83); V/Ic 0.80 (0.73–0.86; referee2 got 0.75, 0.68–0.81, from Mann relations); ASTEP blue/red 0.78 (0.71–0.85). The 97.5% upper end is 0.86–0.89. BT-Settl [M/H] = 0 models; model systematics about ±0.05. |
| White dwarf occulted by the M3V | Teff_WD ≈ 6200 K for R_WD = 0.0125 Rsun (5560–7190 K for 0.015–0.010). g/TESS = 6.0 (4.9–7.4); V/Ic = 4.8. |
| Background K-dwarf EB, unresolved | g/z_s 2.7 (E(B−V) = 0.2) to 1.9 (0.4); V/Ic 2.6–2.1. |

**Two-band precision** (`out/depth_ratio_tests.csv`, `out/two_band_ptp_rejection.csv`):
* **MuSCAT4 (simultaneous):**
  * σ(R) per event 5.0% (r/z, r/i) and 5.7% (g/z, g/i), so **5% needs about 1 event (r/z) or 1.3 (g/z)**;
  * about 4 events if the red noise is twice as large;
  * the PTP posterior is rejected at more than 3σ for 93% after 1 event and 100% after 2 (r/z).
* **LCOGT 1-m (bands on alternate events):** σ(R) 6.4% (r/z) to 8.8% (g/z) per pair of events, so 3–6 events for
  5%; PTP is 74% rejected after 2 events and 98.5% after 4 (r/z).
* **ASTEP+ blue/red (simultaneous):** 12% per event, so about 6 events for 5% (PTP 69% rejected after 4).
* **TRAPPIST V/I+z:** 14 events.
* **Amateur r/Ic:** 18 events.
* The white dwarf and the background K-dwarf EB are rejected by a single g- or V-band event.
* Recommendation: **MuSCAT4 r′ and z_s (with g′ as a free bonus) is the best combination.** A Chilean 1-m pair
  (LCOGT r′/z_s) plus SPECULOOS I+z on the same night is the fallback.

## 3. High-resolution imaging (f07)

**Method.**
* PTP companion posterior: referee2_fpp draws weighted by the r05 PTP evidence. In the sign-fixed REF case:
  * projected separation 0.047 / 0.10 / 0.49″;
  * dI 0.90 / 1.15 / 1.42;
  * dV 1.21 / 1.55 / 1.91;
  * dKs 0.71 / 0.90 / 1.11.
* Curves:
  * **SOAR/HRCam**, verified: Ziegler et al. 2026 (arXiv:2608.26250) Table 2, 5σ ΔI for 2,684 TOIs. Median:
    0.1″ 2.04, 0.2″ 4.04, 1″ 4.75, 1.5″ 5.28, 3″ 6.19. 16th percentile: 1.90, 3.43, 3.95, 4.31, 4.91. The part inside
    0.1″ is assumed.
  * **Zorro**, verified statements, curve shape assumed: inner working angle 0.017″ (832 nm) and 0.026″ (562 nm);
    4–5 mag contrast near the diffraction limit; 5–9 mag out to 1.2″; poorer for V ≳ 15; about −2.5 mag when the
    seeing doubles (Howell & Furlan 2022; Howell et al. 2025).
  * **ERIS**: the NGS mode needs R < 11 (the target has RP 13.07), so use LGS with the target as the tip-tilt star
    (R 7–18 allowed; ESO pages). The curve is assumed.
  * **SPHERE**: SAXO degrades strongly for G > 12–12.5 (Jones et al. 2022). Not recommended.
  * **NaCo** is decommissioned (ESO).

**Fraction of the PTP posterior excluded** (`out/hri_exclusion.csv`; REF and INV, as published and sign-fixed, all
within 0.01 of each other):

| Instrument | PTP excluded | All companions |
|---|---|---|
| Zorro 832 nm nominal | **0.97** | 0.69 |
| Zorro 832 nm, poor seeing | 0.95 | |
| Zorro 562 nm | 0.89 | |
| HRCam, median curve | 0.65 | |
| HRCam, 16th percentile | 0.63 | |
| HRCam, pessimistic | 0.50 | |
| ERIS LGS, nominal | 0.77 | |
| ERIS LGS, pessimistic | 0.51 | |
| SPHERE, degraded | 0.47 | |

* What is left after Zorro is companions projected inside 0.02–0.03″ (1.9–3.0% of the PTP weight).
* **FPP after Zorro alone:** 6.0e-4 (REF) / 1.0e-3 (INV), sign-fixed; 1.4e-3 / 2.5e-3 as published.
* With HRCam (16th percentile) instead: 9.5e-4 / 2.2e-3.
* Files: `out/hri.json`, `figs/hri_ptp.png`.

## 4. Radial velocities (f08)

### Mass and K
For Rp = 4.90 ± 0.28 R_E and M* = 0.370 (`out/f08_rv_atmo.json`, key `mass_K`):

| Relation | Mass (M_E) | K (m/s) |
|---|---|---|
| Chen & Kipping 2017 (verified: C 1.008, S_Neptunian 0.589, 14.6% scatter, T1 2.04 M_E) | 21.3 (16.6–27.3) | 14.8 (11.5–19.0) |
| Otegi et al. 2020, volatile-rich (verified: M = 1.74 R^1.58) | 21.2 (15.7–27.8) | 14.7 (10.9–19.3) |
| Empirical: 10 NASA Exoplanet Archive M-dwarf planets of 4–7 R_E with masses, power-law fit with 0.57 (ln) scatter | 20.0 (11.1–35.8) | 13.9 (7.7–24.9) |

The 95% range of K under Chen & Kipping is 9.0–24.2 m/s.

### Precision per 2400 s on this star
Scaled photon-limited from verified published values:
* **TOI-756** (V 14.6, J 11.1), from arXiv:2510.14927:
  * HARPS: 5.5 m/s photon noise in 2400 s;
  * NIRPS: 8.4 m/s in 3 × 800 s, rising to 15.4 m/s after removing OH-affected lines;
  * fitted jitter about 15 m/s.
* **TOI-674** (V 14.2), from arXiv:2106.01246: HARPS 7.3 m/s in 30 min.

Results:

| Instrument | σ per 2400 s (m/s) |
|---|---|
| HARPS | 8.8 (12.0 if scaled from TOI-674) |
| NIRPS | 9.5 (17.4 on OH-affected dates) |
| NIRPS+HARPS | 6.4 |
| ESPRESSO | 3.4 (scaled from HARPS by area × an **assumed** 1.25 throughput) |

### Epochs for a 5σ K
Taking σ_K = σ_eff √(2/N), at K = 14.8 m/s (K = 11.5 in brackets), for jitter 5 / 10 / 15 m/s:

| Instrument | Epochs |
|---|---|
| ESPRESSO | 9 / 26 / 55 (14 / 43 / 90) |
| NIRPS+HARPS | 16 / 33 / 61 (26 / 54 / 101) |
| HARPS | 24 / 41 / 69 |
| NIRPS | 27 / 44 / 72 |

The jitter of this X-ray-active M3V is unknown; TOI-756 needed about 15 m/s.

### NIRPS scheduling
* The systemic velocity is +16.3 km/s and the BERV spans −26.3 to +26.6 km/s.
* On 40% of observable dates (January–May), |Vsys − BERV| ≤ 10 km/s, which puts the stellar lines on the OH and
  telluric lines, as happened for TOI-756.
* Schedule NIRPS in June–November (`out/berv_season.csv`).

### Instruments that cannot reach Dec −51.3°
Maximum altitude in brackets:
* MAROON-X, KPF, SPIRou, IRD (18.9°, airmass ≥ 3.1);
* HARPS-N, GIANO-B (9.9°);
* NEID (6.7°);
* CARMENES (1.5°);
* SOPHIE (below the horizon);
* HPF (HET's declination limit is −10.3°).

Usable: ESPRESSO, NIRPS/HARPS, CRIRES+ (K band, not competitive here; not quantified), PFS (V = 15.6 is too faint
for iodine RVs), Veloce, SALT-HRS, GHOST, FEROS.

### Low-precision RVs (closing the BD and WD channels)
* K for a BD on the M3V: 2.8 (13 M_J), 4.3 (20), 8.3 (40), 12.0 (60), 15.6 (80) km/s, at e = 0; larger by
  1/√(1 − e²) for eccentric orbits.
* K for a 0.6-Msun white dwarf: about 73 km/s.
* Two epochs at opposite quadratures with σ ≤ 0.5 km/s separate 2K ≥ 5.6 km/s at more than 8σ. Use 3–4 epochs to
  cover eccentric (e ≳ 0.5) orbits and to look for an SB2 line profile.
* Suitable: FEROS, GHOST, SALT-HRS, or any HARPS/NIRPS/ESPRESSO epoch.
* Gaia DR3 already disfavours K ≳ 12 km/s (../bd).

## 5. Atmospheric metrics and Gaia DR4 (f08)

**Planet temperature and metrics.**
* Teq = 472 (461–483) K (A = 0, full redistribution).
* TSM, with Kempton et al. 2018 eq. 1 and scale factor 1.15 (verified): **109 (85–139)** for the Chen & Kipping mass;
  109 (83–148) for Otegi; 116 (65–206) empirical. This is above the Kempton recommended threshold of 90 for 1.5–10 R_E
  (verified).
* ESM = 9.3 (8.2–10.5).

**Ranking.** Among the 19 transiting 4–7 R_E planets with P < 10 d and Teff < 4000 K hosts (NASA Exoplanet Archive,
refreshed 2026-10-17; `out/tsm_esm_rank.csv`; candidate included), it would rank **9th by TSM and 16th by ESM**.
* Planets above it by TSM: AU Mic b, GJ 436 b, GJ 3470 b, TOI-3884 b, TOI-674 b, TOI-3785 b, TOI-1883 b, TOI-1728 b.
* It would be the second Neptune-size planet with a host below 3400 K, after TOI-3884 b.

**Gaia DR4.**
* What is verified (ESA pages):
  * release on 2 December 2026;
  * 66 months of data, from 2014-07-25 to 2020-01-20;
  * epoch astrometry, epoch G/BP/RP photometry, per-CCD G photometry, RVS epoch spectra and per-transit RV parameters,
    XP epoch spectra, NSS orbit and acceleration solutions, and a `vari_planetary_transit` table;
  * DR5 (10.5 years) not before the end of 2030.
* What DR4 gives for this star:
  * no in-transit photometric epoch (ctx/referee_ctx: 96 crossings, none within 6.7 h);
  * epoch RVs, if published for this faint (rv_method 2) star. With DR3's 2.2 km/s combined error from 21 transits,
    single transits are probably good to ≳ 5–10 km/s, so they constrain only K ≳ 10–15 km/s (massive BDs);
  * epoch astrometry to test photocentre wobble. A BD on the M3V gives 0.02–0.12 mas (../bd), which is detectable only
    for the most massive BDs; the precision was not verified here;
  * a possible acceleration solution for the ~17-au PTP companion.
* The in-transit Gaia epoch (2024-10-29) arrives only with DR5.

## 6. Prioritised plan (`out/plan.json`, `out/fpp_ladder.csv`, `figs/fpp_ladder.png`)

**Expected FPP after each step.** Values are r05 numbers, sign-fixed REF / INV (as published REF / INV in brackets).
Two-band steps that include g or V also remove the WD term.

| # | Observation | When | FPP after |
|---|---|---|---|
| 1 | Multi-band seeing-limited transit (MuSCAT4 g′r′i′z_s; or LCOGT 1-m r′/z_s + SPECULOOS I+z). Small apertures (r ≤ 2″) or PSF photometry; Mira bleeds and spikes kept off PA 90°. | 2027-06-01 (SSO); 2027-06-07 and 2027-08-03 (Chile); 2027-06-29 (SAAO) | 1.5e-4 / 4.8e-4 (8.9e-4 / 1.9e-3) |
| 2 | Gemini-S/Zorro 832 + 562 nm speckle (SOAR/HRCam as a cheaper first look) | target season, Feb–Oct | alone 6.0e-4 / 1.0e-3; with 1: 1.4e-4 / 4.5e-4 |
| 3 | 2–4 low-precision RV epochs near quadrature (σ ≤ 0.5 km/s) | any time | 3e-5 (floor = background term carried in r05) |
| 4 | Precise RV mass: ESPRESSO 9–26 epochs, or NIRPS+HARPS 16–33 in June–November | 2027 season | ~4e-5 (bookkeeping), plus the planet's mass and density |

Other ladder entries:
* With the ../bd brown-dwarf prior (f_BD = 1.3e-4) and no follow-up: 1.3e-3 / 3.3e-3.
* With that prior plus steps 1 and 2: 4.6e-5 / 8.8e-5.

## Evidence for / against, as it bears on the design
* **For a planet on the M3V:** every false-positive class that remains has a cheap, decisive test:
  * WD: g/TESS ≈ 6;
  * background EB: g/z ≈ 2–3;
  * PTP: g/z ≈ 0.78, and 97% within Zorro's reach;
  * BD: K ≥ 2.8 km/s.
* The expected planet is an excellent target: K ≈ 15 m/s; TSM ≈ 110.
* **Against, or open:**
  * The FPP numbers assume the TESS signal is astrophysical and periodic. The photometric systematics next to the Mira
    (clean: INV odd-even 2.1σ; late/early depth contrast 1.6σ; contrast moves by 2σ with the regressor set) are not
    in the FPP. Only a ground transit tests them.
  * Ground prospects depend on assumed red noise (0.5–1 ppt per event). If the red noise is doubled, the S/N falls to
    13–15 and 5% depth ratios need about 4 MuSCAT4 events.
  * The PTP chromatic prediction is model-dependent at ±0.05, and spots bias g/z by +1–9%. A measured g/z ≈ 0.9 would
    be ambiguous.
  * The first usable transit is April 2027, and only 2–6 full windows per site per year exist.

## Scripts (`scripts/`)

| Script | What it does | Outputs |
|---|---|---|
| `flib.py` | constants; ephemeris loader (clean / phase-1 / any file); target astrometry; sites; Moffat helpers | — |
| `f01_transits.py [eph.json]` | transit times, σ, barycentric UTC per site, observability windows, Moon, comparison with ctx | `out/transits_all.csv`, `out/transits_obs_all_sites.csv`, `out/obs_<site>.csv`, `out/obs_summary.csv`, `out/f01_transits.json` |
| `f02_gaia_query.py` | Gaia DR3 cone (25″) with XP synthetic photometry; Mira LPV and vari_summary rows | `data/gaia_25arcsec_gspc.csv`, `data/gaia_mira_lpv.csv`, `data/gaia_mira_vari_summary.csv` |
| `f02_mags.py` | cross-match SkyMapper DR4, DECaPS2, 2MASS, VHS (from ../star) with proper motions; adopted magnitudes | `out/mags_table.csv`, `out/mags_adopted.json` |
| `f03_mira.py` | Mira period/phase fit (Gaia + TESS); brightness ranges per band; phase at the observable transits | `out/mira.json`, `out/mira_phase_at_transits.csv`, `figs/mira_fit.png` |
| `f04_chromatic.py` | synthetic photometry (BT-Settl, CK04, blackbody; SVO filters; LD); depth ratios for the planet, PTP, WD, background EB, the 3.1″ star and spots | `out/chromatic.json`, `out/chromatic_ratio_vs_dT.csv`, `out/m3v_band_mags.json` |
| `f05_contam.py` | 2027 geometry; Moffat (+ halo) aperture contamination grid; differential dilution | `out/geometry_2027.csv`, `out/contamination_grid.csv`, `out/contamination_summary.json` |
| `f06_precision.py` | noise model, S/N per event, saturation, depth-ratio tests | `out/precision.csv`, `out/saturation.csv`, `out/depth_ratio_tests.csv`, `out/precision_summary.json` |
| `f07_hri.py` | PTP/BD/TP-diluted draw weights (r05, as published and sign-fixed); contrast curves; excluded fractions; FPP after each step | `out/hri_exclusion.csv`, `out/fpp_after_followup.csv`, `out/hri.json`, `out/ptp_draw_weights.npz` |
| `f07b_ptp_chroma.py` | PTP-posterior distribution of band depth ratios | `out/ptp_chroma.json` |
| `f08_rv_atmo.py` | mass–radius, K, RV precision and epochs, declination reach, BERV/OH season, TSM/ESM and ranking (NASA Exoplanet Archive refresh) | `out/f08_rv_atmo.json`, `out/rv_feasibility.csv`, `out/berv_season.csv`, `out/tsm_esm_rank.csv`, `data/nea_mdwarf_neptunes_refresh.csv` |
| `f09_plan.py` | two-band PTP rejection vs number of events; FPP ladder (WD removed for g/V steps); figures; plan | `out/two_band_ptp_rejection.csv`, `out/fpp_ladder.csv`, `out/plan.json`, `figs/observability.png`, `figs/hri_ptp.png`, `figs/fpp_ladder.png` |
| `f10_sim_check.py` | Monte Carlo of per-event depth errors (white + OU red noise + quadratic baseline, ±2-min timing error) | `out/sim_depth_error.csv` |

**Run order:** f01 → f02_gaia_query → f02_mags → f03 → f04 → f05 → f06 (needs f10's factor 1.42, hard-coded) → f07 →
f07b → f08 → f09. f10 is independent.

**Other folders:**
* `data/filters/`: SVO filter curves.
* `data/claret2011_ld_3300_3400.csv`: limb-darkening coefficients.
* `lit/`: text extracts of the papers used for verified numbers: Ziegler+2026 SOAR III, Howell+2025 Gemini speckle,
  Chen & Kipping 2017, Kempton+2018, TOI-674 (Murgas+2021), TOI-756 (Parc+2025). PDFs were deleted.
* `logs/`: console output.

## Caveats
* **Assumed, not verified:**
  * telescope throughputs, sky brightness, extinction;
  * read noise, full well and readout of SPECULOOS, TRAPPIST, ASTEP+ and the amateur set-up;
  * the red-noise level;
  * C_Y at Siding Spring, Sutherland and Dome C;
  * ESPRESSO's throughput ratio to HARPS;
  * the Zorro and ERIS curve shapes and HRCam inside 0.1″;
  * the Mira's r/z/J/Ks amplitudes;
  * the Mira's cycle-to-cycle phase jitter.
* **Models:**
  * BT-Settl [M/H] = 0 models are about 0.3–0.4 mag too blue in g/V for this metal-rich M3V. Depth ratios are
    differential, but the PTP ratios carry about ±0.05 of model systematics.
  * The ERIS curve is generic. Zorro's 562-nm channel may fail to reach these contrasts at V = 15.6.
* **FPP bookkeeping:**
  * The FPP ladder inherits every prior of `../referee2_fpp` (and its background term of 3e-5, which is not
    re-derived).
  * The sign-fixed version follows `../bd`'s b10 monkey-patch.
  * r05's chroma option tests only V/I for planet-host dilution. The WD removal for g/V two-band steps is done here by
    bookkeeping.
* **Instrument verification:** some instrument facts were not checked against manuals (MuSCAT4 gain, the LCOGT 1-m
  obstruction, ASTEP+ detectors). The Sinistro cameras are due to be replaced by Sophia deep-depletion CCDs in
  2026/2027 (LCO page), which would improve the z band.
