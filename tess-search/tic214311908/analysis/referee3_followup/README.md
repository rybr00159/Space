# Referee 3: adversarial check of the follow-up design (`../followup/`) for TIC 214311908.01

AI-run referee analysis, 2026-10-17. Nothing was posted, uploaded or submitted, and no one was contacted.
Network use was read-only: Gaia archive ADQL, IRSA TAP (NEOWISE-R), the NASA Exoplanet Archive TAP, and two web look-ups
(ESA Gaia DR4 page and news; LCO MuSCAT4 news).

All code is my own (`scripts/`). I read the analyst's scripts but did not import them. I used these data products:
* the clean ephemeris (`../clean/out/ephemeris_final.json`);
* the analyst's posterior-weighted companion draws (`../followup/out/ptp_draw_weights.npz`, built from `../referee2_fpp`);
* the analyst's BT-Settl ratio table and r05 scenario composition, used as inputs to my own bookkeeping;
* DECaPS2 from `../star/data`;
* the M-dwarf list from `../referee_star/data`;
* the Gaia epochs of the Mira from `../ctx/lc`.

Python: `/home/user/work/venv/bin/python`, threads = 1. No jobs are left running. The folder is about 0.6 MB.

## Verdict (short)

The design is sound, and its arithmetic is largely reproduced.

**Reproduced:**
* transit times, the barycentric UTC conversion and the per-site windows (identical);
* magnitudes, Mira-wing contamination, photon and scintillation budgets, K, TSM and ESM, and the Gaia DR4 date and
  content;
* the FPP-ladder bookkeeping.

**Corrections, mostly second order:**
* the PTP chromatic ratios are less chromatic empirically than in BT-Settl and r05;
* the white-dwarf g/TESS ratio is about 9, not 6;
* the Zorro exclusion of 0.97 rests on an assumed inner contrast curve;
* the g-band contamination in large apertures is about 19%, not 13%;
* the Mira phase uncertainty is underestimated.

**Missed points:**
* One multi-band ground transit would also close the dominant brown-dwarf channel (a grazing, eccentric BD on the M3V)
  by its shape alone. Low-precision RV is then not needed to reach an FPP of about 4e-5.
* The X-ray rotation estimate puts P_rot/2 at 4.3–6.4 d, which brackets the orbital period. This matters for the RV
  mass.

## Scripts and outputs

| Script | What it does | Output |
|---|---|---|
| `q01_gaia.py` | Gaia DR3 25″ cone around the M3V, XP synthetic photometry (gspc), the Mira's LPV row | `data/gaia_cone25.csv`, `data/gaia_gspc.csv`, `data/gaia_mira_lpv.csv` |
| `q02_mdwarf_colours.py` | First attempt at an empirical M-dwarf colour sample. Abandoned: the input list is a narrow BP−RP slice. Its output was deleted. | — |
| `q02b_mdwarf_seq.py` | Single M dwarfs within 20 pc (BP−RP 2.2–3.8, RUWE < 1.2) with XP g r i z V Ic and 2MASS Ks | `data/mdw_seq.csv` (312 stars), `data/mdw_seq_step1.csv` |
| `q03_neowise.py` | NEOWISE-R single-exposure W1/W2 photometry of the Mira, 2014–2024 | `data/neowise_mira.csv` |
| `r01_times.py [eph.json]` | Own Roemer, parallax-curvature and Shapiro delay → UTC per site; 2-min observability grid; Moon from exact phase angle; comparison with the analyst; ephemeris spread (INV / REF / phase 1) | `out/r01_transits.csv`, `out/r01_obs_all.csv`, `out/r01_summary.json` |
| `r01b_network.py` | Distinct epochs with a full window anywhere in the network; earliest transit-only coverage; true UTC of the 2024-10-29 Gaia crossing | `out/r01b_network.json` |
| `r03_mira_neowise.py` | Mira period from NEOWISE visit means plus Gaia G, without the TESS fluxes; phase at the 2027 events | `out/r03_mira.json` |
| `r04_contam.py` | Contamination of the M3V aperture by every Gaia source within 25″ (DECaPS2 / XP magnitudes), Moffat β 2.5 and 4, 2-D grid integration | `out/r04_contam.csv`, `out/r04_neighbours.csv` |
| `r05_hri_ladder.py` | Excluded PTP fraction under alternative Zorro and HRCam curves; two-band Bayes factors and 3σ rejection for three colour models and four precisions; FPP ladder by explicit bookkeeping | `out/r05_hri_ladder.json` |
| `r06_rv_tsm.py` | Mass–radius (Chen & Kipping 2017, Otegi+20, empirical NASA Exoplanet Archive fit), K, epochs over the K distribution, TSM, ESM and ranking (own NEA query, transiting only) | `out/r06_rv_tsm.json`, `data/nea_mdw_neptunes.csv` |
| `r07_snr.py` | Own photometric noise model (photon, sky, read noise, Osborn+15 scintillation) with red noise of 0.5, 1.0 and 1.5 ppt | `out/r07_snr.csv`, `out/r07_snr.json` |
| `r08_chroma_empirical.py` | Empirical (non-BT-Settl) PTP depth ratios from the 20-pc sequence and the dKs draws; WD g/TESS from observed colours | `out/r08_chroma.json` |
| `r09_grazing_ground.py` | Asimov Δχ² of the best grazing dark body (k ≥ 0.19, D = 1) against a planet truth, for one MuSCAT4-like r+i+z event | `out/r09_grazing_ground.json` |

## Key results (method, uncertainty, file)

### 1. Transit times and windows (`r01`, `r01b`)
* **Agreement with the analyst** over all 139 epochs (E89–E227):
  * BJD agrees to 4e-5 s;
  * geocentric UTC to ±1 s (rounding);
  * light-time to 13 ms;
  * σ is identical.
* **Windows.** The full windows (Sun < −18° and Sun < −12°) are identical at all 8 sites, with no epoch disagreeing.
  Moon illumination agrees to 0.012 and separation to 0.05°.
* **Ephemeris systematics.**
  * The INV, REF and phase-1 ephemerides differ from the adopted combined one by +0.2 to +1.9, −0.5 to −0.9 and −0.8 to
    −2.2 min over 2026–2028.
  * That is comparable to the quoted 1σ of 1.8–3.1 min. A pipeline systematic of about 1 min should be added in
    quadrature.
  * The ±3σ windows still cover it.
* **Network view.** In 2027, 10 distinct epochs have a full window from Chile, SAAO or Siding Spring; 14 do when Dome C
  is included, and 22 with Sun < −12°. In 2028 the counts are 9 and 14.
  * The analyst's "2–4 per site" is correct, but it understates what the network offers.
* **Earliest coverage.** E114 (2027-03-25 05:51 UTC, Chile) covers the whole transit but has only 0.1 h of pre-ingress
  baseline (Moon 91%). "None before April 2027" is therefore true under the analyst's criteria.
* **Gaia 2024-10-29 crossing.**
  * Relative to mid-transit: +14.9 min, upheld.
  * Its UTC label of 18:20 is BJD_TDB written on the UTC scale. The geocentric UTC is 18:25:32.

### 2. Magnitudes and contamination (`q01`, `r04`)
* **XP synthetic magnitudes reproduced exactly.**
  * M3V: g 16.363, r 14.983, i 13.575, z 12.815, V 15.591, Ic 12.859.
  * Pair A: g 16.442, r 15.312, i 13.833, z 12.883, V 15.947, Ic 13.094.
* **Total contamination** for β 2.5, at (FWHM, r_ap) of (1.5″, 2″), (2″, 3″) and (2.5″, 4″):
  * g: 1.7%, 8.7%, 18.8% (analyst totals 1.7%, 8.9%, 19.7%);
  * z: 0.21%, 1.1%, 2.3%.
  * The analyst's summary quotes the 3.1″ star alone (1.5%, 7.3%, 13%). The faint stars at 4–7″ in 2027.5 add the rest.
* **Differential g−z dilution:**
  * ≤ 1.5% for r ≤ 2″ and FWHM ≤ 1.5″;
  * 15–17% for r = 4″.
* **Mira wings in z at maximum:** 3e-5 (1″, 1.5″) to 3.1e-3 (2.5″, 4″). Upheld.

### 3. Mira period and phase (`r03`; independent of TESS)
* **NEOWISE W1/W2 visit means (22 visits)**, using a sinusoid with extra scatter so that χ²/dof = 1:
  * P = 319.5 ± 22 d;
  * 294 d is disfavoured (Δχ² = 15);
  * 313.5 d is within Δχ² = 3.3.
* **Adding the Gaia G epochs:**
  * P = 318.0 ± 2.6 d;
  * 313.5 d is within Δχ² = 3.5;
  * 294 d (VSX) is excluded (Δχ² = 82).
* **Phase at the 2027 June events:**
  * NEOWISE + Gaia gives 0.37–0.46, against the analyst's 0.54–0.63;
  * the analyst's ±0.06–0.09 is too small; ±0.15 is realistic (the infrared lag is also unmodelled);
  * both put the Mira away from maximum in mid-2027.

### 4. Ground S/N (`r07`)
* **Same red-noise assumption (0.5 ppt):** with independent assumptions I get MuSCAT4 g/r/i/z 20.9/26.4/28.3/28.3, LCOGT
  1-m i 25.0, SPECULOOS 28.0. The analyst's numbers are reproduced.
* **Red noise per event is the controlling, unverified input:**

  | Red noise | MuSCAT4 S/N | σ(r/z) per event | Events for 5% |
  |---|---|---|---|
  | 1.0 ppt | 13–15 | 9.6% | 3.7 |
  | 1.5 ppt | ~10 | 14% | 8 |

  LCOGT 1-m with bands on alternate events needs 3.4, 8.6 and 17 events for red noise of 0.5, 1.0 and 1.5 ppt.

### 5. Chromatic predictions (`r08`, `r05`)
Ratios are median (16–84%) for the sign-fixed REF draws.

| PTP ratio | Empirical, 20-pc M-dwarf sequence | Analyst BT-Settl | referee2 Mann V/I (used inside r05) |
|---|---|---|---|
| r/z | 0.807 (0.751–0.858) | 0.754 (0.675–0.831) | |
| g/z | 0.770 (0.709–0.830) | 0.785 | |
| V/Ic | 0.843 (0.779–0.881) | 0.80 | 0.752 (0.688–0.815) |

* **3σ rejection with one event at σ_R = 5%:**
  * r/z: 78% (empirical) against 93% (BT-Settl);
  * g/z: 97% (empirical) against 83% (BT-Settl).
* **At σ_R = 10%:** 5–13% (empirical) and 9–26% (BT-Settl).
* **Bayes factors at σ_R = 5%:** 7.4e-3 (empirical r/z), 2.2e-3 (BT-Settl r/z) and 6.1e-4 (r05's internal V/I).
  r05 is therefore the most optimistic of the three.
* **Recommendation:** g/z or g/i as the discriminating pair, and r/z only as a secondary check.
* **White dwarf occulted by the M3V:** g/TESS = **9.0** (7.0 at 5500 K to 11.4 at 7200 K), from observed M3V colours
  (g − T = 3.23). The analyst's 6.0 comes from BT-Settl, which is 0.43 mag too bright in g. This changes no conclusion.

### 6. High-resolution imaging (`r05`)
* **PTP weight inside a given separation:**

  | Separation | 0.02″ | 0.03″ | 0.04″ | 0.05″ | 0.06″ | 0.1″ |
  |---|---|---|---|---|---|---|
  | PTP weight inside | 3% | 7% | 12% | 18% | 27% | 50% |

* **Zorro 832 nm, excluded fraction:**
  * 0.97 with the analyst's assumed curve (dI = 1.2 reached at 0.020″);
  * 0.93, 0.88 and 0.81 if dI = 1.2 is reached only at 0.03″, 0.04″ or 0.05″ (plausible for I = 12.9 and V = 15.6).
* **HRCam:** 0.50–0.66. Upheld.

### 7. FPP ladder (`r05`; sign-fixed r05 composition, which referee3_bd reproduces)

| Step | REF | INV |
|---|---|---|
| Now | 1.44e-3 | 3.76e-3 |
| Zorro | 6.0e-4 (6.7–7.3e-4 with conservative curves) | 1.0e-3 (1.26–1.46e-3) |
| Two-band only, σ_R = 5% | 1.65–1.69e-4 | 5.5–5.6e-4 |
| Two-band only, σ_R = 10% | 2.5–3.2e-4 | 0.8–1.0e-3 |
| Zorro (conservative B) + two-band | 1.4–1.6e-4 | 4.6–5.2e-4 |
| + BD closed | 3.0–4.9e-5 | 3.1–8.6e-5 |

The end points are robust. Only the "two-band alone" step depends on the precision and the colour model.

### 8. Missed: the ground light-curve shape closes the grazing-BD channel (`r09`)
* **Method.** Asimov fit of the best grazing dark body (k ≥ 0.19, undiluted) to a planet truth, for one MuSCAT4 event
  (r + i + z at 1.08 / 0.74 / 0.75 ppt per 2 min).
* **Result.**
  * Δχ² = 970 with white noise, 108–242 for β = 2–3.
  * TESS gives Λ_G = 13.6 (REF) and 6.8 (INV).
  * The best grazing solution is k 0.19, b 0.97, a/R* 14.4 (ρ_c ≈ 1.7 g/cm³), as in `../bd`.
* **Consequence.** The BD-on-M3V term (≈ 1.0e-4 REF / 3.7e-4 INV, which dominates after Zorro plus two-band) is closed
  by the same ground event. The FPP then reaches about 4e-5 without RV.
* **Caveats.**
  * Not injection-calibrated.
  * Limb darkening fixed; Asimov data.
  * The β deflation is approximate.
  * BD on a companion (diluted, k ≳ 0.25) is not closed by shape. It is handled by Zorro and the colour test.

### 9. RV, TSM and ESM (`r06`)
* **Mass and K.**
  * Chen & Kipping 2017 (inverted, 14.6% radius scatter): M = 21.3 (16.6–28.3) M_E, K = 14.8 (11.5–19.7) m/s.
  * Otegi: K = 15.3 (13.2–17.7).
  * Empirical (9 transiting M-dwarf planets of 4–7 R_E): K = 14.6 (8.3–25.6).
  * The analyst's values are reproduced.
* **ESPRESSO epochs for a 5σ K.**
  * Median-K values match the analyst: 9 / 26 / 54 at jitter 5 / 10 / 15 m/s.
  * Over the K distribution, the 84th percentile is 14 / 42 / 90 (Chen & Kipping) and 27 / 81 / 170 (empirical).
  * The headline "9–26 epochs" is the median-K case only.
* **Missed:** P_rot from the X-ray activity is 8.5–12.8 d, so P_rot/2 = 4.3–6.4 d brackets P_orb = 5.69 d. Activity can
  mimic or bias K at the planet period. Plan GP or activity-indicator modelling, and use the NIRPS-versus-HARPS
  chromaticity of K.
* **TSM / ESM.** TSM = 108 (81–139) for Chen & Kipping, P(TSM > 90) ≈ 0.7; ESM = 9.1 (7.8–10.7); Teq = 470 (454–486) K.
  * Rank 9 of 19 by TSM and 15 by ESM among transiting 4–7 R_E planets with P < 10 d and Teff < 4000 K. Upheld.
  * Only TOI-3884 b has a host below 3400 K. Upheld.
  * Note: pscomppars contains non-transiting planets with computed radii (GJ 581 b); filter on `tran_flag`.

### 10. Verified on the web
* **Gaia DR4:** 2 December 2026; data 2014-07-25 to 2020-01-20; epoch astrometry, photometry, RVS epoch spectra and
  per-transit RVs, and `vari_planetary_transit`. DR5 not before the end of 2030.
* **MuSCAT4:** on the 2-m FTS at Siding Spring since October 2023.
