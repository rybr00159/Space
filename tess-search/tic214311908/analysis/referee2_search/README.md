# Referee 2: audit of the 'search' analysis of TIC 214311908 (extra planets, aliases, rotation, flares)

Folder: `tic214311908/analysis/referee2_search/`
(AI-run analysis, 2026-10-17. Nothing was posted, uploaded or submitted anywhere, and no one was contacted.)

The analysis under review is in `../search/`. Its scripts were read to find errors; every number below was
re-derived with separate code in `scripts/` (own BLS, own harmonic fits, own trapezoid fits, own cube photometry).
Only data products were read: `../phot/lc_M3V_all.csv`, `../phot/out/per_transit_M3V.csv`,
`../phot/out/timing_table.csv`, `../referee_phot/out/per_transit.csv`, the TESScut cubes in
`/home/user/work/tesscut5/214311906/`, and (to audit them) the analyst's null/injection outputs and
`../search/data/m3v_clean.npz`. Gaia DR3 `vari_rrlyrae` was queried directly (`data/gaia_rrab.csv`).

## Verdict (short)

**Upheld, with independent numbers**

1. **The 0.6485-d signal is the inverted light curve of the Gaia DR3 RRab 5925372615693467904 (200" S).**
   * Coherent 3-harmonic scan (S39-S104): best P = 0.648528 d; Gaia pf = 0.648528 +- 0.000003 d is worse by only
     Delta-chi2 = 0.03. Signal Delta-chi2 = 880 against at most 49 at 119 fake periods. Neighbouring alias teeth worse
     by 106 and 133 (`r01`).
   * Phase: our own 3x3-px aperture photometry of the RRab in the cubes peaks within one 1/30 bin of the maximum
     predicted from the Gaia epoch (propagated phase error 0.013-0.032).
   * The M3V-B waveform is anti-correlated with the RRab's at zero lag: r = -0.84 to -0.92 in S39-S104 (S12 -0.30,
     where the leak is undetected). The leak is 3-13 % of the RRab's own aperture variation (`r02`, `figs/r02_rrab_fold.png`).
   * Per-sector first-harmonic semi-amplitudes are reproduced exactly: 0.73 +- 0.57, 1.84 +- 0.20, 2.46 +- 0.37,
     3.23 +- 0.24, 3.90 +- 0.38, 2.03 +- 0.29 ppt (S12-S104).
   * **New:** the raw pixels around the M3V do *not* carry it. In the raw 3x3 box (Mira core included) the
     first-harmonic amplitude is 0.7 and 0.8 e-/s in S93 and S103 (fake-period p = 0.28, 0.22), where model B shows
     3.2 and 3.8 e-/s. It appears (1.8-2.6 e-/s, inverted, p <= 0.01) only after a field-wide common-mode background
     (median of the distant pixels) is subtracted (`r11`). So it is a pipeline coupling through a field-wide term. The
     analyst's statement that it "does not come through the background" is not supported, and the route remains open.
2. **The leak's effect on per-transit depths is reproduced exactly:** rms change 2.81 ppt, maximum 5.74 ppt (`r01` E).
   * The phase-1 phot per-transit depths regress on our predicted leak bias with slope 1.00 +- 0.38, so they carry
     the leak.
   * The referee's independent depths give slope 0.47 +- 0.28, which is consistent with 0 and with 1 at < 2 sigma.
     Their S103 E55/E56/E59 pattern (19.1/9.2/20.1 e-/s) matches the predicted bias (+5.6/-4.7/+5.7 ppt), so a partial
     leak there cannot be excluded.
3. **Odd/even is not significant.**

   | Light curve and method | Odd - even (ppt) |
   |---|---|
   | Own local fits, cleaned model B | +2.67 +- 1.71 (1.6 sigma) |
   | Own local fits, uncleaned model B | +4.07 +- 1.71 |
   | Referee per-transit depths, raw | +1.57 +- 1.75 |
   | Referee per-transit depths, leak bias subtracted | -0.60 +- 1.75 |
   | Analyst | +3.3 +- 1.8 (fold) / +2.3 +- 1.8 (local) |

4. **The period is unique.**
   * P/2 is excluded: phase-0.5 depth -0.98 +- 0.78 ppt against 14.66 +- 0.78 ppt at phase 0, i.e. 14 sigma for the
     difference. The analyst's ">19 sigma" ignores the error of the primary depth.
   * Cycle count: the next-best assignment of the transit times over 5.55-5.85 d is worse by Delta-chi2 = 1.26e4
     (phot times) and 3.2e4 (referee times), and still >= 3.1e3 with doubled errors (`r04`). The analyst's phot value
     of 1.8e5 came from a +-5-tooth scan; the wider scan finds 5.6038 d at 1.26e4. The conclusion is unchanged.
5. **No second planet.**
   * Own BLS (different grid, durations and code; 1,464,757 frequencies): top peak 26.36 d at Delta-chi2 = 96.7
     (analyst 94.3). It is built from the same single events (S66 ~BTJD 3117, S104 BTJD 4197.97).
   * Real band maxima (0.3-1, 1-3, 3-10, 10-30 d): 39.8, 49.1, 62.3, 96.7. None exceeds any calibrated threshold
     (`r03`, `r06`).
6. **The 5.69-d dip is not a leak from a distant star (new model-free test, `r07`).** The test uses raw cube pixels,
   with a common mode removed and local linear baselines.
   * No pixel more than 2.5 px from the M3V dims at the transit times beyond the fake-epoch distribution in any
     sector (p = 0.11-0.90). Per-pixel errors are 0.3-0.5 e-/s, so a distant source would need eclipses of 100 e-/s
     or more to leak 15 e-/s and would be seen at S/N > 50.
   * The 3x3 box at the M3V (the Mira core is 0.85 px away and is inside it) dims by 13.9 +- 1.4 e-/s, with chi2 4.2/5
     across sectors, without any Mira model or regressor.
7. **No rotational modulation at P_orb.** Per sector, fake-period p = 0.33-0.73 and conservative 95 % upper limits are
   0.46-2.3 % (`r09`); the analyst gives 0.37-2.3 %. 2 P_orb is untestable: there are at most 2.3 cycles per sector,
   degenerate with the trend.
8. **The flare conclusion holds:** the data cannot separate active from inactive stars.
9. **X-ray to rotation:** L_X/L_bol = 4.2e-4 gives P_rot = 8.5 d (Wright+11) or 12.8 d (Wright+18), and 6.5-17 d over
   +-0.3 dex (`r08`). This is consistent with the analyst's ~11 d.

**Revised**

1. **The completeness map is too pessimistic at P < 8 d.** One global threshold (Delta-chi2 = 147) is set by
   long-period single-event systematics and applied at all periods.
   * Per-band thresholds from the analyst's own nulls (FAP 0.25 % per band, 1 % overall) are 57, 61, 92 and 164 for
     0.3-1, 1-3, 3-10 and 10-30 d.
   * With these, the analyst's own 800 injections give the R50 values below (R90 2.3-2.7 R_earth at P < 2 d):

     | P (d) | R50, band thresholds (R_earth) | R50, analyst (R_earth) |
     |---|---|---|
     | 0.5-1 | 2.17 | 2.77 |
     | 1-2 | 2.39 | 3.01 |
     | 2-4 | 3.04 | 3.57 |
     | 4-8 | 3.57 | 4.17 |
     | 8-15 | 4.45 | 4.52 |
     | 15-30 | 5.28 | 5.22 |

   * "Nothing below 2.5 R_earth is recovered" is therefore false for P < 2 d. The real data stay below every band
     threshold, so the null result stands (`r06`).
2. **The analyst's null set is inflated by a gap-split event.** 10 of 16 inverted+rolled maxima, and 7 of the 14
   killed-run maxima, equal 119.5. That is one 45-ppt positive excursion in S12 (BTJD 1634.41; the phase-1 "4-point
   flare") which straddles an 8-h gap in the real data and becomes contiguous once flux values are re-assigned to new
   time stamps.
   * It is not a wrap junction: junction segments give Delta-chi2 <= 7 (`r06b`).
   * Without it, the global 1 % threshold is 137 rather than 147 (Gumbel; 140 from the 30 non-inverted rolls alone).
     Minor.
3. **"2P (an equal-eclipse EB at 11.38 d) cannot be excluded photometrically" is overstated** (`r08`).
   * Odd-even timing offset: -3.3 +- 3.3 min (referee times) and +4.4 +- 4.8 min (phot), so |e cos w| < ~1e-3.
   * A flat-bottomed eclipse needs k <= 0.36 (flat fraction >= 0.47).
   * Every main-sequence pair with k <= 0.36 has an occultation-to-transit depth ratio <= 0.30. It therefore predicts
     odd - even >= 17 ppt at the observed mean depth, excluded at >= 7.7 sigma by every odd/even measurement,
     including the analyst's own +3.3 +- 1.8.
   * 2P survives only for near-twin (k ~ 1) eclipses that are not flat. That is a question for the shape analysis,
     which disfavours V shapes by Delta-chi2 = 17-37.
4. **The TLS 16.79-d "p = 0.42" is circular.**
   * Rolling each sector keeps both events of a two-event candidate, so the rolled copies always contain a similar
     pair.
   * Against the inverted+rolled copies, SDE 12.6 exceeds all 12 (maximum 9.96), so p < 0.08.
   * The candidate is correctly rejected, but by event vetting rather than this p-value. The missing middle event
     gives 2.7 +- 4.3 e-/s against 18.6 expected (~3 sigma), E0 is in scattered light, and E2 is the S104 non-PRF event.
   * Likewise, the per-sector "roll" nulls keep within-orbit periodic signals, so the per-sector p-values are
     conservative, not calibrated.
5. **"The TARS ~8.7-d signal is not on the M3V" is weaker than stated.**
   * Our conservative 95 % upper limits at 8.7 d are 1.6 % (16 e-/s) in S12 and 0.69 % (6.7 e-/s) in S93, against an
     implied 16 and 7.9 e-/s. S12 does not exclude it and S93 only marginally (`r09`).
   * Rotation of the M3V at ~8.7 d is not excluded. It would match the Wright+11 X-ray prediction of 8.5 d.
6. **Flare ED50 depends on the assumed decay times.** The quoted ED50 = 75 s (1.0e33 erg) averages over decay times
   t_1/2 of 2-40 min. For impulsive flares (t_1/2 <= 8 min), the analyst's own injections give ED50 ~ 32-42 s
   (~5e32 erg). Medina+20 active stars would then give ~0.4 flares in 110.8 d, so the conclusion stands.
7. **The single-transit threshold of S/N 11.3** is the median of the maxima of other, noisier light curves, not a
   calibrated false-alarm level. "0 credible single transits" rests on the vetting of the two events above the
   inverted-curve level of 7.1. That vetting is reasonable.

**New checks (missed tests)**

1. Removing the RRab leak does not change the global transit shape on the phot curve (`r10`):

   | Quantity | Before | After |
   |---|---|---|
   | T14 | 1.81 h | 1.79 h |
   | Flat fraction | 0.46 | 0.45 |
   | V-shape Delta-chi2 | 18.7 | 17.0 |
   | Depth | 16.1 | 15.8 ppt |

   The phot-versus-referee shape difference (T14 1.8 h against 1.57 h) is therefore not due to the RRab.
2. Depth by epoch on the cleaned curve: S39+S66 12.6 +- 1.3 ppt against S103+S104 17.1 +- 1.7 ppt, a 2.1 sigma
   difference (2.5 sigma before cleaning). Mild, and worth watching for chromatic or blend effects. The all-transit
   depth chi2 is 17.8/19.

## Key numbers

| Quantity | Value | Method | Uncertainty / calibration | File |
|---|---|---|---|---|
| RRab period in M3V-B curve | 0.648528 d (Gaia 0.648528 +- 0.000003; Delta-chi2 0.03) | coherent 3-harmonic scan S39-S104 | signal Delta-chi2 880 vs fake max 49; teeth +106/+133 | out/r01_rrab_leak.json (B) |
| Leak vs RRab waveform | r = -0.84...-0.92 (S39-S104); leak 3-13 % of RRab aperture amplitude | fold at Gaia epoch, own cube aperture photometry | RRab max within 1 bin (1/30) of Gaia prediction (phase error 0.013-0.032) | out/r02_rrab_phase.json |
| RRab signal in raw 3x3 at M3V | 0.7 (S93), 0.8 (S103) e-/s, p 0.28/0.22; with field common mode 1.8-2.6 e-/s (p <= 0.01) | 3-harmonic fit, raw cube | 100 fake periods | out/r11_rawleak.json |
| Per-transit depth change from leak | rms 2.81, max 5.74 ppt; chi2 29.1/19 -> 17.8/19 | trapezoid (1.6 h, flat 0.7) + linear baseline +-0.3 d | per-sector empirical errors 3.1-7.4 ppt (100 fake epochs) | out/r01_rrab_leak.json (E) |
| Depth vs predicted leak bias (slope) | phot strict 1.00 +- 0.38; referee 0.47 +- 0.28 | weighted regression | phot errors x1.3, referee calibrated | out/r01_rrab_leak.json (D) |
| Odd - even | +2.67 +- 1.71 ppt (cleaned B, own); referee +1.57 +- 1.75 raw, -0.60 +- 1.75 bias-subtracted | per-transit weighted means | empirical per-sector errors | out/r01_rrab_leak.json (E) |
| P/2 exclusion | phase 0.5 depth -0.98 +- 0.78 vs 14.66 ppt -> 14 sigma | fixed trapezoid fold | 400 random phases | out/r04_period.json |
| Next-best cycle assignment | Delta-chi2 1.26e4 (phot), 3.2e4 (referee); >= 3.1e3 with errors x2 | timing chi2 over 5.55-5.85 d | calibrated timing errors | out/r04_period.json |
| Own BLS top peak | 26.36 d, Delta-chi2 96.7, 3 events (max fraction 0.49) | bls2, 0.3-30 d, 1.46e6 freqs | see r06 for null calibration | out/r03_real.json |
| Analyst band thresholds (1 % total) | 57 / 61 / 92 / 164 (0.3-1/1-3/3-10/10-30 d) | Gumbel per band, FAP 0.25 % each | 46 analyst nulls | out/r06_thresholds.json |
| R50 with band thresholds | 2.17, 2.39, 3.04, 3.57, 4.45, 5.28 R_earth (0.5-1 ... 15-30 d) | analyst's 800 injections, logistic | binomial | out/r06_thresholds.json |
| Global threshold without gap-joined event | 137 (vs 147) | Gumbel | 46 + 7 maxima | out/r06_thresholds.json, out/r06b_gap_events.json |
| Raw-pixel dip at M3V (3x3) | 13.9 +- 1.4 e-/s, chi2 4.2/5 over 6 sectors | common-mode-removed cube, local baselines | 80 fake epochs per sector | out/r07_pixdepth.json |
| Far-pixel dimming (> 2.5 px) | max S/N 1.7-5.0, p_far 0.11-0.90 | same | fake-epoch maxima | out/r07_pixdepth.json |
| 2P EB (MS pairs, k <= 0.36) | predicted odd - even >= 17 ppt; excluded >= 7.7 sigma | blackbody 600-1000 nm, dwarf sequence | | out/r08_2p_xray.json |
| Odd-even timing | -3.3 +- 3.3 min (referee), +4.4 +- 4.8 min (phot) | linear ephemeris + parity term | calibrated errors | out/r08_2p_xray.json |
| P_rot from X-rays | 8.5 d (W11), 12.8 d (W18); 6.5-17 d over +-0.3 dex | Rossby relations, M = 0.37 | scatter, counting error | out/r08_2p_xray.json |
| Modulation at P_orb | 0.07-1.07 %, p 0.33-0.73; UL 0.46-2.3 % | sinusoid + quadratic, 1-h bins | 150 fake periods | out/r09_sync.json |
| TARS 8.7 d on M3V | S12 0.55 % (UL 1.6 %), S93 0.30 % (UL 0.69 %) | same | same | out/r09_sync.json |
| Shape before / after RRab removal | T14 1.81 / 1.79 h; flat 0.46 / 0.45; V Delta-chi2 18.7 / 17.0 | global trapezoid | formal | out/r10_shape.json |
| Flare ED50, t_1/2 <= 8 min | ~32-42 s (~5e32 erg) | analyst's injection table, re-binned | binomial | (from ../search/out/flare_injections.csv) |

### Referee-2 BLS null calibration (r03 nulls, r06)

These nulls are independent of the analyst's: 16 within-sector 1-d block shuffles, 10 inverted block shuffles, and
6 copies with whole-sector time shifts of 0-60 d (no wrap).

| Quantity | Value | Note |
|---|---|---|
| Null maxima | median 105.4 | n = 32 |
| Global 1 % threshold (Gumbel) | 149.4 | |
| Global 1 % threshold, gap-split S12 event removed | 140.6 | that event recurs as 119.44 in 6 of the 10 inverted shuffles: the same artefact as in the analyst's nulls |
| Band thresholds (FAP 0.25 % each) | 58.7 / 65.0 / 96.6 / 164.0 | 0.3-1 / 1-3 / 3-10 / 10-30 d; analyst's nulls give 57 / 61 / 92 / 164 |
| Real top peak | 96.7, p_emp = 0.72 | |
| Real band maxima | p_emp = 0.69, 0.31, 0.50, 0.72 | |

**No candidate in any band.**

Own injections (`r05`; 6 per cell) recover dchi2 within ~20 % of the analyst's at the same (P, Rp). This
independently validates their injection chain. Recovery with band thresholds and with the global threshold:

| P (d) | Rp (R_earth) | Band thresholds | Global threshold |
|---|---|---|---|
| 0.7 | 2.5 | 6/6 | 1/6 |
| 1.5 | 2.5 | 3/6 | 0/6 |
| 1.5 | 3.0 | 6/6 | 2/6 |
| 3 | 3.0 | 4/6 | 1/6 |
| 6 | 3.5 | 3/6 | 0/6 |
| 12 | 5.5 | 5/6 | 6/6 |
| 22 | 5.5 | 6/6 | 6/6 |

At P >= 10 d the band threshold is no better than the global one.

## Scripts (`scripts/`, run with /home/user/work/venv/bin/python, OMP_NUM_THREADS <= 2)

| Script | What it does | Outputs |
|---|---|---|
| `rlib2.py` | loading/rescaling, binning, harmonic design matrix, exposure-integrated trapezoid, a/R and T14 | - |
| `bls2.py` | independent numba box search (fine phase bins, cumulative sums, 6 durations 0.35-1.4 x T14c, drift-limited grid) | - |
| `r01_rrab_leak.py` | LS per sector with block-shuffle nulls; coherent period scan; per-sector harmonic fits at the Gaia period; cleaned curve; predicted per-transit leak bias; regressions of phot/referee depths; own local depths, odd/even, epoch groups | out/r01_rrab_leak.json, data/m3v_clean_r2.npz |
| `r02_rrab_phase.py` | own RRab aperture photometry from the cubes; folds at the Gaia epoch; correlation with the M3V leak | out/r02_rrab_phase.json, figs/r02_rrab_fold.png |
| `r03_bls.py` | own coherent BLS (real, inverted) and null realisations (blk, blki, shift) | out/r03_real.json, out/r03_null_*.json, data/r03_real_pgram.npz |
| `run_nulls.sh`, `run_nulls2.sh` | background runners for r03 nulls (the first was stopped after 4 blk realisations) | logs/ |
| `r04_period.py` | cycle-count scan, cycle bootstrap, alias phases | out/r04_period.json |
| `r05_inject.py` | own injection-recovery at 36 (P, Rp) cells x 6 | out/r05_inject_51.csv |
| `r06_thresholds.py` | Gumbel thresholds (global, per band), R50 from the analyst's injections with band thresholds, own nulls/injections | out/r06_thresholds.json |
| `r06b_junction.py` | tests the wrap-junction hypothesis (refuted); with `gap`, finds the gap-split S12 event | out/r06b_junction.json, out/r06b_gap_events.json |
| `r07_pixdepth.py` | model-free whole-cutout per-pixel dip S/N maps | out/r07_pixdepth.json, figs/r07_pixdepth.png |
| `r08_2p_xray.py` | odd/even timing, MS-pair 2P test, X-ray Rossby rotation | out/r08_2p_xray.json |
| `r09_sync.py` | modulation at P_orb, 2 P_orb, 8.7 d with fake-period calibration | out/r09_sync.json |
| `r10_shape.py` | global trapezoid on det vs cleaned flux | out/r10_shape.json |
| `r11_rawleak.py` | RRab signal in raw pixels near the M3V (with/without common mode) | out/r11_rawleak.json |

Logs in `logs/`.
