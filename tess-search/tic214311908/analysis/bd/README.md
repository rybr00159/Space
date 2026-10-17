# Brown-dwarf and companion-hosted channels for the 5.69-d signal on TIC 214311908 (folder bd/)

Folder: `tic214311908/analysis/bd/`

AI-run analysis, 2026-10-17. Nothing was posted, uploaded or submitted, and no one was contacted. Network use was
read-only:
* arXiv API and arXiv e-print sources;
* Crossref metadata (to turn DOIs into bibcodes);
* the NASA Exoplanet Archive TAP service and the Gaia DR3 archive;
* the ATMO 2020 model download;
* one page of the Carmichael brown-dwarf catalogue.

**STATUS: INCOMPLETE.** The workflow asked for the final answer before the full scenario-posterior run finished, so
those jobs were stopped. The table below says what is finished and what is not.

| Part | State |
|---|---|
| Literature priors (b01) | done |
| BD model priors (b03) | done |
| Data sets (b04) | done |
| Real-data shape fits and profiles (b05, b05b) | done |
| Eccentric-secondary scan (b07) | done |
| Gaia RV and RUWE likelihood ratios (b08, b08b) | done |
| Correction of the photo-eccentric sign error in ../referee2_fpp (b10) | done |
| Planet-truth injections (b06) | done: 40 REF, 46 INV |
| Grazing-BD-truth and companion-truth injections (b06) | **not run** |
| Likelihood grid (b09a) | **stopped at row ~36 of 160**; no output file |
| Full posterior (b09), follow-up design (b11), figure (b12), calibration (b06b) | **not run** (scripts written and dry-run tested on a synthetic grid) |

The posterior numbers quoted below come from the corrected referee2_fpp framework (b10), with the new priors applied
by explicit rescaling. They are labelled as such.

## Key results

### 1. Shared sign error in both FPP codes (b10)
`../fpp/scripts/f04_bayes.py` and `../referee2_fpp/scripts/r05_fpp.py` both evaluate the light curve at
rho_circ = rho_host g^-3, with g = (1 + e sin w)/sqrt(1 - e^2).
* The correct relation is rho_circ = rho_host g^+3 (Dawson & Johnson 2012, 2012ApJ...756..122D).
* I checked it with batman. With e = 0.5 and w = 90 deg, T14 is 1.108 h against 1.918 h for a circular orbit, so
  rho_c/rho = 5.19 = g^3. With w = 270 deg, rho_c/rho = 0.192 = g^-3 (`out/b10_r05_eccsign.json`).
* The error attaches the periastron boost in transit probability, (1 + e sin w)/(1 - e^2), to long durations. Long
  durations are exactly what a grazing eccentric BD needs.

Re-running r05 read-only with only this sign fixed (`out/b10_r05_eccsign.json`):

| Light curve (beta noise, Kipping e) | Quantity | As published | Sign fixed |
|---|---|---|---|
| REF | FPP | 0.0031 | **0.0014** |
| REF | PTP | 1.6e-3 | 8.4e-4 |
| REF | BD on the M3V | 7.7e-4 | 9.9e-5 |
| REF | BD total | 9.5e-4 | 1.3e-4 |
| INV | FPP | 0.0074 | **0.0034** |
| INV | PTP | 4.2e-3 | 2.5e-3 |
| INV | BD total | 2.6e-3 | 4.3e-4 |

With circular orbits nothing changes. With white noise every BD term is below 1e-6.

### 2. Light-curve shape (b05, b05b, b06): own fits, LD prior, dT0 free, per-window baselines

**Data.**
* REF: 21 windows, N = 4663.
* INV-clean: phot model-B flux re-detrended here, minus the referee2_search RRab model; 22 windows, N = 3688. It
  matches the delivered cleaned flux to 0.09 ppt rms.
* ../clean/out held no cleaned light curve when this ran.

**Red-noise factors (beta) from fake epochs.** REF 1.11-1.37 by sector; INV 1.00-1.03.

**Planet fit.**

| | k | b | rho_c (g/cm3) | T14 (h) | Flat fraction |
|---|---|---|---|---|---|
| REF | 0.117 | 0.24 | 14.9 | 1.61 | 0.78 |
| INV | 0.123 | 0.73 | 4.8 | 1.84 | 0.57 |

**Grazing BD on the M3V** (D = 1, k 0.19-0.33):
* The best solution sits at the k = 0.19 bound, with b 0.97 and **rho_c = 1.8 g/cm3**.
* Lambda_G = chi2(grazing BD) - chi2(planet) = **13.6 (REF) / 6.8 (INV)** with beta noise; 21.6 / 7.1 with white noise.
* Profile in rho_c: the minimum is at 1.7 g/cm3. At 0.8 g/cm3 the penalty is +31 (REF) / +21 (INV); at 3.5 it is
  +30 (REF) / +22 (INV).
* So a grazing BD needs rho_c/rho* ~ 0.18, which means e >~ 0.55 with the transit near apoastron.

**Companion host with a large occulter** (k >= 0.25, D free), profile in k (REF / INV):

| k | 0.25 | 0.30 | 0.35 | 0.40 | 0.45 | 0.50 | >= 0.6 |
|---|---|---|---|---|---|---|---|
| Penalty | 1.3 / 0.0 | 2.4 / 0.1 | 3.7 / 0.5 | 5.4 / 1.2 | 7.4 / 2.4 | 9.8 / 3.9 | falls back to grazing: 14.4 / 7.4 |

At k = 0.40 the best dilution is D ~ 0.09. That is self-consistent with a ~0.9 R_J body on a ~0.2 Msun companion.

**Planet-truth injections at fake epochs** (`out/b06_inject_*_planet.csv`):
* Lambda_G = 11.6 +- 5.7 (REF, n = 40) and 8.0 +- 5.1 (INV, n = 46), with minima of 1.1 and 0.06.
* So the observed values are typical of a planet.
* No grazing-BD-truth injections were run, so the calibrated likelihood ratio is not available. referee2_shape's
  V-shape-truth calibration (p ~ 0.04 INV, < 0.04 REF) remains the calibrated reference.

### 3. Gaia (b08, b08b)
**RV error.** The target has rv_method_used = 2, the combined CCF of 21 transits (12 deblended). For such stars
rv_amplitude_robust, rv_chisq_pvalue and rv_renormalised_gof are not given (NaN).
* Comparison: 3778 DR3 M dwarfs with rv_method_used = 2 and rv_expected_sig_to_noise within 15% of the target's.
  Their error percentiles are 1.48 (1%), 1.81 (16%), 2.26 (50%) and 3.41 (84%) km/s. The target's 2.217 km/s sits at
  the 47th percentile.
* CCF smearing model: CCF sigma 12/16/20 km/s, peak 0.3/0.6. The likelihood ratio LR(K) for e = 0 is:

  | K (km/s) | <= 6 | 8 | 10 | 12 | 14 | 16 |
  |---|---|---|---|---|---|---|
  | LR | ~1 | 0.66-1.06 | 0.24-1.03 | 0.05-0.92 | 0.01-0.71 | 0.003-0.46 |

  It is weaker at e = 0.5.
* A BD on the M3V has K1 2.9-15.5 km/s (13-78.5 M_J, e = 0). Massive BDs are therefore disfavoured; BDs below ~40 M_J
  are not constrained.
* Companion-hosted BD (companion carries a fraction f_c of the RVS light): LR >= 0.79 for f_c <= 0.2, so no constraint.

**RUWE.** For a BD on the M3V the photocentre orbit is alpha = 0.023-0.117 mas.
* LR = 0.98 (13 M_J) to 0.57 (78.5 M_J), with sigma_CCD 0.22 mas and the clean single-star comparison (median RUWE
  1.14; the target is at the 17th percentile).
* Weak.

### 4. Literature priors (b01, b03)
**Transiting BDs around M dwarfs.** `out/b01_mdwarf_bds.csv` has 15 systems, each with its source.
* Only two are at P = 3-10 d: TOI-2119 b (64 M_J, P 7.20 d, e 0.337) and TOI-5610 b (40 M_J, P 7.95 d, e 0.354).
* None of the M-dwarf BDs at P < 16 d has e > 0.36.
* Two are grazing (NLTT 41135 b, TOI-1278 B).

**Occurrence f_BD at P = 3-10 d.** 1.3e-4 per star (68%: 5e-5 to 2.5e-4; 95% upper 3.5e-4). Ingredients:
* Bryant+23 (2023MNRAS.521.3663B) TESS giant-size rate around M dwarfs, 0.194 +- 0.072%;
* BD share of the confirmed M-dwarf transiting companions at P < 10 d, weighted by 1/p_tr: 0.15 (9 BDs against 45
  giant planets from the NASA Exoplanet Archive);
* share of that BD rate at 3-10 d: 0.46.

This is 8x below the 1e-3 that referee2_fpp assumed. The FGK cross-check is about 0.2% at P < 10 d (Csizmadia+15,
2015A&A...584A..13C).

**Eccentricity.** Ferreira+26 (arXiv:2607.02682) give short-period BDs Beta(0.426, 2.354), mean 0.15, and long-period
BDs Beta(1.88, 2.47). P(e > 0.5):

| Prior | P(e > 0.5) |
|---|---|
| Ferreira short-period | 0.09 |
| Ferreira long-period | 0.41 |
| Kipping planets | 0.10 |
| Hierarchical fit to the M-dwarf BDs at 2.5-16 d | 0.08 |

The transit-probability-weighted probability of rho_c/rho* <= 0.2 is 0.7-0.9% for every prior except the long-period
one (2.4%).

**Tides** (Jackson+08, 2008ApJ...678.1396J). tau_e at a = 0.045 au around the 0.37 Msun M3V is set by the stellar
tide:

| Q'_* | tau_e (65 M_J) |
|---|---|
| 1e5 | 1.8 Gyr |
| 1e6 | 18 Gyr |
| 1e7 | 170 Gyr |

Brown-dwarf tides need Q'_BD <~ 1e5 to matter. With the Ferreira+26 Q values (Q_* = 10^6.0) e ~ 0.5 is not damped in a
few Gyr, so tides do not exclude the eccentric channel.

**BD radius** (ATMO 2020 + 9 +- 10% inflation): R_BD 0.80-0.99 R_J (68%), k = 0.19-0.30 (95%). TESS-band flux ratio
BD/M3V: median 5e-7, > 1 ppt in only 3.6% of draws. So the occultation constrains only young or massive BDs.

**Occultation.** The scan over phases 0.2-0.8 and durations 0.4-2.4 h (b07) finds nothing beyond 3.0 sigma. The 3-sigma
limits are 1.9-4.2 ppt.

### 5. Updated channel estimates
These are rescaled from the sign-fixed r05; the full b09 run was not completed.
* f_BD 1e-3 -> 1.3e-4 (x0.13). This gives a BD total of ~1.7e-5 (REF) / ~5.5e-5 (INV).
* The Gaia RV factor for the mass-weighted BD prior is about 0.3-0.6. That puts BD on the M3V at ~1e-5 (REF) to ~3e-5
  (INV), BD on companions at ~3e-6 / 1e-5, and all BD channels <~ 1e-4 even with the pessimistic long-period e prior.
* PTP becomes the only significant residual: 8.4e-4 (REF) / 2.5e-3 (INV).
* FPP (P(not a planet on the M3V or its twin)): ~0.0013 (REF) / ~0.0030 (INV).
* P(planet on some star of the bound system) >~ 0.999 (REF) / >~ 0.997 (INV).

## Scripts (scripts/, /home/user/work/venv/bin/python)

| Script | Purpose | Outputs |
|---|---|---|
| `blib.py` | constants, helpers, photo-eccentric factor g (correct sign) | - |
| `lclib.py` | fast exposure-integrated transit templates (validated to < 16 ppm against batman), dataset class with profiled per-window baselines | - |
| `fitlib.py` | constrained ML fits: planet / grazing BD / companion-hosted body | - |
| `b01_priors.py` | M-dwarf BD table, f_BD estimate, e distributions, tidal timescales | `out/b01_mdwarf_bds.csv`, `out/b01_priors.json`, `data/ecc_draws.npz` |
| `b03_bd_models.py` | ATMO 2020 radius / Teff priors, inflation, TESS flux ratio | `data/bd_draws.npz`, `out/b03_bd_models.json` |
| `b04_data.py` | REF real and fake windows; INV-clean real and 48 fake ephemerides | `data/{ref,inv}_{real,fake}.npz`, `out/b04_data.json` |
| `b05_fit_real.py` | beta from fake epochs, real fits, profiles in k and rho_c | `out/b05_beta.json`, `out/b05_fits.json` |
| `b05b_comp_profile.py` | companion-host profile in k; k = 0.4 truth | `out/b05b_comp_profile.json`, `out/b05_comp_k40.json` |
| `b06_inject.py` | fake-epoch injection calibration of Lambda_G and Lambda_C | `out/b06_inject_{ref,inv}_planet.csv` (the other truths were not run) |
| `b06b_calibrate.py` | calibrated likelihood ratios | not run |
| `b07_secondary_scan.py` | occultation scan versus phase and duration | `out/b07_secondary_scan.{npz,json}` |
| `b08a_gaia_query.py`, `b08_gaia.py`, `b08b_ruwe.py` | Gaia RV comparison sample, CCF-smearing likelihood ratio, RUWE likelihood ratio | `data/gaia/rv_comparison.csv`, `out/b08_gaia.json`, `out/b08_rv_lr.npz` |
| `b09a_grid.py`, `b09_posterior.py`, `b11_followup.py`, `b12_figs.py` | likelihood grid, full posterior, RV / imaging design, figure | not completed (stopped) |
| `b10_r05_eccsign.py` | batman check of the sign; r05 re-evaluated with the sign fixed | `out/b10_r05_eccsign.json` |
| `run_*.sh` | job lanes | logs in `logs/` |

Other inputs:
* `data/models/` holds the ATMO 2020 CEQ tracks (2020A&A...637A..38P).
* `data/nea/` holds the NASA Exoplanet Archive query.
* `lit/` holds the arXiv abstracts and the Ferreira+26 and Larsen+25 TeX sources used for the priors.
