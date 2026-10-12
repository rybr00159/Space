# TIC 87055315: what eclipses the 9.3" neighbour of TIC 87055304?

Investigator area: `/home/user/work/tic87/neighbour/`. All scripts are in `scripts/` and run with
`/home/user/work/venv/bin/python` (TRICERATOPS with `/home/user/work/venv_tri/bin/python`), niced, OMP_NUM_THREADS=1.
All numbers are collected in `results.json` (built by `scripts/s12_collect.py` from `results/*.json`, `mcmc/*.json`
and `tri/*.json`). Figures are in `fig/`.

Inputs taken as established: three independent localisations put the P = 4.136 d signal on
**TIC 87055315 = Gaia DR3 6720986208895208960** (G = 14.02, T = 13.54), 9.33" at PA 45 deg from TIC 87055304. See
`../transit/`, `../loc_independent/` (target excluded at 8.6-10 sigma) and `../psf_photometry/` (f_target = 0.16 [0.07, 0.25]).

## Bottom line

**The transit-like signal is almost certainly the primary eclipse of an F/G subgiant + M-dwarf eclipsing binary
(an "EBLM"-type system). It is not a planet.**

| quantity | value |
|---|---|
| host TIC 87055315 | Teff 5850 (+110/-160) K, R* = 1.93 +- 0.09 R_sun, M* = 1.17 (+0.17/-0.10) M_sun, rho* = 0.172 (+0.035/-0.027) rho_sun, d = 1294 +- 53 pc, age ~5.6 (+1.5/-2.1) Gyr, turn-off / early subgiant |
| eclipse on the neighbour (undiluted) | 21.4 +- 0.6 ppt; k = 0.1394 +- 0.0031; b = 0.49 (+0.10/-0.19); T14 = 4.35 +- 0.13 h; flat-bottomed (T23 = 2.98 h) |
| companion radius | **0.269 +- 0.013 R_sun = 2.61 +- 0.13 R_Jup** (2-sigma range 2.37-2.89 R_Jup); 2.40 +- 0.12 R_Jup if only 84 % of the decrement is on the neighbour |
| twin EB at 2P = 8.27 d | rejected by the eclipse shape: chi2 worse by 49.5 than the transit model for the same number of free parameters, i.e. ~4.7-5.8 sigma after the red-noise correction (beta = 1.2-1.5) |
| occultation (pair-based) | 0.79 +- 0.77 ppt; 2-sigma upper limit 2.33 ppt, so T_b < 3510 K and L2 < 0.010 L_sun. An M dwarf of this size predicts 1.7 ppt (1.2 sigma above the measured value); a brown dwarf or planet predicts ~0.004 ppt |
| ellipsoidal | 0.87 +- 0.66 ppt (neighbour scale). A 0.26 M_sun companion predicts 0.69 ppt, so the data are consistent with it but do not measure the mass |
| odd/even | 21.65 +- 1.17 vs 21.31 +- 1.05 ppt (0.2 sigma, pair LC); psf_photometry neighbour-only: 1.3 sigma |
| TRICERATOPS (TIC 87055315 as target) | FPP = 1.000 in all 4 runs; NFPP 0.75-0.99 (priors only, no localisation). Conditional on the source being TIC 87055315: EB 0.65, PEB 0.22, SEB 0.09, all planet scenarios 0 (section 4) |
| novelty | no TOI or CTOI, and no ExoFOP follow-up. Not in Kostov+25 tables 2-4, Prsa+22, VSX, ASAS-SN or Gaia DR3 vari_* (all checked within 3" and 21"). Not in Green+23. No SIMBAD object, and no web hits for any identifier |

- **Why not a planet or a brown dwarf.**
  - The companion radius, 2.61 R_Jup, is about 25 % larger than the largest known transiting giant planets
    (~2.0-2.1 R_Jup). Those are all hotter, inflated hot Jupiters; here T_eq is only ~1500 K.
  - Brown dwarfs at the host's age of several Gyr have R ~ 0.8-1.0 R_Jup.
  - The radius is robust:
    - R* comes from the parallax plus the SED (1.93 +- 0.09 R_sun; Gaia FLAME gives 2.03, GSP-Phot 1.84, TIC 1.91);
    - k is the same within 3 % in every extraction and detrending (0.137-0.145);
    - any extra, unrecognised dilution would only make the companion larger.
- **Why an M dwarf.** A 0.27 R_sun main-sequence star has M ~ 0.26 M_sun (Boyajian+12 relation; 0.24 if 5 %
  inflated) and T ~ 3300 K. This agrees with everything measured:
  - the occultation limit (T_b < 3510 K at 2 sigma);
  - the ellipsoidal amplitude (0.87 +- 0.66 vs 0.69 ppt predicted);
  - odd = even, and a flat-bottomed, non-grazing shape.
- **Confidence.**
  - **High** that the eclipser is stellar-sized, and therefore not a planet.
  - **Moderate** that it is specifically a ~0.25 M_sun M dwarf orbiting TIC 87055315 itself. The mass is inferred from
    the radius, not measured.
  - An unresolved companion of the neighbour (within ~1") hosting the eclipse cannot be excluded with TESS data alone.
    Shape arguments, however, require that host to supply at least ~60 % of the neighbour's light (section 3e), and the
    eclipser would still be stellar-sized.
- **What would settle it.**
  - Seeing-limited photometry at one of the predicted times (section 6). The neighbour should drop by ~2.1 % and TIC
    87055304 by < 0.5 ppt.
  - Two or three RVs of TIC 87055315 (V = 14.15) at quadrature. An M-dwarf companion gives K ~ 27 km/s; a planet would
    give < 0.2 km/s.

**For the original question (a planet around TIC 87055304) this supports a false positive.** For TIC 87055315 itself
the verdict is "stellar eclipsing binary (EBLM-type)": supports false positive (as a planet candidate).

## 1. Stellar parameters of TIC 87055315 (`scripts/s01_catalogs.py`, `s04_stellar.py`, `s07_gaia_rv.py` -> `results/stellar.json`, `fig/sed_hrd_nb.png`)

**Gaia DR3 gaia_source.**
- Astrometry: parallax 0.7354 +- 0.0227 mas; Lindegren+21 zero point -0.0432 mas (`zero_point` package), giving a
  corrected parallax of 0.7785 +- 0.031 mas (error x1.3 plus a 0.010 mas floor). The inverse parallax is 1284 +- 51 pc;
  BJ21 r_geo = 1290 (+33/-36) pc.
- Astrometric quality:
  - RUWE 1.019, astrometric_excess_noise 0;
  - ipd_frac_multi_peak 1 %, ipd_gof_harmonic_amplitude 0.031;
  - **non_single_star = 0**, no NSS solution of any kind.
- Variability and EB tables: vari_summary, vari_eclipsing_binary and vari_classifier_result have no rows;
  phot_variable_flag = NOT_AVAILABLE.
- **RV.** -45.97 +- 7.34 km/s from 7 transits, 3 of them deblended, with the faint-star method (rv_method_used = 2);
  template 5250 K / logg 4.0.
  - Test (`s07`): among 20,606 Gaia stars within 12 deg with the same method, G_RVS +- 0.2 mag, template Teff
    5000-6000 K, 5-9 transits and similar expected S/N, the median error is 6.5 km/s. Ours lies at the **64th
    percentile**, so it carries no information on RV variability. A K ~ 27 km/s orbit is neither indicated nor excluded.
  - The RV and the distance (1.29 kpc versus ~0.36 kpc for TIC 87055304) show that the two stars are an unrelated
    chance alignment.

**Catalogue estimates** (all cached in `data/cat/`):

| source | Teff (K) | logg | [M/H] | R (R_sun) | M (M_sun) | other |
|---|---|---|---|---|---|---|
| TIC v8.2 | 5821 +- 122 | 3.90 | - | 1.906 | 1.05 | rho 0.152, d 1270 pc (DR2), E(B-V) 0.103 |
| Gaia GSP-Phot (MARCS) | 5573 (5541-5643) | 4.08 | -0.30 | 1.84 (1.78-1.90) | - | distance_gspphot 1170 pc, A0 0.21 |
| Gaia FLAME | - | - | - | 2.03 (1.97-2.12) | 1.154 (1.108-1.201) | L 3.49, age 6.5 (5.7-7.3) Gyr, evolstage 468 |
| Gaia MSC (binary model) | 5701 + 4854 | 4.76 / 4.56 | -0.29 | - | - | distance_msc 738 pc, which conflicts with the parallax |
| StarHorse 2019 (I/349) | 7131 | 4.09 | -0.30 | - | 1.37 | DR2-based; A_V = 1.14 exceeds the total column, so not used |
| StarHorse 2021 (I/354) | 6299 | 3.97 | -0.16 | - | 1.23 | d 1300 pc, A_V 0.69 |
| SHBoost (V/160) | 5796 | 3.89 | +0.07 | - | 1.21 +- 0.09 | A_V 0.30 |
| Zhang+23 (XP) | 5826 +- 58 | 3.87 | +0.05 | - | - | |
| Huson+25 | 5715 | 4.18 | -0.13 | - | - | |
| Ye+25 | 5857 | 3.62 | +0.12 | - | - | |
| Huang+22 (SMSS) | 5872 +- 67 | - | +0.22 | - | - | dwarf / turn-off, age 5.3 +- 0.9 Gyr |
| **this work** | **5853 (+113/-164)** | **3.95 +- 0.06** | **0.00 +- 0.18** | **1.927 (+0.088/-0.082)** | **1.17 (+0.17/-0.10)** | rho **0.172 (+0.035/-0.027)**, L 3.88 +- 0.40, age 5.6 (+1.5/-2.1) Gyr |

The catalogue median is Teff = 5821 K (sd 225 K) and [M/H] = -0.06 (sd 0.21).

**Method of "this work".**
- **SED fit.**
  - Model: Castelli-Kurucz 2004 atmospheres with Fitzpatrick-99 extinction (local copy of the stellar/ SED library).
  - Data: 12 bands (Gaia G/BP/RP, 2MASS JHKs, SkyMapper g/r/i/z PSF, AllWISE W1/W2), with error floors of
    0.01-0.05 mag.
  - Free parameters: Teff, logg, [M/H], A_V, d, R and a jitter term.
  - Priors: the parallax; [M/H] ~ N(-0.05, 0.20); A_V soft-capped at the total column, 0.36.
  - Sampling: emcee.
  - Result: chi2 = 15.3 for 12 bands without jitter. A_V = 0.34 (+0.07/-0.12) and d = 1294 +- 53 pc.
  - The SED-synthetic T = 13.544 +- 0.010, which agrees with the Stassun relation (13.540) and TIC (13.547).
- **Mass, age and density.**
  - These come from a grid-based MIST v1.2 posterior (read-only `stellar/cache/mist_compact.npz`).
  - Each isochrone point is weighted by (log Teff, log R) from the SED posterior plus [Fe/H], times a Kroupa IMF and a
    uniform-age prior.
  - The posterior is 54 % post-main-sequence (subgiant) and 46 % main sequence, i.e. a turn-off star.
  - Caveat: the MIST [Fe/H] grid is coarse (0.25 dex), and n_eff = 155.

**TESS-band flux ratio** (`s03_fluxratio.py`, `results/fluxratio.json`):
- r = F_nb/F_tg = 0.1897 (Stassun+19 relation applied to DR3 photometry), 0.1868 (TIC v8.2) and 0.1871 (Gaia XP
  sampled spectra integrated over the TESS response).
- Adopted **r = 0.1878 +- 0.0016**, so the dilution of the neighbour in the PSF pair light curve is
  **D = F_t/F_n = 5.32 +- 0.05**. The neighbour supplies 15.8 % of the pair light.

## 2. Depth on the neighbour and transit fit

### 2a. Light curves (`scripts/s02_lc.py` -> `data/lc_*.npy`)
- **pair.** F_target + F_neighbour from the psf_photometry nominal PRF fits (`../psf_photometry/lc/nominal/psf_s*.npz`).
  - Cuts are those of psf_photometry: background < median + 6 MAD, chi2 < 5x median, and S104 BTJD 4178.8-4179.2
    removed.
  - **Why the pair.** The PSF flux split between the two stars, which are 0.45 px apart, is degenerate with
    registration errors of a few 0.01 px. The neighbour-only PSF flux is therefore 6-15x noisier (flux-transfer noise)
    and has a biased normalisation, whereas the pair sum does not depend on the split.
  - **Assumption.** All of the decrement is on the neighbour. The psf_photometry f_target = 0.16 +- 0.09 is treated
    as a systematic (see 2d).
- **qlp.** QLP SAP_FLUX of the target file, with the same cuts as transit/prep.py. Used for cross-checks.

### 2b. Dilution implied by the QLP aperture (`s03`, `fig/fluxratio_dilution.png`)
- **Prediction.** The TESS PRF with Gaia positions predicts that the 2.5-px QLP aperture contains 89-91 % of the
  target's and 84-90 % of the neighbour's PSF.
- **QLP normalisation.** QLP normalises the target light curve to the target's catalogue flux. For a neighbour host the
  effective dilution is therefore D_QLP = (f_t/f_n)/r - 1 = 4.35-4.60 (S13-S93) and 3.94 for S104 (3x3 box).
- **Empirical check.** Comparing the same estimator on QLP SAP and on the PSF pair (same epochs) gives
  D_QLP = **4.39 +- 0.15** (S13-S93, chi2 = 2.4 for 3 dof) and 6.0 +- 1.1 (S104).
- **Consequences.**
  - The first round's choice of D = 4.35 was right.
  - The CONTEXT number "3,580 ppm" is a dilution-scaled depth measured with the wrong (b = 0.85) shape. With the
    correct shape, QLP SAP x (1 + 4.35) gives 3951 x 5.35 = 21.1 ppt, the same as the pair.

### 2c. Transit MCMC on the neighbour (batman + emcee)
- **Nominal fit** (`scripts/nb_mcmc_local.py` -> `mcmc/mcmc_pair_local.json`):
  - data: the pair light curve, 25 transits (S13 6, S39 5, S66 4, S93 4, S104 6), +-0.45 d windows, raw normalised
    flux;
  - baselines: a per-transit quadratic, profiled;
  - model: batman quadratic LD with a Claret-2017 prior, exposure-integrated, circular orbit, free rho, b, k, and a
    free D ~ N(5.32, 0.06);
  - sampling: 48 walkers x 10,000 steps, n/tau = 117.

| parameter | nominal (pair, local baselines) | pair, 1-d biweight (`mcmc_pair_Dfree`) | first round, QLP SAP, D = 4.35, S13-S93 (`../transit/mcmc_nbhost.json`) |
|---|---|---|---|
| k = R2/R1 | **0.1394 +- 0.0031** | 0.1374 +- 0.0027 | 0.1414 (+0.0027/-0.0030) |
| b | **0.49 (+0.10/-0.19)** | 0.42 (+0.13/-0.19) | 0.64 (+0.06/-0.09) |
| a/R1 | **7.50 (+0.67/-0.57)** | 7.94 (+0.59/-0.64) | 6.25 (+0.54/-0.42) |
| rho_transit (rho_sun) | **0.33 (+0.10/-0.07)** | 0.39 +- 0.09 | 0.19 (+0.05/-0.04) |
| T14 (h) | **4.35 +- 0.13** | 4.22 (+0.13/-0.11) | 4.81 +- 0.16 |
| T23 (h) | 2.98 (+0.15/-0.17) | 3.00 (+0.10/-0.15) | 2.89 +- 0.21 |
| depth in pair LC (ppt) | 3.39 +- 0.09 | 3.36 +- 0.09 | 3.95 +- 0.11 (QLP SAP) |
| **undiluted depth on neighbour (ppt)** | **21.4 +- 0.6** | 21.2 +- 0.6 | 21.1 (= 3.95 x 5.35) |
| P (d) | 4.136192 +- 0.000005 | 4.136201 +- 0.00001 | 4.13620 +- 0.00001 |

**Duration and density are extraction-dependent at the 1-2 sigma level** (`scripts/diag_t14.py`, `diag_t14b.py` ->
`results/diag_t14*.json`).
- Maximum-likelihood fits give:
  - pair: T14 = 4.26 h with the biweight and 4.44 h with local baselines;
  - QLP SAP: 4.92 h and 4.65 h respectively;
  - QLP S13-S93 with local baselines: 4.93 h.
- Per-sector T14 ranges over 3.8-5.5 h. QLP S39 and S66 are the longest, and these are also the sectors where QLP is
  noisiest. The biweight-detrended pair shows a +0.8 ppt "shoulder" ~4.5 h before mid-transit (`fig/eb2p_fold.png`).
- Adopted: the nominal (local-baseline) values, with systematic ranges of **rho_transit 0.19-0.39**, **T14 4.2-4.9 h**
  and **b 0.42-0.68**.
- k is stable at 0.137-0.145 across every variant.

### 2d. Depth from psf_photometry's neighbour-only light curve
- The cat-normalised depth is 17.9 +- 1.6 ppt, against the ~21.4 ppt expected if all of the pair decrement is on the
  neighbour. psf_photometry attributes the difference to PRF/registration systematics (f_target = 0.16, consistent
  with 0 within its systematics).
- If 16 % of the decrement really were on TIC 87055304, the neighbour's depth would be ~18.0 ppt and k ~ 0.128. This
  is carried below as the "84 %" variant.

## 3. Nature of the companion (`scripts/s06_pair_tests.py`, `s09_eb2p.py`, `s11_companion.py` -> `results/pair_tests.json`, `eb2p.json`, `companion.json`)

### 3a. Density
- **Comparison.** rho_transit (nominal) = 0.33 (+0.10/-0.07) against rho* = 0.172 (+0.035/-0.027). The ratio is
  1.9 (+0.7/-0.5) and z = 2.3 in log space; only 10 % of the rho_transit posterior lies inside the star's 95 % interval.
  The MCMC errors assume white noise. With the pair light curve's beta = 1.2-1.4 on the ingress/egress timescale they
  should be inflated by about that factor, which brings z to ~1.7-1.9.
- **QLP-based fit.** The first round's QLP-based fit (0.19 +0.05/-0.04) agrees with rho* directly.
- **Reading.** I read this as consistent within the extraction systematics (2c), not as evidence against the neighbour
  hosting.
- **Alternatives considered.**
  - A host denser than the neighbour (an unresolved K/M dwarf) would need rho ~ 1-5, against the measured 0.2-0.4.
  - Matching rho* with an eccentric orbit would need e sin(omega) ~ +0.2. That is unlikely at P = 4.1 d with a
    stellar companion, since such orbits are expected to be circularised.
- **Verdict:** neutral (no contradiction with the neighbour hosting).

### 3b. Companion radius
- **Result.** R2 = k R* = **0.269 +- 0.013 R_sun = 2.61 +- 0.13 R_Jup** (2-sigma: 2.37-2.89 R_Jup), computed from
  20,000 independent draws of the nominal k and the SED radius.
- **Error budget.** The uncertainty is dominated by R* (4.4 %). Inflating the k error for red noise (x1.4) changes the
  total only from 4.9 % to 5.4 %.
- **Cross-checks.**
  - The first-round k (0.1414) with this R* gives 2.65 R_Jup.
  - The 84 % variant gives 2.40 +- 0.12 R_Jup (2-sigma 2.17-2.65).
  - With Gaia FLAME R = 2.03, the result would be 2.75 R_Jup.
- **Verdict.** Too large for any known planet, and far too large for a several-Gyr-old brown dwarf. It is the size of
  a ~0.25 M_sun M dwarf.

### 3c. Secondary eclipse, brightness temperature and luminosity limits
- **Phase 0.5 depth.**
  - Method: occultation template (the nominal shape with uniform disc), the same per-transit quadratic baseline
    estimator, and errors from max(white, epoch scatter, 50 control phases). The pair light curve was re-detrended with
    both eclipse phases masked.
  - **Pair-based, neighbour scale: 0.79 +- 0.77 ppt** (2-sigma upper limit 2.33 ppt). This is 2.4x tighter than
    psf_photometry's neighbour-only 0.30 +- 1.87 ppt (limit 4.04 ppt).
- **Phase scan, 0.30-0.70.** The maximum |S/N| is 1.23 at phase 0.62, and the 2-sigma limit anywhere is 2.73 ppt
  (`fig/pair_secondary_scan.png`).
- **Limits on the companion.**
  - Method: blackbody companion against the CK04 TESS-band surface flux of the star (5850 K, logg 3.95), with
    k^2 = 0.0194; reflected light is ignored, which is conservative.
  - **T_b < 3510 K** (2 sigma, pair); < 3930 K (psf limit); < 3620 K at any phase.
  - **L2 < 0.010 L_sun** (pair) or < 0.016 L_sun (psf), for R2 = 0.269 R_sun.
- **Expectations.**
  - A 0.25-0.30 M_sun M dwarf (3300-3400 K) gives 1.7-2.0 ppt. That is allowed, about 1.2 sigma above the measured
    0.79 +- 0.77.
  - A brown dwarf or planet at T_eq ~ 1500 K gives 0.004 ppt.
  - An M dwarf hotter than ~3500 K (M > 0.4 M_sun, R > 0.38 R_sun) is disfavoured, but its radius would not match k
    anyway.
- **Verdict:** consistent with an M dwarf and with a dark companion. It does not discriminate between them, but it
  excludes any luminous (T > 3.5-3.9 kK) companion.

### 3d. Odd/even and an equal-mass EB at 2P
- **Odd/even** (pair, absolute epoch parity): odd 21.65 +- 1.17, even 21.31 +- 1.05 ppt. The difference is
  0.34 +- 1.57 ppt (0.2 sigma). psf_photometry's neighbour-only light curve gives -4.0 +- 3.2 ppt (1.3 sigma).
- **Twin EB at 2P = 8.27 d** (`s09`). The model has:
  - two identical stars (k = 0.9-1.1, equal surface brightness), circular orbit, with the primary and secondary
    identical;
  - free b, a/R, k, LD, T0, baselines, and a free third-light fraction (f_eb).

  It was fitted by maximum likelihood from 40 starts, on the same data and noise as the transit fit.
  - Best twin: b = 1.83 (grazing), f_eb = 1.0, k = 1.10 (at the prior edge).
  - **chi2 is worse by 49.5** than the best flat-bottomed transit model, which has one free parameter fewer (D is
    fixed by its prior).
  - Red noise of the pair light curve (`scripts/beta_pair.py` -> `results/beta_pair.json`; out-of-eclipse binned rms
    against white): beta = 1.10, 1.21, 1.38 and 1.51 for 0.5-, 1-, 2- and 4.5-h bins. Dividing the white-noise Δchi2 of
    49.5 by beta^2 = 1.5-2.3 gives 22-34, i.e. **4.7-5.8 sigma**.
- **V-shape test.** A grazing transit model (b > 1 - k) at P is also worse, by Δchi2 = 32.7 (3.8-4.7 sigma after the red-noise correction).
- **Verdict.** The eclipse is flat-bottomed and non-grazing. **An equal-mass EB at 2P is excluded by the shape**;
  odd/even equality alone could not exclude it.

### 3e. Ellipsoidal, beaming and reflection (phase curve)
- **Method.**
  - Fit: raw pair flux per segment (25 segments of 2-7 d) times a quadratic, plus cos 4 pi phi, sin 2 pi phi and
    cos 2 pi phi, with both eclipses masked.
  - Errors: the scatter at 60 control periods (2.6-7 d).
  - Injection: a 1-ppt ellipsoidal signal is recovered with a factor of 1.000.
- **Results**, converted to the neighbour scale x (1 + D):
  - **A_ell = 0.87 +- 0.66 ppt** (2-sigma upper limit 2.2 ppt);
  - A_beam = 3.7 +- 1.8 ppt, which is systematics: it would need K ~ 280 km/s;
  - A_refl = -2.4 +- 1.9 ppt.
- **Expectations**, with alpha_ell = 1.33 and a/R = 7.5:
  - 0.25 M_sun gives A_ell = 0.67 ppt and A_beam = 0.35 ppt;
  - 0.075 M_sun gives A_ell = 0.20 ppt;
  - a 1 M_Jup planet gives 0.003 ppt.
- **Mass limit.** M2 < 0.82 M_sun (2 sigma), which is uninformative between an M dwarf and a brown dwarf.
- **Verdict:** inconclusive. It is consistent with the M dwarf, and it excludes only a near-equal-mass companion.

### 3f. Blend within ~1" of the neighbour
- **Shape limit on extra dilution.** The ingress/egress duration limits k to <~ 0.19 (b = 0). The undiluted depth can
  therefore be at most ~3.6 %, so any unresolved host must supply >= ~60 % of the neighbour's TESS flux. It would be
  nearly as bright as TIC 87055315.
- **Gaia constraints.** RUWE 1.02, ipd_frac_multi_peak 1 % and the absence of an NSS solution argue against such a
  companion beyond ~0.1-0.2". The MSC two-star solution (738 pc) is inconsistent with the parallax.
- **What remains open.** A close bound companion with similar brightness (a hierarchical triple) cannot be excluded.
  Even then, the eclipser would be ~0.2-0.3 R_sun, i.e. stellar.

## 4. TRICERATOPS with TIC 87055315 as target (`scripts/s08_tri_prep.py`, `tri/run_tri_queue.sh`, `tri/tri_*.json`)

**Light-curve choice.**
- The fixed runner calls `calc_depths(tdepth)` and `calc_probs(...)` with TRICERATOPS' default
  **dilution_corrected=True**. The light curve must therefore be normalised to the *target's* out-of-transit flux, i.e.
  the on-target, PDCSAP-like convention, and TRI_DEPTH must be the on-target depth.
- TRICERATOPS then computes the dip in the aperture itself, as tdepth x fluxratio_target (TIC 87055304 is inside its
  5x5 aperture), and builds each nearby star's light curve from it.
- I therefore gave it **the neighbour's undiluted light curve**: the PSF pair light curve rescaled by
  f_nb = 1 + (f_pair - 1)(1 + D), D = 5.32 (`data/tri_lc_nb_200s.npz`).
- Feeding the blended SAP or pair light curve (3.4 ppt) to this runner would make TRICERATOPS treat the 3.4-ppt dip as
  the neighbour's own depth, 6x too shallow. That would favour planets and faint nearby hosts.
- Only S66, S93 and S104 (200 s) were used. TRICERATOPS models a 2-min exposure, so the 30-min (S13) and 10-min (S39)
  data would smear ingress and egress.
- **TRI_DEPTH = 0.0212**, the MCMC undiluted depth. The runner's automatic estimate (1 - min of a 3-bin running mean)
  is noise-biased for this light curve.
- Settings: point sigma 0.0149, N = 200,000, TRI_RADIUS = 10 px, sectors 66, 93 and 104, P = 4.1362006 d.

| run | stellar parameters | FPP | NFPP | top scenarios (prob) |
|---|---|---|---|---|
| A (`tri_A_ticparams`) | TIC v8.2: R 1.906, M 1.05, plx 0.761 | **1.000** | 0.992 | NTP on TIC 87055367 (T 13.23, 34") 0.34; NTP on TIC 87055378 (T 12.99, 42") 0.28; NEB on 87055378 0.16; NEB on 87055367 0.14; NEB on 87055318 0.03; NTP on 87055304 0.02; EB on target 0.006 |
| B (`tri_B_adopted`) | this work: R 1.93, M 1.17, plx 0.7785 | **1.000** | 0.749 | NTP on 87055378 0.23; **EB on target 0.16**; NEB on 87055378 0.16; NTP on 87055367 0.15; NEB on 87055367 0.12; **PEB 0.06**; NTP on 87055304 0.05 |
| C (`tri_C_adopted_repeat`) | as B (Monte Carlo repeat) | **1.000** | 0.800 | NEB on 87055378 0.21; NTP on 87055378 0.20; NTP on 87055367 0.20; **EB on target 0.15**; NEB on 87055367 0.11; NTP on 87055304 0.06; **PEB 0.04** |
| D (`tri_D_ontarget_adopted`; `scripts/run_tri_ontarget.py`) | as B, nearby stars removed (signal assumed on TIC 87055315, as the localisation shows) | **1.000** | 0 (by construction) | **EB 0.652** (best M_EB 0.248 M_sun, R_EB 0.273 R_sun, b 0.55); **PEB 0.218**; **SEB 0.094**; DEB 0.035; all planet scenarios (TP, PTP, STP, DTP, BTP) and all x2P scenarios 0.000 |

Summary over runs B and C: FPP = 1.000 +- 0.000 and NFPP = 0.77 +- 0.03.

- **FPP = 1 in every run.** A 2.1 % flat-bottomed transit on a 1.93 R_sun star needs R_p ~ 29 R_Earth. TRICERATOPS'
  planet-radius prior stops at 20 R_Earth, and its best TP fit sits at that edge (R_p = 19.98 R_Earth). No planet
  scenario can reproduce the depth, so the conclusion is driven by the companion's size. This is consistent with 3b.
- **NFPP is high (0.75-0.99) only because of priors.** TRICERATOPS does not use the pixel-level localisation. It
  assigns most of the probability to the bright stars TIC 87055367 and 87055378 (34" and 42", ~1.6-2 TESS pixels) and
  to TIC 87055304. The three localisation analyses exclude all of these: psf_photometry's pixel-map chi2 is >= 164
  for every other star within 60", against 53 for the neighbour; loc_independent puts the source 0.66" from the
  neighbour.
- **Conditional on the localisation** (run D), the eclipse is an EB on TIC 87055315 (65 %), or an EB on an unresolved
  bound companion (PEB 22 %, SEB 9 %). TRICERATOPS' best-fitting EB companion (0.25 M_sun, 0.27 R_sun) agrees with
  section 3.
- The full-run spread (A against B/C: NFPP 0.99 against 0.75-0.80) is driven by the stellar parameters and Monte
  Carlo noise. The runner saves only the top 8 scenarios, so the B/C breakdowns are partial.
- **Verdict:** supports false positive (as a planet). TRICERATOPS adds little beyond the radius argument.

## 5. Novelty (`scripts/s05_novelty.py` -> `results/novelty.json`; `data/novelty/`)
- **ExoFOP-TESS JSON** (target.php?id=87055315&json, fetched 2026-10-12): **no TOIs, no CTOIs**, and no planet
  parameters, imaging, spectroscopy, time series, files or stellar companions. The JSON has no notes field.
- **MAST.** The only TESS products for TIC 87055315 are the S13 QLP and GSFC-eleanor-lite HLSPs. There are no SPOC
  2-min, TESS-SPOC FFI or DV products, so there is no SPOC TCE.
- **VizieR.** Each catalogue's existence was first verified with an unconstrained 1-row query. Cone searches at the
  DR3 position, with radii of 3" and 21":
  - Kostov+25 (J/ApJS/279/50) tables 2, 3 and 4: 0 and 0;
  - VSX (B/vsx/vsx): 0 and 0;
  - ASAS-SN variables (II/366/catv2021): 0 and 0;
  - Gaia DR3 vari EB (I/358/veb) and classifier (I/358/vclassre): 0 and 0;
  - Prsa+22 TESS EBs: 0 and 0.
  - Note: the VizieR ID tried for vari_summary (I/358/vari) does not exist; the Gaia archive query of
    gaiadr3.vari_summary by source_id returned 0 rows.
- **Gaia archive** (by source_id): vari_summary, vari_eclipsing_binary, vari_classifier_result, nss_two_body_orbit,
  nss_acceleration_astro, nss_non_linear_spectro and nss_vim_fl all have **0 rows**.
- **Green+23** (J/MNRAS/522/29/table2) constrained by TIC: **TIC 87055315 is not in the table**. TIC 87055304, whose
  light curve contains the diluted signal, is in it with Score = 0 and Porb = 2.04 d.
- **SIMBAD.** No object for "Gaia DR3 6720986208895208960", "TIC 87055315" or "2MASS J18102616-4359369", and nothing
  within 10". The query works: a 20" cone returns TYC 7912-71-1, the target.
- **Web search** (WebSearch, 2026-10-12), for "TIC 87055315", "6720986208895208960", "2MASS J18102616-4359369", and
  "Gaia DR3 6720986208895208960" OR "TIC 87055304" OR "TYC 7912-71-1": **no page mentions any identifier.**
- **Verdict.** As far as these catalogues go, the 4.136-d eclipsing binary TIC 87055315 is uncatalogued. It is a new
  EBLM-type candidate, not a new planet candidate.

## 6. Follow-up and predicted eclipse times (`scripts/s10_ephem.py` -> `results/ephemeris.json`)

**Decisive test: seeing-limited photometry** (1-2 m class, r'/i'/z', FWHM < 4"; the two stars are 9.3" apart).
- The neighbour (G = 14.0) should show a flat-bottomed ~2.1 % eclipse lasting 4.35 h. In a red band the depth is the
  same within ~1 %, because an M-dwarf companion adds < 1 % of the light.
- TIC 87055304 (G = 12.2) should show < 0.5 ppt.
- One full event decides the host. Two RV epochs (K ~ 27 km/s expected for 0.26 M_sun) would then confirm the M dwarf.

**Ephemeris and uncertainties.**
- Ephemeris A (the brief): P = 4.136195 +- 0.000007 d, T0 = 1658.41325 +- 0.0027 BTJD.
- sigma(T_n) = sqrt(sigma_T0^2 + (n sigma_P)^2), which is ~8 min in 2027.
- This work's local-baseline MCMC, ephemeris B (P = 4.1361919 +- 0.0000053, T0 = 3271.53062 +- 0.0010 BTJD at
  E = 390), agrees with ephemeris A to within 0.6 min in 2027, with sigma ~2.7 min.

**Visibility** (CTIO, SAAO, Siding Spring). An event is counted as fully observable if the target is above 30 deg and
the Sun below -12 deg from T1 - 1 h to T4 + 1 h.
- **No full event is observable in late 2026.** The target sets in the evening, and solar conjunction is around
  21 December.
- Partial chances in late 2026: epoch 647 at 2026-10-21 00:47 UTC (BTJD 4334.5314 +- 0.0053), with 49 % of the eclipse
  observable from CTIO (evening, setting); and epoch 646 at 2026-10-16 21:30 UTC, with 23 % from SAAO.
- The first full windows are in April 2027.

| epoch | BTJD (eph. A) | sigma (min) | UTC mid-eclipse | site with full window |
|---|---|---|---|---|
| 693 | 4524.7964 | 8.0 | 2027-04-29 07:00 | CTIO |
| 698 | 4545.4774 | 8.0 | 2027-05-19 23:19 | SAAO |
| 703 | 4566.1583 | 8.1 | 2027-06-09 15:39 | Siding Spring |
| 707 | 4582.7031 | 8.1 | 2027-06-26 04:43 | CTIO |
| 712 | 4603.3841 | 8.2 | 2027-07-16 21:04 | SAAO |

More events, through 2027-08-23, are in `results/ephemeris.json`. UTC times include the barycentric correction toward
the star.

## 7. Open issues and caveats
- **The mass is not measured.** "M dwarf" rests on the radius (2.6 R_Jup), the occultation limit and plausibility.
  A very young, inflated brown dwarf cannot reach 2.6 R_Jup at the host's age, but only RVs give the mass.
- **Duration and density** depend on extraction and detrending (rho_transit 0.19-0.39). The nominal value is 2.3 sigma
  above rho*. I do not consider this significant given the ~0.1 spread between methods.
- **Neighbour-only depth.** psf_photometry finds 17.9 +- 1.6 ppt, against 21.4 ppt if all of the decrement is on the
  neighbour. The companion radius changes by -8 % in the 84 % variant, which does not change the conclusion.
- **Stellar mass.** The mass of the neighbour (1.17 +0.17/-0.10) comes from a coarse-[Fe/H] MIST grid. FLAME (1.15),
  SHBoost (1.21) and StarHorse21 (1.23) agree.
- **Gaia RV.** The error is unremarkable (64th percentile), so Gaia gives no RV-variability constraint either way.
- **TRICERATOPS.** It does not use the pixel-level localisation. Its NFPP for TIC 87055315 as target reflects only
  priors on the nearby bright stars (section 4).

## Files
- **Scripts** (`scripts/`):
  - `s01_catalogs.py`: TIC, Gaia archive and VizieR queries, cached in `data/cat/`;
  - `s02_lc.py`: light curves;
  - `s03_fluxratio.py`: r and the QLP dilution;
  - `s04_stellar.py`: SED + MIST;
  - `s05_novelty.py`;
  - `s06_pair_tests.py`: secondary, odd/even and phase curve;
  - `s07_gaia_rv.py`;
  - `s08_tri_prep.py`;
  - `s09_eb2p.py`;
  - `s10_ephem.py`;
  - `s11_companion.py`;
  - `s12_collect.py`: writes `results.json`;
  - `nb_mcmc_local.py` (nominal MCMC), `nb_mcmc.py` with `nbmodel.py` (biweight MCMC; copies of transit/ code), and
    `diag_t14*.py`;
  - `run_tri_ontarget.py`: a copy of the fixed runner plus `remove_nearby_stars()`;
  - helpers `nbtools.py`, `psfcommon.py` and `sedlib.py`, which are copies.
- **Results:**
  - `results.json`, the top-level summary;
  - `results/*.json`;
  - `mcmc/mcmc_pair_local.{json,npz}` and `mcmc/mcmc_pair_Dfree.{json,npz}`;
  - `tri/tri_*.json` and `tri/tri_*.log`.
- **Figures:**
  - `fig/sed_hrd_nb.png`: SED and MIST posterior;
  - `fig/fluxratio_dilution.png`;
  - `fig/eb2p_fold.png`: transit versus twin-EB fold and residuals;
  - `fig/pair_secondary_scan.png`;
  - `fig/pair_phasefold.png`.
- **Logs** are in `logs/`.
