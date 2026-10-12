# TIC 87055304 — shared context for the adversarial vetting campaign

Goal: try every reasonable test that could prove the planet hypothesis WRONG for the transit-like signal on
TIC 87055304, with rigour at least equal to the TIC 4206066 manuscript (text: /home/user/work/tic87/reference_manuscript.txt).
People will try to poke holes in this because it was found with AI: every number needs a method, an uncertainty,
and a reproducible script. Report negative results honestly. Do not oversell.

## Target
- TIC 87055304 = TYC 7912-71-1 = Gaia DR3 6720986204578421888. RA 272.606557, Dec -43.995371 (J2000, TIC);
  Gaia DR3 epoch-2016 position RA 272.60647748, Dec -43.99546646; pm (-5.09, -16.20) mas/yr.
- TESS mag T = 11.72; Gaia G = 12.186, BP = 12.513, RP = 11.652, BP-RP = 0.861; phot_bp_rp_excess_factor = 1.2254.
- Gaia DR3 astrometry: parallax 2.807 +- 0.357 mas, RUWE 18.84, astrometric_excess_noise 2.32 mas (sig 5824),
  visibility_periods_used 16, ipd_frac_multi_peak 0, ipd_gof_harmonic_amplitude 0.049, ipd_frac_odd_win 0,
  duplicated_source False, non_single_star 0 (no NSS solution of any kind).
- Gaia DR3 RVS: RV -13.33 +- 0.92 km/s, rv_nb_transits 6, rv_amplitude_robust 5.20 km/s, rv_chisq_pvalue 0.932,
  rv_renormalised_gof -1.56, rv_template_teff 5500, rv_template_logg 3.5. has_rvs False; has_xp_continuous True.
- Stellar parameters disagree:
  - TIC v8.2: Teff 7276 K, R 1.79 R_sun, mass missing, d 464 pc (DR2-based).
  - Gaia GSP-Phot: Teff 5696 K, logg 3.50, [M/H] -0.52, R 3.49, distance_gspphot 998 pc (photometric distance;
    conflicts with the parallax distance of 356 pc).
  - Gaia FLAME: R 1.24, M 0.95, age 11.2 Gyr, evolstage 380 (FLAME uses the parallax, which RUWE 18.8 may corrupt).
  - GSP-Spec: none.
- Contamination: TIC contamination ratio 0.51. Gaia neighbours (data/gaia_dr3_90arcsec.csv):
  - G 14.02 at 8.94", the only bright one (TIC 87055315, T ~ 13.5);
  - G 20.94 at 4.35" and G 18.61 at 5.29", both T ~ 17.7-17.9;
  - many fainter stars (Galactic latitude ~ -12 deg; crowded).

## Signal (from earlier work; to be re-derived independently where relevant)
- P = 4.136195 +- 0.000007 d, T0 = 1658.41325 +- 0.0027 BTJD (BJD_TDB - 2457000).
- Duration ~ 4.5 h (0.1875 d); flat-bottomed.
- Sectors 13 (30-min FFI), 39 (10-min), 66, 93 and 104 (200-s). Roughly 25 transits, 2019-2026.
- QLP SAP flux with a 1-d transit-masked biweight: depth 3,580 +- 120 ppm (combined).
  - Per-sector depths 4101, 4238, 3792, 3629 and 2939 ppm for S13, 39, 66, 93, 104.
  - S104 used a different QLP aperture ('3x3 square'; the other sectors used '2.5:4.0:3.0', i.e. a 2.5 px radius).
- Batman least-squares fit on SAP: Rp/R* 0.0668 +- 0.0022 (not corrected for dilution), b 0.85 +- 0.04,
  a/R* 4.30 +- 0.44, rho* 0.063 rho_sun (circular orbit).
- Odd 3529 +- 172, even 3624 +- 169 ppm. Secondary at phase 0.5: 109 +- 203 ppm.
  - A full-phase scan (0.05-0.95) found one feature at phase 0.39. It is a single scattered-light event at the start
    of S104 (BTJD ~4178.97; background 1.8x normal) and should be masked.
- Phase-curve harmonics at P are not significant: cos2phi -267 +- 32 formal, but the rms at 40 control periods is
  209-241 ppm.
- QLP DET_FLUX (and KSPSAP_FLUX) shows +0.5-1 ppt "shoulders" around each transit and underestimates the depth by
  25-45 %. It is a detrending artefact absent in SAP flux, so use SAP or your own photometry, not DET_FLUX, for depths
  and shapes.
- QLP-pipeline difference-image centroids (our earlier quick method, src/centroid.py): S66 1.1" (0.3 sigma) and
  S104 1.7" (0.3 sigma) from the target; S93 failed (NaN).
- Earlier TRICERATOPS runs (single run each, N = 200,000):
  - QLP DET light curve, TIC star (R 1.79, M 1.3 set by hand): FPP 0.037, NFPP < 1e-4 (TP 0.78, PTP 0.16, STP 0.04).
  - SAP light curve, R 2.7, M 1.3: FPP 0.52 (TP 0.38, EB 0.31, PEB 0.18, PTP 0.08, DEB 0.03, DTP 0.02).
- Not a TOI/CTOI/confirmed planet (nothing within 3'), not a SPOC TCE, not in Kostov+25 eclipse-like list, Green+23
  ellipsoidal score 0, no variability catalogue within 15". SIMBAD: TYC 7912-71-1 type '*'.
- Gaia RVs constant (6 epochs). A transiting companion of stellar mass on the target would give K ~ 10-30 km/s.

## Data (all under /home/user/work/tic87/data)
- `qlp_raw_columns.npy`: dict sector -> all QLP columns (TIME, SAP_FLUX, DET_FLUX, KSPSAP_FLUX if present, QUALITY,
  SAP_BKG, SAP_X, SAP_Y, ORBITID, CADENCENO, ...). Load with np.load(..., allow_pickle=True).item().
- `lc_sap_clean.npy`: dict sector -> {t, f (SAP / transit-masked 1-d biweight), f_raw, trend, exptime, sigma}.
  Quality == 0 only; native cadence.
- Original QLP FITS: `mast_dl/mastDownload/HLSP/hlsp_qlp_tess_ffi_s00XX-0000000087055304_tess_v01_llc/*.fits`.
- Independent S13 pipeline: `mast_dl/.../hlsp_gsfc-eleanor-lite_tess_ffi_s0013-0000000087055304_tess_v1.0_lc.fits`.
- The 9" neighbour TIC 87055315 has its own S13 QLP and eleanor-lite light curves in `mast_dl/`.
  If the signal were on that neighbour, its own light curve would show it deeper.
- TESScut 21x21 cutouts, all 5 sectors (calibrated FFI cutouts; FLUX, FLUX_ERR, QUALITY, TIME; WCS in ext 2):
  `tesscut/tess-s00XX-*_21x21_astrocut.fits`.
- `gaia_dr3_90arcsec.csv`: full Gaia DR3 gaia_source rows within 90" (sep column in arcsec).
- `gaia_dr3_astrophysical_parameters_target.csv`: the target's full Gaia DR3 astrophysical_parameters row.

## Software
- Python: `/home/user/work/venv/bin/python`. Installed: numpy, scipy, astropy, astroquery, batman, emcee, wotan,
  lightkurve, exovetter, leo_vetter, PRF (TESS_PRF package, `import PRF`), tess_stars2px (tess-point), celerite2,
  numba, matplotlib, pandas.
  - Install more with `uv pip install --python /home/user/work/venv/bin/python <pkg>`; there is no pip in the venv.
- TRICERATOPS lives in a separate venv: `/home/user/work/venv_tri/bin/python`. Runner: /home/user/work/src/run_tri.py.
  - Args: tic sectors P t0 dur lcfile.npz out.json; env TRI_MASS, TRI_RAD, TRI_PLX, TRI_RADIUS.
  - The lc npz needs t, f and sigma; the runner folds and bins to 5 min.
- Existing helpers (read them, copy them, do not edit them in place) in /home/user/work/src:
  - vet.py: batman model (`model(t, P, t0, rp, a, b, exptime)`), least-squares `fit`, `rho_star(a, P)`,
    `gaia_neighbours`.
  - centroid.py: TESScut difference-image localisation. offset_test.py: directional test against a neighbour.
  - stage2.py: multi-sector loading. run_tri.py: TRICERATOPS.
- Network: MAST, Gaia archive (slow at times; use timeouts and retries), VizieR, SIMBAD and ExoFOP are reachable
  through the proxy. Do NOT use OpenRouter or any paid API.

## House rules for every investigator
- Work ONLY under /home/user/work/tic87/<your_area>/ (create it). Never modify files in /home/user/work/src, data/,
  or another investigator's folder. Never touch the git repository /home/user/Space.
- Machine: 4 CPU cores and 15 GB RAM, shared by up to 8 investigators.
  - Prefix heavy commands with `nice -n 10`, set OMP_NUM_THREADS=1, keep any single job under ~30 min wall time.
  - Run long jobs in the background, writing logs to your folder. Never run more than 2 heavy processes at once.
- Never kill processes you did not start. Do not use `pgrep -f`/`pkill -f` with patterns that also match your own
  shell command (this has bitten us before); kill by PID.
- Time convention: BTJD = BJD_TDB - 2457000.
- Every result must give: the exact method, the number with uncertainty, the script path that produced it, and an
  honest verdict.
  - Verdict scale: supports planet / neutral / supports false positive / inconclusive.
  - Save figures as PNG in your folder.
- Adversarial mindset: your job is to try to break the planet hypothesis, and to check whether earlier numbers in this
  file are wrong.

## UPDATE (after the first, interrupted round): the signal may be on the 8.94" neighbour
A previous investigator (folder /home/user/work/tic87/transit/) found, before being interrupted, that the transit
signal appears to originate on the 8.94" neighbour.
- Neighbour: Gaia DR3 6720986208895208960 = TIC 87055315; G 14.02, T ~ 13.54, parallax 0.735 mas, RUWE 1.02,
  TIC Teff 5821 K, TIC R 1.906 R_sun.
- Evidence:
  - Two-source PRF amplitude fits put almost all of the signal on the neighbour: f_target = 0.02 +- 0.04 combined.
  - Difference-image centroids sit near the neighbour's sky position in all 5 sectors.
  - Injection tests reproduce the real behaviour only when the injected signal is on the neighbour.
  - Per-pixel depth maps are fit best by the neighbour (delta chi2 75-145 per sector versus the target).
- An MCMC with the neighbour as host (dilution 4.35) gives k = 0.141, b = 0.64, rho = 0.19 rho_sun.
- This contradicts the earlier quick claim (src/centroid.py; README) that the source is on target and the 9"
  neighbour is "excluded at >6 sigma". That earlier claim must be treated as unverified and possibly wrong.
- The decisive question is now: which star hosts the signal? Everything else is secondary.
