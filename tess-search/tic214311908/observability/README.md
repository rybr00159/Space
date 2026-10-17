# Observability tables, TIC 214311908 transits 2026-11 to 2028-12

One plain CSV per site (header in the first line). Each lists the epochs with any usable time at the site (Sun below
-18 deg and the target at airmass < 2). All times come from the adopted ephemeris (`../ephemeris_adopted.json`) and were
computed by `../scripts/w02_predictions.py`. A window counts as full when T14 (1.70 h) + 2 x 3 sigma of timing margin +
2 x 1.5 h of baseline is entirely usable. AI-run analysis; not reviewed by a professional astronomer.

| File | Site | Coordinates (deg, deg, m) |
|---|---|---|
| `obs_CTIO_SOAR_GeminiS.csv` | CTIO / Cerro Pachon (SOAR, Gemini-S), Chile | CTIO: lat -30.169, lon -70.8063, 2207 m; Cerro Pachon: lat -30.2408, lon -70.7367, 2722 m |
| `obs_LCO.csv` | Las Campanas Observatory, Chile | Las Campanas: lat -29.0146, lon -70.6926, 2380 m |
| `obs_La_Silla.csv` | La Silla (TRAPPIST-S, ExTrA, ESO 3.6 m), Chile | La Silla: lat -29.2563, lon -70.738, 2400 m |
| `obs_Paranal.csv` | Paranal (VLT, SPECULOOS-South), Chile | Paranal: lat -24.6272, lon -70.4042, 2635 m |
| `obs_SAAO.csv` | SAAO Sutherland (LCOGT, SALT), South Africa | SAAO: lat -32.3794, lon 20.8107, 1798 m |
| `obs_Siding_Spring.csv` | Siding Spring (FTS/MuSCAT4, LCOGT, AAT), Australia | Siding Spring: lat -31.2733, lon 149.0614, 1165 m |
| `obs_Dome_C.csv` | Dome C (ASTEP+), Antarctica | Dome C: lat -75.0997, lon 123.3325, 3233 m |

Columns:
* `epoch`: transit number from the adopted T0; `bjd_tdb_mid`: mid-transit, target-barycentred BJD_TDB.
* `utc_mid`: geocentric UTC of mid-transit, to the second. The other UTC columns are rounded to the minute:
  `utc_window_start`, `utc_window_end` (full window), `utc_ingress`, `utc_egress` (T14 only, no margin).
* `sigma_min`, `sigma3_min`: 1 sigma and 3 sigma timing uncertainty (min); `window_h`: window length (h).
* `full_window`: the whole window is usable with the Sun below -18 deg (all member sites); `full_window_sun12`: the
  same with the Sun below -12 deg; `transit_covered`: transit +- 3 sigma usable; `good_partial`: transit +- 3 sigma
  usable with at least 0.75 h of baseline on each side, but not a full window.
* `baseline_pre_h`, `baseline_post_h`: usable baseline before and after (h); `frac_window_usable`: usable fraction
  of the window.
* `alt_mid`, `alt_min_window` (deg), `airmass_max_window`, `sun_max_window` (deg): target and Sun during the window.
* `moon_illum` (0-1), `moon_sep_deg`, `moon_alt_mid` (deg): the Moon at mid-transit.
