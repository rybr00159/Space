"""Follow-up feasibility for TIC 214311908.

1. Recompute predicted mid-transit times from the ephemeris (ephemeris_v2.json; its covariance is already scaled by
   chi2/dof = 53.4/11) for 2026-10-17 .. 2027-12-31 and compare with predicted_transits_2027.csv.
2. Observability from CTIO, La Silla, SAAO (Sutherland) and Siding Spring: target altitude > 30 deg and Sun below
   -12 deg (also -18 deg) for the whole window Tc +- (T14/2 + 1 h), T14 = 1.63 h; also the stricter window that adds
   2 sigma of timing uncertainty on each side. Moon illumination and separation reported.
3. Photometric precision of 0.4 m and 1.0 m telescopes on the M3V (SkyMapper DR4 PSF mags g 15.82, r 14.91, i 13.39,
   z 12.77) with a CCD noise model (photon, sky, read noise, scintillation with the Osborn et al. 2015 factor 1.5 on
   Young's formula), and the per-event transit S/N for depths 1.5% and 1.8% with white noise and with a 1 ppt red-noise
   floor per 30-min block (Pont et al. 2006 style).
4. Contamination from the Mira (17.8") and the M-dwarf pair (11.7") for Moffat PSFs, and Mira saturation estimates.
Outputs: out/followup_transits.csv, out/followup.json
"""
import os, sys, json, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun, get_body
import astropy.units as u
import common as C

out = {}
TARGET = SkyCoord(260.9696737 * u.deg, -51.3238212 * u.deg)   # Gaia J2016; PM 0.14"/yr irrelevant here
cov = np.array(C.EPH["cov"])
# ---------------------------------------------------------------- 1. ephemeris
start, stop = Time("2026-10-17").jd - 2457000, Time("2028-01-01").jd - 2457000
ks = np.arange(np.ceil((start - C.T0) / C.P), np.floor((stop - C.T0) / C.P) + 1).astype(int)
tc = C.T0 + ks * C.P
sig = np.sqrt(cov[0, 0] + ks ** 2 * cov[1, 1] + 2 * ks * cov[0, 1]) * 1440
csv = pd.read_csv("/home/user/Space/tess-search/vetting_2026-10-16/results/214311906/predicted_transits_2027.csv")
m = pd.merge(csv, pd.DataFrame(dict(epoch=ks, mine=tc + 2457000, mysig=sig)), on="epoch")
out["csv_check"] = dict(n_csv=len(csv), n_matched=len(m), max_abs_diff_min=float(np.max(np.abs(m.bjd_tdb - m.mine)) * 1440),
                        max_sigma_diff_min=float(np.max(np.abs(m.sigma_min - m.mysig))),
                        csv_epochs=[int(csv.epoch.min()), int(csv.epoch.max())],
                        note="CSV bjd_tdb is rounded to 1e-4 d (0.14 min); sigma already includes the chi2/dof = 4.86 inflation")
# BJD_TDB -> UTC at the observatory: barycentric correction for the target (subtract light-travel time)
SITES = {"CTIO": (-30.1690, -70.8063, 2207), "La Silla": (-29.2563, -70.7380, 2400),
         "SAAO Sutherland": (-32.3794, 20.8107, 1798), "Siding Spring": (-31.2733, 149.0617, 1165)}
half_win = C.T14 / 2 + 1 / 24.0          # d
rows = []
for name, (lat, lon, h) in SITES.items():
    loc = EarthLocation.from_geodetic(lon * u.deg, lat * u.deg, h * u.m)
    for k, t_bjd, s_min in zip(ks, tc, sig):
        # convert BJD_TDB to TDB at the site: t_site = t_bjd - ltt(t_site); iterate once
        tt = Time(t_bjd + 2457000, format="jd", scale="tdb", location=loc)
        ltt = tt.light_travel_time(TARGET, kind="barycentric")
        tsite = Time(t_bjd + 2457000 - ltt.to(u.day).value, format="jd", scale="tdb", location=loc)
        ext = 2 * s_min / 1440
        grid = tsite + np.linspace(-(half_win + ext), half_win + ext, 61) * u.day
        fr = AltAz(obstime=grid, location=loc)
        alt = TARGET.transform_to(fr).alt.deg
        sun = get_sun(grid).transform_to(fr).alt.deg
        core = np.abs((grid - tsite).to(u.day).value) <= half_win + 1e-9
        moon = get_body("moon", tsite, loc)
        sunc = get_sun(tsite)
        elong = sunc.separation(moon).rad
        illum = float((1 - np.cos(elong)) / 2)
        moon_sep = float(moon.separation(TARGET).deg)
        ok12 = bool(np.all(alt[core] > 30) and np.all(sun[core] < -12))
        ok18 = bool(np.all(alt[core] > 30) and np.all(sun[core] < -18))
        ok12_ext = bool(np.all(alt > 30) and np.all(sun < -12))
        dark_core = (alt[core] > 30) & (sun[core] < -12)
        rows.append(dict(site=name, epoch=int(k), bjd_tdb=round(t_bjd + 2457000, 5), utc_mid=tsite.utc.iso[:16], sigma_min=round(s_min, 1),
                         full_window_ok_sun12=ok12, full_window_ok_sun18=ok18, full_window_plus_2sigma_ok=ok12_ext,
                         frac_window_observable=round(float(dark_core.mean()), 2),
                         alt_mid=round(float(alt[30]), 1), alt_min_window=round(float(alt[core].min()), 1),
                         sun_max_window=round(float(sun[core].max()), 1), moon_illum=round(illum, 2), moon_sep_deg=round(moon_sep, 1)))
df = pd.DataFrame(rows)
df.to_csv("out/followup_transits.csv", index=False)
good = df[df.full_window_ok_sun12]
out["observable_full_window_sun_lt_-12"] = {s: good[good.site == s].utc_mid.tolist() for s in SITES}
out["n_observable"] = {s: dict(sun12=int(df[(df.site == s) & df.full_window_ok_sun12].shape[0]),
                               sun18=int(df[(df.site == s) & df.full_window_ok_sun18].shape[0]),
                               plus2sigma=int(df[(df.site == s) & df.full_window_plus_2sigma_ok].shape[0])) for s in SITES}
late26 = df[df.utc_mid < "2027-01-01"]
out["late_2026"] = dict(n_transits=int(late26.epoch.nunique()),
                        observable_full=late26[late26.full_window_ok_sun12][["site", "utc_mid", "alt_min_window", "sun_max_window"]].to_dict("records"),
                        best_partial=late26.sort_values("frac_window_observable", ascending=False).head(6)[["site", "utc_mid", "frac_window_observable", "alt_mid"]].to_dict("records"))
any_good = good.groupby("epoch").site.apply(list)
out["epochs_observable_from_any_site"] = int(len(any_good))
out["first_observable"] = good.sort_values("bjd_tdb").head(5)[["site", "utc_mid", "sigma_min", "moon_illum", "moon_sep_deg"]].to_dict("records")
out["season_months_observable"] = sorted(set(u_[:7] for u_ in good.utc_mid))

# ---------------------------------------------------------------- 3. precision model
mags = {"g": 15.817, "r": 14.908, "i": 13.387, "z": 12.768}         # SkyMapper DR4 PSF (out/mags.json)
mira_mags = {"g": 15.12, "r": 13.95, "i": 11.14, "z": 9.93}          # one SkyMapper epoch; Mira varies (V 12.9-17.0)
lam = {"g": 0.48, "r": 0.62, "i": 0.76, "z": 0.91}; dlam = {"g": 0.15, "r": 0.14, "i": 0.15, "z": 0.13}
thr = {"g": 0.45, "r": 0.50, "i": 0.40, "z": 0.20}                   # atmosphere x optics x QE (typical deep-depletion CCD)
sky_mag = {"g": 21.8, "r": 20.9, "i": 20.0, "z": 19.0}                # dark-sky AB mag/arcsec^2; +1.5 mag brighter at full Moon
TELS = {"0.4m": dict(D=0.40, obsc=0.30, rn=10.0, scale=0.57, alt=2200.0, texp=120.0, overhead=10.0),
        "1.0m": dict(D=1.00, obsc=0.30, rn=8.0, scale=0.39, alt=2200.0, texp=60.0, overhead=20.0)}
def ph_rate_ab0(b):           # photons s^-1 m^-2 for AB mag 0 over the band
    return 3631e-26 / 6.626e-34 * dlam[b] / lam[b]
def prec(tel, b, seeing=2.0, airmass=1.3, moon_bright=False):
    T = TELS[tel]; area = np.pi * (T["D"] / 2) ** 2 * (1 - T["obsc"] ** 2)
    k_ext = {"g": 0.18, "r": 0.10, "i": 0.05, "z": 0.04}[b]
    star = ph_rate_ab0(b) * 10 ** (-0.4 * (mags[b] + k_ext * airmass)) * area * thr[b] * T["texp"]
    rap = max(1.5 * seeing, 2.5); npix = np.pi * rap ** 2 / T["scale"] ** 2
    sky = ph_rate_ab0(b) * 10 ** (-0.4 * (sky_mag[b] - (1.5 if moon_bright else 0))) * area * thr[b] * T["texp"] * np.pi * rap ** 2
    noise_ph = np.sqrt(star + sky + npix * T["rn"] ** 2) / star
    scint = 1.5 * 0.09 * (T["D"] * 100) ** (-2 / 3) * airmass ** 1.75 * np.exp(-T["alt"] / 8000) / np.sqrt(2 * T["texp"])
    tot = np.hypot(noise_ph, scint)
    cadence = T["texp"] + T["overhead"]
    # peak pixel of the M3V and of the Mira (Gaussian core, FWHM = seeing)
    sig_pix = seeing / 2.3548 / T["scale"]; peak_frac = 1 / (2 * np.pi * sig_pix ** 2)
    mira = ph_rate_ab0(b) * 10 ** (-0.4 * (mira_mags[b] + k_ext * airmass)) * area * thr[b] * T["texp"]
    return dict(star_e=star, sky_e=sky, aperture_radius_arcsec=rap, photon_noise_ppt=noise_ph * 1e3, scint_ppt=scint * 1e3,
                per_exposure_ppt=tot * 1e3, cadence_s=cadence, per_min_ppt=tot * 1e3 * np.sqrt(cadence / 60),
                peak_e_M3V=star * peak_frac, peak_e_Mira=mira * peak_frac,
                texp_for_M3V_peak_30ke=T["texp"] * 3e4 / (star * peak_frac))
def snr(tel, b, depth, red_ppt=1.0, **kw):
    p = prec(tel, b, **kw); n_in = C.T14 * 86400 * 0.85 / p["cadence_s"]   # ~85% of T14 at full depth-equivalent
    n_out = 2 * 3600 / p["cadence_s"]                                         # 1 h before + 1 h after
    sw = p["per_exposure_ppt"] * 1e-3 * np.sqrt(1 / n_in + 1 / n_out)
    n_blocks = C.T14 * 24 * 60 / 30
    sr = np.hypot(sw, red_ppt * 1e-3 / np.sqrt(n_blocks))
    return dict(depth=depth, sigma_white=sw, snr_white=depth / sw, sigma_red=sr, snr_red=depth / sr, n_in=n_in, n_out=n_out)
out["precision"] = {}
for tel in TELS:
    for b in mags:
        p = prec(tel, b); p2 = prec(tel, b, seeing=3.0, airmass=1.8, moon_bright=True)
        out["precision"][f"{tel}_{b}"] = dict(nominal={k: round(float(v), 3) for k, v in p.items()},
                                              poor_seeing3_airmass1p8_moon={k: round(float(v), 3) for k, v in p2.items()},
                                              snr={f"{d:.3f}": {k: round(float(v), 4) for k, v in snr(tel, b, d).items()} for d in (0.015, 0.018)},
                                              snr_poor={f"{d:.3f}": {k: round(float(v), 4) for k, v in snr(tel, b, d, red_ppt=2.0, seeing=3.0, airmass=1.8, moon_bright=True).items()} for d in (0.015, 0.018)})
# ---------------------------------------------------------------- 4. contamination (Moffat)
def moffat_frac(r_src, rap, fwhm, beta):
    a = fwhm / (2 * np.sqrt(2 ** (1 / beta) - 1))
    # fraction of the source's total flux inside an aperture of radius rap centred r_src away (numerical)
    x = np.linspace(-rap, rap, 201); X, Y = np.meshgrid(x, x); inside = X ** 2 + Y ** 2 <= rap ** 2
    I = (beta - 1) / (np.pi * a ** 2) * (1 + ((X + r_src) ** 2 + Y ** 2) / a ** 2) ** (-beta)
    return float(np.sum(I[inside]) * (x[1] - x[0]) ** 2)
cont = {}
for fwhm in (1.2, 2.0, 3.0):
    for beta in (2.5, 4.0):
        rap = max(1.5 * fwhm, 2.5)
        fm = moffat_frac(17.8, rap, fwhm, beta); fp = moffat_frac(11.7, rap, fwhm, beta); fs = moffat_frac(0.0, rap, fwhm, beta)
        rec = {}
        for b in mags:
            ratio_mira = 10 ** (-0.4 * (mira_mags[b] - mags[b])); ratio_mira_max = ratio_mira * 10 ** (0.4 * 2.0)   # Mira 2 mag brighter at max
            ratio_pair = 10 ** (-0.4 * ({"g": 16.019, "r": 15.063, "i": 13.363, "z": 12.621}[b] - mags[b]))
            rec[b] = dict(mira_over_target_in_aperture=fm * ratio_mira / fs, mira_at_max=fm * ratio_mira_max / fs, pair_over_target=fp * ratio_pair / fs)
        cont[f"fwhm{fwhm}_beta{beta}"] = dict(aperture_radius=rap, **{b: {k: float(f"{v:.2e}") for k, v in r.items()} for b, r in rec.items()})
out["contamination"] = cont
out["notes"] = ["Precision model: throughputs, sky brightness and read noise are generic assumptions (uncertain by ~30%); scintillation"
                " uses Young (1967) x 1.5 (Osborn et al. 2015) for 2200 m altitude.",
                "Red-noise floor 1 ppt per 30-min block (nominal) or 2 ppt (poor conditions) is an assumption typical of ground-based"
                " transit photometry of a T~13 star; real systematics near a bright saturated Mira could be worse.",
                "The Mira: V 12.9-17.0 (VSX, P 294 d); SkyMapper single epoch i 11.14, z 9.93, i.e. 2.2-2.8 mag brighter than the M3V in"
                " i/z; at maximum it may be ~2 mag brighter still. Its PSF core will saturate in i/z on a 1 m in 60 s (peak_e_Mira)."]
json.dump(out, open("out/followup.json", "w"), indent=1, default=str)
print(json.dumps({k: out[k] for k in ("csv_check", "n_observable", "epochs_observable_from_any_site", "first_observable", "season_months_observable")}, indent=1, default=str))
print("late 2026:", json.dumps(out["late_2026"], indent=1, default=str))
for k, v in out["precision"].items():
    n = v["nominal"]; s = v["snr"]["0.015"]; sp = v["snr_poor"]["0.015"]
    print(f"{k:8s} per-exp {n['per_exposure_ppt']:.2f} ppt (phot {n['photon_noise_ppt']:.2f}, scint {n['scint_ppt']:.2f}), per-min {n['per_min_ppt']:.2f}; "
          f"S/N(1.5%) white {s['snr_white']:.0f} red {s['snr_red']:.1f}; poor {sp['snr_red']:.1f}; peak M3V {n['peak_e_M3V']:.0f} e, Mira {n['peak_e_Mira']:.0f} e")
print(json.dumps(out["contamination"], indent=0)[:2500])
