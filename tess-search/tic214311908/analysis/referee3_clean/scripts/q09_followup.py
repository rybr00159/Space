"""q09: follow-up numbers an observer would use.
(1) Predicted mid-transit times 2026-10 .. 2027-12 from the analyst's adopted ephemeris (clean/out/ephemeris_final.json,
    combined, covariance at E=0) and from the referee3 ephemerides (q05), with 1-sigma errors; difference to the
    phase-1 predictions in ctx/out/followup_transits.csv.
(2) BJD_TDB -> UTC at the observer (geocentre and CTIO): UTC = BJD_TDB - Roemer(Earth/site, star) - (TDB - UTC),
    iterated, astropy ephemeris 'builtin' (accurate to < 0.1 s for the Roemer term); compared with the CSV utc_mid.
(3) Expected ground-based depths and S/N: M3V magnitudes (ctx/out/mags.json), dilution by the 3.1" (G 17.48) star in
    a 3-5" aperture, and the 1-sigma photometric precision needed per 2-min bin for a 5-sigma (and 10-sigma) single
    transit detection with T14 = 1.70 h, T23/T14 = 0.72 (depth fraction of a box: ~0.86).
Output: out/q09_followup.json, out/q09_predictions.csv
"""
import json, numpy as np, pandas as pd
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
R3 = "tic214311908/analysis/referee3_clean"
TOP = "tic214311908/analysis"
star = SkyCoord(260.969674 * u.deg, -51.323821 * u.deg, pm_ra_cosdec=-113.67 * u.mas / u.yr, pm_dec=85.54 * u.mas / u.yr,
                distance=64.31 * u.pc, obstime=Time(2016.0, format="jyear"))
E = json.load(open(f"{TOP}/clean/out/ephemeris_final.json"))["combined"]
T0, P = E["T0_epoch0_BJD_TDB"], E["P_d"]; C = np.array(E["cov_T0epoch0_P"])
q5 = json.load(open(f"{R3}/out/q05_timing.json"))


def bjd_to_utc(bjd, loc=None):
    t = Time(bjd, format="jd", scale="tdb")
    for _ in range(3):
        tu = Time(t.utc.jd, format="jd", scale="utc", location=loc)
        ltt = tu.light_travel_time(star.apply_space_motion(new_obstime=tu) if False else star, kind="barycentric", location=loc)
        t = Time(bjd, format="jd", scale="tdb") - ltt
    return t.utc


def main():
    csv = pd.read_csv(f"{TOP}/ctx/out/followup_transits.csv")
    ep = np.unique(csv.epoch)
    rows = []
    ctio = EarthLocation.of_site("ctio") if False else EarthLocation(lat=-30.1690 * u.deg, lon=-70.8063 * u.deg, height=2207 * u.m)
    for e in ep:
        T = T0 + P * e; s = np.sqrt(C[0, 0] + e ** 2 * C[1, 1] + 2 * e * C[0, 1]) * 1440
        r = dict(epoch=int(e), bjd_tdb_clean=T, sig_min_clean=s)
        for k in ("inv", "ref", "comb_cor", "comb_ind"):
            if k in q5:
                Tk = q5[k]["T0_E0"] + 2457000 + q5[k]["P"] * e
                r[f"dT_{k}_min"] = (Tk - T) * 1440
        old = csv[csv.epoch == e].iloc[0]
        r["dT_phase1_csv_min"] = (old.bjd_tdb - T) * 1440; r["sig_phase1_csv_min"] = old.sigma_min
        tu = bjd_to_utc(T, ctio)
        r["utc_ctio"] = tu.iso[:19]
        r["csv_utc_minus_mine_min"] = (Time(old.utc_mid, scale="utc").jd - tu.jd) * 1440 - (old.bjd_tdb - T) * 1440
        rows.append(r)
    df = pd.DataFrame(rows); df.to_csv(f"{R3}/out/q09_predictions.csv", index=False)
    print(df[["epoch", "bjd_tdb_clean", "sig_min_clean", "dT_phase1_csv_min", "sig_phase1_csv_min", "csv_utc_minus_mine_min"] + [c for c in df.columns if c.startswith("dT_") and "phase1" not in c]].round(2).iloc[::8].to_string())
    m = json.load(open(f"{TOP}/ctx/out/mags.json"))["M3V"]
    # S/N requirement: depth d, flat-bottom fraction ~0.72 of T14, ingress/egress half weight
    T14 = 1.70 * 60; f_eff = 0.72 + 0.5 * 0.28        # effective full-depth duration fraction
    n_in = T14 * f_eff / 2.0                           # 2-min bins in transit (effective)
    out = dict(mags=dict(V=m["gaia_transform"]["V"], I_C=m["gaia_transform"]["I_C"], G=m["gaia"]["G"],
                         SMSS_i=m["smss4"]["iPSF"], SMSS_z=m["smss4"]["zPSF"]),
               predictions_summary=dict(sig_min_2026_10=float(df.sig_min_clean.min()), sig_min_2027_12=float(df.sig_min_clean.max()),
                                        max_abs_shift_vs_phase1_min=float(np.abs(df.dT_phase1_csv_min).max()),
                                        max_abs_csv_utc_error_min=float(np.abs(df.csv_utc_minus_mine_min).max())))
    for d in (15.2, 17.0):
        for dil in (1.0, 0.96):
            dd = d * dil
            out[f"precision_per_2min_for_5sigma_depth{d}_dil{dil}"] = float(dd / 5 * np.sqrt(n_in / 2))   # ppt per 2-min bin (white noise; OOT baseline as long as the transit)
            out[f"precision_per_2min_for_10sigma_depth{d}_dil{dil}"] = float(dd / 10 * np.sqrt(n_in / 2))
    out["n_in_2min_eff"] = n_in
    json.dump(out, open(f"{R3}/out/q09_followup.json", "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
