"""Predicted Gaia field-of-view transits of TIC 214311908 from the commanded scanning law (out/gaia_scanlaw_near.csv,
made by c07), and which of them fall inside predicted TESS transits.

For each 10-s sample: target offset from the FOV centre in the tangent plane, split into along-scan (AL, along the
scan position angle) and across-scan (AC). A crossing is where AL changes sign; time interpolated linearly (spin
60"/s). Observed if |AC| < 0.42 deg (all real DR3 epochs of the Mira had |AC| <= 0.405, see c10b). Times: bjd_fovX column
(barycentric TCB for the FOV direction) converted to TDB (TCB - TDB = L_B (JD_TCB - 2443144.5003725) * 86400 s).
DR4 = data up to 2020-01-20; DR5 = to 2025-01-15 (end of science observations).
Output: out/gaia_fov_transits.csv, out/gaia_in_transit.json
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd
import common as C

RA0, DEC0 = 260.969674, -51.323821
L_B = 1.550519768e-8
d = pd.read_csv("out/gaia_scanlaw_near.csv")
rows = []
for fov in (1, 2):
    ra_c, de_c, th = np.radians(d[f"ra_fov{fov}"].values), np.radians(d[f"dec_fov{fov}"].values), np.radians(d[f"scan_angle_fov{fov}"].values)
    a0, d0 = np.radians(RA0), np.radians(DEC0)
    cosc = np.sin(de_c) * np.sin(d0) + np.cos(de_c) * np.cos(d0) * np.cos(a0 - ra_c)
    xi = np.cos(d0) * np.sin(a0 - ra_c) / cosc
    eta = (np.cos(de_c) * np.sin(d0) - np.sin(de_c) * np.cos(d0) * np.cos(a0 - ra_c)) / cosc
    AL = np.degrees(xi * np.sin(th) + eta * np.cos(th)); AC = np.degrees(xi * np.cos(th) - eta * np.sin(th))
    near = np.degrees(np.arccos(np.clip(cosc, -1, 1))) < 1.5
    t = d[f"bjd_fov{fov}"].values; tj = d["jd_time"].values
    idx = np.where(near[:-1] & near[1:] & (np.sign(AL[:-1]) != np.sign(AL[1:])) & (np.abs(tj[1:] - tj[:-1]) < 30 / 86400))[0]
    for i in idx:
        f = AL[i] / (AL[i] - AL[i + 1])
        tc = t[i] + f * (t[i + 1] - t[i]); ac = AC[i] + f * (AC[i + 1] - AC[i])
        rows.append(dict(fov=fov, bjd_tcb=tc + 2455197.5, ac_deg=ac, scan_angle=np.degrees(th[i]) % 360))
x = pd.DataFrame(rows).sort_values("bjd_tcb")
x = x[np.abs(x.ac_deg) < 0.42].copy()   # 0.42: every Gaia DR3 epoch of the Mira matched a crossing with |AC| <= 0.405 (c10b)
x["bjd_tdb"] = x.bjd_tcb - L_B * (x.bjd_tcb - 2443144.5003725)
x["btjd"] = x.bjd_tdb - 2457000.0
from astropy.time import Time
x["utc"] = Time(x.bjd_tdb.values, format="jd", scale="tdb").utc.iso
x["release"] = np.where(x.bjd_tdb < 2458868.5, "DR4", np.where(x.bjd_tdb < 2460690.5, "DR5", "after_end"))
ep = np.round((x.btjd - C.T0) / C.P); x["epoch"] = ep.astype(int)
x["tc_pred"] = C.T0 + ep * C.P
x["sigma_tc_min"] = np.sqrt(C.T0_ERR**2 + (ep * C.P_ERR)**2 + 2 * ep * C.EPH["cov"][0][1]) * 1440
x["dt_min"] = (x.btjd - x.tc_pred) * 1440
x["frac_in_transit"] = C.template(x.btjd.values, 4.4 / 86400)        # G-band CCD integration ~4.4 s
half = C.T14 / 2 * 1440
# probability that the epoch is inside the transit, marginalised over the ephemeris uncertainty (Gaussian)
from scipy.stats import norm
x["p_in_transit"] = norm.cdf((half - x.dt_min) / x.sigma_tc_min) - norm.cdf((-half - x.dt_min) / x.sigma_tc_min)
x.to_csv("out/gaia_fov_transits.csv", index=False)
res = {}
for rel in ("DR4", "DR5"):
    s = x[x.release == rel]
    res[rel] = dict(n_fov_transits=int(len(s)), expected_in_transit=float(s.p_in_transit.sum()),
                    in_or_near_transit=s[s.p_in_transit > 0.01][["utc", "btjd", "fov", "ac_deg", "epoch", "dt_min", "sigma_tc_min", "p_in_transit", "frac_in_transit"]].round(4).to_dict("records"))
res["note"] = ("Commanded scanning law (cosmos.esa.int full-mission file); realised attitude may differ by <~30 arcsec; outages and "
               "rejected transits are not modelled: in the DR3 window 14 of 48 predicted |AC|<0.35 crossings of the Mira have no DR3 epoch. |AC| < 0.42 deg.")
res["duty_cycle"] = float(C.T14 / C.P)
json.dump(res, open("out/gaia_in_transit.json", "w"), indent=1, default=float)
print("FOV transits:", len(x), "DR4:", (x.release == "DR4").sum(), "DR5:", (x.release == "DR5").sum())
for rel in ("DR4", "DR5"):
    print(rel, "expected in-transit epochs:", round(res[rel]["expected_in_transit"], 2))
    for r in res[rel]["in_or_near_transit"]:
        print("   ", r)
