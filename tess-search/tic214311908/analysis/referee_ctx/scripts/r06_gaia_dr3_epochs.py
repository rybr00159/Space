"""R06: phases of the Mira's Gaia DR3 G epochs (same field-of-view crossings as the M3V, 18" away) on the ephemeris.
g_transit_time = BJD(TCB) - 2455197.5; TCB-TDB ~ 19 s here (neglected, << timing sigma). Also the full ctx crossing list
(gaia_fov_transits.csv): recompute dt_min and sigma independently from the ephemeris covariance.
Output: out/r06_gaia.json
"""
import json, numpy as np, pandas as pd
E = json.load(open("/home/user/Space/tess-search/vetting_2026-10-16/results/214311906/ephemeris_v2.json"))
P, T0, C = E["P"], E["T0"], np.array(E["cov"])
g = pd.read_csv("tic214311908/analysis/ctx/lc/gaia_epoch/mira_dr3_epoch_G.csv")
btjd = g.g_transit_time.dropna().values + 2455197.5 - 2457000.0
k = np.round((btjd - T0) / P); dt = (btjd - (T0 + k * P)) * 1440
sig = np.sqrt(C[0, 0] + k ** 2 * C[1, 1] + 2 * k * C[0, 1]) * 1440
out = dict(n_dr3=len(btjd), min_abs_dt_min=float(np.min(np.abs(dt))), sigma_at_min=float(sig[np.argmin(np.abs(dt))]),
           n_within_T14half_plus_3sigma=int(np.sum(np.abs(dt) < 48.9 + 3 * sig)))
f = pd.read_csv("tic214311908/analysis/ctx/out/gaia_fov_transits.csv")
kk = np.round((f.btjd.values - T0) / P); dd = (f.btjd.values - (T0 + kk * P)) * 1440
ss = np.sqrt(C[0, 0] + kk ** 2 * C[1, 1] + 2 * kk * C[0, 1]) * 1440
out["max_abs_diff_dt_vs_ctx_min"] = float(np.max(np.abs(dd - f.dt_min.values)))
out["max_abs_diff_sigma_vs_ctx_min"] = float(np.max(np.abs(ss - f.sigma_tc_min.values)))
for rel in ("DR4", "DR5"):
    m = f.release.values == rel
    if rel == "DR5": m = np.ones(len(f), bool)   # DR5 contains all crossings
    a = np.abs(dd[m]); i = np.argsort(a)[:4]
    out[rel] = dict(n=int(m.sum()), closest=[dict(utc=f.utc.values[m][j], dt_min=round(float(dd[m][j]), 1), sigma=round(float(ss[m][j]), 1)) for j in i])
print(json.dumps(out, indent=1))
json.dump(out, open("out/r06_gaia.json", "w"), indent=1)
