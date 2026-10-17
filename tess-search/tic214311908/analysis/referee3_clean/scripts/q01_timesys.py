"""q01: independent check of the time system of the cleaned light curves.
For 40 cadences per sector: TESS spacecraft barycentric position from JPL Horizons (id -95, vectors wrt SSB),
Roemer delay for the M3V direction, compared with the TESScut TIMECORR (FFI reference direction).
The analyst's time_btjd_tdb - time_btjd_tesscut must equal Roemer(target) - TIMECORR (to ~0.1 s, Shapiro ignored).
Output: out/q01_timesys.json"""
import numpy as np, glob, json
from astropy.io import fits
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.jplhorizons import Horizons
import pandas as pd
OUT = "tic214311908/analysis/referee3_clean/out"
CL = "tic214311908/analysis/clean"
M3V = SkyCoord(260.969674 * u.deg, -51.323821 * u.deg)
n = np.array(M3V.cartesian.xyz.value)
C_AU_D = 173.144632674  # au per day
inv = pd.read_csv(f"{CL}/lc_M3V_INV_clean.csv", usecols=["time_btjd_tdb", "time_btjd_tesscut", "sector"])
res = {}
ADD = {}
for fn in sorted(glob.glob("/home/user/work/tesscut5/214311906/*.fits")):
    sec = int(fn.split("tess-s")[1][:4])
    with fits.open(fn) as h:
        t = np.array(h[1].data["TIME"], float); tc = np.array(h[1].data["TIMECORR"], float)
    ok = np.isfinite(t)
    idx = np.linspace(0, ok.sum() - 1, 40).astype(int)
    tt = t[ok][idx]; tcc = tc[ok][idx]
    tsc = tt - tcc + 2457000.0   # spacecraft TDB (JD)
    # Horizons with TDB epochs; vectors in ICRF wrt SSB (@0)
    vec = Horizons(id="-95", location="@0", epochs=list(tsc)).vectors(refplane="earth")
    r = np.column_stack([np.array(vec["x"]), np.array(vec["y"]), np.array(vec["z"])])  # au, ICRF equatorial
    roemer = r @ n / C_AU_D   # days
    add = (roemer - tcc) * 86400
    ADD[f"t_{sec}"] = tt; ADD[f"add_{sec}"] = add
    # analyst's shift at the same cadences
    g = inv[inv.sector == sec]
    m = np.isin(np.round(g.time_btjd_tesscut.values, 6), np.round(tt, 6))
    an = (g.time_btjd_tdb.values - g.time_btjd_tesscut.values)[m] * 86400
    an_t = g.time_btjd_tesscut.values[m]
    # match
    mine = np.interp(an_t, tt, add)
    res[sec] = dict(add_target_minus_timecorr_s_mean=float(add.mean()), add_min=float(add.min()), add_max=float(add.max()),
                    roemer_minus_timecorr_rms_vs_analyst_s=float(np.sqrt(np.mean((mine - an) ** 2))) if len(an) else None,
                    max_abs_diff_vs_analyst_s=float(np.max(np.abs(mine - an))) if len(an) else None, n_matched=int(len(an)))
    print(sec, res[sec], flush=True)
json.dump(res, open(f"{OUT}/q01_timesys.json", "w"), indent=1)
np.savez(f"{OUT}/q01_addtab.npz", **ADD)
