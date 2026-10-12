"""Predicted transit (eclipse) mid-times of TIC 87055315 for ground-based seeing-limited follow-up, late 2026 - 2027.

Ephemeris A (task brief): P = 4.136195 +- 0.000007 d, T0 = 1658.41325 +- 0.0027 BTJD;
  sigma(T_n) = sqrt(sigma_T0^2 + (n sigma_P)^2) (no covariance term; conservative for this reference epoch).
Ephemeris B (this work, mcmc/mcmc_pair_local.json): T0 at epoch E_ref = 390 (S66) and P with their covariance from the
  posterior samples.
Visibility: for each event between 2026-10-15 and 2027-12-31, the window [T_c - T14/2 - 1 h, T_c + T14/2 + 1 h] is
checked at 10-min steps for CTIO/LCO (Chile), SAAO (Sutherland) and Siding Spring (Australia): target altitude > 30 deg
and Sun altitude < -12 deg throughout. Times BTJD (BJD_TDB - 2457000) and UTC (BJD-UTC difference ~ 69 s ignored at the
minute level; barycentric correction < 8.3 min is applied via astropy light-travel time).
Output: results/ephemeris.json
"""
import os, json
import numpy as np
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
ROOT = "/home/user/work/tic87/neighbour"
P, eP, T0, eT0 = 4.136195, 0.000007, 1658.41325, 0.0027
mc = json.load(open(os.path.join(ROOT, "mcmc/mcmc_pair_local.json")))
ch = np.load(os.path.join(ROOT, "mcmc/mcmc_pair_local.npz"))
c = ch["chain"]; names = list(ch["names"])
T0b, Pb = c[:, names.index("T0")], c[:, names.index("P")]
E_ref = mc["E_ref"]
T14 = mc["params"]["T14_h"]["med"] / 24
star = SkyCoord(272.609023 * u.deg, -43.993633 * u.deg)
sites = {"CTIO": EarthLocation(lat=-30.1690 * u.deg, lon=-70.8063 * u.deg, height=2200 * u.m),
         "SAAO": EarthLocation(lat=-32.3794 * u.deg, lon=20.8107 * u.deg, height=1798 * u.m),
         "SSO": EarthLocation(lat=-31.2733 * u.deg, lon=149.0644 * u.deg, height=1165 * u.m)}
t_start = Time("2026-10-15T00:00:00", scale="utc"); t_end = Time("2027-12-31T00:00:00", scale="utc")


def bjd_to_utc(btjd, loc):
    """Approximate: subtract the barycentric light-travel time toward the star at that time."""
    t = Time(btjd + 2457000.0, format="jd", scale="tdb", location=loc)
    ltt = t.light_travel_time(star, kind="barycentric")
    return (t - ltt).utc


ev = []
n0 = int(np.ceil((t_start.tdb.jd - 2457000 - T0) / P)); n1 = int(np.floor((t_end.tdb.jd - 2457000 - T0) / P))
for n in range(n0, n1 + 1):
    tA = T0 + n * P; sA = float(np.hypot(eT0, n * eP))
    nb = n - E_ref
    tb = T0b + nb * Pb
    vis = {}
    for name, loc in sites.items():
        tc = bjd_to_utc(tA, loc)
        grid = tc + np.arange(-(T14 / 2 + 1 / 24), T14 / 2 + 1 / 24 + 1e-6, 10 / 1440) * u.day
        fr = AltAz(obstime=grid, location=loc)
        alt = star.transform_to(fr).alt.deg; salt = get_sun(grid).transform_to(fr).alt.deg
        vis[name] = bool(np.all(alt > 30) and np.all(salt < -12))
    if any(vis.values()):
        ev.append(dict(epoch=n, BTJD_A=round(tA, 5), sigma_A_min=round(sA * 1440, 1), BTJD_B=round(float(np.median(tb)), 5),
                       sigma_B_min=round(float(np.std(tb)) * 1440, 1),
                       UTC_mid=bjd_to_utc(tA, sites["CTIO"]).iso[:16], full_window_visible=vis))
# late 2026: no full window is observable (target low in evening twilight, solar conjunction ~ Dec 21);
# list events whose T1-T4 interval is at least partly observable (alt > 30, Sun < -12) from any site
part = []
for n in range(n0, int(np.floor((Time("2026-12-31").tdb.jd - 2457000 - T0) / P)) + 1):
    tA = T0 + n * P; best_site, best_frac = None, 0.0
    for name, loc in sites.items():
        tc = bjd_to_utc(tA, loc)
        grid = tc + np.arange(-T14 / 2, T14 / 2 + 1e-6, 5 / 1440) * u.day
        fr = AltAz(obstime=grid, location=loc)
        ok = (star.transform_to(fr).alt.deg > 30) & (get_sun(grid).transform_to(fr).alt.deg < -12)
        if ok.mean() > best_frac: best_site, best_frac = name, float(ok.mean())
    if best_frac > 0:
        part.append(dict(epoch=n, BTJD_A=round(tA, 5), sigma_A_min=round(float(np.hypot(eT0, n * eP)) * 1440, 1),
                         UTC_mid=bjd_to_utc(tA, sites["CTIO"]).iso[:16], best_site=best_site, frac_T14_observable=round(best_frac, 2)))
res = dict(late2026_partial=part, ephemeris_A=dict(P=P, eP=eP, T0=T0, eT0=eT0), ephemeris_B=dict(P=float(np.median(Pb)), eP=float(np.std(Pb)),
           T0_Eref=float(np.median(T0b)), eT0=float(np.std(T0b)), E_ref=E_ref, corr=float(np.corrcoef(T0b, Pb)[0, 1])),
           T14_h=T14 * 24, window_rule="target alt > 30 deg and Sun alt < -12 deg from T1 - 1 h to T4 + 1 h",
           n_visible=len(ev), events=ev)
json.dump(res, open(os.path.join(ROOT, "results/ephemeris.json"), "w"), indent=1)
for e in ev: print(e)
for e in part: print('partial', e)
