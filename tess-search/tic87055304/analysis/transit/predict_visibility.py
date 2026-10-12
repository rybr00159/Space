"""Transit predictions 2026-11-01 .. 2028-12-31 with ground visibility flags for follow-up.

Input: predictions_2026-11_2028-12.csv (from ttv.py; linear ephemeris from injection-calibrated transit times and
the MCMC ephemeris; adopted sigma = the larger of the two). For each transit and for CTIO (Chile), SAAO
(Sutherland), Siding Spring (Australia) and La Silla, flags whether the full window
[Tc - T14/2 - 1 h - 3 sigma, Tc + T14/2 + 1 h + 3 sigma] has Sun altitude < -12 deg and target altitude > 30 deg
('full'), or whether at least the ingress or egress half is ('partial').
Writes predictions_visibility.csv.
"""
import os, json, csv
import numpy as np
import warnings
warnings.filterwarnings("ignore")
from astropy.utils import iers
iers.conf.auto_download = False
iers.conf.auto_max_age = None
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
import astropy.units as u

HERE = os.path.dirname(os.path.abspath(__file__))
rows = list(csv.DictReader(open(os.path.join(HERE, "predictions_2026-11_2028-12.csv"))))
T14 = json.load(open(os.path.join(HERE, "mcmc_free.json")))["params"]["T14_h"]["med"] / 24
target = SkyCoord(272.60647748 * u.deg, -43.99546646 * u.deg)
sites = {"CTIO": EarthLocation(lat=-30.169 * u.deg, lon=-70.806 * u.deg, height=2207 * u.m),
         "SAAO": EarthLocation(lat=-32.379 * u.deg, lon=20.811 * u.deg, height=1798 * u.m),
         "SSO": EarthLocation(lat=-31.273 * u.deg, lon=149.064 * u.deg, height=1165 * u.m),
         "LaSilla": EarthLocation(lat=-29.257 * u.deg, lon=-70.738 * u.deg, height=2400 * u.m)}
out = []
for r in rows:
    Tc = float(r["BJD_TDB"]); sig = float(r.get("sigma_adopted_min", r["sigma_min"])) / 1440
    half = T14 / 2 + 1 / 24 + 3 * sig
    tt = Time(np.linspace(Tc - half, Tc + half, 40), format="jd", scale="tdb")  # BJD ~ JD within minutes here (sufficient for flags)
    row = dict(epoch=r["epoch"], BJD_TDB=r["BJD_TDB"], sigma_min=r.get("sigma_adopted_min", r["sigma_min"]), UTC=r["UTC_ISO"])
    sun = get_sun(tt)
    for name, loc in sites.items():
        fr = AltAz(obstime=tt, location=loc)
        alt_t = target.transform_to(fr).alt.deg; alt_s = sun.transform_to(fr).alt.deg
        good = (alt_t > 30) & (alt_s < -12)
        n = len(tt)
        if good.all():
            row[name] = "full"
        elif good[: n // 2 + 3].all() or good[n // 2 - 3:].all():
            row[name] = "partial"
        else:
            row[name] = ""
    out.append(row)
with open(os.path.join(HERE, "predictions_visibility.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
nfull = {k: sum(1 for r in out if r[k] == "full") for k in sites}
npart = {k: sum(1 for r in out if r[k] == "partial") for k in sites}
print("full:", nfull, "partial:", npart)
for r in out:
    if any(r[k] == "full" for k in sites):
        print(r["UTC"], r["BJD_TDB"], r["sigma_min"][:4], [k for k in sites if r[k] == "full"])
