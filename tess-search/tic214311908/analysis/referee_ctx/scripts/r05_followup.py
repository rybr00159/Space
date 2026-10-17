"""R05: predicted transit times and ground-based observability, recomputed independently.

1. Ephemeris and covariance from ephemeris_v2.json; epochs from 2026-10-17 to 2027-12-31; sigma from the full
   covariance (already chi2/dof-inflated).
2. predicted_transits_2027.csv: compare bjd_tdb AND utc_mid. UTC from BJD_TDB by removing the barycentric light-travel
   time towards the M3V (astropy, geocentre; site differences < 0.03 s) and TDB->UTC.
3. Observability: CTIO, La Silla, SAAO, Siding Spring; window = mid +- (T14/2 + 1 h) = +-1.815 h, also +- 2 sigma;
   altitude > 30 deg and Sun < -12 deg at every 3-min step.
Output: out/r05_followup.json, out/r05_observable.csv
"""
import json
import numpy as np, pandas as pd
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
import warnings; warnings.filterwarnings("ignore")

E = json.load(open("/home/user/Space/tess-search/vetting_2026-10-16/results/214311906/ephemeris_v2.json"))
P, T0 = E["P"], E["T0"]; C = np.array(E["cov"])
# M3V position at 2027.0 (Gaia DR3 + PM)
ra0, de0, pmra, pmde = 260.96967369918, -51.32382121180, -113.671152, 85.543550
dt = 11.0
tgt = SkyCoord((ra0 + pmra * dt / 3.6e6 / np.cos(np.radians(de0))) * u.deg, (de0 + pmde * dt / 3.6e6) * u.deg)
SITES = {"CTIO": (-30.1690, -70.8063, 2207), "La Silla": (-29.2567, -70.7300, 2400),
         "SAAO": (-32.3794, 20.8107, 1798), "Siding Spring": (-31.2733, 149.0644, 1165)}
k0 = int(np.ceil((Time("2026-10-17").jd - 2457000 - T0) / P)); k1 = int((Time("2027-12-31").jd - 2457000 - T0) / P)
ks = np.arange(k0, k1 + 1)
bjd = T0 + ks * P + 2457000.0
sig = np.sqrt(C[0, 0] + ks ** 2 * C[1, 1] + 2 * ks * C[0, 1]) * 1440


def bjd_to_utc(b):
    t = Time(b, format="jd", scale="tdb")
    for _ in range(2):
        ltt = t.light_travel_time(tgt, kind="barycentric", location=EarthLocation.from_geocentric(0, 0, 0, unit=u.m))
        t = Time(b - ltt.to(u.day).value, format="jd", scale="tdb")
    return t.utc


utc = bjd_to_utc(bjd)
out = {}
csv = pd.read_csv("/home/user/Space/tess-search/vetting_2026-10-16/results/214311906/predicted_transits_2027.csv")
mine = pd.DataFrame(dict(epoch=ks, bjd=bjd, utc=[x.iso for x in utc], sig=sig))
m = csv.merge(mine, on="epoch")
dt_bjd = (m.bjd_tdb - m.bjd) * 1440
dt_utc = np.array([(Time(a, scale="utc") - Time(b, scale="utc")).to(u.min).value for a, b in zip(m.utc_mid, m.utc)])
naive = np.array([(Time(a, scale="utc") - Time(b, format="jd", scale="tdb").utc).to(u.min).value for a, b in zip(m.utc_mid, m.bjd_tdb)])
out["csv_check"] = dict(n=len(m), max_abs_dbjd_min=float(np.max(np.abs(dt_bjd))), max_abs_dsigma_min=float(np.max(np.abs(m.sigma_min - m.sig))),
                        utc_mid_minus_true_utc_min=[float(np.min(dt_utc)), float(np.max(dt_utc))],
                        utc_mid_minus_naive_tdb_as_utc_min=[float(np.min(naive)), float(np.max(naive))],
                        note="CSV utc_mid equals BJD_TDB converted to UTC without removing the barycentric light-travel time "
                             "if utc_mid_minus_naive is ~0 (+-0.5 min rounding)")
print(out["csv_check"])

rows = []
grid = np.linspace(-1.815, 1.815, 73) / 24.0
for name, (lat, lon, h) in SITES.items():
    loc = EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=h * u.m)
    for i, k in enumerate(ks):
        tt = utc[i] + grid * u.day
        fr = AltAz(obstime=tt, location=loc)
        alt = tgt.transform_to(fr).alt.deg; sun = get_sun(tt).transform_to(fr).alt.deg
        okw = np.all((alt > 30) & (sun < -12))
        g2 = np.linspace(-1.815 - 2 * sig[i] / 60, 1.815 + 2 * sig[i] / 60, 73) / 24.0
        tt2 = utc[i] + g2 * u.day; fr2 = AltAz(obstime=tt2, location=loc)
        ok2 = np.all((tgt.transform_to(fr2).alt.deg > 30) & (get_sun(tt2).transform_to(fr2).alt.deg < -12))
        frac = np.mean((alt > 30) & (sun < -12))
        rows.append(dict(site=name, epoch=int(k), utc_mid=utc[i].iso[:16], sigma_min=round(sig[i], 1), full_ok=bool(okw), full_2sig_ok=bool(ok2), frac=round(frac, 2)))
df = pd.DataFrame(rows); df.to_csv("out/r05_observable.csv", index=False)
good = df[df.full_ok]
out["n_site_epochs_full"] = int(len(good)); out["n_epochs_full_any_site"] = int(good.epoch.nunique())
out["n_epochs_full_2sigma_any_site"] = int(df[df.full_2sig_ok].epoch.nunique())
out["per_site"] = good.groupby("site").size().to_dict()
out["months"] = sorted(set(x[:7] for x in good.utc_mid))
out["first"] = good.sort_values("utc_mid").head(3).to_dict("records")
late = df[df.utc_mid < "2027-01-01"]
out["late2026_best"] = late.sort_values("frac", ascending=False).head(4).to_dict("records")
print(json.dumps({k: v for k, v in out.items() if k != "csv_check"}, indent=0, default=str))
json.dump(out, open("out/r05_followup.json", "w"), indent=1, default=str)
