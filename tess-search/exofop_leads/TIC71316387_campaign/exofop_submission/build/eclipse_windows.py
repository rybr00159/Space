"""Upcoming primary eclipses of the TIC 71316387 EB observable from three southern sites.
Ephemeris (BJD_TDB): T0 = 2460165.670416 +/- 0.00018, P = 1.30367228 +/- 0.00000019 d.
A night event: Sun below -12 deg and target above 30 deg altitude from mid-eclipse - 1.3 h to + 1.3 h
(the eclipse lasts 1.61 h; this leaves ~30 min of baseline either side)."""
import json
import numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
import astropy.units as u

T0, ST0 = 2460165.670416, 0.00018
P, SP = 1.30367228, 0.00000019
c = SkyCoord(330.887948 * u.deg, -42.835511 * u.deg)
SITES = {'CTIO': EarthLocation(lat=-30.169 * u.deg, lon=-70.806 * u.deg, height=2207 * u.m),
         'SAAO': EarthLocation(lat=-32.380 * u.deg, lon=20.811 * u.deg, height=1798 * u.m),
         'SSO': EarthLocation(lat=-31.273 * u.deg, lon=149.062 * u.deg, height=1165 * u.m)}
start, stop = Time('2026-10-17T00:00:00', scale='utc'), Time('2027-01-20T00:00:00', scale='utc')


def geo_utc(bjd, loc):
    """BJD_TDB -> UTC at the site (iterate the barycentric light-travel time)."""
    t = Time(bjd, format='jd', scale='tdb', location=loc)
    for _ in range(3):
        t = Time(bjd, format='jd', scale='tdb', location=loc) - t.light_travel_time(c, kind='barycentric')
        t = Time(t.jd, format='jd', scale='tdb', location=loc)
    return t.utc


rows = []
n0 = int(np.ceil((start.tdb.jd - T0) / P)); n1 = int(np.floor((stop.tdb.jd - T0) / P))
for n in range(n0, n1 + 1):
    bjd = T0 + n * P
    sig_min = np.hypot(ST0, n * SP) * 1440
    ok = {}
    for name, loc in SITES.items():
        tm = geo_utc(bjd, loc)
        grid = tm + np.linspace(-1.3, 1.3, 9) * u.h
        fr = AltAz(obstime=grid, location=loc)
        alt = c.transform_to(fr).alt.deg
        sun = get_sun(grid).transform_to(fr).alt.deg
        if alt.min() > 30 and sun.max() < -12:
            ok[name] = float(alt.max())
    if ok:
        tm = geo_utc(bjd, SITES['CTIO'])
        rows.append(dict(epoch=n, bjd_tdb=round(bjd, 5), utc_mid=(tm + 30 * u.s).isot[:16], utc_mid_s=tm.isot[:19], sigma_min=round(float(sig_min), 2), sites=sorted(ok)))   # rounded to the nearest minute
json.dump(rows, open('/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/obsnote/report/eclipse_windows.json', 'w'), indent=1)
for r in rows:
    print(r)
print(len(rows))
