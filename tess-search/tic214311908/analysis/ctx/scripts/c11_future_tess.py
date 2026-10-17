"""Which TESS sectors observe the Mira / TIC 214311908, past and planned (tess-point 0.9.5.post1, the latest on PyPI,
pointings through Sector 134), with sector dates from the HEASARC TESS sector page (fetched 2026-10-17, ext/heasarc_sector.html)
and the MIT orbit-times file (ext/TESS_orbit_times.csv). Also: the closest approach of any planned camera footprint
(how far outside the field the target is), via a grid of offsets.  Output: out/future_tess.json
"""
import re, html, json
import numpy as np
from tess_stars2px import tess_stars2px_function_entry as t2p, TESS_Spacecraft_Pointing_Data
import importlib.metadata as md
out = {"tess_point_version": md.version("tess-point")}
sp = TESS_Spacecraft_Pointing_Data(); out["max_sector_in_tess_point"] = int(max(s for s in sp.sectors if s < 1000))
txt = html.unescape(re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", open("ext/heasarc_sector.html").read()))).replace("\xa0", " ")
dates = {int(m.group(1)): (m.group(2), m.group(3)) for m in re.finditer(r"s(\d{2,3}) [A-Za-z ]+? [\d.]+, -?[\d.]+ -?[\d.]+ (\d{4} \w{3} \d{2}) (\d{4} \w{3} \d{2})", txt)}
import pandas as pd
ot = pd.read_csv("ext/TESS_orbit_times.csv", comment="#")
for sec, gg in ot.groupby("Sector"):
    if int(sec) not in dates:
        dates[int(sec)] = (str(gg["Start of Orbit"].min())[:10], str(gg["End of Orbit"].max())[:10])
out["n_sector_dates_parsed"] = len(dates)
for lab, ra, dec in (("Mira_TIC214311906", 260.977650, -51.323690), ("M3V_TIC214311908", 260.969674, -51.323821)):
    r = t2p([1], [ra], [dec])
    secs = sorted(set(int(s) for s in r[3]))
    out[lab] = [dict(sector=int(s), camera=int(c), ccd=int(k), col=round(float(x), 1), row=round(float(y), 1), dates=dates.get(int(s)))
                for s, c, k, x, y in zip(r[3], r[4], r[5], r[6], r[7])]
    out[lab + "_future_sectors_after_104"] = [s for s in secs if s > 104]
# how close do planned sectors 105-134 come? scan offsets in Dec/RA up to 15 deg
near = {}
for dd in np.arange(-15, 15.1, 1.0):
    for dr in np.arange(-20, 20.1, 1.0):
        r = t2p([1], [260.97765 + dr / np.cos(np.radians(-51.32))], [-51.32369 + dd])
        for s in r[3]:
            if 105 <= s <= 134:
                d = float(np.hypot(dd, dr)); near[int(s)] = min(near.get(int(s), 99), d)
out["closest_planned_sector_footprints_deg"] = dict(sorted(near.items(), key=lambda kv: kv[1])[:8])
out["sector_dates_105_134"] = {s: dates.get(s) for s in range(105, 135)}
json.dump(out, open("out/future_tess.json", "w"), indent=1, default=str)
print(json.dumps({k: v for k, v in out.items() if k != "sector_dates_105_134"}, indent=1, default=str))
