"""Cross-check of c11 (tess-point's internal pointing table): parse the spacecraft pointings (RA, Dec, roll) of
Sectors 97-134 from the HEASARC TESS sector page (ext/heasarc_sector.html, fetched 2026-10-17), compare them with
tess-point's table, and re-run tess-point with the HEASARC pointings as a sector-override file. Also gives the angular
distance from the target to every camera centre (cameras at 12/36 deg from the boresight along the roll axis).
Output: out/future_tess_crosscheck.json
"""
import re, html, json, os, tempfile
import numpy as np
from tess_stars2px import tess_stars2px_function_entry as t2p, TESS_Spacecraft_Pointing_Data
from astropy.coordinates import SkyCoord
import astropy.units as u

txt = html.unescape(re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", open("ext/heasarc_sector.html").read()))).replace("\xa0", " ")
pat = re.compile(r"s(\d{3}) ([A-Za-z ]+?) (-?[\d.]+), ?(-?[\d.]+) (-?[\d.]+) (\d{4} \w{3} \d{2}) (\d{4} \w{3} \d{2})")
pt = {}
for m in pat.finditer(txt):
    s = int(m.group(1))
    if 97 <= s <= 134 and s not in pt:
        pt[s] = dict(field=m.group(2).strip(), ra=float(m.group(3)), dec=float(m.group(4)), roll=float(m.group(5)), start=m.group(6), stop=m.group(7))
sp = TESS_Spacecraft_Pointing_Data()
comp = {}
for s, p in sorted(pt.items()):
    k = np.where(sp.sectors == s)[0]
    if len(k):
        k = k[0]; comp[s] = dict(heasarc=[p["ra"], p["dec"], p["roll"]], tess_point=[float(sp.ras[k]), float(sp.decs[k]), float(sp.rolls[k])])
fn = os.path.join(tempfile.mkdtemp(), "override.txt")
with open(fn, "w") as f:
    for s, p in sorted(pt.items()):
        f.write(f"{s} {p['ra']} {p['dec']} {p['roll']}\n")
res = {}
tgt = SkyCoord(260.969674, -51.323821, unit="deg")
for lab, ra, dec in (("Mira", 260.977650, -51.323690), ("M3V", 260.969674, -51.323821)):
    r = t2p([1], [ra], [dec], sectorOverrideFile=fn)
    res[lab] = [dict(sector=int(s), camera=int(c), ccd=int(k), col=float(x), row=float(y)) for s, c, k, x, y in zip(r[3], r[4], r[5], r[6], r[7])]
# distance to the camera centres for sectors 105-134 (camera centres: offset along the roll axis from the boresight)
spo = TESS_Spacecraft_Pointing_Data(sectorOverrideFile=fn)
dist = {}
for i, s in enumerate(spo.sectors):
    if s < 105:
        continue
    best = 999
    for cam in range(4):
        # camera centres in ecliptic coords were computed by tess-point (camRa/camDec arrays per sector)
        c = SkyCoord(spo.camRa[cam, i], spo.camDec[cam, i], unit="deg")
        best = min(best, c.separation(tgt).deg)
    dist[int(s)] = round(float(best), 1)
out = dict(n_pointings_parsed=len(pt), pointing_comparison_max_diff_deg=max(abs(a - b) for v in comp.values() for a, b in zip(v["heasarc"], v["tess_point"])),
           pointings_105_134={s: pt[s] for s in sorted(pt) if s >= 105}, on_silicon_with_heasarc_pointings=res,
           closest_camera_centre_deg_105_134=dist,
           note="A camera's field is 24x24 deg, so a target more than ~17 deg from every camera centre cannot be on silicon.")
json.dump(out, open("out/future_tess_crosscheck.json", "w"), indent=1)
print("parsed", len(pt), "pointings; max |HEASARC - tess-point| =", out["pointing_comparison_max_diff_deg"])
print("on silicon (HEASARC pointings):", {k: [x["sector"] for x in v] for k, v in res.items()})
print("closest camera centre (deg):", dist)
