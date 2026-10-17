"""c00: pixel positions of the catalogued Gaia DR3 variables (and other bright sources) in each TESScut cube,
and the pixel masks used to exclude them from every field-wide estimator.
Output: out/varpos.json (per sector: pixel x, y of each variable; number of masked pixels)."""
import sys; sys.path.insert(0, "tic214311908/analysis/phot/scripts")
from photlib import *
V = json.load(open("tic214311908/analysis/search/out/gaia_variables.json"))["variables"]
out = {}
for s in SECTORS:
    C = Cutout(s); W = json.load(open(f"{ROOT}/out/wcs_s{s}.json")); dx0, dy0 = W["global_shift"]
    dt = C.epoch - 2016.0
    rows = []
    for v in V:
        ra = v["ra"] + v["pmra"] * dt / 3.6e6 / np.cos(np.radians(v["dec"])); de = v["dec"] + v["pmdec"] * dt / 3.6e6
        x, y = C.wcs.all_world2pix([ra], [de], 0); x = float(x[0]) + dx0; y = float(y[0]) + dy0
        rows.append(dict(source_id=v["source_id"], G=v["phot_g_mean_mag"], cls=v["best_class_name"], sep=v["sep"], x=x, y=y,
                         inside=bool(-1.5 < x < C.nx + 0.5 and -1.5 < y < C.ny + 0.5)))
    out[s] = dict(camera=C.camera, ccd=C.ccd, vars=rows)
    print(f"S{s} cam{C.camera} ccd{C.ccd}: " + "; ".join(f"{r['cls']} G{r['G']:.1f} ({r['x']:.1f},{r['y']:.1f})" for r in rows if r["inside"]))
save_json(out, "tic214311908/analysis/clean/out/varpos.json")
