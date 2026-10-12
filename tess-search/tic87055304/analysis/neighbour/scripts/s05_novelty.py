"""Novelty checks for TIC 87055315 = Gaia DR3 6720986208895208960.

- ExoFOP-TESS JSON (fetched separately with curl into data/novelty/exofop_87055315.json): TOIs, CTOIs, follow-up.
- VizieR: catalogue existence is verified first (an unconstrained 1-row query must return the expected table), then
  cone searches at the neighbour's Gaia DR3 position (epoch 2016) with radius 3" (Kostov+25 tables, VSX, ASAS-SN,
  Gaia DR3 vari_*) and 21" (1 TESS pixel), and Green+23 (J/MNRAS/522/29/table2) constrained by TIC = 87055315.
- Gaia DR3 archive: vari_summary / vari_eclipsing_binary / vari_classifier_result rows by source_id (cached by
  s01_catalogs.py) are re-read here.
- SIMBAD: by identifier (Gaia DR3 / TIC / 2MASS) and a 10" cone.
Output: results/novelty.json
"""
import os, json, time, warnings
import numpy as np, pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
warnings.filterwarnings("ignore")
ROOT = "/home/user/work/tic87/neighbour"
OUT = os.path.join(ROOT, "data/novelty"); os.makedirs(OUT, exist_ok=True)
c16 = SkyCoord(272.609023 * u.deg, -43.993633 * u.deg)
res = {}


def retry(fn, n=3, wait=10):
    for i in range(n):
        try:
            return fn()
        except Exception as e:
            print("retry", i, repr(e)[:200], flush=True); time.sleep(wait)
    return None


# ExoFOP
try:
    d = json.load(open(os.path.join(OUT, "exofop_87055315.json")))
    res["exofop"] = dict(tois=d.get("tois"), ctois=d.get("ctois"), planet_parameters=d.get("planet_parameters"),
                         imaging=d.get("imaging"), spectroscopy=d.get("spectroscopy"), time_series=d.get("time_series"),
                         files=d.get("files"), stellar_companions=d.get("stellar_companions"),
                         keys=list(d.keys()), star_names=d["basic_info"]["star_names"],
                         notes="no notes field present in the JSON")
except Exception as e:
    res["exofop"] = f"failed: {e}"

from astroquery.vizier import Vizier
cats = {"Kostov25_table2": "J/ApJS/279/50/table2", "Kostov25_table3": "J/ApJS/279/50/table3",
        "Kostov25_table4": "J/ApJS/279/50/table4",
        "VSX": "B/vsx/vsx", "ASASSN_var": "II/366/catv2021", "GaiaDR3_vari_EB": "I/358/veb",
        "GaiaDR3_vari_classifier": "I/358/vclassre", "GaiaDR3_vari_summary": "I/358/vari",
        "Prsa22_TESS_EB": "J/ApJS/258/16/tess-ebs", "Green23_table2": "J/MNRAS/522/29/table2"}
V1 = Vizier(columns=["*"], row_limit=1, timeout=120)
V = Vizier(columns=["**", "_r"], row_limit=-1, timeout=180)
viz = {}
for name, cid in cats.items():
    ent = {}
    chk = retry(lambda: V1.query_constraints(catalog=cid))
    ent["catalogue_exists"] = bool(chk is not None and len(chk) > 0)
    ent["tables_seen"] = [t.meta.get("name", "") for t in chk] if chk is not None else []
    ent["columns"] = list(chk[0].colnames)[:25] if chk is not None and len(chk) else []
    for rad in (3, 21):
        r = retry(lambda: V.query_region(c16, radius=rad * u.arcsec, catalog=cid))
        if r is None:
            ent[f"cone{rad}"] = "query failed"
        elif len(r) == 0:
            ent[f"cone{rad}"] = 0
        else:
            rows = []
            for t in r:
                df = t.to_pandas(); df.insert(0, "_table", t.meta.get("name", ""))
                rows.append(df)
            df = pd.concat(rows, ignore_index=True, sort=False)
            df.to_csv(os.path.join(OUT, f"viz_{name}_cone{rad}.csv"), index=False)
            ent[f"cone{rad}"] = int(len(df))
            ent[f"cone{rad}_rows"] = json.loads(df.head(5).to_json(orient="records"))
    viz[name] = ent
    print(name, {k: v for k, v in ent.items() if not k.endswith("rows") and k != "columns"}, flush=True)
# Green+23 by TIC
for tid in (87055315, 87055304):
    r = retry(lambda: V.query_constraints(catalog="J/MNRAS/522/29/table2", TIC=str(tid)))
    viz[f"Green23_TIC{tid}"] = (json.loads(r[0].to_pandas().to_json(orient="records")) if r is not None and len(r)
                                else ("query failed" if r is None else "not in table2"))
    print("Green23", tid, viz[f"Green23_TIC{tid}"], flush=True)
res["vizier"] = viz

# Gaia DR3 archive tables (cached by s01)
ga = {}
for tab in ["vari_summary", "vari_eclipsing_binary", "vari_classifier_result", "nss_two_body_orbit",
            "nss_acceleration_astro", "nss_non_linear_spectro", "nss_vim_fl"]:
    f = os.path.join(ROOT, f"data/cat/gaiadr3_{tab}_nb.csv")
    try:
        ga[tab] = int(len(pd.read_csv(f)))
    except Exception:
        ga[tab] = "empty / not available"
res["gaia_dr3_archive_rows"] = ga

# SIMBAD
from astroquery.simbad import Simbad
sim = {}
for ident in ["Gaia DR3 6720986208895208960", "TIC 87055315", "2MASS J18102616-4359369"]:
    r = retry(lambda: Simbad.query_object(ident))
    sim[ident] = (json.loads(r.to_pandas().to_json(orient="records")) if r is not None and len(r) else "no SIMBAD object")
    print("SIMBAD", ident, sim[ident], flush=True)
r = retry(lambda: Simbad.query_region(c16, radius=10 * u.arcsec))
sim["cone10"] = json.loads(r.to_pandas().to_json(orient="records")) if r is not None and len(r) else "none"
print("SIMBAD cone", sim["cone10"], flush=True)
res["simbad"] = sim
json.dump(res, open(os.path.join(ROOT, "results/novelty.json"), "w"), indent=1, default=str)
print("DONE")
