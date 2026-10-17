"""n05: every VizieR catalogue with an entry for the M3V (cone 4" around both the J2000 and the J2016 position, plus a
4" cone at the Gaia DR2 epoch 2015.5 is covered by these), compared with the earlier ctx/ list; references (ReadMe
header) of the catalogues that are more than all-sky surveys; X-ray/UV/radio catalogues within 30"; the rows of the
catalogues worth citing. Read-only queries.
Outputs: out/vizier_catalogues_m3v.csv, out/vizier_status.json, raw/vizier_rows_<cat>.csv
Usage: python -I n05_vizier.py <novelty_dir> <ctx_vizier_hits_csv>
"""
import sys, os, re, json, time, datetime
import numpy as np, pandas as pd, requests
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord

D, CTXH = sys.argv[1], sys.argv[2]
RAW, OUT = os.path.join(D, "raw"), os.path.join(D, "out")
NOW = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
C00 = SkyCoord(260.970484934795, -51.3242024657333, unit="deg")
C16 = SkyCoord(260.96967369918, -51.3238212118, unit="deg")
v = Vizier(columns=["**", "_r"], row_limit=20, timeout=300)
rows = []
for lab, c in (("J2000", C00), ("J2016", C16)):
    for k in range(3):
        try:
            res = v.query_region(c, radius=4 * u.arcsec, catalog=None); break
        except Exception as e:
            print("retry", e); time.sleep(20)
    for t in res:
        rows.append(dict(cone=lab, catalog=t.meta.get("name", ""), description=t.meta.get("description", ""), n=len(t),
                         min_r=float(np.nanmin(np.array(t["_r"], dtype=float))) if "_r" in t.colnames else np.nan))
df = pd.DataFrame(rows).sort_values("min_r").drop_duplicates("catalog").sort_values("catalog")
ctx = pd.read_csv(CTXH)
ctx_m3v = set(ctx[ctx.star == "M3V"].catalog)
df["in_ctx_list"] = df.catalog.isin(ctx_m3v)
df.to_csv(os.path.join(OUT, "vizier_catalogues_m3v.csv"), index=False)
print(len(df), "catalogues;", (~df.in_ctx_list).sum(), "not in the ctx/ list:")
print(df[~df.in_ctx_list][["catalog", "description", "min_r"]].to_string())
S = {"queried_utc": NOW, "n_catalogues": int(len(df)), "new_vs_ctx": df[~df.in_ctx_list].catalog.tolist(),
     "ctx_not_found_now": sorted(ctx_m3v - set(df.catalog))}

# references for the non-survey catalogues (ReadMe header)
survey = re.compile(r"^(I/|II/|IV/|V/|VI/|B/)")
refs = {}
for cat in sorted(set(c.rsplit("/", 1)[0] for c in df.catalog if not survey.match(c)) |
                  {"J/A+A/649/A6", "J/MNRAS/506/2269", "J/A+A/657/A7", "J/A+A/682/A34", "J/A+A/684/A121"}):
    try:
        r = requests.get(f"https://cdsarc.cds.unistra.fr/ftp/{cat}/ReadMe", timeout=60)
        head = r.text.split("\n================================================================================", 1)[0]
        head = "\n".join(head.splitlines()[:12])
        refs[cat] = head
    except Exception as e:
        refs[cat] = f"error {e}"
    time.sleep(0.5)
S["references"] = refs
for k, h in refs.items():
    print("----", k); print(h[:600])

# rows of the catalogues worth citing
cite = ["J/A+A/649/A6/table1c", "J/MNRAS/506/2269/catalog", "J/A+A/657/A7/tablea3", "J/A+A/682/A34/erass1-m",
        "J/A+A/684/A121/like_ctp", "J/AJ/165/267/table1", "J/ApJ/984/58/table1", "J/ApJ/891/58/table6",
        "J/A+A/704/A344/dr1ge3", "J/A+A/704/A237/80pcbpok", "J/A+A/692/A150/allskytg", "J/A+A/653/A98/aspic1_1"]
cite += [c for c in df.catalog if not survey.match(c) and c not in cite]
S["rows"] = {}
for cat in cite:
    try:
        res = Vizier(columns=["**", "_r"], row_limit=5).query_region(C16 if "A34" not in cat and "A121" not in cat else C16,
                                                                      radius=(10 if ("A34" in cat or "A121" in cat or "A344" in cat) else 4) * u.arcsec,
                                                                      catalog=cat)
        if len(res):
            t = res[0].to_pandas()
            t.to_csv(os.path.join(RAW, f"vizier_rows_{re.sub('[^A-Za-z0-9]+', '_', cat)}.csv"), index=False)
            S["rows"][cat] = [{k: str(x) for k, x in rr.items() if str(x) not in ("", "nan", "--", "None")} for rr in t.to_dict("records")]
        else:
            S["rows"][cat] = []
    except Exception as e:
        S["rows"][cat] = f"error {e}"
    print("rows", cat, str(S["rows"][cat])[:500])

# X-ray / UV / radio catalogues within 30" (X-ray positions are uncertain by several arcsec)
xr = []
for k in range(3):
    try:
        res = Vizier(columns=["**", "_r"], row_limit=10, timeout=300).query_region(C16, radius=30 * u.arcsec, catalog=None); break
    except Exception as e:
        print("retry", e); time.sleep(20)
kw = re.compile(r"X-ray|xray|ROSAT|RASS|XMM|Chandra|Swift|eROSITA|eRASS|GALEX|ultraviolet|\bUV\b|radio|NVSS|SUMSS|RACS|VLASS|LoTSS|flare", re.I)
for t in res:
    name, desc = t.meta.get("name", ""), t.meta.get("description", "")
    if kw.search(desc) or kw.search(name) or name.startswith(("IX/", "J/A+A/682/A34", "J/A+A/684", "J/A+A/704/A344", "II/335", "II/312", "VIII/")):
        xr.append(dict(catalog=name, description=desc, n=len(t), min_r=float(np.nanmin(np.array(t["_r"], dtype=float)))))
S["xray_uv_radio_within_30arcsec"] = xr
print("X-ray/UV/radio within 30\":", json.dumps(xr, indent=0)[:3000])
json.dump(S, open(os.path.join(OUT, "vizier_status.json"), "w"), indent=1)
