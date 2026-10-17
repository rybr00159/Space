"""n07: details of the X-ray entries near the M3V (eRASS1 main, Freund+24 coronal identifications, Salvato+25 eRASS1
counterparts, and the eRASS:3 catalogue J/A+A/712/A171, which the earlier ctx/ check did not have), and explicit checks of
M-dwarf / nearby-star catalogues that a 4" all-VizieR cone might miss. Read-only.
Outputs: out/xray_mdwarf.json, raw/vizier_xray_<cat>.csv
Usage: python -I n07_xray_mdwarf.py <novelty_dir>
"""
import sys, os, re, json, time, datetime
import numpy as np, pandas as pd, requests
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord

D = sys.argv[1]
RAW, OUT = os.path.join(D, "raw"), os.path.join(D, "out")
NOW = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
C16 = SkyCoord(260.96967369918, -51.3238212118, unit="deg")
C00 = SkyCoord(260.970484934795, -51.3242024657333, unit="deg")
S = {"queried_utc": NOW, "readme": {}, "xray": {}, "mdwarf_catalogues": {}}
for cat in ("J/A+A/712/A171", "J/A+A/704/A344"):
    r = requests.get(f"https://cdsarc.cds.unistra.fr/ftp/{cat}/ReadMe", timeout=60)
    txt = r.text
    S["readme"][cat] = txt.split("\n================================================================================", 1)[0][:1500]
    # keep the abstract and the table list (small)
    open(os.path.join(RAW, f"readme_{re.sub('[^A-Za-z0-9]+', '_', cat)}.txt"), "w").write(txt[:12000])
    print("====", cat); print(S["readme"][cat][:1200])
for cat in ("J/A+A/682/A34/erass1-m", "J/A+A/684/A121/like_ctp", "J/A+A/704/A344/dr1ge3", "J/A+A/712/A171/erass3-m",
            "J/A+A/712/A171/dr2mg", "J/A+A/712/A171/dr2mc"):
    res = Vizier(columns=["**", "_r"], row_limit=20).query_region(C16, radius=30 * u.arcsec, catalog=cat)
    if len(res):
        t = res[0].to_pandas()
        t.to_csv(os.path.join(RAW, f"vizier_xray_{re.sub('[^A-Za-z0-9]+', '_', cat)}.csv"), index=False)
        S["xray"][cat] = [{k: str(x) for k, x in rr.items() if str(x) not in ("", "nan", "--", "None")} for rr in t.to_dict("records")]
    else:
        S["xray"][cat] = []
    print("====", cat, len(S["xray"][cat]))
    for rr in S["xray"][cat]:
        print({k: v for k, v in rr.items() if re.search(r"_r|Name|ERO|IAU|ID|Id|Gaia|Ctp|p_|p-|prob|Prob|Flux|FX|flux|Sep|sep|ePos|pos|RA|DE|Ext|LIKE|Like|DET|Det|class|Class|MJD|cts|Rate", k)})
# explicit M-dwarf / nearby-star catalogues (cone 10" around both epochs)
mcats = {"J/AJ/155/180": "TESS Cool Dwarf Catalog (Muirhead+ 2018)", "J/AJ/158/138": "TIC 8 / CTL (Stassun+ 2019)",
         "J/AJ/142/138": "bright M dwarfs (Lepine & Gaidos 2011)", "J/MNRAS/435/2161": "PPMXL M dwarfs (Frith+ 2013)",
         "J/A+A/650/A201": "10-pc sample (Reyle+ 2021)", "J/AJ/161/234": "K and M dwarfs (Medan+ 2021)",
         "J/A+A/649/A6": "GCNS (Gaia Collab., Smart+ 2021)", "J/ApJS/267/8": "XP metallicities (Andrae+ 2023)",
         "J/AJ/165/267": "stars within 120 pc (Hardegree-Ullman+ 2023)", "J/ApJ/891/58": "TIC exposure times (Suissa+ 2020)",
         "J/A+A/657/A7": "Gaia EDR3 PMa and CPM companions (Kervella+ 2022)", "J/MNRAS/506/2269": "wide binaries (El-Badry+ 2021)",
         "J/ApJS/258/16": "TESS EBs (Prsa+ 2022)", "I/358/vlpv": "Gaia DR3 LPVs (Lebzelter+ 2023)",
         "I/358/vclassre": "Gaia DR3 variability classification (Rimoldini+ 2023)", "J/A+A/677/A137": "Gaia DR3 photometric dispersions (Maiz Apellaniz+ 2023)", "J/MNRAS/486/1907": "ASAS-SN variables (Jayasinghe+ 2019)",
         "J/MNRAS/491/13": "ASAS-SN variables V (Jayasinghe+ 2020)", "II/366": "ASAS-SN variables (Christy+ 2023)"}
for cat, lab in mcats.items():
    hits = []
    for c in (C16, C00):
        res = None
        for k in range(4):
            try:
                res = Vizier(columns=["**", "_r"], row_limit=10, timeout=300).query_region(c, radius=10 * u.arcsec, catalog=cat); break
            except Exception as e:
                err = e; time.sleep(30)
        if res is None:
            hits.append(f"error {str(err)[:120]}"); continue
        for t in res:
            hits.append(dict(table=t.meta.get("name"), n=len(t), min_r=float(np.nanmin(np.array(t["_r"], dtype=float)))))
    # dedupe
    seen, h2 = set(), []
    for h in hits:
        key = h["table"] if isinstance(h, dict) else h
        if key not in seen:
            seen.add(key); h2.append(h)
    S["mdwarf_catalogues"][cat] = {"label": lab, "hits_within_10arcsec": h2}
    print(cat, lab, h2)
    time.sleep(0.5)
json.dump(S, open(os.path.join(OUT, "xray_mdwarf.json"), "w"), indent=1)
