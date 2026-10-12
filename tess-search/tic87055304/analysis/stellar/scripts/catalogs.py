"""Collect every independent catalogue estimate of Teff, logg, [Fe/H], R, M, distance, extinction for
TIC 87055304 = Gaia DR3 6720986204578421888 (Gaia DR2 id looked up via gaiadr3.dr2_neighbourhood).

Output: out/catalog_raw/<cat>.csv (full matched rows), out/catalog_compilation.csv (harmonised table).
Usage: python catalogs.py
"""
import os, sys, json, time, warnings
import numpy as np, pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
warnings.filterwarnings("ignore")

OUT = "/home/user/work/tic87/stellar/out"
RAW = os.path.join(OUT, "catalog_raw"); os.makedirs(RAW, exist_ok=True)
DR3 = 6720986204578421888
TIC = 87055304
c0 = SkyCoord(272.606557 * u.deg, -43.995371 * u.deg)   # TIC (J2000, epoch 2000)

def retry(fn, n=4, wait=10):
    for i in range(n):
        try:
            return fn()
        except Exception as e:
            print("  retry", i, repr(e)[:200]); time.sleep(wait)
    return None

# ---------------------------------------------------------------- TIC v8.2 (MAST)
from astroquery.mast import Catalogs
tic = retry(lambda: Catalogs.query_criteria(catalog="Tic", ID=TIC))
if tic is not None:
    tic.to_pandas().to_csv(os.path.join(RAW, "tic82_mast.csv"), index=False)

# ---------------------------------------------------------------- Gaia DR2 via Gaia archive
from astroquery.gaia import Gaia
Gaia.ROW_LIMIT = -1
def gq(q):
    return Gaia.launch_job(q).get_results().to_pandas()
nb = retry(lambda: gq(f"select * from gaiadr3.dr2_neighbourhood where dr3_source_id={DR3}"))
if nb is not None:
    nb.to_csv(os.path.join(RAW, "dr2_neighbourhood.csv"), index=False); print(nb)
    dr2id = int(nb.sort_values("angular_distance").iloc[0]["dr2_source_id"])
    dr2 = retry(lambda: gq(f"select * from gaiadr2.gaia_source where source_id={dr2id}"))
    if dr2 is not None: dr2.to_csv(os.path.join(RAW, "gaiadr2_source.csv"), index=False)
    ap2 = retry(lambda: gq(f"select * from gaiadr2.astrophysical_parameters where source_id={dr2id}"))
    if ap2 is not None and len(ap2): ap2.to_csv(os.path.join(RAW, "gaiadr2_ap.csv"), index=False)
else:
    dr2id = None
# DR3 extras: Bailer-Jones distances are on VizieR; also gaiadr3 astrophysical_parameters_supp (other GSP-Phot libraries)
sup = retry(lambda: gq(f"select * from gaiadr3.astrophysical_parameters_supp where source_id={DR3}"))
if sup is not None and len(sup): sup.to_csv(os.path.join(RAW, "gaiadr3_ap_supp.csv"), index=False)
for tab in ["gaiadr3.nss_two_body_orbit", "gaiadr3.nss_acceleration_astro", "gaiadr3.nss_non_linear_spectro",
            "gaiadr3.nss_vim_fl", "gaiadr3.vari_summary", "gaiadr3.synthetic_photometry_gspc",
            "gaiadr3.xp_summary", "gaiadr3.galaxy_candidates"]:
    r = retry(lambda: gq(f"select * from {tab} where source_id={DR3}"), n=2)
    n = 0 if r is None else len(r)
    print(tab, "rows:", n)
    if n: r.to_csv(os.path.join(RAW, tab.replace(".", "_") + ".csv"), index=False)

# ---------------------------------------------------------------- VizieR catalogues
cats = {
    "TIC82_viz": "IV/39/tic82",
    "BJ21": "I/352/gedr3dis",
    "BJ18": "I/347/gaia2dis",
    "StarHorse19": "I/349/starhorse",
    "StarHorse21": "I/354/starhorse2021",
    "SHBoost24": "V/160",
    "Zhang23": "J/MNRAS/524/1855",
    "Huson25": "J/ApJ/984/58",
    "Ye25": "J/A+A/695/A75",
    "Kordopatis23": "J/A+A/669/A104",
    "Huang22": "J/ApJ/925/164",
    "Bai19": "J/AJ/158/93",
    "Andrae23": "J/ApJS/267/8",
    "RAVE6": "III/283",
    "GALAH3": "J/MNRAS/506/150",
    "GALAH4": "J/MNRAS/539/1520",   # may not exist
    "GSPC": "I/360/syntphot",
    "Hardegree25_none": "J/AJ/170/1",  # placeholder (expected empty)
}
V = Vizier(columns=["**", "_r"], row_limit=-1, timeout=120)
summary = {}
for name, cid in cats.items():
    res = retry(lambda: V.query_region(c0, radius=5 * u.arcsec, catalog=cid), n=3)
    if res is None:
        summary[name] = "query failed"; continue
    if len(res) == 0:
        summary[name] = "no match within 5 arcsec"; print(name, cid, "-> none"); continue
    rows = []
    for t in res:
        df = t.to_pandas(); df.insert(0, "_table", t.meta.get("name", ""))
        rows.append(df)
    df = pd.concat(rows, ignore_index=True, sort=False)
    df.to_csv(os.path.join(RAW, f"{name}.csv"), index=False)
    summary[name] = f"{len(df)} rows in {sorted(set(df['_table']))}"
    print(name, cid, "->", summary[name])
json.dump(dict(summary=summary, dr2_source_id=dr2id), open(os.path.join(RAW, "summary.json"), "w"), indent=1)
