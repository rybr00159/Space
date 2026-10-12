"""Catalogue compilation for the 9.3" neighbour TIC 87055315 = Gaia DR3 6720986208895208960.

Queries (with retries, all results cached as CSV in ../data/cat/):
  * TIC v8.2 row (MAST Catalogs, catalog='Tic') for TIC 87055315 and the target TIC 87055304
  * Gaia DR3 archive: astrophysical_parameters (GSP-Phot, FLAME, ESP-*, MSC), astrophysical_parameters_supp,
    NSS tables, vari_summary / vari_eclipsing_binary / vari_classifier_result, xp_summary, dr2_neighbourhood
  * VizieR stellar-parameter catalogues (Bailer-Jones+21, StarHorse 2019/2021, SHBoost, Zhang+23, Andrae+23,
    Huson+25, Ye+25, Huang+22, Kordopatis+23), photometry (2MASS, SkyMapper DR4, AllWISE, CatWISE, GSPC, APASS),
    and novelty / variability catalogues (VSX, ASAS-SN, Gaia DR3 EB/classifier, Kostov+25 tables 2-4, Prsa+22, ...)
Usage: python s01_catalogs.py
"""
import os, sys, json, time, warnings
import numpy as np, pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
warnings.filterwarnings("ignore")

ROOT = "/home/user/work/tic87/neighbour"
OUT = os.path.join(ROOT, "data/cat"); os.makedirs(OUT, exist_ok=True)
NB_DR3, NB_TIC = 6720986208895208960, 87055315
TG_DR3, TG_TIC = 6720986204578421888, 87055304
# Gaia DR3 epoch-2016 position of the neighbour; propagate to J2000 for catalogues at epoch ~2000
RA16, DE16, PMRA, PMDE = 272.609023, -43.993633, 0.483269, -5.794434
c16 = SkyCoord(RA16 * u.deg, DE16 * u.deg)
cosd = np.cos(np.radians(DE16))
c00 = SkyCoord((RA16 - PMRA * 16 / 3.6e6 / cosd) * u.deg, (DE16 - PMDE * 16 / 3.6e6) * u.deg)


def retry(fn, n=4, wait=10):
    for i in range(n):
        try:
            return fn()
        except Exception as e:
            print("  retry", i, repr(e)[:200], flush=True); time.sleep(wait)
    return None


summary = {}
# ------------------------------------------------------------------ TIC v8.2
from astroquery.mast import Catalogs
for tic in (NB_TIC, TG_TIC):
    f = os.path.join(OUT, f"tic82_{tic}.csv")
    if not os.path.exists(f):
        r = retry(lambda: Catalogs.query_criteria(catalog="Tic", ID=tic))
        if r is not None:
            r.to_pandas().to_csv(f, index=False)
    summary[f"tic82_{tic}"] = os.path.exists(f)
# TIC stars within 60" (for the flux budget)
f = os.path.join(OUT, "tic82_cone60.csv")
if not os.path.exists(f):
    r = retry(lambda: Catalogs.query_region(c16, radius=60 * u.arcsec, catalog="Tic"))
    if r is not None:
        r.to_pandas().to_csv(f, index=False)
print("TIC done", flush=True)

# ------------------------------------------------------------------ Gaia archive
from astroquery.gaia import Gaia
Gaia.ROW_LIMIT = -1


def gq(q):
    return Gaia.launch_job(q).get_results().to_pandas()


tabs = ["gaiadr3.astrophysical_parameters", "gaiadr3.astrophysical_parameters_supp", "gaiadr3.nss_two_body_orbit",
        "gaiadr3.nss_acceleration_astro", "gaiadr3.nss_non_linear_spectro", "gaiadr3.nss_vim_fl", "gaiadr3.vari_summary",
        "gaiadr3.vari_eclipsing_binary", "gaiadr3.vari_classifier_result", "gaiadr3.xp_summary",
        "gaiadr3.synthetic_photometry_gspc", "gaiadr3.dr2_neighbourhood"]
for tab in tabs:
    for sid, lab in ((NB_DR3, "nb"), (TG_DR3, "tg")):
        f = os.path.join(OUT, f"{tab.replace('.', '_')}_{lab}.csv")
        if os.path.exists(f):
            continue
        key = "dr3_source_id" if "neighbourhood" in tab else "source_id"
        r = retry(lambda: gq(f"select * from {tab} where {key}={sid}"), n=3)
        n = -1 if r is None else len(r)
        summary[f"{tab}_{lab}"] = n
        print(tab, lab, "rows:", n, flush=True)
        if r is not None:
            r.to_csv(f, index=False)
# DR2 counterpart of the neighbour (TIC v8 is DR2-based)
f = os.path.join(OUT, "gaiadr2_source_nb.csv")
if not os.path.exists(f):
    try:
        nb = pd.read_csv(os.path.join(OUT, "gaiadr3_dr2_neighbourhood_nb.csv"))
        dr2id = int(nb.sort_values("angular_distance").iloc[0]["dr2_source_id"])
        r = retry(lambda: gq(f"select * from gaiadr2.gaia_source where source_id={dr2id}"), n=3)
        if r is not None: r.to_csv(f, index=False)
    except Exception as e:
        print("DR2 lookup failed", e)
print("Gaia done", flush=True)

# ------------------------------------------------------------------ VizieR
from astroquery.vizier import Vizier
V = Vizier(columns=["**", "_r"], row_limit=-1, timeout=180)
cats = {
    # stellar parameters
    "TIC82_viz": ("IV/39/tic82", c00, 3), "BJ21": ("I/352/gedr3dis", c16, 2), "StarHorse19": ("I/349/starhorse", c16, 2),
    "StarHorse21": ("I/354/starhorse2021", c16, 2), "SHBoost24": ("V/160", c16, 2), "Zhang23": ("J/MNRAS/524/1855", c16, 2),
    "Huson25": ("J/ApJ/984/58", c16, 2), "Ye25": ("J/A+A/695/A75", c16, 2), "Kordopatis23": ("J/A+A/669/A104", c16, 2),
    "Huang22": ("J/ApJ/925/164", c16, 2), "Andrae23": ("J/ApJS/267/8", c16, 2), "GaiaDR3_paramp": ("I/355/paramp", c16, 2),
    "GaiaDR3_main": ("I/355/gaiadr3", c16, 2), "RAVE6": ("III/283", c16, 3), "GALAH": ("J/MNRAS/506/150", c16, 3),
    "Gaia_ESO": ("J/A+A/666/A121", c16, 3), "Fouesneau22_gspspec": ("I/355/spectra", c16, 2),
    # photometry
    "2MASS": ("II/246/out", c00, 3), "SMSS_DR4": ("II/379/smssdr4", c16, 2), "AllWISE": ("II/328/allwise", c00, 4),
    "CatWISE": ("II/365/catwise", c16, 4), "GSPC": ("I/360/syntphot", c16, 2), "APASS9": ("II/336/apass9", c00, 3),
    "UKIDSS_VVV": ("II/348/vvv2", c16, 2), "DENIS": ("B/denis/denis", c00, 3), "GALEX": ("II/335/galex_ais", c16, 4),
    # novelty / variability (within 21" = 1 TESS pixel)
    "VSX": ("B/vsx/vsx", c16, 21), "ASASSN_var": ("II/366/catv2021", c16, 21), "GaiaDR3_veb": ("I/358/veb", c16, 21),
    "GaiaDR3_vclass": ("I/358/vclassre", c16, 21), "GaiaDR3_vari_sum": ("I/358/vari", c16, 21),
    "Kostov25_t2": ("J/ApJS/279/50/table2", c16, 21), "Kostov25_t3": ("J/ApJS/279/50/table3", c16, 21),
    "Kostov25_t4": ("J/ApJS/279/50/table4", c16, 21), "Prsa22": ("J/ApJS/258/16", c16, 21),
    "Rowan22_ASASSN_EB": ("J/MNRAS/517/2190", c16, 21), "Gao25_TESSvar": ("J/ApJS/276/57", c16, 21),
    "Heinze18_ATLAS": ("J/AJ/156/241", c16, 21), "IJspeert24": ("J/A+A/691/A242", c16, 21),
    "Petrosky21_WISE": ("J/MNRAS/503/3975", c16, 21), "Thiemann21_SWASP": ("J/MNRAS/502/1299", c16, 21),
    "Green23": ("J/MNRAS/522/29", c16, 21), "Fetherolf23_TESSvar": ("J/ApJS/268/4", c16, 21),
    "Howard22_EvryFlare": ("J/ApJ/938/103", c16, 21), "Justesen21_EB": ("J/ApJ/912/123", c16, 21),
}
for name, (cid, cc, rad) in cats.items():
    f = os.path.join(OUT, f"viz_{name}.csv")
    if os.path.exists(f):
        summary[name] = "cached"; continue
    res = retry(lambda: V.query_region(cc, radius=rad * u.arcsec, catalog=cid), n=3)
    if res is None:
        summary[name] = "query failed"; print(name, "FAILED", flush=True); continue
    if len(res) == 0:
        summary[name] = f"no match within {rad} arcsec"; print(name, cid, "-> none", flush=True)
        pd.DataFrame().to_csv(f, index=False); continue
    rows = []
    for t in res:
        df = t.to_pandas(); df.insert(0, "_table", t.meta.get("name", ""))
        rows.append(df)
    df = pd.concat(rows, ignore_index=True, sort=False)
    df.to_csv(f, index=False)
    summary[name] = f"{len(df)} rows in {sorted(set(df['_table']))}"
    print(name, cid, "->", summary[name], flush=True)
json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=1, default=str)
print("ALL DONE")
