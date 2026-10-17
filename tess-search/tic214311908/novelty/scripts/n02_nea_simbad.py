"""n02: NASA Exoplanet Archive (TAP: pscomppars, ps, toi, stellarhosts) and SIMBAD (TAP: identifiers, object types,
bibliography, measurements, neighbours) for TIC 214311908 (M3V) and its neighbours. Read-only GET requests only.
Outputs: out/nea_status.json, out/simbad_status.json, raw/nea_*.csv, raw/simbad_*.csv
Usage: python -I n02_nea_simbad.py <novelty_dir>
"""
import sys, os, io, json, time, datetime, urllib.parse
import requests, pandas as pd

D = sys.argv[1]
RAW, OUT = os.path.join(D, "raw"), os.path.join(D, "out")
UA = {"User-Agent": "tic214-novelty-check (read-only)"}
NOW = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
RA, DEC = 260.970484934795, -51.3242024657333          # J2000, TIC
RA16, DEC16 = 260.96967369918, -51.3238212118           # J2016, Gaia DR3


def nea(q, tag):
    url = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?" + urllib.parse.urlencode({"query": q, "format": "csv"})
    for k in range(4):
        try:
            r = requests.get(url, timeout=180, headers=UA); break
        except Exception as e:
            print("retry", tag, e); time.sleep(10)
    open(os.path.join(RAW, f"nea_{tag}.csv"), "w").write(r.text[:200000])
    try:
        df = pd.read_csv(io.StringIO(r.text))
    except Exception:
        df = pd.DataFrame()
    return {"query": q, "http_status": r.status_code, "n_rows": int(len(df)), "rows": df.head(20).astype(str).to_dict("records")}


out = {"queried_utc": NOW}
import numpy as np
# box search (the archive's web firewall rejects ADQL CONTAINS/POINT queries sent through this proxy); box half-size r_deg
cone = lambda r_deg: (f"ra between {RA - r_deg / np.cos(np.radians(DEC)):.5f} and {RA + r_deg / np.cos(np.radians(DEC)):.5f} "
                      f"and dec between {DEC - r_deg:.5f} and {DEC + r_deg:.5f}")
out["pscomppars_cone_10arcmin"] = nea(f"select pl_name,hostname,tic_id,gaia_dr3_id,ra,dec,disc_year,pl_orbper from pscomppars where {cone(10/60)}", "pscomppars_cone")
out["pscomppars_by_id"] = nea("select pl_name,hostname,tic_id,gaia_dr3_id from pscomppars where tic_id in ('TIC 214311908','TIC 214311906','TIC 214311910') "
                              "or gaia_dr3_id like '%5925373272850784768%'", "pscomppars_id")
out["ps_cone_10arcmin"] = nea(f"select pl_name,hostname,tic_id,default_flag,disc_refname from ps where {cone(10/60)}", "ps_cone")
out["stellarhosts_cone_10arcmin"] = nea(f"select sy_name,hostname,tic_id,ra,dec from stellarhosts where {cone(10/60)}", "stellarhosts_cone")
out["toi_by_tid"] = nea("select toi,tid,ra,dec,tfopwg_disp,pl_orbper from toi where tid in (214311908,214311906,214311910,1515709322,1515709367,1515709327)", "toi_tid")
out["toi_cone_30arcmin"] = nea(f"select toi,tid,ra,dec,tfopwg_disp,pl_orbper from toi where {cone(0.5)}", "toi_cone")
out["toi_table_size"] = nea("select count(*) as n, max(toipfx) as max_toi, max(rowupdate) as last_update from toi", "toi_size")
out["pscomppars_size"] = nea("select count(*) as n from pscomppars", "pscomppars_size")
json.dump(out, open(os.path.join(OUT, "nea_status.json"), "w"), indent=1)
for k, v in out.items():
    if isinstance(v, dict):
        print("NEA", k, v["http_status"], v["n_rows"], v["rows"][:3])

# ---------------- SIMBAD TAP ----------------
def simbad(q, tag):
    url = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync?" + urllib.parse.urlencode(
        {"request": "doQuery", "lang": "adql", "format": "csv", "query": q})
    for k in range(4):
        try:
            r = requests.get(url, timeout=180, headers=UA); break
        except Exception as e:
            print("retry", tag, e); time.sleep(10)
    open(os.path.join(RAW, f"simbad_{tag}.csv"), "w").write(r.text[:200000])
    try:
        df = pd.read_csv(io.StringIO(r.text))
    except Exception:
        df = pd.DataFrame()
    time.sleep(1)
    return df, r.status_code


S = {"queried_utc": NOW, "objects": {}}
NAMES = {"M3V": "Gaia DR3 5925373272850784768", "Mira": "Gaia DR3 5925373268528684288",
         "pairA": "Gaia DR3 5925373272850782208", "star_3p1arcsec": "Gaia DR3 5925373272827260544"}
for lab, name in NAMES.items():
    rec = {"query_name": name}
    b, st = simbad("SELECT b.oid, b.main_id, b.otype, b.sp_type, b.ra, b.dec, b.pmra, b.pmdec, b.plx_value, b.rvz_radvel, "
                   "b.nbref FROM basic AS b JOIN ident AS i ON i.oidref = b.oid "
                   f"WHERE i.id = '{name}'", f"basic_{lab}")
    rec["http_status"] = st
    if len(b) == 0:
        rec["found"] = False; S["objects"][lab] = rec; continue
    rec["found"] = True
    rec["basic"] = b.iloc[0].astype(str).to_dict()
    oid = int(b.iloc[0]["oid"])
    ids, _ = simbad(f"SELECT id FROM ident WHERE oidref = {oid}", f"ident_{lab}")
    rec["identifiers"] = ids["id"].astype(str).tolist() if len(ids) else []
    ot, _ = simbad(f"SELECT otype, origin FROM otypes WHERE oidref = {oid}", f"otypes_{lab}")
    rec["otypes"] = ot.astype(str).to_dict("records")
    bib, _ = simbad("SELECT r.bibcode, r.\"year\", r.journal, r.title FROM has_ref AS h JOIN ref AS r ON r.oidbib = h.oidbibref "
                    f"WHERE h.oidref = {oid}", f"bib_{lab}")
    rec["bibliography"] = bib.astype(str).to_dict("records")
    for tab in ("mesVar", "mesSpT", "mesFe_h", "mesRot", "mesVelocities", "mesDistance", "mesXmm"):
        m, st = simbad(f"SELECT * FROM {tab} WHERE oidref = {oid}", f"{tab}_{lab}")
        rec[tab] = {"http_status": st, "n": int(len(m)), "rows": m.head(10).astype(str).to_dict("records")}
    S["objects"][lab] = rec
    print("SIMBAD", lab, rec["basic"].get("main_id"), rec["basic"].get("otype"), "nbref", rec["basic"].get("nbref"),
          "| ids:", rec["identifiers"], "| bib:", [x["bibcode"] for x in rec["bibliography"]],
          "| mes:", {t: rec[t]["n"] for t in ("mesVar", "mesSpT", "mesFe_h", "mesRot", "mesVelocities", "mesDistance", "mesXmm")})
# everything SIMBAD has within 60 arcsec of the M3V (J2000 position)
nb, st = simbad("SELECT b.main_id, b.otype, b.sp_type, b.ra, b.dec, b.nbref, "
                f"DISTANCE(POINT('ICRS', b.ra, b.dec), POINT('ICRS', {RA}, {DEC}))*3600 AS sep_arcsec FROM basic AS b "
                f"WHERE CONTAINS(POINT('ICRS', b.ra, b.dec), CIRCLE('ICRS', {RA}, {DEC}, {60/3600}))=1 ORDER BY sep_arcsec", "cone60")
S["cone_60arcsec"] = nb.astype(str).to_dict("records")
print("SIMBAD cone 60\":", nb[["main_id", "otype", "nbref", "sep_arcsec"]].to_string() if len(nb) else "empty")
# any SIMBAD object with a planet/candidate type within 30 arcmin
pl, st = simbad("SELECT b.main_id, b.otype, b.ra, b.dec, "
                f"DISTANCE(POINT('ICRS', b.ra, b.dec), POINT('ICRS', {RA}, {DEC}))*60 AS sep_arcmin FROM basic AS b "
                "WHERE b.otype IN ('Pl', 'Pl?', 'EB*', 'EB?') AND "
                f"CONTAINS(POINT('ICRS', b.ra, b.dec), CIRCLE('ICRS', {RA}, {DEC}, 0.5))=1 ORDER BY sep_arcmin", "planets_EB_30arcmin")
S["planet_or_EB_types_within_30arcmin"] = {"n": int(len(pl)), "n_planet_types": int(pl.otype.isin(["Pl", "Pl?"]).sum()) if len(pl) else 0,
                                           "within_3arcmin": pl[pl.sep_arcmin < 3].astype(str).to_dict("records") if len(pl) else []}
print("SIMBAD Pl/EB within 30':", S["planet_or_EB_types_within_30arcmin"])
json.dump(S, open(os.path.join(OUT, "simbad_status.json"), "w"), indent=1)
