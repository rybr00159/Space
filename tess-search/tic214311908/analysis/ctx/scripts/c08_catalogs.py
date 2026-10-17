"""Literature/catalogue status of TIC 214311906 / 214311908 / 214311910: ExoFOP (TOI/CTOI/notes/files), the full TOI and
CTOI tables, SIMBAD (identifiers, object types, bibliography), and every VizieR catalogue with an entry within 20" of
each star (variable-star and eclipsing-binary catalogues highlighted). Output: out/catalog_status.json, out/vizier_hits.csv
"""
import json, io, warnings, re
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, requests
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
import astropy.units as u
from astropy.coordinates import SkyCoord

stars = {"214311906": ("Mira", 260.9776003, -51.3237201), "214311908": ("M3V", 260.9696737, -51.3238212),
         "214311910": ("pairA", 260.9748734, -51.3240552)}
out = {"exofop": {}, "toi_table": {}, "ctoi_table": {}, "simbad": {}, "vizier": {}}
# ExoFOP target JSON (already fetched in ext/exofop) + obsnotes page
for tic in stars:
    d = json.load(open(f"ext/exofop/target_{tic}.json"))
    out["exofop"][tic] = {k: len(d[k]) for k in ("tois", "ctois", "planet_parameters", "stellar_companions", "imaging", "spectroscopy", "time_series", "files")}
    try:
        r = requests.get(f"https://exofop.ipac.caltech.edu/tess/view_obsnotes.php?id={tic}", timeout=60)
        txt = re.sub(r"<[^>]+>", " ", r.text); txt = re.sub(r"\s+", " ", txt)
        out["exofop"][tic]["obsnotes_page_status"] = r.status_code
        out["exofop"][tic]["obsnotes_mentions_tic"] = txt.count(tic)
        out["exofop"][tic]["obsnotes_excerpt"] = txt[:300]
    except Exception as e:
        out["exofop"][tic]["obsnotes_error"] = str(e)
# Full TOI and CTOI tables
for name, url in (("toi_table", "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv"),
                  ("ctoi_table", "https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv")):
    r = requests.get(url, timeout=180); df = pd.read_csv(io.StringIO(r.text), low_memory=False)
    ticcol = [c for c in df.columns if c.lower().startswith("tic")][0]
    out[name]["n_rows"] = len(df); out[name]["tic_column"] = ticcol
    tics = df[ticcol].astype(str).str.replace(r"\.0$", "", regex=True)
    for tic in list(stars) + ["1515709322", "1515709367"]:
        out[name][tic] = int((tics == tic).sum())
    # also anything within 2 arcmin of the Mira
    rac = [c for c in df.columns if c.lower().startswith("ra")][0]; decc = [c for c in df.columns if c.lower().startswith("dec")][0]
    try:
        cc = SkyCoord(df[rac].astype(str).values, df[decc].astype(str).values, unit=(u.hourangle, u.deg))
    except Exception:
        cc = SkyCoord(pd.to_numeric(df[rac], errors="coerce").values, pd.to_numeric(df[decc], errors="coerce").values, unit="deg")
    sep = cc.separation(SkyCoord(260.97765, -51.32369, unit="deg")).arcmin
    out[name]["within_10arcmin"] = df[sep < 10][[ticcol]].astype(str).values.ravel().tolist()
    out[name]["min_sep_arcmin"] = float(np.nanmin(sep))
# SIMBAD
Simbad.add_votable_fields("otype", "otypes", "ids", "sp_type", "flux(V)", "flux(G)") if hasattr(Simbad, "add_votable_fields") else None
for tic, (lab, ra, dec) in stars.items():
    try:
        r = Simbad.query_region(SkyCoord(ra, dec, unit="deg"), radius=8 * u.arcsec)
        if r is None:
            out["simbad"][tic] = "no object within 8 arcsec"; continue
        rows = []
        for x in r:
            row = {k: str(x[k]) for k in r.colnames if k.lower() in ("main_id", "otype", "otypes", "ids", "sp_type", "flux_v", "flux_g", "v", "g", "ra", "dec")}
            rows.append(row)
        out["simbad"][tic] = rows
        # bibliography count
        for row in rows:
            mid = row.get("main_id") or row.get("MAIN_ID")
            try:
                q = f"SELECT bibcode, journal, title FROM ref JOIN has_ref ON oidref = ref.oidbib JOIN ident ON ident.oidref = has_ref.oidref WHERE id = '{mid}'"
                b = Simbad.query_tap(q)
                row["n_refs"] = len(b); row["refs"] = [f"{y['bibcode']}: {str(y['title'])[:120]}" for y in b]
            except Exception as e:
                row["refs_error"] = str(e)[:200]
    except Exception as e:
        out["simbad"][tic] = f"error {e}"
# VizieR: all catalogues with an entry within 20" (5" for the M3V), flagging variability catalogues
v = Vizier(columns=["*"], row_limit=50)
hits = []
for tic, (lab, ra, dec) in stars.items():
    res = v.query_region(SkyCoord(ra, dec, unit="deg"), radius=(5 if lab != "Mira" else 3) * u.arcsec, catalog=None)
    for t in res:
        desc = t.meta.get("description", "")
        hits.append(dict(tic=tic, star=lab, catalog=t.meta.get("name", ""), description=desc, n=len(t)))
hv = pd.DataFrame(hits)
hv.to_csv("out/vizier_hits.csv", index=False)
kw = re.compile(r"variab|eclips|binar|transit|planet|period|rotation|flare|TOI|exoplanet|light curve|lightcurve|VSX|Mira|long-period", re.I)
out["vizier"] = {lab: hv[(hv.star == lab)].catalog.tolist() for lab in ("Mira", "M3V", "pairA")}
out["vizier_variability_related"] = hv[hv.description.str.contains(kw) | hv.catalog.str.contains(r"vsx|I/358|ASAS|ATLAS|J/ApJS/258/16|J/A\+A/652", case=False)].to_dict("records")
json.dump(out, open("out/catalog_status.json", "w"), indent=1, default=str)
print(json.dumps({k: out[k] for k in ("exofop", "toi_table", "ctoi_table")}, indent=1, default=str)[:4000])
print("SIMBAD:", json.dumps(out["simbad"], indent=1, default=str)[:4000])
print("VizieR variability-related:", json.dumps(out["vizier_variability_related"], indent=1)[:6000])
print("n VizieR catalogues:", {k: len(v) for k, v in out["vizier"].items()})
