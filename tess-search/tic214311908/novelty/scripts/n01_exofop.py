"""n01: ExoFOP-TESS status of TIC 214311908 (M3V), 214311906 (Mira), 214311910 (pair A/B) and the pair's split TIC
entries (1515709322, 1515709367), plus the 3.1" neighbour TIC 1515709327.
Read-only GET requests only.
  - full TOI and CTOI tables (CSV) -> raw/toi_table.csv.gz, raw/ctoi_table.csv.gz (gzipped copies) and searches by
    TIC ID and by position
  - target JSON, user-upload counts, file list, imaging/spectroscopy/time-series summaries, stellar companions, observing
    notes page -> raw/exofop_<TIC>_*.{json,txt}
Output: out/exofop_status.json
Usage: python -I n01_exofop.py <novelty_dir>
"""
import sys, os, io, re, json, gzip, hashlib, time, datetime
import numpy as np, pandas as pd, requests
import astropy.units as u
from astropy.coordinates import SkyCoord

D = sys.argv[1]
RAW, OUT = os.path.join(D, "raw"), os.path.join(D, "out")
UA = {"User-Agent": "tic214-novelty-check (read-only)"}
NOW = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
# J2000 (TIC epoch) positions
M3V = SkyCoord(260.970484934795, -51.3242024657333, unit="deg")
TICS = {"214311908": "M3V (target)", "214311906": "Mira", "214311910": "pair A/B (TIC 8 'SPLIT' parent)",
        "1515709322": "pair, split entry 1", "1515709367": "pair, split entry 2 (duplicate of 214311910)",
        "1515709327": "3.1 arcsec neighbour (G=17.5)"}


def get(url, timeout=180):
    for k in range(3):
        try:
            r = requests.get(url, timeout=timeout, headers=UA)
            return r
        except Exception as e:
            err = e; time.sleep(5)
    raise err


out = {"queried_utc": NOW, "tables": {}, "targets": {}}
for name, url in (("toi", "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv"),
                  ("ctoi", "https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv")):
    r = get(url)
    txt = r.text
    df = pd.read_csv(io.StringIO(txt), low_memory=False, comment=None)
    ticcol = [c for c in df.columns if c.lower().startswith("tic")][0]
    tics = df[ticcol].astype(str).str.replace(r"\.0$", "", regex=True).str.strip()
    rac = [c for c in df.columns if c.lower().startswith("ra")][0]
    decc = [c for c in df.columns if c.lower().startswith("dec")][0]
    try:
        cc = SkyCoord(df[rac].astype(str).values, df[decc].astype(str).values, unit=(u.hourangle, u.deg))
    except Exception:
        cc = SkyCoord(pd.to_numeric(df[rac], errors="coerce").values, pd.to_numeric(df[decc], errors="coerce").values,
                      unit="deg")
    sep = cc.separation(M3V).arcmin
    rec = {"url": url, "http_status": r.status_code, "n_rows": int(len(df)), "sha256": hashlib.sha256(txt.encode()).hexdigest(),
           "columns": list(df.columns), "tic_column": ticcol,
           "matches_by_tic": {t: int((tics == t).sum()) for t in TICS},
           "n_within_10arcmin": int((sep < 10).sum()), "n_within_60arcmin": int((sep < 60).sum()),
           "min_sep_arcmin": float(np.nanmin(sep))}
    j = int(np.nanargmin(sep))
    idcol = name.upper()
    rec["nearest"] = {str(k): str(v) for k, v in df.iloc[j].items() if k in (idcol, ticcol, rac, decc, "TFOPWG Disposition",
                                                                          "TESS Disposition", "User Disposition", "Period (days)",
                                                                          "Comments", "Notes")}
    rec["within_60arcmin"] = [dict(id=str(df.iloc[i][idcol]), tic=str(tics.iloc[i]), sep_arcmin=round(float(sep[i]), 2))
                              for i in np.argsort(sep) if sep[i] < 60]
    # compact extract (header + rows within 60 arcmin) instead of the full multi-MB table
    o = np.argsort(sep); o = o[sep[o] < 60]
    df.iloc[o].assign(sep_arcmin=np.round(sep[o], 2)).to_csv(os.path.join(RAW, f"{name}_rows_within_60arcmin.csv"), index=False)
    # freshness: highest ID and latest dates
    num = pd.to_numeric(df[idcol], errors="coerce")
    rec["max_id"] = float(np.nanmax(num))
    for c in df.columns:
        if re.search(r"date|lastmod", c, re.I):
            vals = pd.to_datetime(df[c], errors="coerce")
            if vals.notna().any():
                rec.setdefault("latest_dates", {})[c] = str(vals.max().date())
    # free-text mentions of our identifiers anywhere in the table
    rec["text_mentions"] = {s: int(txt.count(s)) for s in ("214311908", "214311906", "214311910", "5925373272850784768",
                                                          "J17235292", "J17235464")}
    out["tables"][name] = rec
    print(name, rec["n_rows"], rec["matches_by_tic"], "min sep", round(rec["min_sep_arcmin"], 2), "arcmin", rec["nearest"])

for tic, lab in TICS.items():
    t = {"label": lab}
    r = get(f"https://exofop.ipac.caltech.edu/tess/target.php?id={tic}&json")
    open(os.path.join(RAW, f"exofop_{tic}_target.json"), "w").write(r.text)
    try:
        d = r.json()
        t["target_json_keys"] = {k: (len(v) if isinstance(v, (list, dict)) else v) for k, v in d.items()
                                 if k in ("tois", "ctois", "planet_parameters", "stellar_companions", "imaging", "spectroscopy",
                                          "time_series", "files", "tags", "notes", "obs_notes", "nearby_targets")}
        t["all_list_lengths"] = {k: len(v) for k, v in d.items() if isinstance(v, list)}
        bi = d.get("basic_info", {})
        t["basic_info"] = {k: bi.get(k) for k in ("tic_id", "star_names", "ra", "dec", "pm_ra", "pm_dec", "confirmed_planets",
                                                  "kepler_name", "k2_name", "tic_contamination_ratio") if k in bi}
    except Exception as e:
        t["target_json_error"] = f"{r.status_code} {str(e)[:200]} {r.text[:200]}"
    for key, url in (("uploads", f"https://exofop.ipac.caltech.edu/tess/download_uploads.php?target={tic}"),
                     ("filelist", f"https://exofop.ipac.caltech.edu/tess/download_filelist.php?id={tic}&output=csv"),
                     ("imaging", f"https://exofop.ipac.caltech.edu/tess/download_imaging.php?output=csv&target={tic}"),
                     ("spect", f"https://exofop.ipac.caltech.edu/tess/download_spect.php?output=csv&target={tic}"),
                     ("tseries", f"https://exofop.ipac.caltech.edu/tess/download_tseries.php?output=csv&target={tic}"),
                     ("stellarcomp", f"https://exofop.ipac.caltech.edu/tess/download_stellarcomp.php?output=csv&target={tic}"),
                     ("obsnotes_page", f"https://exofop.ipac.caltech.edu/tess/view_obsnotes.php?id={tic}")):
        r = get(url)
        txt = r.text
        if key == "obsnotes_page":
            plain = re.sub(r"<script.*?</script>", " ", txt, flags=re.S)
            plain = re.sub(r"<style.*?</style>", " ", plain, flags=re.S)
            plain = re.sub(r"<[^>]+>", " ", plain); plain = re.sub(r"\s+", " ", plain).strip()
            open(os.path.join(RAW, f"exofop_{tic}_{key}.txt"), "w").write(plain[:5000])
            t[key] = {"http_status": r.status_code, "text_len": len(plain), "excerpt": plain[:600]}
        else:
            open(os.path.join(RAW, f"exofop_{tic}_{key}.txt"), "w").write(txt[:20000])
            lines = [x for x in txt.strip().splitlines() if x.strip()]
            t[key] = {"http_status": r.status_code, "n_lines": len(lines), "content": txt.strip()[:800]}
        time.sleep(1)
    out["targets"][tic] = t
    print(tic, lab, json.dumps({k: (v if k not in ("obsnotes_page",) else v["excerpt"][:200]) for k, v in t.items()}, default=str)[:1500])

json.dump(out, open(os.path.join(OUT, "exofop_status.json"), "w"), indent=1, default=str)
