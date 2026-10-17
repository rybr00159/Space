"""n09: NASA ADS full-text search (read-only GET requests to api.adsabs.harvard.edu/v1/search/query).
  - Targets: TIC 214311908, Gaia DR3 5925373272850784768 and 2MASS J17235292-5119271 (the M3V), and TIC 214311906 (the
    Mira), each in every written form that ADS indexes as a different token (no-space, hyphen, ".01", no "J" prefix).
  - Field `full` = body text + acknowledgements + abstract + title + keywords. No year or collection filter.
  - Every hit is checked in its highlight snippets: ADS expands synonyms (e.g. "TIC" -> "Tycho Input Catalog"), so a
    hit counts only if the identifier itself is in the returned text.
  - Positive controls of the same written forms (TIC 172900988, TRAPPIST-1, Proxima Cen, Barnard's star), body-only vs
    abstract counts, index freshness and the share of records with indexed full text.
Needs an ADS API token: set ADS_API_TOKEN, or run where a proxy adds it (as here).
Outputs: out/ads_fulltext.json, raw/ads_<query>.json (only target queries with hits)
Usage: python -I n09_ads_fulltext.py <novelty_dir>
"""
import sys, os, json, time, datetime, re
import requests

D = sys.argv[1]
RAW, OUT = os.path.join(D, "raw"), os.path.join(D, "out")
API = "https://api.adsabs.harvard.edu/v1/search/query"
HDR = {"User-Agent": "tic214-novelty-check (read-only)"}
if os.environ.get("ADS_API_TOKEN"):
    HDR["Authorization"] = "Bearer " + os.environ["ADS_API_TOKEN"]
NOW = lambda: datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
TARGETS = {
    "M3V TIC 214311908": ["214311908", "TIC 214311908", "TIC214311908", "TIC-214311908", "214311908.01"],
    "M3V Gaia DR3 5925373272850784768": ["5925373272850784768", "Gaia DR3 5925373272850784768"],
    "M3V 2MASS J17235292-5119271": ["J17235292-5119271", "2MASS J17235292-5119271", "2MASSJ17235292-5119271",
                                    "17235292-5119271", "J17235292", "5119271"],
    "Mira TIC 214311906": ["214311906", "TIC 214311906", "TIC214311906", "TIC-214311906", "214311906.01"],
}
CONTROLS = {   # same written forms as the targets, for stars that papers do name
    "TIC": ["172900988", "TIC 172900988", "TIC172900988"],
    "Gaia DR3": ["2635476908753563008", "5853498713190525696", "4472832130942575872"],
    "2MASS": ["J23062928-0502285", "J23062928", "23062928-0502285"],
}
VENUES = ["AJ", "ApJ", "ApJL", "ApJS", "MNRAS", "A&A", "PASP", "RNAAS", "AN", "JAVSO", "OEJV", "IBVS", "PZ", "ATel"]


def ads(q, rows=0, fl="bibcode", **extra):
    """One ADS query; retries on 429/5xx. Returns the decoded JSON (or {'error': ...})."""
    params = {"q": q, "rows": rows, "fl": fl, **extra}
    for k in range(4):
        try:
            r = requests.get(API, params=params, headers=HDR, timeout=90)
        except Exception as e:
            err = {"error": str(e)[:200]}; time.sleep(10 * (k + 1)); continue
        if r.status_code == 429 or r.status_code >= 500:
            err = {"error": f"HTTP {r.status_code}"}; time.sleep(30 * (k + 1)); continue
        try:
            j = r.json()
        except Exception:
            return {"error": f"HTTP {r.status_code}: {r.text[:200]}"}
        j["_http_status"], j["_ratelimit_remaining"] = r.status_code, r.headers.get("X-RateLimit-Remaining")
        time.sleep(1)
        return j
    return err


def n_found(q):
    j = ads(q)
    return j["response"]["numFound"] if "response" in j else None


def norm(s):
    """Text with tags and every non-alphanumeric character removed, upper case (for the identifier-in-snippet test)."""
    return re.sub(r"[^0-9A-Z]", "", re.sub(r"<[^>]+>", "", s).upper())


def search(field, ident, rows=100, keep=None):
    q = f'{field}:"{ident}"'
    j = ads(q, rows=rows, fl="id,bibcode,title,year,doctype,pub", hl="true", **{"hl.fl": "body,abstract,title,ack",
            "hl.snippets": 3, "hl.fragsize": 200})
    if "response" not in j:
        return {"q": q, "error": j.get("error", str(j)[:200])}
    hl, key = j.get("highlighting", {}), norm(ident)
    hits = []
    for d in j["response"]["docs"]:
        snips = [s for fld in ("body", "abstract", "title", "ack") for s in hl.get(d["id"], {}).get(fld, [])]
        hits.append({"bibcode": d["bibcode"], "year": d.get("year"), "doctype": d.get("doctype"),
                     "title": (d.get("title") or [""])[0][:150], "identifier_in_text": any(key in norm(s) for s in snips),
                     "snippets": [re.sub(r"\s+", " ", s)[:300] for s in snips[:2]]})
    rec = {"q": q, "numFound": j["response"]["numFound"], "n_identifier_in_text": sum(h["identifier_in_text"] for h in hits),
           "hits": hits if keep is None else hits[:keep]}
    if keep is None and hits:
        open(os.path.join(RAW, "ads_" + re.sub("[^A-Za-z0-9]+", "_", q).strip("_") + ".json"), "w").write(
            json.dumps(j, indent=1))
    return rec


out = {"queried_utc": NOW(), "api": API, "field_full": "body + ack + abstract + title + keywords",
       "filters": "none (all ADS collections, all years)", "targets": {}, "controls": {}}
for grp, lst in TARGETS.items():
    out["targets"][grp] = {}
    for ident in lst:
        rec = search("full", ident)
        out["targets"][grp][ident] = rec
        print("target", repr(rec["q"]), rec.get("numFound"), "identifier in text:", rec.get("n_identifier_in_text"),
              [h["bibcode"] for h in rec.get("hits", [])][:10], rec.get("error", ""), flush=True)

for grp, lst in CONTROLS.items():
    for ident in lst:
        rec = {f: n_found(f'{f}:"{ident}"') for f in ("full", "body", "abs")}
        rec["body_not_abs"] = n_found(f'body:"{ident}" -abs:"{ident}"')
        rec["sample"] = search("body", ident, rows=5, keep=3)["hits"]
        out["controls"][f"{grp}: {ident}"] = rec
        print("control", grp, repr(ident), {k: v for k, v in rec.items() if k != "sample"},
              [(h["bibcode"], h["identifier_in_text"]) for h in rec["sample"]], flush=True)
# A decimal number is one token: "700.01" does not match "700" (why the ".01" forms are searched separately)
out["controls"]["tokens: 700.01 without 700"] = n_found('body:"700.01" -body:"700"')
# The J-prefixed half of a 2MASS name is one token ("23062928" does not match "J23062928"); a Unicode minus is a hyphen
out["controls"]["tokens: 23062928 (2MASS half without the J)"] = n_found('body:"23062928"')
out["controls"]["tokens: J23062928−0502285 (Unicode minus)"] = n_found('body:"J23062928−0502285"')
print("tokens", {k: v for k, v in out["controls"].items() if k.startswith("tokens")}, flush=True)

j = ads("bibstem:arXiv entdate:[NOW-7DAYS TO *]", rows=1, fl="bibcode,entry_date", sort="entry_date desc")
out["freshness"] = {"arxiv_entered_last_7d": j["response"]["numFound"],
                    "arxiv_entered_last_7d_with_body": n_found("has:body bibstem:arXiv entdate:[NOW-7DAYS TO *]"),
                    "latest_entry": [(d["bibcode"], d.get("entry_date")) for d in j["response"]["docs"]]}
print("freshness", out["freshness"], flush=True)

cov = {}
for label, base in [("refereed astronomy 2019-2026", "database:astronomy property:refereed year:2019-2026"),
                    ("arXiv 2019-2026", "bibstem:arXiv year:2019-2026")] + \
                   [(v + " 2019-2026", f'bibstem:"{v}" year:2019-2026') for v in VENUES]:
    n, nb = n_found(base), n_found("has:body " + base)
    cov[label] = {"records": n, "with_body": nb, "fraction": round(nb / n, 3) if n and nb is not None else None}
    print("coverage", label, cov[label], flush=True)
out["coverage_has_body"] = cov
out["queried_utc_end"] = NOW()
json.dump(out, open(os.path.join(OUT, "ads_fulltext.json"), "w"), indent=1, ensure_ascii=False)
