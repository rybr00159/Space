"""n04: literature searches by identifier (read-only GET requests).
  - arXiv API (export.arxiv.org): searches all metadata fields (title, abstract, comments, journal-ref) for each identifier.
  - Semantic Scholar snippet search (full-text passages of open-access papers): moved to n04c (its public API rate-limits
    hard; the keyword paper search ignores quotes and returned only irrelevant "TiC"/"tic" papers, so it was dropped).
ADS is not queried here (no ADS API token was available when this ran); see n09_ads_fulltext.py.
Outputs: out/literature_search.json, raw/arxiv_*.xml (only non-empty responses), raw/s2_*.json
Usage: python -I n04_literature.py <novelty_dir>
"""
import sys, os, json, time, datetime, re, urllib.parse
import requests

D = sys.argv[1]
RAW, OUT = os.path.join(D, "raw"), os.path.join(D, "out")
UA = {"User-Agent": "tic214-novelty-check (read-only)"}
NOW = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
IDENTS = {
    "M3V": ["TIC 214311908", "214311908", "5925373272850784768", "2MASS J17235292-5119271", "J17235292-5119271",
            "J172352.86-511926.1", "UCAC4 194-162330", "APASS 24986318", "1eRASS J172352.9-511928", "J172352.9-511928",
            "PIC DR1 9844052"],
    "Mira": ["TIC 214311906", "214311906", "5925373268528684288", "J17235464-5119254", "ASASSN-V J172354.64-511925.3",
             "J172354.64-511925.3", "UCAC4 194-162350", "J172354.61-511925.3"],
    "pair": ["TIC 214311910", "214311910", "5925373272850782208", "5925373272827253632", "J17235419-5119283",
             "UCAC4 194-162344"],
    "coords": ["17 23 52.9 -51 19 27", "172352-511927", "J1723-5119"],
}
out = {"queried_utc": NOW, "arxiv": {}}
for grp, lst in IDENTS.items():
    for ident in lst:
        q = f'all:"{ident}"'
        url = "http://export.arxiv.org/api/query?" + urllib.parse.urlencode({"search_query": q, "max_results": 20})
        try:
            txt = requests.get(url, timeout=90, headers=UA).text
            m = re.search(r"<opensearch:totalResults[^>]*>(\d+)<", txt)
            n = int(m.group(1)) if m else txt.count("<entry>")
            titles = re.findall(r"<entry>.*?<title>(.*?)</title>", txt, flags=re.S)
            ids_ = re.findall(r"<entry>\s*<id>(.*?)</id>", txt, flags=re.S)
            out["arxiv"][ident] = {"n": n, "hits": [f"{i}: {re.sub(chr(10), ' ', t)}" for i, t in zip(ids_, titles)]}
            if n:
                open(os.path.join(RAW, f"arxiv_{re.sub('[^A-Za-z0-9]+', '_', ident)}.xml"), "w").write(txt)
        except Exception as e:
            out["arxiv"][ident] = {"error": str(e)[:200]}
        print("arXiv", repr(ident), out["arxiv"][ident]); time.sleep(3.5)
for grp, lst in ({} if True else IDENTS).items():   # Semantic Scholar part disabled here; see n04c
    for ident in lst:
        if grp == "coords":
            continue
        for kind, base in (("s2_paper_search", "https://api.semanticscholar.org/graph/v1/paper/search"),
                           ("s2_snippet_search", "https://api.semanticscholar.org/graph/v1/snippet/search")):
            params = {"query": f'"{ident}"', "limit": 10}
            if kind == "s2_paper_search":
                params["fields"] = "title,year,externalIds"
            rec = None
            for k in range(4):
                try:
                    r = requests.get(base, params=params, timeout=90, headers=UA)
                    if r.status_code == 429:
                        time.sleep(15 * (k + 1)); continue
                    rec = {"http_status": r.status_code}
                    try:
                        j = r.json()
                    except Exception:
                        j = {"raw": r.text[:300]}
                    data = j.get("data", []) if isinstance(j, dict) else []
                    rec["n_returned"] = len(data)
                    rec["total"] = j.get("total") if isinstance(j, dict) else None
                    if kind == "s2_paper_search":
                        rec["hits"] = [f'{d.get("year")}: {d.get("title")}' for d in data]
                    else:
                        hits = []
                        for d in data:
                            p = d.get("paper", {}); s = d.get("snippet", {})
                            txt = s.get("text", "")
                            hits.append({"paper": p.get("title"), "corpusId": p.get("corpusId"),
                                         "identifier_in_snippet": ident.replace('"', '') in txt, "snippet": txt[:300]})
                        rec["hits"] = hits
                    if not data and isinstance(j, dict) and "data" not in j:
                        rec["response"] = str(j)[:300]
                    break
                except Exception as e:
                    rec = {"error": str(e)[:200]}; time.sleep(10)
            out[kind][ident] = rec
            print(kind, repr(ident), json.dumps(rec)[:600]); time.sleep(4)
json.dump(out, open(os.path.join(OUT, "literature_search.json"), "w"), indent=1)
