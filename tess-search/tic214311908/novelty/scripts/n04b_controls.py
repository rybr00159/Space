"""n04b: (1) positive controls for the arXiv identifier searches in n04 (the same query form must find known objects);
(2) Semantic Scholar snippet search (full-text passages of open-access papers) for the main identifiers of the M3V and
the Mira, with positive controls, run slowly because the public API rate-limits hard (HTTP 429).
Output: out/literature_controls.json
Usage: python -I n04b_controls.py <novelty_dir>
"""
import sys, os, re, json, time, datetime, urllib.parse
import requests

D = sys.argv[1]
UA = {"User-Agent": "tic214-novelty-check (read-only)"}
out = {"queried_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"), "arxiv_controls": {},
       "s2_snippet": {}}
for ident in ("TIC 172900988", "172900988", "TOI-3884", "TOI 3884"):
    url = "http://export.arxiv.org/api/query?" + urllib.parse.urlencode({"search_query": f'all:"{ident}"', "max_results": 5})
    txt = requests.get(url, timeout=90, headers=UA).text
    m = re.search(r"<opensearch:totalResults[^>]*>(\d+)<", txt)
    out["arxiv_controls"][ident] = int(m.group(1)) if m else -1
    print("arXiv control", ident, out["arxiv_controls"][ident], flush=True); time.sleep(3.5)
time.sleep(60)
for ident in ("TOI-3884", "TIC 172900988", "TIC 214311908", "5925373272850784768", "J17235292-5119271",
              "1eRASS J172352.9-511928", "TIC 214311906", "5925373268528684288", "ASASSN-V J172354.64-511925.3"):
    rec = {"http_status": None}
    for k in range(5):
        try:
            r = requests.get("https://api.semanticscholar.org/graph/v1/snippet/search",
                             params={"query": f'"{ident}"', "limit": 10}, timeout=90, headers=UA)
        except Exception as e:
            rec = {"error": str(e)[:200]}; time.sleep(30); continue
        rec["http_status"] = r.status_code
        if r.status_code == 429:
            time.sleep(30 * (k + 1)); continue
        try:
            data = r.json().get("data", [])
        except Exception:
            data = []
        rec.update({"n": len(data),
                    "with_identifier_in_text": sum(ident in d.get("snippet", {}).get("text", "") for d in data),
                    "hits": [{"title": d.get("paper", {}).get("title"),
                              "identifier_in_text": ident in d.get("snippet", {}).get("text", ""),
                              "snippet": d.get("snippet", {}).get("text", "")[:200]} for d in data]})
        break
    out["s2_snippet"][ident] = rec
    print("S2 snippet", ident, json.dumps(rec)[:400], flush=True); time.sleep(20)
json.dump(out, open(os.path.join(D, "out", "literature_controls.json"), "w"), indent=1)
