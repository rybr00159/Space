"""Second round of arXiv API queries (missing titles from c14 and Neptune-size occurrence around M dwarfs).
Output: ext/lit/arxiv_abstracts_2.json"""
import json, time, re, urllib.parse, requests
import xml.etree.ElementTree as ET
Q = {
 "HU19": 'au:Hardegree_Ullman AND ti:"M Dwarfs" AND ti:Kepler',
 "Clark12": 'au:Clark AND au:Blake AND ti:eclipsing',
 "Shan15": 'au:Shan AND au:Johnson AND ti:binary',
 "Libby23_TOI3884": 'au:Libby_Roberts AND ti:TOI',
 "Murgas21_TOI674": 'au:Murgas AND ti:Neptunes',
 "Kanodia20_TOI1728": 'au:Kanodia AND ti:"super Neptune"',
 "Dreizler20_LP714": 'au:Dreizler AND ti:"Neptune desert"',
 "Bonfils12_GJ3470": 'au:Bonfils AND ti:"Uranus"',
 "MoeKratter21": 'au:Moe AND au:Kratter AND ti:binaries',
 "Nep_occ_M": 'abs:"Neptune" AND abs:"occurrence" AND abs:"M dwarf" AND abs:"TESS"',
 "Nep_occ_M2": 'ti:"Neptune" AND ti:"M dwarf" AND abs:"occurrence"',
 "Mdwarf_EB_occ": 'abs:"eclipsing binaries" AND abs:"M dwarfs" AND abs:"occurrence" AND abs:"period"',
 "Kepler_Mdwarf_giants": 'abs:"M dwarfs" AND abs:"occurrence" AND abs:"Kepler" AND abs:"Neptune-size"',
 "Bryson_Mdwarf": 'au:Hardegree_Ullman AND ti:"Scaling K2"',
 "Sabotta21": 'au:Sabotta AND ti:CARMENES AND ti:occurrence',
 "Pinamonti22": 'au:Pinamonti AND ti:occurrence AND ti:HADES',
 "Ribas23": 'au:Ribas AND ti:"CARMENES" AND ti:"Data Release 1"',
}
out = {}
ns = {"a": "http://www.w3.org/2005/Atom", "ar": "http://arxiv.org/schemas/atom"}
for k, q in Q.items():
    url = "http://export.arxiv.org/api/query?" + urllib.parse.urlencode({"search_query": q, "max_results": 8})
    for attempt in range(3):
        try:
            r = requests.get(url, timeout=60); r.raise_for_status(); break
        except Exception as e:
            print(k, "retry", e); time.sleep(5)
    root = ET.fromstring(r.text); ents = []
    for e in root.findall("a:entry", ns):
        ents.append(dict(id=e.find("a:id", ns).text.split("/abs/")[-1], title=re.sub(r"\s+", " ", e.find("a:title", ns).text),
                         authors=[a.find("a:name", ns).text for a in e.findall("a:author", ns)][:6],
                         journal=(e.find("ar:journal_ref", ns).text if e.find("ar:journal_ref", ns) is not None else None),
                         doi=(e.find("ar:doi", ns).text if e.find("ar:doi", ns) is not None else None),
                         published=e.find("a:published", ns).text[:10], abstract=re.sub(r"\s+", " ", e.find("a:summary", ns).text).strip()))
    out[k] = ents
    print(k, [(x["id"], x["title"][:80]) for x in ents], flush=True)
    time.sleep(3.5)
json.dump(out, open("ext/lit/arxiv_abstracts_2.json", "w"), indent=1)
