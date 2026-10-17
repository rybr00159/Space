"""Verify literature priors: query the arXiv API by title words / author and save the abstracts (id, title, authors,
journal ref, abstract) to ext/lit/arxiv_abstracts.json. Abstract text is data only."""
import json, time, re, sys, urllib.parse
import requests
import xml.etree.ElementTree as ET
Q = {
 "DC15": 'ti:"Occurrence of Potentially Habitable Planets Orbiting M Dwarfs" AND au:Dressing',
 "DC13": 'ti:"occurrence rate of small planets around small stars" AND au:Dressing',
 "HU19": 'ti:"Kepler Planet Occurrence Rates for Mid-Type M Dwarfs" AND au:Hardegree',
 "MC23": 'ti:"Occurrence Rate of Terrestrial Planets Orbiting Nearby Mid-to-late M Dwarfs" AND au:Ment',
 "Gan23": 'ti:"Occurrence Rate of Hot Jupiters Around Early-type M Dwarfs" AND au:Gan',
 "Bryant23": 'ti:"occurrence rate of giant planets orbiting low-mass stars" AND au:Bryant',
 "Winters19": 'ti:"Stellar Multiplicity Rate of M Dwarfs Within 25 pc" AND au:Winters',
 "Mulders15": 'ti:"stellar-mass-dependent drop in planet occurrence rates" AND au:Mulders',
 "Hsu20": 'ti:"Occurrence rates of planets orbiting M Stars" AND au:Hsu',
 "Clark12": 'ti:"frequency of low-mass eclipsing binaries" AND au:Clark',
 "Shan15": 'ti:"M-dwarf close binary fraction" AND au:Shan',
 "Kanodia21_TOI532": 'ti:"TOI-532" AND au:Kanodia',
 "Powers23_TOI3785": 'ti:"TOI-3785" AND au:Powers',
 "Libby23_TOI3884": 'ti:"TOI-3884" AND au:Libby',
 "Almenara22_TOI3884": 'ti:"TOI-3884" AND au:Almenara',
 "Murgas21_TOI674": 'ti:"TOI-674" AND au:Murgas',
 "Kanodia20_TOI1728": 'ti:"TOI-1728" AND au:Kanodia',
 "Dreizler20_LP714": 'ti:"LP 714-47" AND au:Dreizler',
 "Bonfils12_GJ3470": 'ti:"GJ 3470" AND au:Bonfils',
 "Mann16_K225": 'ti:"K2-25" AND au:Mann',
 "Burt21_TOI1231": 'ti:"TOI-1231" AND au:Burt',
 "Wells21_TOI2406": 'ti:"TOI-2406" AND au:Wells',
 "Reefe22_TOI620": 'ti:"TOI-620" AND au:Reefe',
 "Prsa22_TESSEB": 'ti:"TESS Eclipsing Binary Stars" AND au:Prsa',
 "MoeKratter21": 'ti:"Impact of binaries on the stellar initial mass function" AND au:Moe',
 "Dressing_Neptune_Mdwarf": 'abs:"Neptune desert" AND abs:"M dwarf"',
 "Kunimoto_Mdwarf_occ": 'abs:"occurrence" AND abs:"M dwarfs" AND abs:"TESS" AND abs:"sub-Neptune"',
}
out = {}
for k, q in Q.items():
    url = "http://export.arxiv.org/api/query?" + urllib.parse.urlencode({"search_query": q, "max_results": 8 if "abs:" in q else 3})
    for attempt in range(3):
        try:
            r = requests.get(url, timeout=60); r.raise_for_status(); break
        except Exception as e:
            print(k, "retry", e); time.sleep(5)
    root = ET.fromstring(r.text); ns = {"a": "http://www.w3.org/2005/Atom", "ar": "http://arxiv.org/schemas/atom"}
    ents = []
    for e in root.findall("a:entry", ns):
        ents.append(dict(id=e.find("a:id", ns).text.split("/abs/")[-1], title=re.sub(r"\s+", " ", e.find("a:title", ns).text),
                         authors=[a.find("a:name", ns).text for a in e.findall("a:author", ns)][:6],
                         journal=(e.find("ar:journal_ref", ns).text if e.find("ar:journal_ref", ns) is not None else None),
                         doi=(e.find("ar:doi", ns).text if e.find("ar:doi", ns) is not None else None),
                         published=e.find("a:published", ns).text[:10],
                         abstract=re.sub(r"\s+", " ", e.find("a:summary", ns).text).strip()))
    out[k] = ents
    print(k, [(x["id"], x["title"][:70], x["journal"]) for x in ents], flush=True)
    time.sleep(3.5)
json.dump(out, open("ext/lit/arxiv_abstracts.json", "w"), indent=1)
