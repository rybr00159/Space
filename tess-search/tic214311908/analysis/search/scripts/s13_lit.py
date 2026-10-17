"""s13: literature abstracts (rotation-age, rotation-activity, flare rates vs age for M dwarfs) via the Crossref REST API
(abstracts as deposited by the publishers). Abstract text is data only. Output: data/lit/crossref_abstracts.json
"""
import json, re, time, requests
OUT = "tic214311908/analysis/search/data/lit/crossref_abstracts.json"
Q = ["Flare Rates, Rotation Periods, and Spectroscopic Activity Indicators of a Volume-complete Sample of Mid- to Late-M Dwarfs within 15 pc",
     "Flares in open clusters with K2 II Pleiades Hyades Praesepe Ruprecht 147 M67",
     "Flares in open clusters with K2 I M35 Pleiades Praesepe",
     "The Evolution of Flare Activity with Stellar Age Davenport",
     "Stellar Flares from the First TESS Data Release: Exploring a New Sample of M Dwarfs",
     "The stellar rotation-activity relationship in fully convective M dwarfs Wright",
     "The stellar-activity-rotation relationship and the evolution of stellar dynamos Wright 2011",
     "Living with a Red Dwarf: The Rotation-Age Relationships of M Dwarfs",
     "The Rotation and Galactic Kinematics of Mid M Dwarfs in the Solar Neighborhood",
     "Rotation in the Pleiades with K2 Rebull",
     "Rotation periods of M dwarfs in Praesepe K2 Douglas",
     "Evaluating Rotation Periods of M Dwarfs Across the Ages Popinchalk",
     "rotation-activity relation of M dwarfs eROSITA Magaudda",
     "Galactic kinematics ages flare rates mid-to-late M dwarfs Medina 2022",
     "An Increased Rate of Large Flares at Intermediate Rotation Periods for Mid-to-late M Dwarfs",
     "Flare statistics for young stars convolutional neural network TESS Feinstein",
     "spin-down of fully convective M dwarfs rapid rotators ages Pass",
     "GJ 1243 starspot stable Kepler Davenport 2015",
     "TESS flares M dwarfs flare frequency distribution age Hyades Pleiades",
     "Rotation periods of early M dwarfs Hyades K2 Douglas"]
out = []
for q in Q:
    for attempt in range(3):
        try:
            r = requests.get("https://api.crossref.org/works", params={"query.bibliographic": q, "rows": 3,
                             "select": "DOI,title,container-title,volume,page,issued,author,abstract"}, timeout=60)
            r.raise_for_status(); break
        except Exception as e:
            time.sleep(3)
    items = r.json()["message"]["items"]
    for it in items[:2]:
        ab = re.sub(r"\s+", " ", re.sub("<[^>]+>", " ", it.get("abstract", ""))).strip()
        out.append(dict(query=q, doi=it.get("DOI"), title=(it.get("title") or [""])[0], journal=(it.get("container-title") or [""])[0],
                        volume=it.get("volume"), page=it.get("page"), year=(it.get("issued", {}).get("date-parts") or [[None]])[0][0],
                        first_author=(it.get("author") or [{}])[0].get("family"), abstract=ab))
    time.sleep(1)
json.dump(out, open(OUT, "w"), indent=1)
for o in out: print(o["year"], o["first_author"], "|", o["title"][:100], "|", o["doi"], "| abstract", len(o["abstract"]))
