import warnings, glob, re
warnings.filterwarnings("ignore")
from astroquery.mast import Observations
import xml.etree.ElementTree as ET
obs = Observations.query_criteria(target_name=["61141412"], obs_collection=["TESS", "HLSP"])
obs = obs[[p in ("SPOC", "TESS-SPOC") for p in obs["provenance_name"]]]
prods = Observations.get_product_list(obs)
x = prods[[str(f).endswith("dvr.xml") for f in prods["productFilename"]]]
Observations.download_products(x, download_dir="dv")
for f in sorted(glob.glob("dv/**/*dvr.xml", recursive=True)):
    root = ET.parse(f).getroot()
    ns = {"dv": root.tag.split("}")[0].strip("{")} if root.tag.startswith("{") else {}
    print("\n" + f.split("/")[-1])
    for pr in root.iter():
        if pr.tag.endswith("planetResults"):
            num = pr.attrib.get("planetNumber")
            fit = [c for c in pr.iter() if c.tag.endswith("allTransitsFit")]
            vals = {}
            if fit:
                for p in fit[0].iter():
                    if p.tag.endswith("modelParameter") and p.attrib.get("name") in ("orbitalPeriodDays", "transitEpochBtjd", "transitDepthPpm", "transitDurationHours", "planetRadiusEarthRadii", "impactParameter", "ratioSemiMajorAxisToStarRadius"):
                        vals[p.attrib["name"]] = round(float(p.attrib["value"]), 5)
            mes = [c.attrib for c in pr.iter() if c.tag.endswith("planetCandidate")]
            ms = {k: mes[0].get(k) for k in ("maxMultipleEventSigma", "maxSingleEventSigma", "observedTransitCount", "orbitalPeriodInDays", "epochTjd") if mes and k in mes[0]}
            print(f"  TCE {num}: {vals} {ms}")
