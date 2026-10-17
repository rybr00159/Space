"""Download the light-curve files (ELEANOR-LITE, TARS, TEQUILA, QLP) for the three stars and the Mira's close neighbours."""
import warnings, os, shutil; warnings.filterwarnings("ignore")
from astroquery.mast import Observations
from astropy.table import Table
t = Table.read("out/mast_timeseries_60arcsec.ecsv")
keep = {"214311906", "214311908", "214311910", "1515709367", "1515709322", "s0012-cam1-ccd2-x1515-y0269"}
t = t[[str(n) in keep for n in t["target_name"]]]
os.makedirs("lc/raw", exist_ok=True)
for r in t:
    url = r["dataURL"]; fn = os.path.join("lc", os.path.basename(url))
    if r["provenance_name"] == "QLP" and str(r["target_name"]) in ("214311906", "214311908", "214311910"):
        src = "/home/user/work/lc_cache/" + os.path.basename(url)
        if os.path.exists(src): shutil.copy(src, fn); print("copied", fn); continue
    if not os.path.exists(fn):
        st = Observations.download_file(url, local_path=fn)
        print(st, fn, os.path.getsize(fn) if os.path.exists(fn) else None)
