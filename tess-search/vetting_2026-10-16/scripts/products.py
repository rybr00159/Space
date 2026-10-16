import sys, warnings
sys.path.insert(0, "/home/user/work/src"); warnings.filterwarnings("ignore")
import pandas as pd
from astroquery.mast import Observations
for tic in [90086266, 61332759, 214311906, 55829034, 61141412]:
    obs = Observations.query_criteria(target_name=[str(tic)], obs_collection=["TESS", "HLSP"])
    df = obs.to_pandas()[["provenance_name", "sequence_number", "t_exptime", "dataURL"]]
    df = df[df.dataURL.astype(str).str.contains("lc.fits|llc.fits|dvr|dvt", regex=True) | df.provenance_name.isin(["SPOC", "TESS-SPOC", "QLP"])]
    g = df.groupby("provenance_name").sequence_number.apply(lambda s: sorted(set(int(x) for x in s)))
    print(tic, dict(g))
    # data validation products (TCEs) from SPOC / TESS-SPOC
    try:
        prods = Observations.get_product_list(obs[[p in ("SPOC", "TESS-SPOC") for p in obs["provenance_name"]]])
        dv = [str(p) for p in prods["productFilename"] if ("dvr" in str(p) or "dvt" in str(p) or "dvs" in str(p))]
        print("   DV products:", dv[:12] if dv else "none")
    except Exception as e:
        print("   DV lookup failed:", e)
