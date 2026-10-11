"""Fetch TIC v8.2 parameters for every TIC ID in a QLP target list (MAST API, threaded batches)."""
import sys, time
from concurrent.futures import ThreadPoolExecutor
import numpy as np, pandas as pd
from tic_mast import fetch

KEEP = ["ID", "Tmag", "Teff", "e_Teff", "logg", "rad", "e_rad", "mass", "rho", "d", "plx", "pmRA", "pmDEC",
        "disposition", "ra", "dec", "GAIA", "GAIAmag", "gaiabp", "gaiarp", "Vmag", "Jmag", "Hmag", "Kmag",
        "contratio", "lumclass", "objType", "splists", "gaiaqflag", "duplicate_id", "numcont"]

src, out, batch = sys.argv[1], sys.argv[2], int(sys.argv[3])
ids = pd.read_csv(src, dtype={"#TIC_ID": str})["#TIC_ID"].astype(np.int64).values
chunks = [ids[i:i + batch] for i in range(0, len(ids), batch)]
t0 = time.time(); parts = []
def job(c):
    df = fetch(c)
    return df[[k for k in KEEP if k in df.columns]]
with ThreadPoolExecutor(4) as ex:
    for n, df in enumerate(ex.map(job, chunks)):
        parts.append(df)
        if n % 50 == 0: print(f"{n}/{len(chunks)} {time.time()-t0:.0f}s", flush=True)
res = pd.concat(parts, ignore_index=True)
res.to_parquet(out); print("done", len(res), "of", len(ids), f"{time.time()-t0:.0f}s")
