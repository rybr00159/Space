"""Fetch TIC v8.2 parameters by ID from the MAST catalog API, in batches."""
import sys, json, time
import numpy as np, pandas as pd, requests
URL = "https://mast.stsci.edu/api/v0/invoke"
COLS = "ID,Tmag,Teff,logg,rad,mass,d,disposition,ra,dec,GAIA,GAIAmag,gaiabp,gaiarp,contratio,lumclass,Kmag,objType"

def fetch(ids):
    req = {"service": "Mast.Catalogs.Filtered.Tic", "format": "json", "pagesize": len(ids) + 10, "page": 1,
           "params": {"columns": COLS, "filters": [{"paramName": "ID", "values": [str(i) for i in ids]}]}}
    for attempt in range(5):
        try:
            r = requests.post(URL, data={"request": json.dumps(req)}, timeout=300)
            r.raise_for_status()
            return pd.DataFrame(r.json()["data"])
        except Exception as e:
            print("retry", attempt, e, file=sys.stderr); time.sleep(3 * 2 ** attempt)
    raise RuntimeError("MAST failed")

if __name__ == "__main__":
    t0 = time.time(); ids = pd.read_csv(sys.argv[1], dtype={"#TIC_ID": str})["#TIC_ID"].astype(np.int64).values[:int(sys.argv[2])]
    df = fetch(ids); print(len(df), f"{time.time()-t0:.1f}s"); print(df.head(3).T)
