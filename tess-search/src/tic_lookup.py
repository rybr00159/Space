"""Fetch TIC v8.2 parameters for a list of TIC IDs via VizieR TAP upload (IV/39/tic82)."""
import sys, io, time
import numpy as np, pandas as pd, requests
from astropy.table import Table

TAP = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
COLS = 'c."TIC", c."Tmag", c."Teff", c."logg", c."Rad", c."Mass", c."Dist", c."Disp", c."RAJ2000", c."DEJ2000", c."GAIA", c."Gmag", c."Bpmag", c."Rpmag", c."Contrat", c."Lclass", c."Kmag"'

def lookup(tics, tmag_max=None):
    t = Table({"tic": np.asarray(tics, dtype=np.int64)})
    buf = io.BytesIO(); t.write(buf, format="votable"); buf.seek(0)
    where = f' WHERE c."Tmag" <= {tmag_max}' if tmag_max else ""
    q = f'SELECT {COLS} FROM "IV/39/tic82" AS c JOIN TAP_UPLOAD.t AS u ON c."TIC" = u.tic{where}'
    for attempt in range(4):
        try:
            r = requests.post(TAP, data={"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv",
                                         "MAXREC": 2000000, "QUERY": q, "UPLOAD": "t,param:t"},
                              files={"t": ("t.xml", buf.getvalue())}, timeout=900)
            r.raise_for_status()
            return pd.read_csv(io.StringIO(r.text))
        except Exception as e:
            print("retry", attempt, e, file=sys.stderr); time.sleep(5 * 2 ** attempt)
    raise RuntimeError("TAP failed")

if __name__ == "__main__":
    src, out, tmax, chunk = sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4])
    ids = pd.read_csv(src, dtype={"#TIC_ID": str})["#TIC_ID"].astype(np.int64).values
    parts = []
    for i in range(0, len(ids), chunk):
        t0 = time.time(); df = lookup(ids[i:i + chunk], tmax); parts.append(df)
        print(f"{i + chunk}/{len(ids)} rows={len(df)} {time.time() - t0:.0f}s", flush=True)
    pd.concat(parts).to_parquet(out)
