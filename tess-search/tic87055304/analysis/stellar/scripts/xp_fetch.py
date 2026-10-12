"""Fetch Gaia DR3 XP continuous coefficients for (a) the target and (b) a training set of stars on the same
sightline with good parallaxes, for the empirical XP k-NN absolute-magnitude test (xp_knn.py).
Training set: within 1.5 deg, has_xp_continuous, 10.5 < G < 14.0, 0.60 < BP-RP < 1.15, parallax_over_error > 10,
RUWE < 1.4.  Output: out/xp/train_meta.csv, out/xp/xp_coeffs.npz
"""
import os, time, warnings
import numpy as np, pandas as pd
from astroquery.gaia import Gaia
warnings.filterwarnings("ignore")
OUT = "/home/user/work/tic87/stellar/out/xp"; os.makedirs(OUT, exist_ok=True)
TID = 6720986204578421888
Gaia.ROW_LIMIT = -1
fmeta = os.path.join(OUT, "train_meta.csv")
COLS = """source_id, ra, dec, parallax, parallax_error, ruwe, phot_g_mean_mag, phot_g_mean_flux, bp_rp,
           phot_bp_rp_excess_factor, teff_gspphot, logg_gspphot, mh_gspphot, distance_gspphot, azero_gspphot, ag_gspphot"""
COND = """has_xp_continuous='t' and phot_g_mean_mag between 10.5 and 14.0 and bp_rp between 0.60 and 1.15
           and parallax_over_error > 10 and ruwe < 1.4"""
if not os.path.exists(fmeta):
    parts = []
    # synchronous queries are capped at 2000 rows -> split into annuli (each well below the cap)
    for r1, r2 in [(0.0, 1.0), (1.0, 1.4), (1.4, 1.7)]:
        q = f"""select {COLS} from gaiadr3.gaia_source where 1=contains(point('ICRS', ra, dec), circle('ICRS', 272.606477, -43.995466, {r2}))
               and 0=contains(point('ICRS', ra, dec), circle('ICRS', 272.606477, -43.995466, {r1 if r1 > 0 else 1e-6})) and {COND}"""
        for i in range(5):
            try:
                p = Gaia.launch_job(q).get_results().to_pandas(); break
            except Exception as e:
                print("retry", e); time.sleep(30)
        print("annulus", r1, r2, len(p)); assert len(p) < 2000
        parts.append(p)
    p = Gaia.launch_job(f"select {COLS} from gaiadr3.gaia_source where source_id={TID}").get_results().to_pandas()
    parts.append(p)
    meta = pd.concat(parts, ignore_index=True).drop_duplicates("source_id")
    meta.to_csv(fmeta, index=False)
meta = pd.read_csv(fmeta)
print("training+target rows", len(meta), "target present", (meta.source_id == TID).any())
fco = os.path.join(OUT, "xp_coeffs.npz")
ids = list(meta.source_id.astype(np.int64))
have = {}
if os.path.exists(fco):
    z = np.load(fco); have = {int(s): (b, r, eb, er) for s, b, r, eb, er in zip(z["sid"], z["bp"], z["rp"], z["bpe"], z["rpe"])}
todo = [s for s in ids if s not in have]
print("to fetch", len(todo))
for i in range(0, len(todo), 100):
    chunk = todo[i:i + 100]
    for k in range(5):
        try:
            res = Gaia.load_data(ids=chunk, data_release="Gaia DR3", retrieval_type="XP_CONTINUOUS",
                                 data_structure="RAW", format="csv", verbose=False)
            break
        except Exception as e:
            print("retry load_data", e); time.sleep(30); res = None
    if res is None: continue
    for key, tabs in res.items():
        for t in (tabs if isinstance(tabs, list) else [tabs]):
            df = t.to_pandas()
            for _, r in df.iterrows():
                def arr(x):
                    if isinstance(x, str): return np.array([float(v) for v in x.strip("()[]").split(",")])
                    return np.array(x, float)
                try:
                    have[int(r["source_id"])] = (arr(r["bp_coefficients"]), arr(r["rp_coefficients"]),
                                                 arr(r["bp_coefficient_errors"]), arr(r["rp_coefficient_errors"]))
                except Exception as e:
                    print("parse fail", e)
    sid = np.array(list(have.keys()), dtype=np.int64)
    np.savez_compressed(fco, sid=sid, bp=np.array([have[s][0] for s in sid]), rp=np.array([have[s][1] for s in sid]),
                        bpe=np.array([have[s][2] for s in sid]), rpe=np.array([have[s][3] for s in sid]))
    print("fetched", len(have), flush=True)
