"""Collect broadband photometry of TIC 87055304 from VizieR (GALEX, APASS, SkyMapper DR1.1/DR4, Tycho-2, 2MASS,
AllWISE, CatWISE2020, DENIS, Gaia DR3), keeping every source within 15" so that blending can be judged.
Output: out/phot_raw/<cat>.csv and out/photometry_raw_summary.txt
"""
import os, time, warnings
import numpy as np, pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
warnings.filterwarnings("ignore")
OUT = "/home/user/work/tic87/stellar/out/phot_raw"; os.makedirs(OUT, exist_ok=True)
c0 = SkyCoord(272.606557 * u.deg, -43.995371 * u.deg)
cats = {"GALEX_GUVcat": "II/335/gr6xy", "GALEX_DR5": "II/312/ais", "APASS9": "II/336/apass9",
        "SMSS_DR1": "II/358/smss", "SMSS_DR4": "II/379/smssdr4", "Tycho2": "I/259/tyc2",
        "2MASS": "II/246/out", "AllWISE": "II/328/allwise", "CatWISE": "II/365/catwise",
        "DENIS": "B/denis/denis", "UCAC4": "I/322A/out", "GaiaDR3": "I/355/gaiadr3", "VHS": "II/367/vhs_dr5",
        "unWISE": "II/363/unwise", "VPHAS": "II/341/vphasp", "XMM_OM": "II/370/xmmom6s"}
V = Vizier(columns=["**", "_r"], row_limit=-1, timeout=120)
lines = []
for name, cid in cats.items():
    res = None
    for i in range(3):
        try:
            res = V.query_region(c0, radius=15 * u.arcsec, catalog=cid); break
        except Exception as e:
            print(name, "retry", e); time.sleep(10)
    if res is None or len(res) == 0:
        lines.append(f"{name} ({cid}): no source within 15\""); print(lines[-1]); continue
    df = pd.concat([t.to_pandas().assign(_table=t.meta.get("name", "")) for t in res], ignore_index=True, sort=False)
    df = df.sort_values("_r")
    df.to_csv(os.path.join(OUT, f"{name}.csv"), index=False)
    lines.append(f"{name} ({cid}): {len(df)} rows, r = {', '.join(f'{x:.2f}' for x in df['_r'].values[:8])}")
    print(lines[-1])
open("/home/user/work/tic87/stellar/out/photometry_raw_summary.txt", "w").write("\n".join(lines) + "\n")
