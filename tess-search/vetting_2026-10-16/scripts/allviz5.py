import sys, re, pickle, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
T = {90086266: (276.27770, -42.64403, 2.339469), 61332759: (287.61974, -44.73357, 0.740299),
     214311906: (260.97765, -51.32369, 5.69149), 55829034: (272.03163, -35.69140, 0.616146), 61141412: (286.96474, -43.00969, 6.5285)}
v = Vizier(columns=["*", "_r"], row_limit=200, timeout=300)
for tic, (ra, dec, P) in T.items():
    res = v.query_region(SkyCoord(ra * u.deg, dec * u.deg), radius=30 * u.arcsec)
    out = {k: (res[k].meta.get("description", ""), res[k].to_pandas()) for k in res.keys()}
    pickle.dump(out, open(f"out/{tic}/allviz30.pkl", "wb"))
    print(f"\n=== TIC {tic} (P {P}): {len(out)} tables within 30\"")
    for k, (desc, df) in out.items():
        if not re.search(r"variab|eclips|binar|period|rotat|ellipsoid|transit|planet|TESS|candidate", desc, re.I): continue
        # print any period-like column values and whether they match P, 2P, P/2
        pcols = [c for c in df.columns if re.fullmatch(r"(?i)(per|period|p|porb|prot|p_?rot|p_?orb|perIn|Pday|P1|Per1|period_?\w*)", c)]
        rmin = df["_r"].min() if "_r" in df else np.nan
        vals = []
        for c in pcols:
            for x in pd.to_numeric(df[c], errors="coerce").dropna().values[:5]:
                tag = ""
                for kk, lab in ((1, "=P"), (2, "=2P"), (0.5, "=P/2"), (3, "=3P"), (1/3, "=P/3")):
                    if abs(x / P - kk) < 0.01 * kk: tag = lab
                vals.append(f"{c}={x:.5g}{tag}")
        print(f"  {k:28s} r={rmin:5.1f}\" n={len(df):3d} {desc[:70]:70s} {' '.join(vals[:6])}")
