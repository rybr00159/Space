import warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from astroquery.gaia import Gaia
from nbr import tmag
T = {214311906: (260.97765, -51.32369, 0.0024, 5925373268528684288), 55829034: (272.03163, -35.69140, 0.0016, 4038825960551504768)}
for tic, (ra, dec, depth, gid) in T.items():
    q = f"""SELECT g.source_id, g.ra, g.dec, g.phot_g_mean_mag AS G, g.bp_rp, g.parallax, g.pmra, g.pmdec, g.ruwe, g.ref_epoch,
           g.phot_variable_flag, g.ipd_frac_multi_peak, DISTANCE(POINT({ra},{dec}), POINT(g.ra,g.dec))*3600 AS sep
           FROM gaiadr3.gaia_source g WHERE CONTAINS(POINT(g.ra,g.dec), CIRCLE({ra},{dec},{150/3600}))=1"""
    g = Gaia.launch_job_async(q).get_results().to_pandas()
    g["T"] = tmag(g.G, g.bp_rp); Tt = float(g.loc[g.source_id == gid, "T"].iloc[0])
    g["dT"] = g["T"] - Tt; g["frac"] = 10 ** (-0.4 * g.dT)
    g["need_depth"] = depth * (1 + g.frac) / g.frac   # intrinsic eclipse depth needed if this star alone is diluted by the target
    g.to_csv(f"out/{tic}/gaia_150arcsec_all.csv", index=False)
    m = g[(g.need_depth < 0.5) & (g.source_id != gid)].sort_values("sep")
    print(f"\nTIC {tic}: {len(g)} Gaia sources within 150\"; target T {Tt:.2f}; {len(m)} could cause the depth with an eclipse < 50%")
    print(m[["source_id", "sep", "G", "bp_rp", "T", "dT", "parallax", "ruwe", "need_depth"]].head(20).round(3).to_string(index=False))
