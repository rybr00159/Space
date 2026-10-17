"""s04b: like s04, but saves the full 21x21 difference-image cube D (float32) per sector, plus all Gaia stars
(T < 19) with their cutout pixel positions, for full-cutout searches of periodic signals (s05b).
Output: data/pixfull_s{sector}.npz (large; deleted at the end of the analysis, see README).
"""
import sys
sys.path.insert(0, "tic214311908/analysis/phot/scripts")
import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
from photlib import *
from p05_phot import setup, photometry
OUT = "tic214311908/analysis/search/data"
for s in ([int(a) for a in sys.argv[1:]] or SECTORS):
    S = setup(s); C, W, J, g, prf, rBA = S
    R = photometry(s, S=S, verbose=False, return_D=True)
    L = np.load(f"{ROOT}/lc_allstars_s{s:03d}.npz")
    gg = g[g["T"] < 19.0]
    np.savez_compressed(f"{OUT}/pixfull_s{s}.npz", t_tdb=L["t_tdb"], ok=L["quality_ok"], intr=L["intransit_window"],
                        D=R["D"].astype(np.float32), b=R["b_img"], exptime=R["exptime"],
                        star_id=gg.source_id.values, star_x=gg.x.values, star_y=gg.y.values, star_T=gg["T"].values,
                        star_name=gg.name.values.astype(str), jac=C.jacobian(*R["pos"]["M3V"]),
                        m3v_xy=np.array(R["pos"]["M3V"]), mira_xy=np.array(R["pos"]["Mira"]), pair_xy=np.array(R["pos"]["pair"]))
    print(s, "saved", flush=True)
