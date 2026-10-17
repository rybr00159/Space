"""s04: extract and cache, per sector, the difference-image cube D of the phot pipeline (p05, read-only import) for the
pixels around the M3V/Mira, plus the per-cadence quantities needed to analyse periodic signals at the pixel level.

D(t) = pixel(t) - a - c dB - d dx - e dy (static scene, background and pointing response removed; keeps the Mira's
variation and every star's short-term variability). Saved: data/pix_s{sector}.npz with
  t_tdb (target-barycentred, from phot/lc_allstars), ok (phot quality_ok), intr (5.69-d windows), D (nt x 11 x 11 crop
  centred on the M3V), b (Mira empirical image, crop), F_mira_B / F_M3V_B / F_pair_A / F_M3V_A (e-/s), x0, y0 (crop
  origin), positions of the M3V, pair, Mira and the two faint neighbours (crop pixel coordinates), PRF images of each
  (crop), and the WCS Jacobian (arcsec per pixel) at the M3V.
Nothing is written into ../phot.
"""
import sys
sys.path.insert(0, "tic214311908/analysis/phot/scripts")
import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
from photlib import *
from p05_phot import setup, photometry

OUT = "tic214311908/analysis/search/data"
NEIGH = dict(bg3=5925373272827260544, bg5=5925373272827262208)
HALF = 5


def main(sectors):
    for s in sectors:
        S = setup(s); C, W, J, g, prf, rBA = S
        R = photometry(s, S=S, verbose=False, return_D=True)
        L = np.load(f"{ROOT}/lc_allstars_s{s:03d}.npz")
        assert len(L["t_tess"]) == len(R["t"])
        m3 = R["pos"]["M3V"]
        xc, yc = int(round(m3[0])), int(round(m3[1]))
        x0, y0 = xc - HALF, yc - HALF
        sl = (slice(None), slice(y0, y0 + 2 * HALF + 1), slice(x0, x0 + 2 * HALF + 1))
        D = R["D"][sl].astype(np.float32)
        pos = {k: (v[0] - x0, v[1] - y0) for k, v in R["pos"].items()}
        for nm, sid in NEIGH.items():
            r = g[g.source_id == sid].iloc[0]
            pos[nm] = (float(r.x) - x0, float(r.y) - y0)
        prfs = {}
        for k, (x, y) in pos.items():
            if k == "Mira": continue
            prfs[k] = prf.image(x + x0, y + y0)[sl[1:]]
        jac = C.jacobian(*m3)
        nB = R["FIT"]["B"]["names"]; nA = R["FIT"]["A"]["names"]
        np.savez_compressed(f"{OUT}/pix_s{s}.npz", t_tdb=L["t_tdb"], ok=L["quality_ok"], intr=L["intransit_window"], D=D,
                            b=R["b_img"][sl[1:]], F_mira_B=R["FIT"]["B"]["F"][:, nB.index("Mira")],
                            F_M3V_B=R["FIT"]["B"]["F"][:, nB.index("M3V")], F_M3V_A=R["FIT"]["A"]["F"][:, nA.index("M3V")],
                            F_pair_A=R["FIT"]["A"]["F"][:, nA.index("pair")], F_mira_A=R["FIT"]["A"]["F"][:, nA.index("Mira")],
                            x0=x0, y0=y0, pos_names=np.array(list(pos.keys())), pos_xy=np.array(list(pos.values())),
                            prf_names=np.array(list(prfs.keys())), prf_imgs=np.array(list(prfs.values())), jac=jac,
                            bkg=R["bkg"], dx=R["dx"], dy=R["dy"], exptime=R["exptime"], sig=R["sig"])
        print(f"S{s}: crop x0 {x0} y0 {y0}; nt {len(R['t'])}; positions " + ", ".join(f"{k} ({v[0]:.2f},{v[1]:.2f})" for k, v in pos.items()), flush=True)
        del R, D


if __name__ == "__main__":
    main([int(a) for a in sys.argv[1:]] or list(SECTORS))
