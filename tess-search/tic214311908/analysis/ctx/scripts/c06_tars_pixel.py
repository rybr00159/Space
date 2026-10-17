"""TARS (TESS All-sky Rotation Survey, Boyle et al.; 1-pixel aperture on TIC 214311908) depths converted to flux lost (e-/s)
using the TESScut flux of the pixel containing the M3V, and compared with the PRF-predicted fraction of a source on the
M3V / pairA / Mira falling in that pixel. Mira's brightness in that pixel changes between sectors, so a constant dF
argues for a constant-brightness host. Output: out/tars_pixel.json
"""
import os, sys, glob, json, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
import PRF
PRFDIR = "<local>/tic713/localise/prf_files/s0004/"
GAIA = "/home/user/Space/tess-search/vetting_2026-10-16/results/214311906/gaia_150arcsec_all.csv"
fold = pd.read_csv("out/fold_products.csv")
out = {}
for sec in (12, 39, 66, 93):
    f = glob.glob(f"/home/user/work/tesscut5/214311906/tess-s{sec:04d}-*_astrocut.fits")[0]
    with fits.open(f) as h:
        w = WCS(h[2].header); d = h[1].data; h0 = h[0].header; h1 = h[1].header
        ok = (d["QUALITY"] == 0) & np.isfinite(d["TIME"])
        tmid = np.nanmedian(d["TIME"][ok])
        g = pd.read_csv(GAIA); dt = 2000.0 + (tmid + 2457000 - 2451545) / 365.25 - 2016.0
        pos = {}
        for k, sid in (("M3V", 5925373272850784768), ("pairA", 5925373272850782208), ("Mira", 5925373268528684288)):
            r = g[g.source_id == sid].iloc[0]
            pos[k] = w.all_world2pix([[r.ra + np.nan_to_num(r.pmra) * dt / 3.6e6 / np.cos(np.radians(r.dec)), r.dec + np.nan_to_num(r.pmdec) * dt / 3.6e6]], 0)[0]
        ix, iy = int(round(pos["M3V"][0])), int(round(pos["M3V"][1]))
        pix = np.nanmedian(d["FLUX"][ok][:, iy, ix])
    prf = PRF.TESS_PRF(int(h0["CAMERA"]), int(h0["CCD"]), sec, int(h1["1CRV4P"]) + ix, int(h1["2CRV4P"]) + iy, localdatadir=PRFDIR)
    frac = {k: float(prf.locate(p[0], p[1], (21, 21))[iy, ix]) for k, p in pos.items()}
    row = fold[(fold["product"] == "TARS") & (fold.sector == str(sec))].iloc[0]
    dF = row.depth * pix; edF = row.err_emp_std * pix
    out[f"S{sec}"] = dict(pixel=[ix, iy], pixel_median_flux=float(pix), tars_depth=row.depth, tars_err_emp=row.err_emp_std,
                         dF_pixel=dF, dF_pixel_err=edF, prf_frac_in_pixel=frac,
                         implied_total_dF_if={k: dF / v for k, v in frac.items()},
                         implied_total_dF_err_if={k: edF / v for k, v in frac.items()})
    print(sec, "pixel", ix, iy, "flux %.0f" % pix, "depth %.0f ppm" % (row.depth * 1e6), "dF_pix %.2f+-%.2f" % (dF, edF),
          "PRF frac", {k: round(v, 3) for k, v in frac.items()}, "implied total dF", {k: round(dF / v, 1) for k, v in frac.items()})
json.dump(out, open("out/tars_pixel.json", "w"), indent=1, default=float)
