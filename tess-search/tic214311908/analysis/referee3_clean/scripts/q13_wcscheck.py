"""q13: sanity check of the pixel positions used in q02 (TESScut WCS + Gaia, no shift): flux-weighted centroids of
isolated bright stars (T < 11, > 5 px from the Mira) in the median image vs their WCS positions, per sector.
Also the EE of the q02 apertures if the M3V position is shifted by the measured mean offset.
Output: out/q13_wcscheck.json"""
import json, glob, sys, numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
sys.path.insert(0, "tic214311908/analysis/phot/scripts")
from photlib import Cutout, PRFModel
R3 = "tic214311908/analysis/referee3_clean"
PH = "tic214311908/analysis/phot"
gaia = pd.read_csv(f"{PH}/data/gaia/gaia_360arcsec.csv")
S2 = json.load(open(f"{R3}/out/q02_bigap.json"))
out = {}
for sec in (12, 39, 66, 93, 103, 104):
    fn = glob.glob(f"/home/user/work/tesscut5/214311906/tess-s{sec:04d}-*_astrocut.fits")[0]
    with fits.open(fn) as h:
        F = np.array(h[1].data["FLUX"], float); q = np.array(h[1].data["QUALITY"]); w = WCS(h[2].header); t = h[1].data["TIME"]
    ok = (q == 0) & np.all(np.isfinite(F.reshape(len(F), -1)), 1)
    med = np.median(F[ok], 0); med -= np.percentile(med, 20)
    ep = 2000 + (np.nanmedian(t) + 2457000 - 2451545) / 365.25
    g = gaia.copy(); dt = ep - 2016
    g["x"], g["y"] = w.all_world2pix((g.ra + np.nan_to_num(g.pmra) * dt / 3.6e6 / np.cos(np.radians(g.dec))).values, (g.dec + np.nan_to_num(g.pmdec) * dt / 3.6e6).values, 0)
    mira = g[g.source_id == 5925373268528684288].iloc[0]
    offs = []
    for _, s in g[g["T"] < 11.0].iterrows():
        if not (2.5 < s.x < 17.5 and 2.5 < s.y < 17.5) or np.hypot(s.x - mira.x, s.y - mira.y) < 5: continue
        nb = g[(np.hypot(g.x - s.x, g.y - s.y) < 2.5) & (g["T"] < s["T"] + 2.5) & (g.source_id != s.source_id)]
        if len(nb): continue
        ix, iy = int(round(s.x)), int(round(s.y)); sub = med[iy - 2:iy + 3, ix - 2:ix + 3]
        # compare with PRF-model centroid at the WCS position (the PRF centroid is not exactly the star position)
        C = Cutout(sec); prf = PRFModel(C); im = prf.image(s.x, s.y)[iy - 2:iy + 3, ix - 2:ix + 3]
        k = np.arange(-2, 3)
        cx = (sub * k[None, :]).sum() / sub.sum() - (im * k[None, :]).sum() / im.sum()
        cy = (sub * k[:, None]).sum() / sub.sum() - (im * k[:, None]).sum() / im.sum()
        offs.append((float(cx), float(cy), float(s["T"])))
    out[sec] = dict(n=len(offs), offsets=offs, mean=[float(np.mean([o[0] for o in offs])) if offs else None, float(np.mean([o[1] for o in offs])) if offs else None])
    print(sec, out[sec]["n"], np.round(out[sec]["mean"], 3) if offs else None, [np.round(o[:2], 2).tolist() for o in offs])
json.dump(out, open(f"{R3}/out/q13_wcscheck.json", "w"), indent=1)
