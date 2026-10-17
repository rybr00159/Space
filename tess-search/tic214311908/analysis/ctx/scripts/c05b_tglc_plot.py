"""Plot the raw TGLC psf_flux light curves (S12, S39) of the six stars, with transit times marked (diagnostic)."""
import glob, sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from astropy.io import fits
import common as C
STARS = {"M3V": 5925373272850784768, "pairA": 5925373272850782208, "pairB": 5925373272827253632,
         "Mira": 5925373268528684288, "ref1": 5925373204131295872, "ref2": 5925373204131294336}
fig, ax = plt.subplots(6, 2, figsize=(14, 14), sharex="col")
for j, sec in enumerate((12, 39)):
    for i, (k, sid) in enumerate(STARS.items()):
        d = fits.getdata(glob.glob(f"lc/tglc/*{sid}-s{sec:04d}*")[0], 1)
        ok = (d["TGLC_flags"] == 0) & ((d["TESS_flags"] & ~2048) == 0)
        a = ax[i, j]; a.plot(d["time"][ok], d["psf_flux"][ok], ",k"); a.plot(d["time"][~ok], d["psf_flux"][~ok], ",r")
        a.set_ylabel(f"{k} S{sec}"); lo, hi = np.nanpercentile(d["psf_flux"][ok], [0.5, 99.5]); a.set_ylim(lo - 0.2 * (hi - lo), hi + 0.2 * (hi - lo))
        for e in range(-400, 0):
            tc = C.T0 + e * C.P
            if d["time"].min() < tc < d["time"].max(): a.axvline(tc, color="b", alpha=0.3)
        print(sec, k, "n", ok.sum(), "neg", np.sum(d['psf_flux'][ok] <= 0), "median", np.nanmedian(d["psf_flux"][ok]), "p1,p99", np.nanpercentile(d["psf_flux"][ok], [1, 99]))
plt.tight_layout(); plt.savefig("fig/tglc_raw_psf_flux.png", dpi=80)
