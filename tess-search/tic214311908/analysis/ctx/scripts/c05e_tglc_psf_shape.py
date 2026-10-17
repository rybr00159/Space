"""Sensitivity of the TGLC S39 psf_flux leak-model result (c05) to the PSF shape. TGLC fits its own empirical ePSF
(11x11 px, oversampled) to each frame, which can differ from the lab-based TESS PRF (PRF package) used in c05,
especially in the wings that set the cross-talk between stars 0.3-0.9 px apart. Variants: the PRF convolved with a
Gaussian of sigma 0.2/0.35/0.5 px (broader core and wings), and the PRF raised to a power 0.9/1.1 and renormalised
(heavier/lighter wings relative to the core). For each: chi2 of the source on the M3V, the pair, the Mira and the two
faint Gaia stars N of the M3V, and the free-position best fit. Output: out/tglc_psf_shape.json
"""
import os, sys, json, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd
from scipy.ndimage import gaussian_filter
import c05_tglc_leak as L

sector = 39
pos, prf0, J, w = L.positions(sector); sxy = L.tglc_star_xy(sector); Ji = np.linalg.inv(J); names = list(L.STARS)
obs, rnd, lvl = L.observed(sector, "psf_flux")
ov = np.array([obs[k] for k in names]); R = np.array([rnd[k] for k in names]); cov = np.cov(R)
g = pd.read_csv(L.GAIA); dt = 2021.42 - 2016.0
g["ra_ep"] = g.ra + np.nan_to_num(g.pmra) * dt / 3.6e6 / np.cos(np.radians(g.dec)); g["dec_ep"] = g.dec + np.nan_to_num(g.pmdec) * dt / 3.6e6
def xy_of(sid):
    r = g[g.source_id == sid].iloc[0]; return tuple(w.all_world2pix([[r.ra_ep, r.dec_ep]], 0)[0])
hyp = {"M3V": pos["M3V"], "pairA": pos["pairA"], "Mira": pos["Mira"],
       "G17.5_3.1asec": xy_of(5925373272827260544), "G17.7_5.5asec": xy_of(5925373272827262208)}

class Wrapped:
    def __init__(self, base, sigma=0.0, power=1.0):
        self.base, self.sigma, self.power = base, sigma, power
    def locate(self, x, y, shape):
        # evaluate on a padded stamp so the convolution has room, then crop
        pad = 6; big = self.base.locate(x + pad, y + pad, (shape[0] + 2 * pad, shape[1] + 2 * pad))
        tot = big.sum()
        if self.power != 1.0:
            big = np.clip(big, 0, None) ** self.power; big *= tot / big.sum()
        if self.sigma > 0:
            big = gaussian_filter(big, self.sigma, mode="constant")
        return big[pad:pad + shape[0], pad:pad + shape[1]]

res = {}
for tag, s, p in (("baseline", 0.0, 1.0), ("gauss0.20", 0.2, 1.0), ("gauss0.35", 0.35, 1.0), ("gauss0.50", 0.5, 1.0),
                  ("power0.90_heavier_wings", 0.0, 0.9), ("power1.10_lighter_wings", 0.0, 1.1)):
    prf = Wrapped(prf0, s, p)
    def model(xy):
        cp, ca = L.leak_coeffs(pos, prf, sxy, xy); return np.array([cp[k] for k in names])
    out = {}
    for k, xy in hyp.items():
        dF, e, c2 = L.gls(ov, cov, model(xy)); out[k] = dict(chi2=round(c2, 2), dF=round(dF, 2))
    best = (1e9, None)
    for dn in np.arange(-15, 20.1, 1.5):
        for de in np.arange(-35, 5.1, 1.5):
            d = Ji @ np.array([de, dn]); c2 = L.gls(ov, cov, model((pos["Mira"][0] + d[0], pos["Mira"][1] + d[1])))[2]
            if c2 < best[0]: best = (c2, (float(de), float(dn)))
    out["free_best_chi2"] = round(best[0], 2); out["free_best_EN_from_Mira"] = best[1]
    out["coeffs_if_M3V"] = dict(zip(names, np.round(model(hyp["M3V"]), 3).tolist()))
    res[tag] = out
    print(tag, {k: v for k, v in out.items()}, flush=True)
res["observed_dF"] = {k: round(obs[k], 2) for k in names}
res["note"] = "dof = 5 for fixed-position hypotheses (6 stars, 1 amplitude); M3V at (E,N) = (-18.4, +0.1) arcsec from the Mira"
json.dump(res, open("out/tglc_psf_shape.json", "w"), indent=1, default=float)
