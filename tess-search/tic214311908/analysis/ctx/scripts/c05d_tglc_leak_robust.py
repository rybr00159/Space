"""Robustness of the TGLC S39 psf_flux leak-model result (c05/c05c) to modelling choices:
  V2  the M3V's TESS image displaced by the TESS lateral-colour effect (2-4" away from the S39 camera-1 boresight,
      PA 165 deg; 713 campaign: ~1.5"/mag of BP-RP, steeper for the reddest stars);
  V3  emulation of TGLC's edge-pixel outlier mask: the Mira's TESS image is displaced ~16" (PA 165) from the Gaia
      position where TGLC subtracts it, leaving a static dipole; edge pixels where that residual deviates by > 1 sd
      from the edge median are masked as TGLC does;
  V4  WCS/PRF jitter: all star/source positions shifted together by +-0.15 px relative to TGLC's 5x5 cut grids.
For each variant: chi2 of the source-on-star hypotheses (M3V, pairA, Mira, the two faint stars 3" and 5.5" N of the
M3V) and the free-position best fit. Output: out/tglc_leak_robust.json
"""
import os, sys, json, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd
import c05_tglc_leak as L

sector = 39
pos, prf, J, w = L.positions(sector); sxy = L.tglc_star_xy(sector); Ji = np.linalg.inv(J); names = list(L.STARS)
xm, ym = pos["Mira"]
obs, rnd, lvl = L.observed(sector, "psf_flux")
ov = np.array([obs[k] for k in names]); R = np.array([rnd[k] for k in names]); cov = np.cov(R)
g = pd.read_csv(L.GAIA); dt = 2021.42 - 2016.0
g["ra_ep"] = g.ra + np.nan_to_num(g.pmra) * dt / 3.6e6 / np.cos(np.radians(g.dec)); g["dec_ep"] = g.dec + np.nan_to_num(g.pmdec) * dt / 3.6e6
def xy_of(sid):
    r = g[g.source_id == sid].iloc[0]; return tuple(w.all_world2pix([[r.ra_ep, r.dec_ep]], 0)[0])
faint = {"G17.5_3.1asec": xy_of(5925373272827260544), "G17.7_5.5asec": xy_of(5925373272827262208)}
pa = np.radians(165.0); unit_EN = np.array([np.sin(pa), np.cos(pa)])   # direction away from the S39 boresight
def shift_xy(xy, arcsec):
    d = Ji @ (unit_EN * arcsec); return (xy[0] + d[0], xy[1] + d[1])

def coeffs(src, pos_=pos, mask=None):
    yy, xx = np.mgrid[0:5, 0:5] - 2.0; c = {}
    for k, (x, y) in pos_.items():
        ox, oy = pos[k][0] - sxy[k][0], pos[k][1] - sxy[k][1]     # cut grid fixed by TGLC's own placement
        Rr = prf.locate(src[0] - ox + 10, src[1] - oy + 10, (25, 25))[10:15, 10:15]
        big = prf.locate(x - ox + 10, y - oy + 10, (25, 25)); PX = big[10:15, 10:15]; PX11 = big[7:18, 7:18].sum()
        A = np.column_stack([PX.ravel() / PX11, xx.ravel(), yy.ravel(), np.ones(25)]); keep = np.ones(25, bool)
        if mask is not None: keep = ~mask[k]
        c[k] = float(np.linalg.lstsq(A[keep], Rr.ravel()[keep], rcond=None)[0][0])
    return np.array([c[k] for k in names])

def evaluate(tag, pos_=pos, mask=None, m3v_shift=0.0):
    hyp = {"M3V": shift_xy(pos_["M3V"], m3v_shift), "pairA": pos_["pairA"], "Mira_gaia": pos_["Mira"],
           "Mira_shift16": shift_xy(pos_["Mira"], 16.0), **{k: (v[0] + pos_["Mira"][0] - pos["Mira"][0], v[1] + pos_["Mira"][1] - pos["Mira"][1]) for k, v in faint.items()}}
    out = {k: L.gls(ov, cov, coeffs(v, pos_, mask))[2] for k, v in hyp.items()}
    grid = np.arange(-40, 20.1, 1.5); best = (1e9, None)
    for dn in np.arange(-25, 25.1, 1.5):
        for de in grid:
            d = Ji @ np.array([de, dn]); c2 = L.gls(ov, cov, coeffs((pos_["Mira"][0] + d[0], pos_["Mira"][1] + d[1]), pos_, mask))[2]
            if c2 < best[0]: best = (c2, (float(de), float(dn)))
    out["free_best_chi2"] = best[0]; out["free_best_EN"] = best[1]
    print(tag, {k: (round(v, 2) if isinstance(v, float) else v) for k, v in out.items()}, flush=True)
    return out

res = {"V1_baseline": evaluate("V1")}
for s in (2.0, 3.0, 4.0):
    res[f"V2_M3V_colour_shift_{s}"] = evaluate(f"V2 shift {s}", m3v_shift=s)
# V3 edge mask emulation
edge = np.array([0, 1, 2, 3, 4, 5, 9, 10, 14, 15, 19, 20, 21, 22, 23, 24])
for ratio in (1.0, 2.0, 4.0):
    mask = {}
    mshift = shift_xy(pos["Mira"], 16.0)
    for k, (x, y) in pos.items():
        ox, oy = x - sxy[k][0], y - sxy[k][1]
        true = ratio * prf.locate(mshift[0] - ox + 10, mshift[1] - oy + 10, (25, 25))[10:15, 10:15]
        model = prf.locate(pos["Mira"][0] - ox + 10, pos["Mira"][1] - oy + 10, (25, 25))[10:15, 10:15]
        r = (true - model).ravel()[edge]; m = np.zeros(25, bool); m[edge[np.abs(r - np.median(r)) > np.std(r)]] = True
        mask[k] = m
    res[f"V3_edge_mask_ratio{ratio}"] = evaluate(f"V3 ratio {ratio}", mask=mask)
for dx, dy in ((0.15, 0), (-0.15, 0), (0, 0.15), (0, -0.15)):
    pos_ = {k: (v[0] + dx, v[1] + dy) for k, v in pos.items()}
    res[f"V4_jitter_{dx}_{dy}"] = evaluate(f"V4 {dx},{dy}", pos_=pos_)
json.dump(res, open("out/tglc_leak_robust.json", "w"), indent=1, default=float)
