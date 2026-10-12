"""Robustness checks of the PRF localisation (prf_localize.py):
 A. flux-weighted centroids (3x3 and 5x5 boxes around the target pixel) of the dF map and of the out-of-transit
    image O, compared with the centroids the TESS PRF predicts for a dF source at the target and at the 9" neighbour
    (same boxes, same pixel grid) -> chi^2 of each hypothesis using bootstrap (over transits) errors;
 B. PRF source position from the dF map for different fit regions (5x5, 7x7, 9x9, 11x11 pixels);
 C. classic difference image (mean out-of-transit within 0.75-1.5 T14 of mid-transit minus mean in-transit within
    0.35 T14) instead of the template dF fit;
 D. dF map from odd and even transits separately (S66+S93+S104 200-s data together are not combined here; per sector).
Writes localize_checks.json.
"""
import os, json, glob
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.optimize import minimize
import prf_localize as pl
from prf_localize import gaia, TARGET_ID, NB_ID

HERE = os.path.dirname(os.path.abspath(__file__))
tgt = gaia[gaia.source_id == TARGET_ID].iloc[0]; nbr = gaia[gaia.source_id == NB_ID].iloc[0]
rng = np.random.default_rng(5)


def centroid(img, box, H):
    n = img.shape[0]; yy, xx = np.mgrid[0:n, 0:n]
    m = (np.abs(xx - H) <= box) & (np.abs(yy - H) <= box)
    w = img * m
    return (w * xx).sum() / w.sum(), (w * yy).sum() / w.sum()


def diff_image(S, blocks, T14):
    t = S["t"]; ins, outs = [], []
    for e, ii, Tc in blocks:
        dt = np.abs(t[ii] - Tc)
        a = ii[dt < 0.35 * T14]; b = ii[(dt > 0.75 * T14) & (dt < 1.5 * T14)]
        if len(a) < 2 or len(b) < 4: continue
        ins.append(S["flux"][a][:, S["ys"], S["xs"]].mean(0)); outs.append(S["flux"][b][:, S["ys"], S["xs"]].mean(0))
    D = np.array(outs) - np.array(ins)
    return D.mean(0), D.std(0) / np.sqrt(len(D)), len(D)


def fit_pos(M, img, err, x0, y0, H, half):
    sl = slice(H - half, H + half + 1)
    wgt = 1 / err[sl, sl].ravel() ** 2
    def chi(p):
        m = M.img(*p)[sl, sl].ravel()
        a = np.sum(m * img[sl, sl].ravel() * wgt) / np.sum(m * m * wgt)
        return np.sum((img[sl, sl].ravel() - a * m) ** 2 * wgt)
    r = minimize(chi, [x0, y0], method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-5))
    return r.x


def main():
    out = {}
    for fpath in sorted(glob.glob("/home/user/work/tic87/data/tesscut/*.fits")):
        S = pl.load_sector(fpath); M = pl.PRFModel(S); H = pl.HALF
        dF, edF, O = pl.dF_map(S, S["blocks"])
        xt, yt = map(float, S["w"].world_to_pixel_values(tgt.ra, tgt.dec)); xn, yn = map(float, S["w"].world_to_pixel_values(nbr.ra, nbr.dec))
        Pt = M.img(xt, yt); Pn = M.img(xn, yn)
        r = dict()
        # A. centroids
        bootc = {1: [], 2: []}
        nb = len(S["blocks"])
        for b in range(200):
            pick = [S["blocks"][i] for i in rng.integers(0, nb, nb)]
            d2, _, _ = pl.dF_map(S, pick)
            for box in (1, 2):
                bootc[box].append(centroid(d2, box, H))
        for box in (1, 2):
            cd = np.array(centroid(dF, box, H)); co = np.array(centroid(O, box, H))
            ct = np.array(centroid(Pt, box, H)); cn = np.array(centroid(Pn, box, H))
            bc = np.array(bootc[box]); cov = np.cov(bc.T)
            ci = np.linalg.inv(cov)
            chi_t = float((cd - ct) @ ci @ (cd - ct)); chi_n = float((cd - cn) @ ci @ (cd - cn))
            r[f"centroid_box{2 * box + 1}"] = dict(dF=cd.tolist(), O=co.tolist(), pred_target=ct.tolist(), pred_neighbour=cn.tolist(),
                                                  err=np.sqrt(np.diag(cov)).tolist(), chi2_target=chi_t, chi2_neighbour=chi_n,
                                                  dF_minus_O_px=(cd - co).tolist())
        # B. region size
        reg = {}
        for half in (2, 3, 4, 5):
            x, y = fit_pos(M, dF, edF, xt, yt, H, half)
            E_, N_ = pl.to_EN(S["w"], x, y, xt, yt)
            reg[f"{2 * half + 1}x{2 * half + 1}"] = [E_, N_]
        r["region_size_EN_arcsec"] = reg
        # C. classic difference image
        D, eD, nd = diff_image(S, S["blocks"], pl.tmpl.T14)
        x, y = fit_pos(M, D, np.where(eD > 0, eD, np.median(eD)), xt, yt, H, 4)
        r["diff_image_EN_arcsec"] = list(pl.to_EN(S["w"], x, y, xt, yt)); r["diff_image_n"] = nd
        # D. odd / even
        for par in (0, 1):
            bl = [b for b in S["blocks"] if (b[0] + 390 - 390) % 2 == par]
            if len(bl) >= 2:
                d2, e2, _ = pl.dF_map(S, bl)
                x, y = fit_pos(M, d2, e2, xt, yt, H, 4)
                r[f"parity{par}_EN_arcsec"] = list(pl.to_EN(S["w"], x, y, xt, yt))
        r["neighbour_EN_arcsec"] = list(pl.to_EN(S["w"], xn, yn, xt, yt))
        out[S["s"]] = r
        print(S["s"], json.dumps({k: (np.round(v, 2).tolist() if isinstance(v, list) else v) for k, v in r.items() if not k.startswith("centroid")}),
              {k: (round(v["chi2_target"], 1), round(v["chi2_neighbour"], 1), np.round(v["dF_minus_O_px"], 3).tolist()) for k, v in r.items() if k.startswith("centroid")}, flush=True)
    json.dump(out, open(os.path.join(HERE, "localize_checks.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
