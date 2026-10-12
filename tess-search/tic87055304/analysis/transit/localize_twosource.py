"""Further localisation checks of the flux-drop source (complements prf_localize.py / localize_checks.py):
 1. Two-source fit of each sector's dF map: dF = a_T PRF_T + a_N PRF_N (positions fixed at Gaia, after the
    PRF/WCS registration offset measured on the out-of-transit image with the target/neighbour flux ratio fixed
    to Gaia); fraction of the flux drop on the target f_T = a_T / (a_T + a_N), errors from bootstrap over transits.
 2. Model-independent PSF: elliptical-free circular Gaussian PSF (width fitted to the out-of-transit image core),
    dF source position fitted with it.
 3. Per-TESS-orbit dF source positions (PRF) for the 200-s sectors.
Writes localize_twosource.json.
"""
import os, json, glob
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.optimize import minimize
import prf_localize as pl
from prf_localize import gaia, TARGET_ID, NB_ID

HERE = os.path.dirname(os.path.abspath(__file__))
tgt = gaia[gaia.source_id == TARGET_ID].iloc[0]; nbr = gaia[gaia.source_id == NB_ID].iloc[0]
rng = np.random.default_rng(3)
H = pl.HALF


def register(M, O, xt, yt, xn, yn):
    """PRF/WCS registration offset from O with all Gaia stars at fixed relative (Gaia) fluxes except a free common
    scale for the target+neighbour pair (ratio fixed), a free scale for the rest and a constant."""
    st = M.stars
    def design(dx, dy):
        pair = np.zeros(O.shape); rest = np.zeros(O.shape)
        for r in st.itertuples():
            im = M.img(r.px + dx, r.py + dy) * 10 ** (-0.4 * (r.Tmag - 10))
            if r.source_id in (TARGET_ID, NB_ID): pair += im
            else: rest += im
        return np.vstack([pair.ravel(), rest.ravel(), np.ones(O.size)]).T
    sl = slice(H - 3, H + 4)
    Oc = O[sl, sl].ravel()
    def chi(p):
        X = design(*p).reshape(O.shape[0], O.shape[1], 3)[sl, sl].reshape(-1, 3)
        c, *_ = np.linalg.lstsq(X, Oc, rcond=None)
        return np.sum((Oc - X @ c) ** 2)
    best = min(((chi((a, b)), a, b) for a in np.arange(-0.5, 0.51, 0.1) for b in np.arange(-0.5, 0.51, 0.1)))
    x0 = np.array([best[1], best[2]])
    r = minimize(chi, x0, method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-3, initial_simplex=np.array([x0, x0 + [0.05, 0], x0 + [0, 0.05]])))
    return r.x


def two_source(M, dF, edF, xt, yt, xn, yn, dx, dy, half=4):
    sl = slice(H - half, H + half + 1)
    Pt = M.img(xt + dx, yt + dy)[sl, sl].ravel(); Pn = M.img(xn + dx, yn + dy)[sl, sl].ravel()
    w = 1 / edF[sl, sl].ravel() ** 2
    X = np.vstack([Pt, Pn]).T
    Xw = X * w[:, None]
    cov = np.linalg.inv(X.T @ Xw); c = cov @ (Xw.T @ dF[sl, sl].ravel())
    chi2 = float(np.sum((dF[sl, sl].ravel() - X @ c) ** 2 * w))
    return c, cov, chi2


def gauss_img(x, y, s, shape):
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    # pixel-integrated circular Gaussian (erf), normalised
    from scipy.special import erf
    def ax(c, g):
        return 0.5 * (erf((g + 0.5 - c) / (np.sqrt(2) * s)) - erf((g - 0.5 - c) / (np.sqrt(2) * s)))
    return ax(x, xx) * ax(y, yy)


def main():
    out = {}
    for fpath in sorted(glob.glob("/home/user/work/tic87/data/tesscut/*.fits")):
        S = pl.load_sector(fpath); M = pl.PRFModel(S)
        xt, yt = map(float, S["w"].world_to_pixel_values(tgt.ra, tgt.dec)); xn, yn = map(float, S["w"].world_to_pixel_values(nbr.ra, nbr.dec))
        dF, edF, O = pl.dF_map(S, S["blocks"])
        dx, dy = register(M, O, xt, yt, xn, yn)
        c, cov, chi2 = two_source(M, dF, edF, xt, yt, xn, yn, dx, dy)
        c0, cov0, chi20 = two_source(M, dF, edF, xt, yt, xn, yn, 0.0, 0.0)
        fb = []
        nb = len(S["blocks"])
        for b in range(300):
            pick = [S["blocks"][i] for i in rng.integers(0, nb, nb)]
            d2, e2, _ = pl.dF_map(S, pick)
            cb, _, _ = two_source(M, d2, e2, xt, yt, xn, yn, dx, dy)
            fb.append(cb[0] / (cb[0] + cb[1]))
        fb = np.array(fb)
        r = dict(registration_px=[float(dx), float(dy)], a_target=float(c[0]), a_target_err=float(np.sqrt(cov[0, 0])),
                 a_neighbour=float(c[1]), a_neighbour_err=float(np.sqrt(cov[1, 1])), corr=float(cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])),
                 f_target=float(c[0] / (c[0] + c[1])), f_target_boot_median=float(np.median(fb)), f_target_boot_sd=float(np.std(fb)),
                 f_target_boot_p97=float(np.percentile(fb, 97.7)), chi2=chi2,
                 f_target_unregistered=float(c0[0] / (c0[0] + c0[1])), a_target_unreg=float(c0[0]), a_neighbour_unreg=float(c0[1]))
        # Gaussian PSF: fit width on O core (target+neighbour pair with Gaia ratio + constant)
        def chiO(p):
            s_ = p[0]
            if not 0.3 < s_ < 2.0: return 1e30
            img = gauss_img(xt + p[1], yt + p[2], s_, S["shape"]) * 10 ** (-0.4 * (tgt.Tmag - 10)) + gauss_img(xn + p[1], yn + p[2], s_, S["shape"]) * 10 ** (-0.4 * (nbr.Tmag - 10))
            sl = slice(H - 1, H + 2)
            X = np.vstack([img[S["ys"], S["xs"]][sl, sl].ravel(), np.ones(9)]).T
            cc, *_ = np.linalg.lstsq(X, O[sl, sl].ravel(), rcond=None)
            return np.sum((O[sl, sl].ravel() - X @ cc) ** 2)
        g = minimize(chiO, [0.8, 0.0, 0.0], method="Nelder-Mead")
        sg, gdx, gdy = g.x
        sl = slice(H - 3, H + 4)
        wgt = 1 / edF[sl, sl].ravel() ** 2
        def chiD(p):
            m = gauss_img(p[0], p[1], sg, S["shape"])[S["ys"], S["xs"]][sl, sl].ravel()
            a = np.sum(m * dF[sl, sl].ravel() * wgt) / np.sum(m * m * wgt)
            return np.sum((dF[sl, sl].ravel() - a * m) ** 2 * wgt)
        gd = minimize(chiD, [xt, yt], method="Nelder-Mead", options=dict(xatol=1e-4))
        r["gauss_sigma_px"] = float(sg); r["gauss_registration_px"] = [float(gdx), float(gdy)]
        r["gauss_source_EN_arcsec"] = list(pl.to_EN(S["w"], gd.x[0] - gdx, gd.x[1] - gdy, xt, yt))
        r["neighbour_EN_arcsec"] = list(pl.to_EN(S["w"], xn, yn, xt, yt))
        # per orbit (200-s sectors): split blocks at the mid-sector gap
        if S["E"] == 200.0:
            tmid = np.median(S["t"])
            gaps = np.where(np.diff(S["t"]) > 0.5)[0]
            tg = S["t"][gaps[np.argmax(np.diff(S["t"])[gaps])]] if len(gaps) else tmid
            orb = {}
            for lab_, bl in (("orbit1", [b for b in S["blocks"] if b[2] <= tg]), ("orbit2", [b for b in S["blocks"] if b[2] > tg])):
                if len(bl) >= 2:
                    d2, e2, _ = pl.dF_map(S, bl)
                    x2, y2, _, _, _ = M.fit_source(d2, e2, xt, yt)
                    orb[lab_] = dict(n=len(bl), EN=list(pl.to_EN(S["w"], x2 - dx, y2 - dy, xt, yt)))
            r["per_orbit"] = orb
        out[S["s"]] = r
        print(S["s"], json.dumps({k: (np.round(v, 3).tolist() if isinstance(v, (list, float)) else v) for k, v in r.items()}), flush=True)
    # combined fraction on target
    a = np.array([out[s]["f_target_boot_median"] for s in out]); e = np.array([out[s]["f_target_boot_sd"] for s in out])
    w = 1 / e ** 2
    out["combined_f_target"] = dict(mean=float(np.sum(a * w) / w.sum()), err=float(1 / np.sqrt(w.sum())), chi2=float(np.sum((a - np.sum(a * w) / w.sum()) ** 2 * w)))
    print(out["combined_f_target"])
    json.dump(out, open(os.path.join(HERE, "localize_twosource.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
