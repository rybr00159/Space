"""Empirical (variability-based) image of the Mira per sector.

The Mira is the only strongly variable star near the M dwarfs (5-15% within a sector). Regressing every pixel's
time series on the Mira's flux proxy F_ap(t) (background-subtracted sum in r = 2.5 px around it) together with the
background level B(t), the pointing offsets dx(t), dy(t) (p04_jitter, from field stars) and a constant gives, per
pixel, b_j = d(pixel_j)/d(F_ap): the image of the *varying* light, i.e. the Mira's own TESS image (normalised so
that the sum inside the proxy aperture is ~1). Its PRF-fitted position, measured relative to the WCS calibrated on
normal-colour field stars, is the Mira's colour-dependent image offset, free of the static degeneracy with the pair.
Data binned to 1 h to keep it fast; cadences with quality != 0 excluded.
Fits to the b image: (1) PRF at free position (with the p02 blur), (2) free position and extra blur,
(3) free position + extra blur + an additional free-amplitude PRF at the M3V (tests whether the M3V varies with the
proxy, which would indicate leakage of the M3V's light into the proxy, expected small).
Outputs: out/mirapsf_s{sector}.json, data/mirapsf_s{sector}.npz (b image and its error), figs/mirapsf_s{sector}.png
"""
import sys; sys.path.insert(0, "tic214311908/analysis/phot/scripts")
from photlib import *
from p03_scene import build_stars
from scipy.optimize import least_squares
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt


def binned(t, Y, width):
    k = np.floor((t - t[0]) / width).astype(int)
    u, inv = np.unique(k, return_inverse=True)
    n = np.bincount(inv)
    out = np.zeros((len(u),) + Y.shape[1:])
    np.add.at(out, inv, Y)
    return out / n.reshape((-1,) + (1,) * (Y.ndim - 1)), n


def run(sector, width=1 / 24):
    C = Cutout(sector); W = json.load(open(f"{ROOT}/out/wcs_s{sector}.json"))
    J = np.load(f"{ROOT}/data/jitter_s{sector}.npz")
    g, rBA = build_stars(C, W)
    prf = PRFModel(C); prf.blur = W["blur"]
    mira = g[g.name == "Mira"].iloc[0]; m3v = g[g.name == "M3V"].iloc[0]; pair = g[g.name == "pair"].iloc[0]
    good = C.fin & (C.q == 0) & np.isfinite(J["dx"]) & np.isfinite(J["bkg"])
    # drop cadences with extreme background (scattered light) or jitter
    b_ = J["bkg"][good]; lim = np.median(b_) + 5 * robust_std(b_)
    good &= (J["bkg"] < lim)
    for k in ("dx", "dy"):
        v = J[k]; good &= np.abs(v - np.nanmedian(v[good])) < 6 * robust_std(v[good]) + 0.02
    t = C.t[good]; cube = C.flux[good].reshape(good.sum(), -1)
    yy, xx = np.mgrid[0:C.ny, 0:C.nx]
    ap = (np.hypot(xx - mira.x, yy - mira.y) <= 2.5).ravel()
    B = J["bkg"][good]
    Fap = cube[:, ap].sum(1) - B * ap.sum()
    X = np.column_stack([np.ones(len(t)), Fap - Fap.mean(), B - B.mean(), J["dx"][good] - J["dx"][good].mean(), J["dy"][good] - J["dy"][good].mean()])
    Xb, n = binned(t, X, width); Yb, _ = binned(t, cube, width)
    keep = n >= max(2, 0.5 * np.median(n))
    Xb, Yb = Xb[keep], Yb[keep]
    coef, *_ = np.linalg.lstsq(Xb, Yb, rcond=None)
    res = Yb - Xb @ coef
    dof = len(Yb) - Xb.shape[1]
    s2 = (res ** 2).sum(0) / dof
    cov_unit = np.linalg.pinv(Xb.T @ Xb)
    bimg = coef[1].reshape(C.ny, C.nx); berr = np.sqrt(s2 * cov_unit[1, 1]).reshape(C.ny, C.nx)
    # red-noise check: correlation of residuals between neighbouring bins
    frac_var = float(np.std(Fap) / np.mean(Fap))
    # ---- PRF fits to the b image within r = 3.5 px of the Mira's Gaia position
    mask = np.hypot(xx - mira.x, yy - mira.y) <= 3.5
    def make(p, mode):
        bl = prf.blur if mode == 1 else np.hypot(prf.blur, p[2])
        cols = [prf.image(mira.x + p[0], mira.y + p[1], blur=bl)[mask], np.ones(mask.sum())]
        if mode == 3: cols.append(prf.image(m3v.x, m3v.y)[mask]); cols.append(prf.image(pair.x, pair.y)[mask])
        return np.column_stack(cols)
    fits_ = {}
    for mode in (1, 2, 3):
        w = 1 / berr[mask]
        def rr(p):
            A = make(p, mode); c, *_ = np.linalg.lstsq(A * w[:, None], bimg[mask] * w, rcond=None); return (A @ c - bimg[mask]) * w
        best = None
        for gx in np.linspace(-1.2, 1.2, 7):
            for gy in np.linspace(-1.2, 1.2, 7):
                p0 = [gx, gy] + ([0.2] if mode > 1 else [])
                c2 = np.sum(rr(p0) ** 2)
                if best is None or c2 < best[0]: best = (c2, p0)
        q = least_squares(rr, best[1], diff_step=1e-3, x_scale=0.05)
        A = make(q.x, mode); c, *_ = np.linalg.lstsq(A * w[:, None], bimg[mask] * w, rcond=None)
        chi2 = float(np.sum(q.fun ** 2)); nd = mask.sum() - len(q.x) - A.shape[1]
        cov = np.linalg.pinv(q.jac.T @ q.jac) * max(1.0, chi2 / nd)
        Jm = C.jacobian(mira.x, mira.y)
        en = Jm @ q.x[:2]; cen = Jm @ cov[:2, :2] @ Jm.T
        model = np.full(bimg.shape, np.nan); model[mask] = A @ c
        fits_[mode] = dict(p=q.x.tolist(), perr=np.sqrt(np.diag(cov)).tolist(), chi2=chi2, dof=int(nd), amp=c.tolist(),
                           offset_EN_arcsec=en.tolist(), offset_EN_err=np.sqrt(np.diag(cen)).tolist(), model=model)
    # positions of the M dwarfs relative to the Mira's Gaia position in arcsec
    Jm = C.jacobian(mira.x, mira.y)
    rel = {nm: (Jm @ np.array([r.x - mira.x, r.y - mira.y])).tolist() for nm, r in (("M3V", m3v), ("pair", pair))}
    out = dict(sector=sector, n_bins=int(len(Yb)), frac_var_proxy=frac_var, proxy_range=[float(np.percentile(Fap, 2)), float(np.percentile(Fap, 98))],
               fits={k: {kk: vv for kk, vv in v.items() if kk != "model"} for k, v in fits_.items()}, rel_EN=rel)
    save_json(out, f"{ROOT}/out/mirapsf_s{sector}.json")
    np.savez(f"{ROOT}/data/mirapsf_s{sector}.npz", b=bimg, berr=berr, mira_xy=[mira.x, mira.y], m3v_xy=[m3v.x, m3v.y], pair_xy=[pair.x, pair.y])
    fig, ax = plt.subplots(1, 4, figsize=(18, 4.5))
    ax[0].imshow(bimg, origin="lower"); ax[0].set_title(f"S{sector}: d pixel / d F_Mira (var {frac_var*100:.1f}%)")
    for k, mode in enumerate((1, 2, 3)):
        R = (bimg - fits_[mode]["model"]) / berr
        ax[k + 1].imshow(R, origin="lower", cmap="RdBu_r", vmin=-5, vmax=5)
        ax[k + 1].set_title(f"mode {mode}: chi2 {fits_[mode]['chi2']:.0f}/{fits_[mode]['dof']}, off EN {np.round(fits_[mode]['offset_EN_arcsec'],1).tolist()}\"")
    for a in ax:
        for r_, c_ in ((mira, "r"), (m3v, "y"), (pair, "m")): a.plot(r_.x, r_.y, "+", color=c_, ms=10, mew=2)
        p = fits_[2]["p"]; a.plot(mira.x + p[0], mira.y + p[1], "kx", ms=10, mew=2)
        a.set_xlim(mira.x - 5.5, mira.x + 5.5); a.set_ylim(mira.y - 5.5, mira.y + 5.5)
    fig.tight_layout(); fig.savefig(f"{ROOT}/figs/mirapsf_s{sector}.png", dpi=75); plt.close(fig)
    f2 = fits_[2]; f3 = fits_[3]
    print(f"S{sector}: proxy var {frac_var*100:.1f}%, {len(Yb)} bins; Mira image offset (px) {np.round(f2['p'][:2],3).tolist()} +- {np.round(f2['perr'][:2],3).tolist()}"
          f" = EN {np.round(f2['offset_EN_arcsec'],1).tolist()} +- {np.round(f2['offset_EN_err'],1).tolist()}\", extra blur {abs(f2['p'][2]):.3f}; chi2 m1 {fits_[1]['chi2']:.0f} m2 {f2['chi2']:.0f} "
          f"m3 {f3['chi2']:.0f} (dof {f2['dof']}); m3 M3V/pair amp rel to Mira {f3['amp'][2]/f3['amp'][0]:+.3f}/{f3['amp'][3]/f3['amp'][0]:+.3f}; "
          f"M3V at EN {np.round(rel['M3V'],1).tolist()}, pair at {np.round(rel['pair'],1).tolist()}", flush=True)


if __name__ == "__main__":
    for s in ([int(a) for a in sys.argv[1:]] or SECTORS):
        run(s)
