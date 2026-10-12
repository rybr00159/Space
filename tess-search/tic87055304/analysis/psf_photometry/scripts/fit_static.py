"""Static PRF model of each sector's median TESScut image.

Model: sum over Gaia DR3 stars (G<19, positions propagated to the sector epoch, WCS of the cutout) of flux * PRF(x+dx0, y+dy0)
+ planar background.  PRF = TESS_PRF smeared by a Gaussian of width sig (fitted).  Free fluxes for stars with T<15.5
inside the stamp (weak ridge prior, 30%, towards the Gaia-predicted flux; NO prior on the target or the neighbour); fainter
stars share one flux zero-point.  Global registration (dx0, dy0) and sig are fitted on pixels > R_EXCL px from the
target, so that the registration is independent of the target/neighbour pair.  Then:
  (1) the pair fluxes are fitted with the registration fixed -> F_nb/F_tgt compared with the Gaia/TIC prediction;
  (2) the pair-only registration is fitted with F_nb/F_tgt fixed to the Gaia prediction -> compared with the field one;
  (3) bright isolated field stars are fitted with free positions -> scatter of (measured - Gaia) = local astrometric
      registration error at a random location; this sets the systematic error budget for the target/neighbour split.
Output: static_s{sector}.npz/json in ../static/, PNG residual maps.
"""
import os, sys, json
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.optimize import least_squares
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import *

OUT = os.path.join(ROOT, "static"); os.makedirs(OUT, exist_ok=True)
R_EXCL = 3.0
TFREE = 15.5


def zp_flux(T):
    return 10 ** (-0.4 * (T - 10.0))   # relative units; scaled by fitted zero points


def build(g, prf, dx, dy, free):
    shape = (21, 21)
    imgs = prf.imgs(g.x.values + dx, g.y.values + dy, shape)          # (n, 21, 21)
    A_free = imgs[free].reshape(free.sum(), -1).T
    faint = (imgs[~free] * zp_flux(g["T"].values[~free])[:, None, None]).sum(0).ravel()
    yy, xx = np.mgrid[0:21, 0:21]
    bk = np.vstack([np.ones(441), (xx.ravel() - 10) / 10, (yy.ravel() - 10) / 10]).T
    return np.hstack([A_free, faint[:, None], bk]), imgs


def solve(A, b, w, prior_mu, prior_sig):
    """weighted LSQ with independent Gaussian priors (sig=inf -> none)."""
    Aw = A * w[:, None]; bw = b * w
    m = np.isfinite(prior_sig)
    if m.any():
        P = np.zeros((m.sum(), A.shape[1])); P[np.arange(m.sum()), np.where(m)[0]] = 1 / prior_sig[m]
        Aw = np.vstack([Aw, P]); bw = np.concatenate([bw, prior_mu[m] / prior_sig[m]])
    x, *_ = np.linalg.lstsq(Aw, bw, rcond=None)
    return x


def static_fit(s, sig0=0.0, fit_sig=True, verbose=True):
    sec = load_sector(s)
    g = gaia_stars(sec)
    prf = PRFModel(sec, sig_px=sig0, xt=g.x[g.Source == TGT].values[0], yt=g.y[g.Source == TGT].values[0])
    M = np.median(sec["F"], axis=0)
    it = np.where(g.Source == TGT)[0][0]; ib = np.where(g.Source == NB)[0][0]
    xt, yt = g.x[it], g.y[it]
    inside = (g.x > -1.5) & (g.x < 21.5) & (g.y > -1.5) & (g.y < 21.5)
    free = np.array((inside & (g["T"] < TFREE)).values, copy=True); free[it] = free[ib] = True
    ifree = np.where(free)[0]
    yy, xx = np.mgrid[0:21, 0:21]
    r = np.hypot(xx - xt, yy - yt).ravel()
    sig_pix = np.sqrt(np.clip(M.ravel(), 1, None) + (0.02 * M.ravel()) ** 2)   # Poisson-like + 2% model error
    w = 1 / sig_pix
    field = r > R_EXCL
    # prior: Gaia-predicted flux * global zero point; zp from a quick first pass
    nfree = free.sum()

    def lin(dx, dy, sig, mask, zp_guess=None):
        if sig != prf.sig: prf.set_sigma(sig)
        A, imgs = build(g, prf, dx, dy, free)
        mu = np.zeros(A.shape[1]); ps = np.full(A.shape[1], np.inf)
        if zp_guess is not None:
            mu[:nfree] = zp_guess * zp_flux(g["T"].values[ifree]); ps[:nfree] = 0.3 * mu[:nfree]
            ps[list(ifree).index(it)] = np.inf; ps[list(ifree).index(ib)] = np.inf
        x = solve(A[mask], M.ravel()[mask], w[mask], mu, ps)
        return x, A, imgs

    # zero point guess from bright free stars with no priors, field pixels
    x, A, _ = lin(0, 0, sig0, np.ones(441, bool))
    br = (g["T"].values[ifree] < 13.5)
    zp = np.median(x[:nfree][br] / zp_flux(g["T"].values[ifree][br]))

    def resid(p, mask):
        dx, dy, sig = p if fit_sig else (p[0], p[1], sig0)
        x, A, _ = lin(dx, dy, max(sig, 0.0), mask, zp)
        return ((M.ravel() - A @ x) * w)[mask]

    p0 = [0.0, 0.0, max(sig0, 0.2)] if fit_sig else [0.0, 0.0]
    lo = [-1, -1, 0.0] if fit_sig else [-1, -1]; hi = [1, 1, 1.5] if fit_sig else [1, 1]
    res_field = least_squares(resid, p0, args=(field,), bounds=(lo, hi), diff_step=1e-3)
    pf = list(res_field.x) + ([] if fit_sig else [sig0])
    J = res_field.jac; cov = np.linalg.pinv(J.T @ J) * (2 * res_field.cost / max(1, field.sum() - len(p0)))
    perr = np.sqrt(np.diag(cov))
    res_all = least_squares(resid, p0, args=(np.ones(441, bool),), bounds=(lo, hi), diff_step=1e-3)
    pa = list(res_all.x) + ([] if fit_sig else [sig0])

    # (1) pair fluxes with field registration, fit on all pixels
    x, A, imgs = lin(pf[0], pf[1], pf[2], np.ones(441, bool), zp)
    Ft, Fn = x[list(ifree).index(it)], x[list(ifree).index(ib)]
    model = (A @ x).reshape(21, 21)
    resid_img = M - model
    # formal covariance of the pair fluxes from pixel noise in the static image is negligible; the systematic error is the
    # registration uncertainty: compute dF/d(offset) numerically
    def pair_ratio(dx, dy, sig):
        x_, _, _ = lin(dx, dy, sig, np.ones(441, bool), zp)
        return x_[list(ifree).index(ib)] / x_[list(ifree).index(it)], x_[list(ifree).index(it)], x_[list(ifree).index(ib)]

    # direction of the neighbour (unit vector) and the derivative of the ratio along it
    u = np.array([g.x[ib] - xt, g.y[ib] - yt]); sep = np.hypot(*u); u /= sep
    rp = pair_ratio(pf[0] + 0.05 * u[0], pf[1] + 0.05 * u[1], pf[2])[0]
    rm = pair_ratio(pf[0] - 0.05 * u[0], pf[1] - 0.05 * u[1], pf[2])[0]
    dratio_ds = (rp - rm) / 0.1
    ratio_gaia = 10 ** (-0.4 * (g["T"][ib] - g["T"][it]))

    # (2) pair-only registration with ratio fixed to Gaia: replace the two free pair fluxes by one combined template
    def resid_pair(p):
        dx, dy = p
        A, imgs_ = build(g, prf, pf[0], pf[1], free)   # field stars at field registration
        # recompute pair columns at the pair-specific offset
        pt = prf.img(xt + dx, yt + dy); pn = prf.img(g.x[ib] + dx, g.y[ib] + dy)
        cols = [list(ifree).index(it), list(ifree).index(ib)]
        keep = [c for c in range(A.shape[1]) if c not in cols]
        A2 = np.hstack([A[:, keep], (pt + ratio_gaia * pn).ravel()[:, None]])
        mu = np.zeros(A2.shape[1]); ps = np.full(A2.shape[1], np.inf)
        kk = [k for k in range(nfree) if k not in cols]
        mu[:len(kk)] = zp * zp_flux(g["T"].values[ifree][kk]); ps[:len(kk)] = 0.3 * mu[:len(kk)]
        xx_ = solve(A2, M.ravel(), w, mu, ps)
        return ((M.ravel() - A2 @ xx_) * w)[r < 2.5]
    rpair = least_squares(resid_pair, [pf[0], pf[1]], diff_step=1e-3)
    Jp = rpair.jac; covp = np.linalg.pinv(Jp.T @ Jp) * (2 * rpair.cost / max(1, (r < 2.5).sum() - 2))

    # (3) free-position fits of bright, isolated field stars (> R_EXCL+1 px from target, no star within 2 mag and 1.5 px)
    star_off = []
    Tv = g["T"].values
    for k in np.where(free & (g["T"].values < 13.8))[0]:
        if k in (it, ib): continue
        if not (1.5 < g.x[k] < 19.5 and 1.5 < g.y[k] < 19.5): continue
        if np.hypot(g.x[k] - xt, g.y[k] - yt) < R_EXCL + 1: continue
        d = np.hypot(g.x - g.x[k], g.y - g.y[k])
        if ((d < 1.5) & (d > 0) & (Tv < Tv[k] + 2.0)).any(): continue
        # model of everything else at the field registration
        others = model - x[list(ifree).index(k)] * imgs[k]
        stamp = (np.hypot(xx - g.x[k] - pf[0], yy - g.y[k] - pf[1]) < 2.5).ravel()
        def rs(p):
            im = prf.img(g.x[k] + pf[0] + p[0], g.y[k] + pf[1] + p[1])
            AA = np.vstack([im.ravel(), np.ones(441)]).T[stamp]
            bb = (M - others).ravel()[stamp] + x[-3]  # add back constant bkg so the column of ones absorbs it
            cc, *_ = np.linalg.lstsq(AA * w[stamp, None], bb * w[stamp], rcond=None)
            return (bb - AA @ cc) * w[stamp]
        rr = least_squares(rs, [0.0, 0.0], diff_step=1e-3)
        star_off.append(dict(source=int(g.Source[k]), T=float(Tv[k]), x=float(g.x[k]), y=float(g.y[k]),
                             ddx=float(rr.x[0]), ddy=float(rr.x[1])))
    so = np.array([[d["ddx"], d["ddy"]] for d in star_off]) if star_off else np.zeros((0, 2))
    # projection onto the neighbour direction
    proj = so @ u if len(so) else np.array([])

    out = dict(sector=s, n_cad=len(sec["t"]), sig_jit=pf[2], dx0=pf[0], dy0=pf[1], dx0_err=float(perr[0]), dy0_err=float(perr[1]),
               dx0_allpix=pa[0], dy0_allpix=pa[1], sig_allpix=pa[2],
               F_tgt=float(Ft), F_nb=float(Fn), ratio_fit=float(Fn / Ft), ratio_gaia=float(ratio_gaia),
               dratio_dshift_along_nb=float(dratio_ds), sep_px=float(sep), nb_unit=u.tolist(),
               pair_offset_with_gaia_ratio=[float(rpair.x[0]), float(rpair.x[1])], pair_offset_err=np.sqrt(np.diag(covp)).tolist(),
               pair_minus_field_along_nb=float((rpair.x[0] - pf[0]) * u[0] + (rpair.x[1] - pf[1]) * u[1]),
               field_star_offsets=star_off,
               field_star_offset_rms=[float(np.std(so[:, 0])) if len(so) else None, float(np.std(so[:, 1])) if len(so) else None],
               field_star_offset_mean=[float(np.mean(so[:, 0])) if len(so) else None, float(np.mean(so[:, 1])) if len(so) else None],
               field_star_proj_nb_rms=float(np.std(proj)) if len(proj) else None,
               resid_rms_pair_region=float(np.sqrt(np.mean((resid_img.ravel()[r < 2.5]) ** 2))),
               resid_frac_pair_peak=float(np.max(np.abs(resid_img.ravel()[r < 2.5])) / M.ravel()[r < 2.5].max()),
               chi2_field=float(2 * res_field.cost), npix_field=int(field.sum()), zp=float(zp),
               bkg=[float(v) for v in x[-3:]], target_xy=[float(xt), float(yt)], nb_xy=[float(g.x[ib]), float(g.y[ib])])
    flux_all = np.where(free, 0.0, x[nfree] * zp_flux(g["T"].values))
    flux_all[ifree] = x[:nfree]
    sc = g[["Source", "x", "y", "T", "Gmag", "r_t"]].copy(); sc["flux"] = flux_all; sc["free_static"] = free
    sc.to_csv(os.path.join(OUT, f"scene_s{s}.csv"), index=False)
    np.savez(os.path.join(OUT, f"static_s{s}.npz"), M=M, model=model, resid=resid_img, x=x, free=free,
             ifree=ifree, imgs_pair=np.array([imgs[it], imgs[ib]]))
    json.dump(out, open(os.path.join(OUT, f"static_s{s}.json"), "w"), indent=1)

    fig, ax = plt.subplots(1, 3, figsize=(13, 4.3))
    for a, im, ti in zip(ax, [M, model, resid_img / np.sqrt(np.clip(M, 1, None))], ["median image", "PRF model", "(data-model)/sqrt(data)"]):
        if "data-model" in ti:
            v = np.percentile(np.abs(im), 99); h = a.imshow(im, origin="lower", cmap="RdBu_r", vmin=-v, vmax=v)
        else:
            h = a.imshow(np.log10(np.clip(im, 1, None)), origin="lower", cmap="viridis")
        plt.colorbar(h, ax=a, fraction=0.046)
        a.plot(xt + pf[0], yt + pf[1], "r+", ms=10); a.plot(g.x[ib] + pf[0], g.y[ib] + pf[1], "mx", ms=8)
        a.set_title(f"S{s} {ti}")
    plt.tight_layout(); plt.savefig(os.path.join(OUT, f"static_s{s}.png"), dpi=90); plt.close()
    if verbose:
        print(json.dumps({k: v for k, v in out.items() if k != "field_star_offsets"}, indent=None))
    return out


if __name__ == "__main__":
    secs = [int(a) for a in sys.argv[1:]] or SECTORS
    for s in secs:
        static_fit(s)
