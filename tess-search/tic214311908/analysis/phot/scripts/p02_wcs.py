"""Per-sector WCS calibration on normal-colour field stars (the Mira and its 3-px surroundings are excluded),
and a measurement of the PRF blur (extra Gaussian width) that best matches the data.

Model of the sector median image (quality==0 cadences): sum over Gaia stars of flux_i * PRF(x_i + dx, y_i + dy)
with (dx, dy) a global shift and the PRF optionally blurred by a Gaussian of sigma `blur` px; stars with T < TFREE
have free fluxes, fainter ones (T < 18.5) enter as one template with Gaia-predicted relative fluxes and a free scale;
background = 2-D quadratic polynomial. Nonlinear parameters (dx, dy, blur) by least squares.
Check: individual PRF fits (free x, y, flux, local constant) to bright, isolated, normal-colour stars; offsets of
each star vs its BP-RP; robust mean offset and its scatter.
Outputs: out/wcs_s{sector}.json, figs/wcs_s{sector}.png
"""
import sys; sys.path.insert(0, "tic214311908/analysis/phot/scripts")
from photlib import *
from scipy.optimize import least_squares
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

TFREE = 14.5


def design(C, prf, g, dx, dy, blur, mask, free_idx, faint_idx):
    cols = []
    for i in free_idx:
        cols.append(prf.image(g.x[i] + dx, g.y[i] + dy, blur=blur)[mask])
    fa = np.zeros(mask.sum())
    for i in faint_idx:
        fa += tmag_flux(g["T"][i]) * prf.image(g.x[i] + dx, g.y[i] + dy, blur=blur)[mask]
    cols.append(fa)
    yy, xx = np.mgrid[0:C.ny, 0:C.nx]; u = (xx[mask] - 10) / 10; v = (yy[mask] - 10) / 10
    cols += [np.ones(mask.sum()), u, v, u * u, v * v, u * v]
    return np.column_stack(cols)


def run(sector):
    C = Cutout(sector); g = C.gaia(18.5); prf = PRFModel(C)
    good = C.fin & (C.q == 0)
    img = np.median(C.flux[good], 0)
    mira = g[g.source_id == IDS["Mira"]].iloc[0]
    yy, xx = np.mgrid[0:C.ny, 0:C.nx]
    mask = np.hypot(xx - mira.x, yy - mira.y) > 3.0
    # stars that matter: inside cutout +- 3 px
    ins = (g.x > -3) & (g.x < C.nx + 2) & (g.y > -3) & (g.y < C.ny + 2)
    g = g[ins].reset_index(drop=True)
    # exclude the Mira and the M-dwarf trio from the global fit (inside the masked region anyway)
    free_idx = [i for i in range(len(g)) if g["T"][i] < TFREE]
    faint_idx = [i for i in range(len(g)) if g["T"][i] >= TFREE]
    b = img[mask]
    def solve(p):
        A = design(C, prf, g, p[0], p[1], abs(p[2]), mask, free_idx, faint_idx)
        sig = 0.02 * np.abs(b) + 5.0
        c, *_ = np.linalg.lstsq(A / sig[:, None], b / sig, rcond=None)
        return A, c, sig
    def res(p):
        A, c, sig = solve(p); return (A @ c - b) / sig
    best = None
    for dx in (-0.3, 0, 0.3):
        for dy in (-0.3, 0, 0.3):
            for bl in (0.05, 0.25):
                c2 = np.sum(res([dx, dy, bl]) ** 2)
                if best is None or c2 < best[0]: best = (c2, [dx, dy, bl])
    r = least_squares(res, best[1], diff_step=[1e-3, 1e-3, 1e-2], x_scale=[0.05, 0.05, 0.05])
    dx, dy, blur = r.x[0], r.x[1], abs(r.x[2])
    J = r.jac; cov = np.linalg.pinv(J.T @ J) * (2 * r.cost / (mask.sum() - J.shape[1] - len(free_idx) - 7))
    A, c, sig = solve(r.x)
    model = np.full(img.shape, np.nan); model[mask] = A @ c
    resid = img - model
    chi2_noblur = np.sum(res([dx, dy, 0.0]) ** 2)
    fl = dict(zip(g.source_id[free_idx].tolist(), c[:len(free_idx)].tolist()))
    faint_scale = float(c[len(free_idx)])
    # ---- individual fits to isolated normal-colour stars (independent of the global solution's shift)
    prf.blur = blur
    indiv = []
    for i in free_idx:
        s_ = g.iloc[i]
        if not (2.0 < s_.x < C.nx - 3.0 and 2.0 < s_.y < C.ny - 3.0): continue
        if np.hypot(s_.x - mira.x, s_.y - mira.y) < 3.5: continue
        if not (np.isfinite(s_.bp_rp)): continue
        d = np.hypot(g.x - s_.x, g.y - s_.y)
        nb = g[(d < 1.5) & (d > 0) & (g["T"] < s_["T"] + 1.5)]
        if len(nb): continue  # isolated only
        m = np.hypot(xx - s_.x, yy - s_.y) <= 2.5
        # subtract every other star with its global-fit flux (or Gaia flux for faint ones)
        oth = np.zeros(img.shape)
        for j in range(len(g)):
            if j == i or np.hypot(g.x[j] - s_.x, g.y[j] - s_.y) > 6: continue
            f_ = fl.get(g.source_id[j], tmag_flux(g["T"][j]) * faint_scale)
            oth += f_ * prf.image(g.x[j] + dx, g.y[j] + dy)
        dat = (img - oth)[m]; w = 1 / (0.02 * np.abs(img[m]) + 5)
        def rr(p):
            M = np.column_stack([prf.image(p[0], p[1])[m], np.ones(m.sum())])
            cc, *_ = np.linalg.lstsq(M * w[:, None], dat * w, rcond=None); return (M @ cc - dat) * w
        q = least_squares(rr, [s_.x + dx, s_.y + dy], diff_step=1e-3, x_scale=0.05)
        Jq = q.jac; cq = np.linalg.pinv(Jq.T @ Jq) * max(1, 2 * q.cost / max(1, m.sum() - 4))
        indiv.append(dict(source_id=int(s_.source_id), T=float(s_["T"]), bp_rp=float(s_.bp_rp), ruwe=float(s_.ruwe),
                          x=float(s_.x), y=float(s_.y), ox=float(q.x[0] - s_.x), oy=float(q.x[1] - s_.y),
                          ex=float(np.sqrt(cq[0, 0])), ey=float(np.sqrt(cq[1, 1]))))
    ind = pd.DataFrame(indiv)
    nrm = ind[(ind.bp_rp < 2.0) & (ind.ruwe < 1.4)]
    med_off = [float(np.median(nrm.ox)), float(np.median(nrm.oy))]
    sc = [float(robust_std(nrm.ox)), float(robust_std(nrm.oy))]
    # colour slope of individual offsets (all colours) - lateral colour within the field
    Acol = np.column_stack([np.ones(len(ind)), ind.bp_rp - 1.0])
    slope = []
    for k in ("ox", "oy"):
        cc, *_ = np.linalg.lstsq(Acol, ind[k].values, rcond=None)
        rs = ind[k].values - Acol @ cc; cv = np.linalg.pinv(Acol.T @ Acol) * robust_std(rs) ** 2
        slope.append((float(cc[1]), float(np.sqrt(cv[1, 1]))))
    ref = ind[ind.source_id == IDS["ref"]]
    out = dict(sector=sector, global_shift=[float(dx), float(dy)], global_shift_err=[float(np.sqrt(cov[0, 0])), float(np.sqrt(cov[1, 1]))],
               blur=float(blur), blur_err=float(np.sqrt(cov[2, 2])), chi2=float(2 * r.cost), chi2_noblur=float(chi2_noblur), npix=int(mask.sum()),
               n_free=len(free_idx), faint_scale=faint_scale,
               indiv_normal_median=med_off, indiv_normal_rstd=sc, n_indiv_normal=int(len(nrm)),
               colour_slope_px_per_mag={"x": slope[0], "y": slope[1]},
               ref_star_offset=ref[["ox", "oy", "ex", "ey"]].values.tolist(),
               resid_rms_rel=float(robust_std((resid[mask]) / np.clip(img[mask], 1, None))),
               jacobian=C.jacobian(10, 10).tolist(), indiv=indiv, free_fluxes={str(k): v for k, v in fl.items()})
    save_json(out, f"{ROOT}/out/wcs_s{sector}.json")
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    ax[0].imshow(img, origin="lower", vmin=np.percentile(img, 5), vmax=np.percentile(img, 99)); ax[0].set_title(f"S{sector} median image")
    v = np.nanpercentile(np.abs(resid[mask] / np.sqrt(np.abs(img[mask]))), 98)
    ax[1].imshow(resid / np.sqrt(np.abs(img)), origin="lower", cmap="RdBu_r", vmin=-v, vmax=v); ax[1].set_title("(data-model)/sqrt(data), Mira region masked")
    ax[2].errorbar(ind.bp_rp, ind.ox, ind.ex, fmt="o", label="x offset"); ax[2].errorbar(ind.bp_rp + 0.03, ind.oy, ind.ey, fmt="s", label="y offset")
    ax[2].axhline(dx, color="C0", ls=":"); ax[2].axhline(dy, color="C1", ls=":"); ax[2].set_xlabel("BP-RP"); ax[2].set_ylabel("offset (px)")
    ax[2].legend(); ax[2].set_ylim(-1, 1)
    fig.tight_layout(); fig.savefig(f"{ROOT}/figs/wcs_s{sector}.png", dpi=80); plt.close(fig)
    print(f"S{sector}: global shift ({dx:+.3f}, {dy:+.3f}) +- ({np.sqrt(cov[0,0]):.3f}, {np.sqrt(cov[1,1]):.3f}) px, blur {blur:.3f} px "
          f"(chi2 {2*r.cost:.0f} vs {chi2_noblur:.0f} unblurred, {mask.sum()} px); individual normal-colour stars (n={len(nrm)}): "
          f"median ({med_off[0]:+.3f}, {med_off[1]:+.3f}) rstd ({sc[0]:.3f}, {sc[1]:.3f}); colour slope x {slope[0][0]:+.3f}+-{slope[0][1]:.3f}, "
          f"y {slope[1][0]:+.3f}+-{slope[1][1]:.3f} px/mag; ref star {np.round(ref[['ox','oy']].values, 3).tolist()}; rel resid rms {out['resid_rms_rel']:.3f}", flush=True)


if __name__ == "__main__":
    for s in ([int(a) for a in sys.argv[1:]] or SECTORS):
        run(s)
