"""Per-cadence simultaneous PRF photometry of the target (TIC 87055304) and the 9" neighbour (TIC 87055315).

Pipeline (per sector):
 1. scene: static PRF fit of the median image (fit_static.py): Gaia DR3 positions at the sector epoch + global registration
    (dx0, dy0) fitted on field pixels > 3 px from the target, PRF = TESS_PRF (sigma smear fitted, ~0), all fluxes.
 2. registration(): per-cadence shift d(t) = (dx, dy) of the whole field from pixels > R_EXCL px from the target
    (target/neighbour pair excluded), by a second-order Taylor model of the static scene:
        I_t = a (M + dx Gx + dy Gy + dx^2/2 Gxx + dx dy Gxy + dy^2/2 Gyy) + b0 + b1 X + b2 Y.
    The same model fitted on the pair region only (r < 2.5 px) gives the pair's own shift (a classical centroid series).
 3. photometry(): in pixels within R_FIT of the target, fit
        I_t = sum_free F_k(t) PRF_k(x_k + dx0 + dx(t) + ox, y_k + dy0 + dy(t) + oy) + fixed faint stars + plane background
    with free fluxes for the target, the neighbour, every star with T < T_FREE within R_FIT + 0.3 px and T < 13.5 within
    R_FIT + 1 px; other stars fixed at their static fluxes (moved with d(t)).
    Registration modes: 'field' (default, from step 2), 'none' (d = 0), 'free' (d fitted per cadence from the local
    region together with the fluxes, Gauss-Newton; shows the pair-centroid / flux-ratio degeneracy).
    Weights 1/FLUX_ERR^2.  Returns fluxes, their per-cadence covariance, background and chi2.
"""
import os, sys, json, time
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import pandas as pd
from common import *

STATIC = os.path.join(ROOT, "static")
R_EXCL = 3.0


def load_scene(s):
    sc = pd.read_csv(os.path.join(STATIC, f"scene_s{s}.csv"))
    js = json.load(open(os.path.join(STATIC, f"static_s{s}.json")))
    return sc, js


def scene_model(prf, X, Y, F, ddx=0.0, ddy=0.0):
    m = np.zeros((21, 21))
    for x, y, f in zip(X, Y, F):
        if f != 0: m += f * prf.img(x + ddx, y + ddy)
    return m


def scene_derivs(prf, X, Y, F, h=0.02):
    m0 = scene_model(prf, X, Y, F)
    mxp = scene_model(prf, X, Y, F, h, 0); mxm = scene_model(prf, X, Y, F, -h, 0)
    myp = scene_model(prf, X, Y, F, 0, h); mym = scene_model(prf, X, Y, F, 0, -h)
    mpp = scene_model(prf, X, Y, F, h, h); mmm = scene_model(prf, X, Y, F, -h, -h)
    mpm = scene_model(prf, X, Y, F, h, -h); mmp = scene_model(prf, X, Y, F, -h, h)
    Gx = (mxp - mxm) / (2 * h); Gy = (myp - mym) / (2 * h)
    Gxx = (mxp - 2 * m0 + mxm) / h ** 2; Gyy = (myp - 2 * m0 + mym) / h ** 2
    Gxy = (mpp - mpm - mmp + mmm) / (4 * h ** 2)
    return m0, Gx, Gy, Gxx, Gxy, Gyy


def registration(sec, prf, sc, js, mask_fn=None):
    X = sc.x.values + js["dx0"]; Y = sc.y.values + js["dy0"]; F = sc.flux.values
    m0, Gx, Gy, Gxx, Gxy, Gyy = scene_derivs(prf, X, Y, F)
    M = np.median(sec["F"], axis=0)
    yy, xx = np.mgrid[0:21, 0:21]
    xt, yt = js["target_xy"][0] + js["dx0"], js["target_xy"][1] + js["dy0"]
    r = np.hypot(xx - xt, yy - yt)
    out = {}
    for name, msk in [("field", r > R_EXCL), ("pair", r < 2.5)]:
        m = msk.ravel()
        w = 1 / np.sqrt(np.clip(M.ravel(), 1, None))[m]
        base = [M.ravel()[m], Gx.ravel()[m], Gy.ravel()[m], np.ones(m.sum())]
        if name == "field":
            base += [(xx.ravel()[m] - 10) / 10, (yy.ravel()[m] - 10) / 10]
        A = np.vstack(base).T
        Aw = A * w[:, None]
        pinv = np.linalg.pinv(Aw)
        D = sec["F"].reshape(len(sec["t"]), -1)[:, m]
        c = (pinv @ (D * w[None, :]).T).T                       # first pass
        a = c[:, 0]; dx = c[:, 1] / a; dy = c[:, 2] / a
        for _ in range(2):                                        # second-order correction
            corr = a[:, None] * (0.5 * dx[:, None] ** 2 * Gxx.ravel()[m] + dx[:, None] * dy[:, None] * Gxy.ravel()[m]
                                 + 0.5 * dy[:, None] ** 2 * Gyy.ravel()[m])
            c = (pinv @ ((D - corr) * w[None, :]).T).T
            a = c[:, 0]; dx = c[:, 1] / a; dy = c[:, 2] / a
        res = (D - (A @ c.T).T) * w[None, :]
        out[name] = dict(dx=dx, dy=dy, a=a, bkg=c[:, 3], chi2=np.sum(res ** 2, 1), npix=int(m.sum()))
    return out


def photometry(sec, prf, sc, js, d, R_FIT=3.5, T_FREE=15.5, mode="field", offset=(0.0, 0.0), extra_free=(),
               data=None, return_design=False):
    """d: (N,2) per-cadence shifts (ignored for mode 'none'). offset: extra constant registration offset (px).
    data: optional replacement pixel cube (for injections)."""
    Fcube = sec["F"] if data is None else data
    N = len(sec["t"])
    X = sc.x.values + js["dx0"] + offset[0]; Y = sc.y.values + js["dy0"] + offset[1]; Fl = sc.flux.values
    src = sc.Source.values
    xt, yt = X[src == TGT][0], Y[src == TGT][0]
    yy, xx = np.mgrid[0:21, 0:21]
    reg = (np.hypot(xx - xt, yy - yt) <= R_FIT).ravel()
    dist = np.hypot(X - xt, Y - yt)
    free = (((sc["T"].values < T_FREE) & (dist < R_FIT + 0.3)) | ((sc["T"].values < 13.5) & (dist < R_FIT + 1.0))
            | (src == TGT) | (src == NB) | np.isin(src, list(extra_free)))
    order = np.argsort(sc["T"].values[free])
    ifree = np.where(free)[0][order]
    ifix = np.where(~free & (dist < R_FIT + 7))[0]
    fx0, fGx, fGy, _, _, _ = scene_derivs(prf, X[ifix], Y[ifix], Fl[ifix])
    fx0, fGx, fGy = fx0.ravel()[reg], fGx.ravel()[reg], fGy.ravel()[reg]
    xr = (xx.ravel()[reg] - xt) / 3; yr = (yy.ravel()[reg] - yt) / 3
    D = Fcube.reshape(N, -1)[:, reg]
    E = sec["E"].reshape(N, -1)[:, reg]
    nf = len(ifree)
    flux = np.zeros((N, nf)); cov2 = np.zeros((N, 2, 2)); bkg = np.zeros((N, 3)); chi2 = np.zeros(N)
    dd = np.zeros((N, 2))
    if mode == "field": dd = np.array(d, float)
    jt = list(src[ifree]).index(TGT); jn = list(src[ifree]).index(NB)
    designs = []
    for i in range(N):
        dx, dy = dd[i]
        w = 1 / E[i]
        for it_ in range(4 if mode == "free" else 1):
            cols = [prf.img(X[k] + dx, Y[k] + dy).ravel()[reg] for k in ifree]
            A = np.vstack(cols + [np.ones(reg.sum()), xr, yr]).T
            b = D[i] - fx0 - dx * fGx - dy * fGy
            if mode == "free" and it_ > 0:
                # add shift columns: sum_k F_k dPRF_k/dx, dy (+ fixed stars' gradient)
                h = 0.02
                gx = sum(f_ * (prf.img(X[k] + dx + h, Y[k] + dy).ravel()[reg] - prf.img(X[k] + dx - h, Y[k] + dy).ravel()[reg]) / (2 * h)
                         for k, f_ in zip(ifree, sol[:nf])) + fGx
                gy = sum(f_ * (prf.img(X[k] + dx, Y[k] + dy + h).ravel()[reg] - prf.img(X[k] + dx, Y[k] + dy - h).ravel()[reg]) / (2 * h)
                         for k, f_ in zip(ifree, sol[:nf])) + fGy
                A2 = np.hstack([A, gx[:, None], gy[:, None]])
                s2, *_ = np.linalg.lstsq(A2 * w[:, None], b * w, rcond=None)
                dx += s2[-2]; dy += s2[-1]
                cols = [prf.img(X[k] + dx, Y[k] + dy).ravel()[reg] for k in ifree]
                A = np.vstack(cols + [np.ones(reg.sum()), xr, yr]).T
                b = D[i] - fx0 - dx * fGx - dy * fGy
            Aw = A * w[:, None]
            sol, *_ = np.linalg.lstsq(Aw, b * w, rcond=None)
        if mode == "free":
            # final covariance including the shift parameters (marginalised)
            A2w = np.hstack([A, gx[:, None], gy[:, None]]) * w[:, None]
            C = np.linalg.pinv(A2w.T @ A2w)
        else:
            C = np.linalg.pinv(Aw.T @ Aw)
        flux[i] = sol[:nf]; bkg[i] = sol[nf:nf + 3]
        cov2[i] = C[np.ix_([jt, jn], [jt, jn])]
        chi2[i] = np.sum(((b - A @ sol) * w) ** 2)
        dd[i] = (dx, dy)
        if return_design and i == 0: designs.append((A, reg))
    res = dict(t=sec["t"], flux=flux, ids=src[ifree], T=sc["T"].values[ifree], dist=dist[ifree], cov_tn=cov2, bkg=bkg,
               chi2=chi2, npix=int(reg.sum()), d=dd, jt=jt, jn=jn, exptime=sec["exptime"])
    if return_design: res["design0"] = designs
    return res


def inject(sec, prf_inj, sc, js, d, source, depth, P_inj, t0_inj, offset=(0.0, 0.0)):
    """Return a copy of the cube with a transit (template shape, ephemeris P_inj/t0_inj) removed from star `source`:
    I -= depth * tmpl(t) * F_static * PRF_inj(x + d(t)).  prf_inj may differ from the extraction PRF."""
    k = np.where(sc.Source.values == source)[0][0]
    X = sc.x.values[k] + js["dx0"] + offset[0]; Y = sc.y.values[k] + js["dy0"] + offset[1]; F0 = sc.flux.values[k]
    tm = template(sec["t"], sec["exptime"], P_inj, t0_inj)
    cube = sec["F"].copy()
    for i in np.where(tm > 0)[0]:
        cube[i] -= depth * tm[i] * F0 * prf_inj.img(X + d[i, 0], Y + d[i, 1])
    return cube


def gauss_sigma(sec, sc, js):
    """Best-fit circular Gaussian sigma for the median image within 3 px of the target (fluxes of stars within 4 px
    free, others from the static scene at their TESS-PRF fluxes, constant background)."""
    from scipy.optimize import minimize_scalar
    M = np.median(sec["F"], axis=0)
    X = sc.x.values + js["dx0"]; Y = sc.y.values + js["dy0"]
    xt, yt = js["target_xy"][0] + js["dx0"], js["target_xy"][1] + js["dy0"]
    yy, xx = np.mgrid[0:21, 0:21]; reg = (np.hypot(xx - xt, yy - yt) < 3).ravel()
    dist = np.hypot(X - xt, Y - yt); fr = (dist < 4) & (sc["T"].values < 16)
    def cost(sg):
        g = GaussPSF(sg)
        cols = [g.img(X[k], Y[k]).ravel()[reg] for k in np.where(fr)[0]]
        fixed = sum(sc.flux.values[k] * g.img(X[k], Y[k]).ravel()[reg] for k in np.where(~fr & (dist < 8))[0])
        A = np.vstack(cols + [np.ones(reg.sum())]).T
        b = M.ravel()[reg] - fixed
        c, *_ = np.linalg.lstsq(A, b, rcond=None)
        return np.sum((b - A @ c) ** 2 / np.clip(M.ravel()[reg], 1, None))
    r = minimize_scalar(cost, bounds=(0.3, 1.5), method="bounded")
    return float(r.x)


def run_sector(s, tag="nominal", mode="field", R_FIT=3.5, T_FREE=15.5, prf_kind="tess", offset=(0.0, 0.0), sig=None,
               ccdpos=None):
    sec = load_sector(s)
    sc, js = load_scene(s)
    if prf_kind == "tess":
        prf = PRFModel(sec, sig_px=js["sig_jit"] if sig is None else sig, xt=js["target_xy"][0], yt=js["target_xy"][1],
                       ccdpos=ccdpos)
    elif prf_kind == "gauss":
        sg = sig if sig is not None else gauss_sigma(sec, sc, js)
        print(f"  S{s}: Gaussian PSF sigma = {sg:.3f} px", flush=True)
        prf = GaussPSF(sg)
    reg = registration(sec, prf, sc, js)
    d = np.vstack([reg["field"]["dx"], reg["field"]["dy"]]).T
    ph = photometry(sec, prf, sc, js, d, R_FIT=R_FIT, T_FREE=T_FREE, mode=mode, offset=offset)
    return sec, sc, js, reg, ph


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--sectors", default="13,39,66,93,104")
    ap.add_argument("--tag", default="nominal")
    ap.add_argument("--mode", default="field")
    ap.add_argument("--rfit", type=float, default=3.5)
    ap.add_argument("--tfree", type=float, default=15.5)
    ap.add_argument("--prf", default="tess")
    ap.add_argument("--sig", type=float, default=None)
    ap.add_argument("--offset_nb", type=float, default=0.0, help="extra registration offset along the target->nb direction (px)")
    ap.add_argument("--pairreg", action="store_true", help="use the pair-only static registration (F_nb/F_t fixed to Gaia)")
    ap.add_argument("--ccdpos", default=None, help="col,row of the CCD position at which to take the PRF (mismatch test)")
    a = ap.parse_args()
    od = os.path.join(ROOT, "lc", a.tag); os.makedirs(od, exist_ok=True)
    for s in [int(v) for v in a.sectors.split(",")]:
        t0 = time.time()
        js = json.load(open(os.path.join(STATIC, f"static_s{s}.json")))
        u = np.array(js["nb_unit"])
        off = a.offset_nb * u
        if a.pairreg:
            off = off + np.array(js["pair_offset_with_gaia_ratio"]) - np.array([js["dx0"], js["dy0"]])
        ccdpos = None if a.ccdpos is None else tuple(float(v) for v in a.ccdpos.split(","))
        sec, sc, js, reg, ph = run_sector(s, a.tag, a.mode, a.rfit, a.tfree, a.prf, tuple(off), a.sig, ccdpos)
        np.savez(os.path.join(od, f"psf_s{s}.npz"), t=ph["t"], flux=ph["flux"], ids=ph["ids"], T=ph["T"], dist=ph["dist"],
                 cov_tn=ph["cov_tn"], bkg=ph["bkg"], chi2=ph["chi2"], npix=ph["npix"], d=ph["d"], jt=ph["jt"], jn=ph["jn"],
                 exptime=ph["exptime"], reg_field=np.vstack([reg["field"][k] for k in ("dx", "dy", "a", "bkg", "chi2")]),
                 reg_pair=np.vstack([reg["pair"][k] for k in ("dx", "dy", "a", "bkg", "chi2")]), cad=sec["cad"])
        ft, fn = ph["flux"][:, ph["jt"]], ph["flux"][:, ph["jn"]]
        print(f"S{s} {a.tag}: N={len(ph['t'])} free={len(ph['ids'])} Ft med {np.median(ft):.1f} Fn med {np.median(fn):.1f} "
              f"corr(Ft,Fn) formal {np.median(ph['cov_tn'][:,0,1]/np.sqrt(ph['cov_tn'][:,0,0]*ph['cov_tn'][:,1,1])):.3f} "
              f"reg rms dx {np.std(reg['field']['dx']):.4f} dy {np.std(reg['field']['dy']):.4f} "
              f"chi2/pix {np.median(ph['chi2'])/ph['npix']:.2f} ({time.time()-t0:.0f}s)", flush=True)
