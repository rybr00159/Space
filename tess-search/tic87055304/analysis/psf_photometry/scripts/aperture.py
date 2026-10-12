"""Aperture photometry versus aperture position, and per-pixel depth map, compared with PRF-model predictions.

Apertures (per sector): single pixels and 2x2 boxes containing / nearest to (a) the target, (b) the neighbour,
(c) the point opposite the neighbour (2 x_t - x_n), plus every single pixel within 2.5 px of the target (pixel map).
Measured quantity: flux decrement dF_ap (e-/s) from the depth estimator (depthfit.py, per-transit quadratic baseline),
and the relative depth dF_ap / F_ap, F_ap = median aperture flux minus the static sky plane.
Prediction under host k: dF_ap = S * frac_k(ap), frac_k = fraction of star k's PRF in the aperture at the in-transit
registration (TESS PRF, field registration), S free (one per sector and hypothesis).  chi2 uses the control-phase
covariance of all apertures (red noise, inter-aperture correlation) plus an optional PRF-model error term
(eps * predicted dF_ap, eps = 0, 0.05, 0.10).  Hypotheses: target, neighbour, any other Gaia star within 60",
and a free-position point source (grid search) for the pixel map.
"""
import os, sys, json
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import *
from psfphot import load_scene
from depthfit import depth_fit, control, summarize

REGS = os.environ.get("AP_REGS", "0") == "1"   # decorrelate aperture fluxes against field registration dx, dy
BKGSUB = os.environ.get("AP_BKG", "1") == "1"  # subtract the per-cadence PSF-fit background plane (default)
OUT = os.path.join(ROOT, "aperture" + ("_regs" if REGS else "") + ("" if BKGSUB else "_nobkg")); os.makedirs(OUT, exist_ok=True)


def chi2_host(y, C, pred):
    Ci = np.linalg.inv(C)
    S = (pred @ Ci @ y) / (pred @ Ci @ pred)
    r = y - S * pred
    return float(r @ Ci @ r), float(S)


def chi2_host_eps(y, C, pred, eps):
    S = chi2_host(y, C, pred)[1]
    for _ in range(3):
        Ce = C + np.diag((eps * S * pred) ** 2)
        c, S = chi2_host(y, Ce, pred)
    return c, S


LAM = 0.5


def run(s):
    sec = load_sector(s); sc, js = load_scene(s)
    prf = PRFModel(sec, sig_px=js["sig_jit"], xt=js["target_xy"][0], yt=js["target_xy"][1])
    nom = np.load(os.path.join(ROOT, "lc", "nominal", f"psf_s{s}.npz"))
    assert np.allclose(nom["t"], sec["t"])
    d = nom["d"]
    t = sec["t"]; ex = sec["exptime"]
    # masks as in analyze_lc
    bk = nom["bkg"][:, 0]; med = np.median(bk); mad = 1.4826 * np.median(np.abs(bk - med))
    m = (bk < med + 6 * mad) & (nom["chi2"] < 5 * np.median(nom["chi2"]))
    dt, ep = phase_dt(t)
    intr = (np.abs(dt) < DUR / 2) & m
    dmean = d[intr].mean(0)
    X = sc.x.values + js["dx0"] + dmean[0]; Y = sc.y.values + js["dy0"] + dmean[1]
    src = sc.Source.values
    it = np.where(src == TGT)[0][0]; ib = np.where(src == NB)[0][0]
    xt, yt, xn, yn = X[it], Y[it], X[ib], Y[ib]
    xo, yo = 2 * xt - xn, 2 * yt - yn
    stat = np.load(os.path.join(ROOT, "static", f"static_s{s}.npz"))
    yy, xx = np.mgrid[0:21, 0:21]
    bkx = stat["x"][-3:]
    sky = bkx[0] + bkx[1] * (xx - 10) / 10 + bkx[2] * (yy - 10) / 10
    Mimg = np.median(sec["F"][m], axis=0) - sky

    def pix_of(x, y): return (int(np.floor(y + 0.5)), int(np.floor(x + 0.5)))
    def box2_of(x, y):
        # 2x2 box whose centre (corner point between 4 pixels) is nearest to (x, y)
        cx, cy = int(np.floor(x)), int(np.floor(y))   # lower-left pixel of the box: pixels cx..cx+1
        return [(cy + a, cx + b) for a in (0, 1) for b in (0, 1)]
    aps = {}
    aps["1px_target"] = [pix_of(xt, yt)]; aps["1px_neighbour"] = [pix_of(xn, yn)]; aps["1px_opposite"] = [pix_of(xo, yo)]
    aps["2x2_target"] = box2_of(xt, yt); aps["2x2_neighbour"] = box2_of(xn, yn); aps["2x2_opposite"] = box2_of(xo, yo)
    # 3x3 around target pixel for reference
    py, px = pix_of(xt, yt)
    aps["3x3_target"] = [(py + a, px + b) for a in (-1, 0, 1) for b in (-1, 0, 1)]
    # drop duplicate apertures (e.g. target and neighbour in the same pixel) - keep the first, record the alias
    aliases = {}
    for k in list(aps.keys()):
        for k2 in list(aps.keys()):
            if k2 == k: break
            if sorted(aps[k2]) == sorted(aps[k]):
                aliases[k] = k2; del aps[k]; break
    named = list(aps.keys())
    # pixel map
    rpix = np.hypot(xx - xt, yy - yt)
    pixlist = [(int(a), int(b)) for a, b in zip(*np.where(rpix < 2.5))]
    for (a, b) in pixlist: aps[f"pix_{a}_{b}"] = [(a, b)]
    names = list(aps.keys())
    series = []
    Fap = []
    # per-cadence local background plane from the nominal PSF fit (b0 + b1 (x-xt0)/3 + b2 (y-yt0)/3), subtracted per pixel
    xt0 = sc.x.values[it] + js["dx0"]; yt0 = sc.y.values[it] + js["dy0"]
    bk3 = nom["bkg"]
    for k in names:
        idx = tuple(np.array(aps[k]).T)
        raw = sec["F"][:, idx[0], idx[1]].sum(1)
        if BKGSUB:
            pxx = idx[1].astype(float); pyy = idx[0].astype(float)
            bsum = len(pxx) * bk3[:, 0] + bk3[:, 1] * np.sum((pxx - xt0) / 3) + bk3[:, 2] * np.sum((pyy - yt0) / 3)
            raw = raw - bsum
            Fap.append(float(np.median(raw[m])))
        else:
            Fap.append(Mimg[idx].sum())
        series.append(raw)
    series = np.array(series)
    # host fractions at in-transit registration, averaged over in-transit cadences (subsample for speed)
    sub = np.where(intr)[0][::max(1, intr.sum() // 40)]
    def frac_for(kstar, xoff=0.0, yoff=0.0, xy=None):
        acc = np.zeros(len(names))
        for i in sub:
            if xy is None:
                im = prf.img(sc.x.values[kstar] + js["dx0"] + d[i, 0] + xoff, sc.y.values[kstar] + js["dy0"] + d[i, 1] + yoff)
            else:
                im = prf.img(xy[0] + d[i, 0] - dmean[0], xy[1] + d[i, 1] - dmean[1])
            acc += np.array([im[tuple(np.array(aps[k]).T)].sum() for k in names])
        return acc / len(sub)
    fr_t = frac_for(it); fr_n = frac_for(ib)
    # depths
    tm_ = t[m]
    regs = None
    if REGS:
        dxm, dym = d[m, 0], d[m, 1]
        regs = [dxm, dym, dxm ** 2, dym ** 2, dxm * dym]
    r = summarize(tm_, list(series[:, m]), ex, names, regs=regs)
    y = np.array([r[k]["A"] for k in names]); C = np.array(r["_cov_control"])
    e_ad = np.array([r[k]["err"] for k in names])
    # rescale covariance so its diagonal equals the adopted (max) errors
    sd = np.sqrt(np.diag(C)); Cn = C / np.outer(sd, sd)
    # shrink the sample correlation matrix towards the identity (the control-phase sample is small and its windows
    # overlap); LAM = 0.5 adopted, 1.0 = diagonal only
    Cn = (1 - LAM) * Cn + LAM * np.eye(len(Cn))
    C_ad = Cn * np.outer(e_ad, e_ad)
    out = dict(sector=s, apertures={}, target_xy=[xt, yt], nb_xy=[xn, yn], opp_xy=[xo, yo], aliases=aliases, shrinkage=LAM)
    for j, k in enumerate(names):
        out["apertures"][k] = dict(pixels=aps[k], dF=float(y[j]), dF_err=float(e_ad[j]), F_ap=float(Fap[j]),
                                    depth=float(y[j] / Fap[j]), depth_err=float(e_ad[j] / Fap[j]),
                                    frac_target=float(fr_t[j]), frac_nb=float(fr_n[j]))
    # chi2 tests
    tests = {}
    for setname, sel in [("named", [names.index(k) for k in named if k != "3x3_target"]),
                         ("named_diag", [names.index(k) for k in named if k != "3x3_target"]),
                         ("pixelmap", [names.index(k) for k in names if k.startswith("pix_")])]:
        sel = np.array(sel)
        ys, Cs = y[sel], C_ad[np.ix_(sel, sel)]
        if setname == "named_diag": Cs = np.diag(np.diag(Cs))
        tt = {}
        for eps in (0.0, 0.05, 0.10):
            ct, St = chi2_host_eps(ys, Cs, fr_t[sel], eps); cn, Sn = chi2_host_eps(ys, Cs, fr_n[sel], eps)
            tt[f"eps{eps}"] = dict(chi2_target=ct, S_target=St, chi2_nb=cn, S_nb=Sn, dof=len(sel) - 1, dchi2=ct - cn)
        tests[setname] = tt
    # other Gaia stars as hosts (pixel map, eps 0.05)
    sel = np.array([names.index(k) for k in names if k.startswith("pix_")])
    others = []
    dist = np.hypot(X - xt, Y - yt)
    for kstar in np.where(dist < 60 / 21.0)[0]:
        fr = frac_for(kstar)
        if fr[sel].sum() < 1e-3: continue
        c, S_ = chi2_host_eps(y[sel], C_ad[np.ix_(sel, sel)], fr[sel], 0.05)
        F_star = sc.flux.values[kstar]
        others.append(dict(source=int(src[kstar]), T=float(sc["T"].values[kstar]), sep_arcsec=float(dist[kstar] * 21.0),
                           chi2=c, S=S_, flux_static=float(F_star), implied_depth=float(S_ / max(F_star, 1e-3))))
    others.sort(key=lambda z: z["chi2"])
    out["other_hosts_pixelmap_eps0.05"] = others[:12]
    # free-position point source on the pixel map (grid around target, step 0.02 px)
    gx = np.arange(-1.0, 1.0001, 0.02); best = None; grid = np.zeros((len(gx), len(gx)))
    Cs = C_ad[np.ix_(sel, sel)]
    for a_, ox in enumerate(gx):
        for b_, oy in enumerate(gx):
            im = prf.img(xt + ox, yt + oy)
            fr = np.array([im[tuple(np.array(aps[names[j]]).T)].sum() for j in sel])
            c, S_ = chi2_host(y[sel], Cs, fr)
            grid[b_, a_] = c
            if best is None or c < best[0]: best = (c, ox, oy, S_)
    # 1-sigma region of the free position (delta chi2 < 2.30)
    okg = grid < best[0] + 2.30
    GX, GY = np.meshgrid(gx, gx)
    out["free_position"] = dict(chi2_min=best[0], dx_px=best[1], dy_px=best[2], S=best[3],
                                dx_range=[float(GX[okg].min()), float(GX[okg].max())], dy_range=[float(GY[okg].min()), float(GY[okg].max())],
                                nb_rel_px=[float(xn - xt), float(yn - yt)],
                                dist_to_nb_px=float(np.hypot(best[1] - (xn - xt), best[2] - (yn - yt))),
                                dist_to_target_px=float(np.hypot(best[1], best[2])))
    # sky offset of the best position (arcsec E, N)
    wcs = sec["wcs"]
    ra0, de0 = wcs.pixel_to_world_values(xt - dmean[0] - js["dx0"], yt - dmean[1] - js["dy0"])
    ra1, de1 = wcs.pixel_to_world_values(xt + best[1] - dmean[0] - js["dx0"], yt + best[2] - dmean[1] - js["dy0"])
    out["free_position"]["offset_arcsec_EN"] = [float((ra1 - ra0) * np.cos(np.radians(de0)) * 3600), float((de1 - de0) * 3600)]
    out["tests"] = tests
    np.savez(os.path.join(OUT, f"grid_s{s}.npz"), grid=grid, gx=gx, nb=np.array([xn - xt, yn - yt]))
    json.dump(out, open(os.path.join(OUT, f"aperture_s{s}.json"), "w"), indent=1, default=float)
    np.savez(os.path.join(OUT, f"pixmap_s{s}.npz"), names=np.array(names), y=y, err=e_ad, Fap=np.array(Fap), fr_t=fr_t, fr_n=fr_n,
             C=C_ad, pix=np.array([aps[k][0] if len(aps[k]) == 1 else (-1, -1) for k in names]))
    return out


if __name__ == "__main__":
    secs = [int(v) for v in sys.argv[1].split(",")] if len(sys.argv) > 1 else SECTORS
    for s in secs:
        o = run(s)
        print(f"S{s}:", {k: (round(v['depth'] * 1e3, 2), round(v['depth_err'] * 1e3, 2)) for k, v in o["apertures"].items() if not k.startswith("pix_")})
        for k, v in o["tests"].items():
            print("   ", k, {e: (round(q["chi2_target"], 1), round(q["chi2_nb"], 1), q["dof"]) for e, q in v.items()})
        print("    free position:", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in o["free_position"].items()})
        print("    best other hosts:", [(h["source"] % 100000, round(h["T"], 1), round(h["sep_arcsec"], 1), round(h["chi2"], 1), round(h["implied_depth"], 3)) for h in o["other_hosts_pixelmap_eps0.05"][:5]], flush=True)
