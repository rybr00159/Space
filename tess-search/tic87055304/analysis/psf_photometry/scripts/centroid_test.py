"""Semi-independent check: in-transit motion of the target/neighbour pair image relative to the field.

registration() fits every cadence with I = a M + dx Gx + dy Gy + c on (i) field pixels (> 3 px from the target) and
(ii) the pair region (< 2.5 px).  The series (pair - field) is the pair's own apparent motion (a centroid-like series,
free of the per-star flux split).  Its in-transit change is measured with the depth estimator (template shape, local
quadratic baselines, control-phase errors) and compared with the change predicted by applying the same linear fit to the
model difference image  -dF_pair * PRF_host  for host = target, neighbour (dF_pair from the nominal PSF photometry).
"""
import os, sys, json
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from common import *
from psfphot import load_scene, scene_derivs
from depthfit import summarize
from analyze_lc import masks


def main():
    nom = json.load(open(os.path.join(ROOT, "lc", "nominal", "depths.json")))
    out = {}
    tot = {"t": 0.0, "n": 0.0}
    for s in SECTORS:
        sec = load_sector(s); sc, js = load_scene(s)
        prf = PRFModel(sec, sig_px=js["sig_jit"], xt=js["target_xy"][0], yt=js["target_xy"][1])
        d = np.load(os.path.join(ROOT, "lc", "nominal", f"psf_s{s}.npz"))
        m = masks(d)
        rf, rp = d["reg_field"], d["reg_pair"]
        u = np.array(js["nb_unit"]); v = np.array([-u[1], u[0]])
        ddx = rp[0] - rf[0]; ddy = rp[1] - rf[1]
        along = ddx * u[0] + ddy * u[1]; perp = ddx * v[0] + ddy * v[1]
        r = summarize(d["t"][m], [along[m], perp[m]], float(d["exptime"]), ["along", "perp"])
        # predictions: linear pair-region fit applied to model difference images
        X = sc.x.values + js["dx0"]; Y = sc.y.values + js["dy0"]; F = sc.flux.values
        m0, Gx, Gy, _, _, _ = scene_derivs(prf, X, Y, F)
        M = np.median(sec["F"], axis=0)
        yy, xx = np.mgrid[0:21, 0:21]
        xt, yt = js["target_xy"][0] + js["dx0"], js["target_xy"][1] + js["dy0"]
        reg = (np.hypot(xx - xt, yy - yt) < 2.5).ravel()
        w = 1 / np.sqrt(np.clip(M.ravel(), 1, None))[reg]
        A = np.vstack([M.ravel()[reg], Gx.ravel()[reg], Gy.ravel()[reg], np.ones(reg.sum())]).T
        pinv = np.linalg.pinv(A * w[:, None])
        dF = nom["sectors"][str(s)]["dF_t"] + nom["sectors"][str(s)]["dF_n"]
        src = sc.Source.values
        pred = {}
        for host, sid in (("target", TGT), ("neighbour", NB)):
            k = np.where(src == sid)[0][0]
            dimg = -dF * prf.img(X[k], Y[k])
            c = pinv @ (dimg.ravel()[reg] * w)
            # in-transit change of (a, dx*a, dy*a): shift = c[1]/a, a~1; depthfit measures A with y = base - A*tmpl
            pdx, pdy = c[1], c[2]
            pred[host] = dict(along=-(pdx * u[0] + pdy * u[1]), perp=-(pdx * v[0] + pdy * v[1]))
        meas_a, err_a = r["along"]["A"], r["along"]["err"]
        meas_p, err_p = r["perp"]["A"], r["perp"]["err"]
        chi_t = ((meas_a - pred["target"]["along"]) / err_a) ** 2 + ((meas_p - pred["target"]["perp"]) / err_p) ** 2
        chi_n = ((meas_a - pred["neighbour"]["along"]) / err_a) ** 2 + ((meas_p - pred["neighbour"]["perp"]) / err_p) ** 2
        tot["t"] += chi_t; tot["n"] += chi_n
        out[str(s)] = dict(meas_along_mpx=meas_a * 1e3, err_along_mpx=err_a * 1e3, meas_perp_mpx=meas_p * 1e3, err_perp_mpx=err_p * 1e3,
                           pred_target_along_mpx=pred["target"]["along"] * 1e3, pred_nb_along_mpx=pred["neighbour"]["along"] * 1e3,
                           pred_target_perp_mpx=pred["target"]["perp"] * 1e3, pred_nb_perp_mpx=pred["neighbour"]["perp"] * 1e3,
                           chi2_target=chi_t, chi2_nb=chi_n,
                           frac_along=(meas_a - pred["target"]["along"]) / (pred["neighbour"]["along"] - pred["target"]["along"]))
        q = out[str(s)]
        print(f"S{s}: pair-minus-field in-transit shift along nb {q['meas_along_mpx']:+.2f}+-{q['err_along_mpx']:.2f} mpx "
              f"(pred target {q['pred_target_along_mpx']:+.2f}, nb {q['pred_nb_along_mpx']:+.2f}); perp {q['meas_perp_mpx']:+.2f}+-{q['err_perp_mpx']:.2f} "
              f"(pred {q['pred_target_perp_mpx']:+.2f}/{q['pred_nb_perp_mpx']:+.2f}); chi2 t {chi_t:.1f} n {chi_n:.1f}; "
              f"fraction of the way target->nb {q['frac_along']:.2f}", flush=True)
    out["total"] = dict(chi2_target=tot["t"], chi2_nb=tot["n"], dof=2 * len(SECTORS))
    print("TOTAL chi2 target %.1f, neighbour %.1f (dof %d)" % (tot["t"], tot["n"], 2 * len(SECTORS)))
    os.makedirs(os.path.join(ROOT, "centroid"), exist_ok=True)
    json.dump(out, open(os.path.join(ROOT, "centroid", "centroid_test.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
