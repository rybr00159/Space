"""Injection-recovery / leakage test for the simultaneous PRF photometry.

For each sector a transit with the template shape and the real period is injected at a fake phase (0.27 or 0.73, away
from the real transit and its phase 0.5) into the PIXEL data of ONE star: I -= tmpl(t) * dF_inj * PRF_inj(x_star + d(t)),
with dF_inj = 0.0036 * (F_t + F_n) (i.e. the observed pair depth).  The whole pipeline (registration + photometry) is rerun
on the injected cube and the depth estimator is applied at the injected ephemeris to both PSF light curves.
Recovery = dF(host)/dF_inj, leakage = dF(other)/dF_inj, measured on (injected - uninjected) light curves (noise-free
response) and on the injected light curve itself (noisy; error from the nominal control-phase analysis).
Injection PRFs: the extraction PRF itself ('same'), and deliberately mismatched ones: TESS PRF smeared by 0.4 px,
a Gaussian PSF (sigma fitted), the TESS PRF taken at a distant CCD position, and the star's position offset by
+-0.05 px along the target->neighbour direction (registration / astrometric error).
"""
import os, sys, json, time
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from common import *
from psfphot import load_scene, registration, photometry, inject, gauss_sigma
from depthfit import depth_fit

OUT = os.path.join(ROOT, "injection"); os.makedirs(OUT, exist_ok=True)


def measure(ph, t0i, m=None):
    ft = ph["flux"][:, ph["jt"]]; fn = ph["flux"][:, ph["jn"]]
    t = ph["t"]; ex = ph["exptime"]
    a = depth_fit(t, ft, ex, P, t0i, mask=m); b = depth_fit(t, fn, ex, P, t0i, mask=m)
    return a["A"], b["A"], a["n"]


def main(sectors):
    nominal = json.load(open(os.path.join(ROOT, "lc", "nominal", "depths.json")))
    results = []
    for s in sectors:
        t_start = time.time()
        sec = load_sector(s); sc, js = load_scene(s)
        prf = PRFModel(sec, sig_px=js["sig_jit"], xt=js["target_xy"][0], yt=js["target_xy"][1])
        reg = registration(sec, prf, sc, js)
        d = np.vstack([reg["field"]["dx"], reg["field"]["dy"]]).T
        base = photometry(sec, prf, sc, js, d)
        Ft0 = sc.flux[sc.Source == TGT].values[0]; Fn0 = sc.flux[sc.Source == NB].values[0]
        dF_inj = 0.0036 * (Ft0 + Fn0)
        u = np.array(js["nb_unit"])
        prf_smear = PRFModel(sec, sig_px=0.4, xt=js["target_xy"][0], yt=js["target_xy"][1])
        prf_gauss = GaussPSF(gauss_sigma(sec, sc, js))
        prf_far = PRFModel(sec, sig_px=0.0, ccdpos=(300.0, 300.0))
        configs = [("same", prf, (0, 0)), ("smear0.4", prf_smear, (0, 0)), ("gauss", prf_gauss, (0, 0)),
                   ("prf_ccd300", prf_far, (0, 0)), ("pos+0.05", prf, tuple(0.05 * u)), ("pos-0.05", prf, tuple(-0.05 * u))]
        nom = nominal["sectors"][str(s)]
        for phase in (0.27, 0.73):
            t0i = T0 + phase * P
            bt, bn, _ = measure(base, t0i)
            for host, src, F0 in (("target", TGT, Ft0), ("neighbour", NB, Fn0)):
                for name, pinj, off in configs:
                    cube = inject(sec, pinj, sc, js, d, src, dF_inj / F0, P, t0i, offset=off)
                    reg_i = registration(dict(sec, F=cube), prf, sc, js)
                    d_i = np.vstack([reg_i["field"]["dx"], reg_i["field"]["dy"]]).T
                    ph = photometry(sec, prf, sc, js, d_i, data=cube)
                    at, an, ntr = measure(ph, t0i)
                    rt, rn = at - bt, an - bn            # noise-free response
                    rec = dict(sector=s, phase=phase, host=host, inj_prf=name, dF_inj=float(dF_inj), n_tr=ntr,
                               resp_t=float(rt / dF_inj), resp_n=float(rn / dF_inj),
                               noisy_dFt=float(at), noisy_dFn=float(an), err_t=nom["dF_t_err"], err_n=nom["dF_n_err"],
                               f_target_noisy=float(at / (at + an)), f_target_resp=float(rt / (rt + rn)),
                               max_reg_change=float(np.max(np.abs(d_i - d))))
                    results.append(rec)
                    print(f"S{s} ph{phase} inj->{host:9s} [{name:10s}] response: target {rec['resp_t']:+.3f}, nb {rec['resp_n']:+.3f} "
                          f"(f_t {rec['f_target_resp']:+.3f}); noisy dFt {at:.1f}+-{nom['dF_t_err']:.1f} dFn {an:.1f}+-{nom['dF_n_err']:.1f} "
                          f"inj {dF_inj:.1f}; reg change {rec['max_reg_change']:.1e}", flush=True)
        print(f"S{s} done in {time.time()-t_start:.0f}s", flush=True)
        json.dump(results, open(os.path.join(OUT, "injection_results.json"), "w"), indent=1)
    return results


if __name__ == "__main__":
    secs = [int(v) for v in sys.argv[1].split(",")] if len(sys.argv) > 1 else SECTORS
    main(secs)
