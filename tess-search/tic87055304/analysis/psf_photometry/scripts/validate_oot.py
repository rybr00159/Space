"""Out-of-transit validation of the separated PSF light curves (nominal photometry).

Per sector, out of transit (|dt| > 0.75 * duration from any transit), after removing a 1-d running median:
  - point-to-point scatter (std of first differences / sqrt 2) of Ft, Fn, Ft+Fn (e-/s and relative);
  - median formal (per-cadence LSQ) errors and correlation coefficient of Ft, Fn;
  - empirical correlation of the detrended Ft and Fn (full cadence and 4.5-h binned);
  - rms in 4.5-h bins versus the white-noise expectation (beta factor);
  - correlation of the detrended Ft with the field registration (dx, dy) and background.
"""
import os, json
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.ndimage import median_filter
from common import *
from analyze_lc import masks


def runmed(t, y, w=1.0):
    out = np.empty_like(y)
    for i in range(len(t)):
        m = np.abs(t - t[i]) < w / 2
        out[i] = np.median(y[m])
    return out


def main(tag="nominal"):
    res = {}
    for s in SECTORS:
        d = np.load(os.path.join(ROOT, "lc", tag, f"psf_s{s}.npz"))
        m = masks(d)
        t = d["t"]; dt, _ = phase_dt(t)
        oot = m & (np.abs(dt) > 0.75 * DUR)
        ft = d["flux"][:, int(d["jt"])]; fn = d["flux"][:, int(d["jn"])]
        t_, a, b = t[oot], ft[oot], fn[oot]
        ra = a - runmed(t_, a); rb = b - runmed(t_, b)
        p2p = lambda x: float(np.std(np.diff(x)) / np.sqrt(2))
        cov = d["cov_tn"][oot]
        form_t = float(np.median(np.sqrt(cov[:, 0, 0]))); form_n = float(np.median(np.sqrt(cov[:, 1, 1])))
        form_rho = float(np.median(cov[:, 0, 1] / np.sqrt(cov[:, 0, 0] * cov[:, 1, 1])))
        ex = float(d["exptime"]); nb = max(1, int(round(DUR * 86400 / ex)))
        k = len(ra) // nb
        ba = ra[:k * nb].reshape(k, nb).mean(1); bb = rb[:k * nb].reshape(k, nb).mean(1)
        beta_t = float(np.std(ba) / (np.std(ra) / np.sqrt(nb))); beta_n = float(np.std(bb) / (np.std(rb) / np.sqrt(nb)))
        rf = d["reg_field"][:, oot]
        res[str(s)] = dict(n_oot=int(oot.sum()), cadence_s=ex,
                           p2p_t=p2p(ra), p2p_n=p2p(rb), p2p_sum=p2p(ra + rb),
                           p2p_t_rel_ppt=p2p(ra) / np.median(a) * 1e3, p2p_n_rel_ppt=p2p(rb) / np.median(b) * 1e3,
                           p2p_sum_rel_ppt=p2p(ra + rb) / np.median(a + b) * 1e3,
                           formal_err_t=form_t, formal_err_n=form_n, formal_corr=form_rho,
                           emp_corr=float(np.corrcoef(ra, rb)[0, 1]), emp_corr_binned=float(np.corrcoef(ba, bb)[0, 1]),
                           rms_bin_t=float(np.std(ba)), rms_bin_n=float(np.std(bb)), beta_t=beta_t, beta_n=beta_n,
                           corr_t_dx=float(np.corrcoef(ra, rf[0] - runmed(t_, rf[0]))[0, 1]),
                           corr_t_dy=float(np.corrcoef(ra, rf[1] - runmed(t_, rf[1]))[0, 1]),
                           corr_t_bkg=float(np.corrcoef(ra, d["bkg"][oot, 0] - runmed(t_, d["bkg"][oot, 0]))[0, 1]))
        q = res[str(s)]
        print(f"S{s}: p2p Ft {q['p2p_t']:.1f} (formal {form_t:.1f}) Fn {q['p2p_n']:.1f} (formal {form_n:.1f}) sum {q['p2p_sum']:.1f} e-/s; "
              f"rel: t {q['p2p_t_rel_ppt']:.2f} n {q['p2p_n_rel_ppt']:.1f} sum {q['p2p_sum_rel_ppt']:.2f} ppt; corr formal {form_rho:.2f} emp {q['emp_corr']:.2f} "
              f"binned {q['emp_corr_binned']:.2f}; beta t {beta_t:.2f} n {beta_n:.2f}; corr(Ft, dx,dy,bkg) {q['corr_t_dx']:.2f} {q['corr_t_dy']:.2f} {q['corr_t_bkg']:.2f}", flush=True)
    json.dump(res, open(os.path.join(ROOT, "lc", tag, "oot_validation.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
