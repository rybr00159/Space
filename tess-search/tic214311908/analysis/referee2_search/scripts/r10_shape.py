"""r10: does removing the RRab leak change the transit shape inferred from the phot model-B light curve?

Global exposure-integrated trapezoid (T0, P fixed at the phot ephemeris; depth, T14, flat fraction and a small time
offset free) fitted to all cadences within +-0.3 d of the 20 transits, on the phot detrended flux (det) and on the
referee-2 RRab-cleaned flux (r01), per-sector robust white weights. Also the V-shape (flat = 0) refit and its
Delta-chi2. Errors: profile Delta-chi2 = 1 with white weights (formal only; shapes are compared, not calibrated).
Output: out/r10_shape.json
"""
import sys
sys.path.insert(0, "tic214311908/analysis/referee2_search/scripts")
from rlib2 import *
from scipy.optimize import least_squares


def main():
    z = np.load(f"{ME}/data/m3v_clean_r2.npz")
    t, s, cad = z["t"], z["s"], z["cad"]
    x = phase_dist(t)
    k = np.abs(x) < 0.3
    out = {}
    for nm, f in (("det", z["det"]), ("clean", z["clean"])):
        y = f[k] - 1; tt = x[k]; cc = cad[k]; ss = s[k]
        w = np.zeros(len(y))
        for q in np.unique(ss):
            m = ss == q
            oot = m & (np.abs(tt) > 0.06)
            y[m] -= np.median(y[oot]); w[m] = 1 / rstd(y[oot]) ** 2
        sw = np.sqrt(w)

        def model(p, flat=None):
            dep, T14, fl, dt = p if flat is None else (p[0], p[1], flat, p[2])
            return -dep * trapz_model(tt, dt, T14, np.clip(fl, 0, 0.999), cc)

        def res(p, flat=None):
            return (y - model(p, flat)) * sw
        r = least_squares(res, [0.015, 0.07, 0.6, 0.0], bounds=([0, 0.02, 0, -0.02], [0.05, 0.2, 0.999, 0.02]), x_scale=[0.003, 0.01, 0.1, 0.002])
        chi = float(np.sum(r.fun ** 2))
        rv = least_squares(lambda p: res(p, 0.0), [0.015, 0.08, 0.0], bounds=([0, 0.02, -0.02], [0.08, 0.25, 0.02]))
        # profile in flat fraction
        prof = []
        for fl in np.linspace(0, 0.95, 20):
            rr = least_squares(lambda p: res(p, fl), [r.x[0], r.x[1], r.x[3]], bounds=([0, 0.02, -0.02], [0.08, 0.25, 0.02]))
            prof.append((float(fl), float(np.sum(rr.fun ** 2) - chi), float(rr.x[1] * 24)))
        J = r.jac; cov = np.linalg.inv(J.T @ J)
        out[nm] = dict(depth_ppt=float(r.x[0] * 1e3), depth_err=float(np.sqrt(cov[0, 0]) * 1e3), T14_h=float(r.x[1] * 24),
                       T14_err_h=float(np.sqrt(cov[1, 1]) * 24), flat=float(r.x[2]), flat_err=float(np.sqrt(cov[2, 2])),
                       dt_min=float(r.x[3] * 1440), chi2=chi, n=int(len(y)), dchi2_V=float(np.sum(rv.fun ** 2) - chi), flat_profile=prof)
        print(nm, {kk: (round(v, 3) if isinstance(v, float) else v) for kk, v in out[nm].items() if kk != "flat_profile"}, flush=True)
        print("   profile:", [(round(a, 2), round(b, 1), round(c, 2)) for a, b, c in prof], flush=True)
    jdump(out, f"{ME}/out/r10_shape.json")


if __name__ == "__main__":
    main()
