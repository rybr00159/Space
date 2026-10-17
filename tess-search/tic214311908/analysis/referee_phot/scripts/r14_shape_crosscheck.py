"""Referee r14: transit shape (T14, flat fraction) cross-check. The referee's global trapezoid fit (r06 model,
least squares with free T0, P, T14, flat, depth; exposure-integrated) is applied to (a) the investigator's
lc_M3V_all.csv (model B, quality_ok & no flare, within +-0.35 d of transits, flux_err weights) and (b) the referee's
local-window LC (r05) — so differences between data and between fitting codes can be separated. Also a profile in
flat fraction. Output: out/shape_crosscheck.json"""
import sys; sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
from rlib import *
from r06_transit import shape_eph, load, window_sigma
from scipy.optimize import least_squares
G = json.load(open(f"{R}/out/transit_fit.json"))["global"]
T0r, Pr = G["T0"], G["P"]


def fit(t, f, e, ex, sign):
    def res(p):
        dT0, dP, T14, fl, dep = p; y = np.zeros(len(t))
        for x_ in np.unique(ex):
            m = ex == x_; y[m] = shape_eph(t[m], x_, T0r + dT0, Pr + dP, T14, fl)
        return (f - sign * dep * y) / e
    best = None
    for T14 in (0.06, 0.068, 0.076, 0.084):
        for fl in (0.3, 0.5, 0.7, 0.85):
            q = least_squares(res, [0, 0, T14, fl, 0.015 if sign < 0 else 15], bounds=([-0.02, -5e-4, 0.03, 0, 0], [0.02, 5e-4, 0.12, 0.98, 60]),
                              x_scale=[1e-3, 1e-5, 3e-3, 0.05, 0.002 if sign < 0 else 1], diff_step=[1e-4, 1e-7, 1e-3, 1e-2, 1e-3])
            if best is None or q.cost < best.cost: best = q
    q = best; cov = np.linalg.pinv(q.jac.T @ q.jac) * max(1, 2 * q.cost / (len(t) - 5)); pe = np.sqrt(np.diag(cov))
    prof = []
    for flv in (0.0, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9):
        qq = least_squares(lambda p: res([p[0], p[1], p[2], flv, p[3]]), [q.x[0], q.x[1], q.x[2], q.x[4]], x_scale=[1e-3, 1e-5, 3e-3, 0.002 if sign < 0 else 1], diff_step=[1e-4, 1e-7, 1e-3, 1e-3])
        prof.append((flv, float(2 * qq.cost - 2 * q.cost), float(qq.x[2] * 24)))
    return dict(T0=T0r + q.x[0], P=Pr + q.x[1], T14_h=q.x[2] * 24, T14_err_h=pe[2] * 24, flat=q.x[3], flat_err=pe[3], depth=q.x[4], depth_err=pe[4], chi2=2 * q.cost, n=len(t), profile=prof)


out = {}
M = pd.read_csv(f"{R}/../phot/lc_M3V_all.csv")
ok = (M.quality_ok == 1) & (M.flare_flag == 0) & np.isfinite(M.flux_norm)
ph = np.abs(((M.time_btjd_tdb - T0r) / Pr + 0.5) % 1 - 0.5) * Pr
d = M[ok & (ph < 0.35)]
out["investigator_lc"] = fit(d.time_btjd_tdb.values, d.flux_norm.values - 1, d.flux_err.values, d.cadence_s.values / 86400, -1)
# without S66 E-128 (their timing outlier)
ep = np.round((d.time_btjd_tdb - T0r) / Pr).astype(int)
d2 = d[ep != -128]
out["investigator_lc_noE-128"] = fit(d2.time_btjd_tdb.values, d2.flux_norm.values - 1, d2.flux_err.values, d2.cadence_s.values / 86400, -1)
D = {s: load(s) for s in SECTORS}
t = np.concatenate([D[s]["t"] for s in SECTORS]); f = np.concatenate([D[s]["aM"] for s in SECTORS])
e = np.concatenate([window_sigma(D[s]) for s in SECTORS]); ex = np.concatenate([np.full(len(D[s]["t"]), D[s]["exptime"]) for s in SECTORS])
out["referee_lc"] = fit(t, f, e, ex, -1.0 if False else -1.0) if False else None
def fit_e(t, f, e, ex):
    return fit(t, f, e, ex, +1)  # referee flux change is negative in transit: f = -dep*y -> use sign handling below
# referee flux change aM: f = -dep*y; pass -f so that the same 'f - sign*dep*y' form applies with sign=+1
out["referee_lc"] = fit(t, -f, e, ex, +1)
for k, v in out.items():
    print(f"{k}: T14 {v['T14_h']:.2f}+-{v['T14_err_h']:.2f} h, flat {v['flat']:.2f}+-{v['flat_err']:.2f}, depth {v['depth']:.4g}, T0 {v['T0']:.5f}, P {v['P']:.7f}; flat profile (flat, dchi2, T14): {[(a, round(b,1), round(c,2)) for a,b,c in v['profile']]}")
save_json(out, f"{R}/out/shape_crosscheck.json")
