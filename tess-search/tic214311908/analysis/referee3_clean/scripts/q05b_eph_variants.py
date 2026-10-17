"""q05b: sensitivity of the linear ephemeris to the timing-error model (same per-transit times as the analyst, q05
reproduces them to 0.1 min): (a) analyst's per-event calibrated errors; (b) sector-uniform std from q05 LC-level
injections; (c) sector-uniform robust std; (d) analyst errors with a 5-min floor; (e) jackknife over sectors;
(f) Monte Carlo prediction spread for 2027 (E=130) across variants.  Output: out/q05b_eph_variants.json"""
import json, numpy as np, pandas as pd
R3 = "tic214311908/analysis/referee3_clean"
CL = "tic214311908/analysis/clean"
a = pd.read_csv(f"{CL}/out/per_transit_clean.csv"); cal = json.load(open(f"{R3}/out/q05_timing.json"))["calibration"]
def wls(E, T, s):
    A = np.column_stack([np.ones(len(E)), E]) / s[:, None]; c, *_ = np.linalg.lstsq(A, T / s, rcond=None)
    cov = np.linalg.inv(A.T @ A); r = (T - c[0] - c[1] * E) / s; return c, cov, float(r @ r)
def pred(c, cov, e=130):
    return c[0] + c[1] * e, np.sqrt(cov[0, 0] + e * e * cov[1, 1] + 2 * e * cov[0, 1]) * 1440
out = {}
for pl in ("inv", "ref"):
    g = a[a.pipeline == pl]; E = g.epoch.values.astype(float); T = g.tc.values
    errs = {"analyst": g.tc_err_cal.values, "q05_std": np.array([cal[f"{pl}_{s}"]["std_min"] for s in g.sector]) / 1440,
            "q05_rstd": np.array([cal[f"{pl}_{s}"]["rstd_min"] for s in g.sector]) / 1440,
            "analyst_floor5": np.maximum(g.tc_err_cal.values, 5 / 1440)}
    for k, s in errs.items():
        c, cov, chi2 = wls(E, T, s); Tp, sp = pred(c, cov)
        out[f"{pl}_{k}"] = dict(T0=c[0], T0_err=np.sqrt(cov[0, 0]), P=c[1], P_err=np.sqrt(cov[1, 1]), chi2=chi2, dof=len(E) - 2, T130=Tp, sig130_min=sp)
    # jackknife over sectors with analyst errors
    Ps = []
    for s in np.unique(g.sector):
        m = g.sector.values != s; c, cov, _ = wls(E[m], T[m], errs["analyst"][m]); Ps.append(c[1])
    Ps = np.array(Ps); n = len(Ps)
    out[f"{pl}_jackknife_P_err"] = float(np.sqrt((n - 1) / n * np.sum((Ps - Ps.mean()) ** 2)))
    out[f"{pl}_jackknife_P_values"] = Ps.tolist()
ref130 = out["ref_analyst"]["T130"]
for k, v in out.items():
    if isinstance(v, dict):
        print(f"{k:22s} T0 {v['T0']:.5f}+-{v['T0_err']:.5f} P {v['P']:.7f}+-{v['P_err']:.7f} chi2 {v['chi2']:.1f}/{v['dof']} T(E130)-REFanalyst {(v['T130']-ref130)*1440:+.1f} min, sig {v['sig130_min']:.1f}")
    else:
        print(k, np.round(v, 8) if isinstance(v, float) else v)
json.dump(out, open(f"{R3}/out/q05b_eph_variants.json", "w"), indent=1, default=float)
