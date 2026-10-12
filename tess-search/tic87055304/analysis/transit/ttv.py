"""Linear ephemeris, TTV statistics, period-derivative limit, O-C plot and follow-up predictions for TIC 87055304.

Input: transit_fits_<tag>.json (per-transit mid-times with formal errors and per-sector injection calibration).
Calibrated error per transit = formal * max(1, s_sector), where s_sector is the injection scale factor
(std of (recovered - injected)/formal; robust MAD version reported too).
Ephemeris T(E) = T0 + P (E - E_ref) by weighted least squares; chi^2/dof and p-value of residuals; also a
quadratic ephemeris T = T0 + P m + 0.5 dP/dE m^2 for the period-derivative limit (dP/dt = (dP/dE)/P).
Odd/even timing offset and per-sector mean O-C are reported.
Predictions for 2026-11-01 .. 2028-12-31 (BJD_TDB), with sigma from the full T0-P covariance (also from the MCMC
posterior, if available; the larger sigma is adopted).
Usage: python ttv.py [tag]
"""
import os, sys, json
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.time import Time

HERE = os.path.dirname(os.path.abspath(__file__))
tag = sys.argv[1] if len(sys.argv) > 1 else "raw"
js = json.load(open(os.path.join(HERE, f"transit_fits_{tag}.json")))
st = json.load(open(os.path.join(HERE, "ml_start.json")))
E_REF = st["E_ref"]
T0_tmpl, P_tmpl = js["template"][0], js["template"][1]
# epochs in the template are relative to the template T0 (which is at E_REF); convert to absolute epoch numbers
rows = [r for r in js["transits"] if r.get("used") and "time_tc" in r and np.isfinite(r.get("time_err", np.nan)) and not r.get("time_at_edge")]
scale_key = "time_scale_sd"
inj = js["injection"]
E = np.array([r["epoch"] + E_REF for r in rows], float)
sec = np.array([r["sector"] for r in rows])
tc = np.array([r["time_tc"] for r in rows])
ef = np.array([r["time_err"] for r in rows])
sc = np.array([max(1.0, inj[str(s)][scale_key]) for s in sec])
e = ef * sc


def linfit(E, tc, e, Eref):
    m = E - Eref
    A = np.vstack([np.ones_like(m), m]).T
    W = np.diag(1 / e ** 2)
    cov = np.linalg.inv(A.T @ W @ A)
    p = cov @ A.T @ W @ tc
    r = tc - A @ p
    chi2 = float(np.sum((r / e) ** 2))
    return p, cov, r, chi2


p, cov, r, chi2 = linfit(E, tc, e, E_REF)
dof = len(tc) - 2
pval = float(stats.chi2.sf(chi2, dof))
# formal-error-only version
pf, covf, rf, chi2f = linfit(E, tc, ef, E_REF)
# quadratic ephemeris
m = E - E_REF
A = np.vstack([np.ones_like(m), m, 0.5 * m ** 2]).T
Wm = np.diag(1 / e ** 2)
covq = np.linalg.inv(A.T @ Wm @ A); pq = covq @ A.T @ Wm @ tc
chi2q = float(np.sum(((tc - A @ pq) / e) ** 2))
dPdE, sdPdE = pq[2], np.sqrt(covq[2, 2])
dPdt = dPdE / p[1]; sdPdt = sdPdE / p[1]  # dimensionless (d/d)
ms_per_yr = 86400e3 * 365.25
# odd/even timing
odd = (E % 2 == 1)
def wmean(x, w):
    return np.sum(x * w) / np.sum(w), 1 / np.sqrt(np.sum(w))
mo, so = wmean(r[odd], 1 / e[odd] ** 2); me, se = wmean(r[~odd], 1 / e[~odd] ** 2)
# per-sector mean O-C
persec = {}
for s in np.unique(sec):
    k = sec == s
    mm, ss = wmean(r[k], 1 / e[k] ** 2)
    persec[int(s)] = dict(n=int(k.sum()), mean_OC_min=mm * 1440, err_min=ss * 1440,
                         chi2=float(np.sum((r[k] / e[k]) ** 2)), rms_min=float(np.std(r[k]) * 1440))
chi2_sec = float(sum(((v["mean_OC_min"]) / v["err_min"]) ** 2 for v in persec.values()))

out = dict(tag=tag, E_ref=E_REF, n_transits=int(len(tc)), error_scale_key=scale_key,
           sector_scale={s: max(1.0, inj[s][scale_key]) for s in inj},
           P=float(p[1]), P_err=float(np.sqrt(cov[1, 1])), T0_BTJD=float(p[0]), T0_err=float(np.sqrt(cov[0, 0])),
           T0_BJD=float(p[0] + 2457000), cov_T0_P=float(cov[0, 1]),
           chi2=chi2, dof=dof, chi2_red=chi2 / dof, p_value=pval, rms_OC_min=float(np.std(r) * 1440),
           chi2_formal_errors=chi2f, chi2_red_formal=chi2f / dof,
           max_abs_OC_sigma=float(np.max(np.abs(r / e))),
           quadratic=dict(dPdE_d=float(dPdE), dPdE_err_d=float(sdPdE), dPdt=float(dPdt), dPdt_err=float(sdPdt),
                          dPdt_ms_per_yr=float(dPdt * ms_per_yr), dPdt_err_ms_per_yr=float(sdPdt * ms_per_yr),
                          abs_dPdt_95_ms_per_yr=float((abs(dPdt) + 1.96 * sdPdt) * ms_per_yr),
                          chi2=chi2q, dchi2_vs_linear=chi2 - chi2q),
           odd_even_timing=dict(odd_mean_OC_min=mo * 1440, odd_err_min=so * 1440, even_mean_OC_min=me * 1440, even_err_min=se * 1440,
                                diff_min=(mo - me) * 1440, diff_sigma=float((mo - me) / np.hypot(so, se))),
           per_sector=persec, chi2_sector_means=chi2_sec,
           transits=[dict(epoch=int(Ei), sector=int(s), tc_BTJD=float(t), err_formal_min=float(a * 1440), err_cal_min=float(b * 1440),
                          OC_min=float(x * 1440), OC_sigma=float(x / b)) for Ei, s, t, a, b, x in zip(E, sec, tc, ef, e, r)])
# MCMC ephemeris for comparison and predictions
mc_path = os.path.join(HERE, "mcmc_free.json")
if os.path.exists(mc_path):
    mcj = json.load(open(mc_path)); ch = np.load(os.path.join(HERE, "mcmc_free.npz"))["chain"]
    out["mcmc_ephemeris"] = dict(T0=float(np.median(ch[:, 0])), T0_err=float(np.std(ch[:, 0])), P=float(np.median(ch[:, 1])),
                                 P_err=float(np.std(ch[:, 1])), cov=float(np.cov(ch[:, 0], ch[:, 1])[0, 1]))
    out["timing_vs_mcmc"] = dict(dT0_sigma=float((p[0] - np.median(ch[:, 0])) / np.hypot(np.sqrt(cov[0, 0]), np.std(ch[:, 0]))),
                                 dP_sigma=float((p[1] - np.median(ch[:, 1])) / np.hypot(np.sqrt(cov[1, 1]), np.std(ch[:, 1]))))
# predictions
t1 = Time("2026-11-01T00:00:00", scale="tdb").jd - 2457000; t2 = Time("2028-12-31T23:59:59", scale="tdb").jd - 2457000
n1 = int(np.ceil((t1 - p[0]) / p[1])); n2 = int(np.floor((t2 - p[0]) / p[1]))
pred = []
for n in range(n1, n2 + 1):
    T = p[0] + n * p[1]
    sT = np.sqrt(cov[0, 0] + n ** 2 * cov[1, 1] + 2 * n * cov[0, 1])
    row = dict(epoch=int(n + E_REF), BJD_TDB=float(T + 2457000), sigma_min=float(sT * 1440))
    if "mcmc_ephemeris" in out:
        me_ = out["mcmc_ephemeris"]; sM = np.sqrt(me_["T0_err"] ** 2 + n ** 2 * me_["P_err"] ** 2 + 2 * n * me_["cov"])
        row["BJD_TDB_mcmc"] = float(me_["T0"] + n * me_["P"] + 2457000); row["sigma_mcmc_min"] = float(sM * 1440)
        row["sigma_adopted_min"] = float(max(sT, sM) * 1440)
    row["UTC_ISO"] = Time(T + 2457000, format="jd", scale="tdb").utc.iso[:16]
    pred.append(row)
import csv
with open(os.path.join(HERE, f"predictions_2026-11_2028-12.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(pred[0].keys())); w.writeheader(); w.writerows(pred)
out["predictions_summary"] = dict(n=len(pred), first=pred[0], last=pred[-1],
                                  sigma_range_min=[pred[0].get("sigma_adopted_min", pred[0]["sigma_min"]), pred[-1].get("sigma_adopted_min", pred[-1]["sigma_min"])],
                                  T14_h=js["template_T14_h"])
json.dump(out, open(os.path.join(HERE, f"ttv_{tag}.json"), "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k not in ("transits",)}, indent=1))

# O-C figure
fig, ax = plt.subplots(2, 1, figsize=(10, 6.5), gridspec_kw=dict(height_ratios=[2, 1]))
cols = {13: "C0", 39: "C1", 66: "C2", 93: "C3", 104: "C4"}
for s in np.unique(sec):
    k = sec == s
    ax[0].errorbar(tc[k], r[k] * 1440, e[k] * 1440, fmt="o", color=cols[int(s)], ms=4, label=f"S{s} (n={k.sum()})")
ax[0].axhline(0, color="k", lw=0.6)
tt = np.linspace(tc.min() - 50, tc.max() + 50, 300); mm = (tt - p[0]) / p[1]
band = np.sqrt(cov[0, 0] + mm ** 2 * cov[1, 1] + 2 * mm * cov[0, 1]) * 1440
ax[0].fill_between(tt, -band, band, color="0.8", alpha=0.6, label="linear ephemeris 1-sigma")
ax[0].set_ylabel("O - C (min)"); ax[0].legend(fontsize=7, ncol=3)
ax[0].set_title(f"TIC 87055304: P = {p[1]:.7f} +- {np.sqrt(cov[1, 1]):.7f} d, T0 = {p[0] + 2457000:.5f} +- {np.sqrt(cov[0, 0]):.5f} BJD_TDB; "
                f"chi2/dof = {chi2:.1f}/{dof}", fontsize=9)
ax[0].set_xlabel("BTJD")
ax[1].plot(tc, r / e, "ko", ms=3); ax[1].axhline(0, color="k", lw=0.6)
for y in (-2, 2): ax[1].axhline(y, color="0.6", ls=":")
ax[1].set_ylabel("O - C / sigma"); ax[1].set_xlabel("BTJD")
fig.tight_layout(); fig.savefig(os.path.join(HERE, f"oc_{tag}.png"), dpi=130)
