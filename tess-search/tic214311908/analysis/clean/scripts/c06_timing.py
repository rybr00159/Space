"""c06: per-transit depths and mid-times on the cleaned light curves of both pipelines, injection-calibrated errors,
TTV chi2, odd/even, and the FINAL linear ephemeris (each pipeline and their combination).

Per-transit fit (identical to the one used on the pixel-level injections, c05.fit_window): fixed batman shape of the
adopted solution (k 0.1195, b 0.52, a/R* 25.4, u 0.125/0.533), free depth (mid-transit units) and mid-time (grid
+-0.05 d, step 0.001 d, parabolic refinement), linear baseline, sigma = window OOT robust std.
Calibration (c05 pixel injections, matched templates: invprf for INV, refprf for REF; every fake transit of the sector):
  timing:  sigma_cal = tc_err x std(pulls) of the sector (pull = (tc_fit - tc_true) / tc_err); robust std and the
           fraction of > 3-sigma outliers are reported.  For S12 (30-min) the per-sector pull std is used likewise.
  depth:   sigma_cal = std of the no-injection depths at the same fake epochs (empirical single-transit error).
TTV: chi2 of the times about the best linear ephemeris with calibrated errors; p from chi2(n-2) and from 20000 Monte
Carlo draws of the empirical pull distribution.
Combination: per event, inverse-variance average of the two pipelines' times (variance including the INV-REF
correlation of the timing errors at the same fake transits, c05), then a weighted linear fit (adopted); a full GLS is
reported as a variant.
Outputs: out/per_transit_clean.csv, out/timing_clean.json, out/ephemeris_final.json
"""
import sys
sys.path.insert(0, "tic214311908/analysis/clean/scripts")
from cleanlib import CL, P1, T01, F_M3V, SECTORS, robust_std, jd
from c05_inject import fit_window, accept
from c03_windows import load_windows
from scipy import stats
import numpy as np, pandas as pd, json

rng = np.random.default_rng(606)


def calib():
    """per sector, per pipeline: timing pull stats and empirical depth error from c05; INV-REF timing correlation."""
    out = {}
    for s in SECTORS:
        d = pd.read_csv(f"{CL}/out/inject_pix_s{s}.csv")
        c = {}
        for pipe, var, pre in (("inv", "invprf", "inv_raw_"), ("ref", "refprf", "ref_")):
            x = d[d.variant == var]
            pull = (x[f"{pre}dt_min"] / (x[f"{pre}tc_err"] * 1440)).values
            pull = pull[np.isfinite(pull)]
            dn = x[f"{pre}depth_fix_noinj"].values; dn = dn[np.isfinite(dn)]
            rec = x[f"{pre}depth_fix"].values; rec = rec[np.isfinite(rec)]
            c[pipe] = dict(n=int(len(pull)), pull_std=float(np.std(pull)), pull_rstd=float(robust_std(pull)),
                           pull_out3=float(np.mean(np.abs(pull) > 3)), pull_mean=float(np.mean(pull)), pulls=pull.tolist(),
                           depth_err_emp=float(np.std(dn)), depth_noinj_mean=float(np.mean(dn)),
                           depth_rec_mean=float(np.mean(rec)), depth_rec_std=float(np.std(rec)),
                           dt_std_min=float(np.std(x[f"{pre}dt_min"].dropna())))
        a = d[d.variant == "invprf"].set_index("tc"); b = d[d.variant == "refprf"].set_index("tc")
        j = a.join(b, lsuffix="_a", rsuffix="_b")
        dta, dtb = j["inv_raw_dt_min_a"], j["ref_dt_min_b"]; m = np.isfinite(dta) & np.isfinite(dtb)
        c["rho_dt"] = float(np.corrcoef(dta[m], dtb[m])[0, 1]) if m.sum() > 4 else 0.0
        da, db = j["inv_raw_depth_fix_noinj_a"], j["ref_depth_fix_noinj_b"]; m2 = np.isfinite(da) & np.isfinite(db)
        c["rho_depth"] = float(np.corrcoef(da[m2], db[m2])[0, 1]) if m2.sum() > 4 else 0.0
        c["n_pair"] = int(m.sum())
        out[str(s)] = c
    return out


def per_transit(W, pipe, CAL):
    rows = []
    for i, ii in enumerate(W.idx):
        s = int(W.sector[i]); tc0 = float(W.tc0[i])
        r = fit_window(W.t[ii], W.f[ii], tc0, float(W.ex[ii][0]))
        if r is None: continue
        c = CAL[str(s)][pipe]
        rows.append(dict(pipeline=pipe, sector=s, epoch=int(W.epoch[i]), tc_pred=tc0, tc=r["tc_free"], tc_err_formal=r["tc_err"],
                         tc_err_cal=r["tc_err"] * max(c["pull_std"], 1.0), depth_ppt=r["depth_fix"], depth_err_formal=r["depth_fix_err"],
                         depth_err_cal=max(c["depth_err_emp"], r["depth_fix_err"]), depth_free_ppt=r["depth_free"], sig_ppt=r["sig"] * 1e3,
                         n_pts=int(len(ii))))
    return pd.DataFrame(rows)


def lin_eph(E, t, e, Eref=0):
    A = np.column_stack([np.ones(len(E)), E - Eref]) / e[:, None]
    c, *_ = np.linalg.lstsq(A, t / e, rcond=None)
    cov = np.linalg.inv(A.T @ A); r = t / e - A @ c
    return c, cov, float(r @ r), r


def ttv_mc(E, sectors, e, CAL, pipe, n=20000):
    """chi2 distribution under the null using the empirical pulls (resampled per sector) scaled to the calibrated errors."""
    out = np.empty(n)
    pulls = {s: np.asarray(CAL[str(s)][pipe]["pulls"]) / max(CAL[str(s)][pipe]["pull_std"], 1.0) for s in set(sectors)}
    for k in range(n):
        x = np.array([rng.choice(pulls[s]) for s in sectors]) * e
        _, _, c2, _ = lin_eph(E, x, e)
        out[k] = c2
    return out


def gls(E, t, C, Eref):
    A = np.column_stack([np.ones(len(E)), E - Eref])
    Ci = np.linalg.inv(C)
    cov = np.linalg.inv(A.T @ Ci @ A); c = cov @ A.T @ Ci @ t
    r = t - A @ c
    return c, cov, float(r @ Ci @ r)


def main():
    CAL = calib()
    Wi = load_windows(f"{CL}/data/win_inv_raw.npz"); Wr = load_windows(f"{CL}/data/win_ref.npz")
    D = pd.concat([per_transit(Wi, "inv", CAL), per_transit(Wr, "ref", CAL)], ignore_index=True)
    out = dict(calibration={s: {p: {k: v for k, v in CAL[s][p].items() if k != "pulls"} if isinstance(CAL[s][p], dict) else CAL[s][p]
                                for p in CAL[s]} for s in CAL})
    res = {}
    for pipe in ("inv", "ref"):
        d = D[(D.pipeline == pipe) & np.isfinite(D.tc) & np.isfinite(D.tc_err_cal)].copy()
        E = d.epoch.values.astype(float); t = d.tc.values; e = d.tc_err_cal.values; ef = d.tc_err_formal.values
        w = 1 / e ** 2; Eref = float(np.round(np.sum(w * E) / np.sum(w)))
        c, cov, chi2, r = lin_eph(E, t, e, Eref)
        c0, cov0, _, _ = lin_eph(E, t, e, 0)
        _, _, chi2f, _ = lin_eph(E, t, ef, 0)
        er = np.array([ef[i] * max(CAL[str(sv)][pipe]["pull_rstd"], 1.0) for i, sv in enumerate(d.sector.values)])
        _, _, chi2r, rr_ = lin_eph(E, t, er, 0)
        n_out3 = int(np.sum(np.abs(rr_) > 3))
        n = len(E); dof = n - 2
        mc = ttv_mc(E, d.sector.values, e, CAL, pipe)
        scale = max(1.0, np.sqrt(chi2 / dof))
        d["oc_min"] = (t - (c0[0] + c0[1] * E)) * 1440
        D.loc[d.index, "oc_min"] = d["oc_min"]
        # depths: weighted mean, chi2, odd/even
        dep = d.depth_ppt.values; de = d.depth_err_cal.values
        wm = np.sum(dep / de ** 2) / np.sum(1 / de ** 2); wme = 1 / np.sqrt(np.sum(1 / de ** 2))
        chi2d = float(np.sum(((dep - wm) / de) ** 2))
        odd = (d.epoch.values % 2) != 0
        def wmean(m): return np.sum(dep[m] / de[m] ** 2) / np.sum(1 / de[m] ** 2), 1 / np.sqrt(np.sum(1 / de[m] ** 2))
        do, doe = wmean(odd); dv, dve = wmean(~odd)
        res[pipe] = dict(n=n, Eref=Eref, T0_ref=float(c[0]), P=float(c[1]), cov_ref=cov.tolist(), T0_ref_err=float(np.sqrt(cov[0, 0])),
                         P_err=float(np.sqrt(cov[1, 1])), T0_E0=float(c0[0]), T0_E0_err=float(np.sqrt(cov0[0, 0])), cov_E0=cov0.tolist(),
                         chi2_cal=chi2, chi2_formal=chi2f, chi2_robust_scaled=chi2r, n_gt3sigma_robust=n_out3, dof=dof, p_chi2=float(stats.chi2.sf(chi2, dof)), p_mc=float(np.mean(mc >= chi2)),
                         oc_rms_min=float(np.std(d.oc_min)), err_scale_applied=float(scale),
                         depth_wmean_ppt=float(wm), depth_wmean_err=float(wme), depth_chi2=chi2d, depth_dof=n - 1,
                         p_depth_chi2=float(stats.chi2.sf(chi2d, n - 1)), odd_ppt=[float(do), float(doe)], even_ppt=[float(dv), float(dve)],
                         odd_minus_even=[float(do - dv), float(np.hypot(doe, dve))])
        print(f"{pipe}: n {n}; P {c[1]:.7f}+-{np.sqrt(cov[1,1]):.7f}; T0(E{Eref:.0f}) {c[0]:.5f}+-{np.sqrt(cov[0,0]):.5f}; TTV chi2 {chi2:.1f}/{dof} "
              f"(p {stats.chi2.sf(chi2, dof):.2f}, MC {np.mean(mc >= chi2):.2f}; formal {chi2f:.1f}; robust-scaled {chi2r:.1f}, {n_out3} > 3 sigma); O-C rms {np.std(d.oc_min):.1f} min; "
              f"depth {wm:.2f}+-{wme:.2f} chi2 {chi2d:.1f}/{n-1}; odd-even {do-dv:+.2f}+-{np.hypot(doe,dve):.2f}", flush=True)
    # ---- combination. (i) ADOPTED: per event, inverse-variance average of the INV and REF times (non-negative
    # weights), with the variance of the average including the INV-REF error correlation rho (per sector, c05);
    # then a weighted linear fit. (ii) Variant: full GLS of all 39 measurements with the 2x2 per-event covariance
    # (optimal under the Gaussian model, but with rho > sigma_REF / sigma_INV it gives INV negative weight, which makes
    # it fragile to mis-estimated rho; reported only).
    a = D[D.pipeline == "inv"].set_index("epoch"); b = D[D.pipeline == "ref"].set_index("epoch")
    rows = []; ev = []
    for ep in sorted(set(a.index) | set(b.index)):
        m = []
        for pipe, tab in (("inv", a), ("ref", b)):
            if ep in tab.index and np.isfinite(tab.loc[ep, "tc"]):
                m.append((pipe, int(tab.loc[ep, "sector"]), float(tab.loc[ep, "tc"]), float(tab.loc[ep, "tc_err_cal"])))
                rows.append((ep, pipe, m[-1][1], m[-1][2], m[-1][3]))
        if len(m) == 1:
            ev.append((ep, m[0][2], m[0][3], 1))
        else:
            (_, sct, t1, e1), (_, _, t2, e2) = m
            w1, w2 = 1 / e1 ** 2, 1 / e2 ** 2; rho = CAL[str(sct)]["rho_dt"]
            tc = (w1 * t1 + w2 * t2) / (w1 + w2)
            var = (w1 ** 2 * e1 ** 2 + w2 ** 2 * e2 ** 2 + 2 * w1 * w2 * rho * e1 * e2) / (w1 + w2) ** 2
            ev.append((ep, tc, np.sqrt(var), 2))
    Ee = np.array([x[0] for x in ev], float); te = np.array([x[1] for x in ev]); ee = np.array([x[2] for x in ev])
    we = 1 / ee ** 2; Eref = float(np.round(np.sum(we * Ee) / np.sum(we)))
    c, cov, chi2, _ = lin_eph(Ee, te, ee, Eref); c0, cov0, _, _ = lin_eph(Ee, te, ee, 0)
    dof = len(Ee) - 2
    res["combined"] = dict(method="per-event inverse-variance average (non-negative weights) + WLS", n_events=len(Ee),
                           n_both=int(sum(x[3] == 2 for x in ev)), Eref=Eref, T0_ref=float(c[0]), P=float(c[1]), cov_ref=cov.tolist(),
                           T0_ref_err=float(np.sqrt(cov[0, 0])), P_err=float(np.sqrt(cov[1, 1])), T0_E0=float(c0[0]),
                           T0_E0_err=float(np.sqrt(cov0[0, 0])), cov_E0=cov0.tolist(), chi2=chi2, dof=dof, p_chi2=float(stats.chi2.sf(chi2, dof)),
                           rho_dt_by_sector={s: CAL[s]["rho_dt"] for s in CAL},
                           events=[dict(epoch=int(x[0]), tc=float(x[1]), err=float(x[2]), n_pipelines=int(x[3])) for x in ev])
    print(f"combined (per-event average): {len(Ee)} events; P {c[1]:.7f}+-{np.sqrt(cov[1,1]):.7f}; T0(E{Eref:.0f}) {c[0]:.5f}+-{np.sqrt(cov[0,0]):.5f}; "
          f"T0(E0) {c0[0]:.5f}+-{np.sqrt(cov0[0,0]):.5f}; chi2 {chi2:.1f}/{dof}", flush=True)
    E = np.array([r[0] for r in rows], float); t = np.array([r[3] for r in rows]); e = np.array([r[4] for r in rows])
    C = np.diag(e ** 2)
    for i in range(len(rows)):
        for j in range(i + 1, len(rows)):
            if rows[i][0] == rows[j][0] and rows[i][1] != rows[j][1]:
                C[i, j] = C[j, i] = CAL[str(rows[i][2])]["rho_dt"] * e[i] * e[j]
    cg, covg, chi2g = gls(E, t, C, Eref); cg0, covg0, _ = gls(E, t, C, 0.0)
    res["combined_gls_variant"] = dict(Eref=Eref, T0_ref=float(cg[0]), P=float(cg[1]), T0_ref_err=float(np.sqrt(covg[0, 0])),
                                       P_err=float(np.sqrt(covg[1, 1])), T0_E0=float(cg0[0]), T0_E0_err=float(np.sqrt(covg0[0, 0])),
                                       chi2=chi2g, dof=len(E) - 2)
    print(f"combined GLS variant: P {cg[1]:.7f}+-{np.sqrt(covg[1,1]):.7f}; T0(E0) {cg0[0]:.5f}+-{np.sqrt(covg0[0,0]):.5f}", flush=True)
    out.update(res)
    D.to_csv(f"{CL}/out/per_transit_clean.csv", index=False, float_format="%.6f")
    jd(out, f"{CL}/out/timing_clean.json")
    # ---- final ephemeris file
    fin = {}
    for k in ("inv", "ref", "combined"):
        r = res[k]
        sc = max(1.0, np.sqrt(r.get("chi2_cal", r.get("chi2")) / r["dof"]))
        covr = np.array(r["cov_ref"]) * sc ** 2; cov0 = np.array(r["cov_E0"]) * sc ** 2
        fin[k] = dict(reference_epoch=int(r["Eref"]),
                      T0_ref_BTJD_TDB=r["T0_ref"], T0_ref_BJD_TDB=r["T0_ref"] + 2457000.0, T0_ref_err_d=float(np.sqrt(covr[0, 0])),
                      P_d=r["P"], P_err_d=float(np.sqrt(covr[1, 1])), cov_T0ref_P=covr.tolist(),
                      corr_T0ref_P=float(covr[0, 1] / np.sqrt(covr[0, 0] * covr[1, 1])),
                      T0_epoch0_BTJD_TDB=r["T0_E0"], T0_epoch0_BJD_TDB=r["T0_E0"] + 2457000.0, T0_epoch0_err_d=float(np.sqrt(cov0[0, 0])),
                      cov_T0epoch0_P=cov0.tolist(), chi2=float(r.get("chi2_cal", r.get("chi2"))), dof=int(r["dof"]), error_scale_applied=float(sc))
    fin["adopted"] = "combined"
    fin["notes"] = ("Linear ephemeris T(E) = T0 + P (E - E_ref) from per-transit mid-times of the cleaned light curves (variables excluded "
                    "from field-wide estimators), errors calibrated by pixel-level injections (c05); BTJD_TDB = BJD_TDB - 2457000, "
                    "target-barycentred. Epoch 0 = the phase-1 numbering (T ~ 3840.914 BTJD). Errors scaled by sqrt(chi2/dof) when > 1. "
                    "Combined = per event the inverse-variance average of the INV and REF times, its variance including the per-sector INV-REF "
                    "timing-error correlation measured at the same fake transits (c05), then a weighted linear fit.")
    jd(fin, f"{CL}/out/ephemeris_final.json")


if __name__ == "__main__":
    main()
