"""Transit timing from the deblended M3V light curve, error calibration by light-curve-level injections, and
comparison with the QLP event times (ephemeris_v2.json, S103/S104 corrected) that gave chi2 = 53.4 / 11.

1. Per-transit times (p08, model B, strict quality, coverage >= 50%): linear ephemeris fit, chi2, O-C.
2. Error calibration: in every sector, 60 synthetic transits (global depth/shape of p08) are multiplied into the
   detrended M3V light curve at random epochs away from the real transits and data gaps; each is fitted exactly like
   a real transit; pulls (t_fit - t_true) / sigma give the factor by which the per-transit errors must be scaled.
3. QLP: the same epochs measured on QLP SAP flux of TIC 214311908 (its own detrending, corrected times), and the
   published per-event QLP times from ephemeris_v2.json; differences to our times.
Outputs: out/timing.json, out/timing_table.csv, figs/oc.png
"""
import sys; sys.path.insert(0, "tic214311908/analysis/phot/scripts")
from photlib import *
from tlib import trap_frac, fit_single
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
RNG = np.random.default_rng(11)
TR = json.load(open(f"{ROOT}/out/transits.json"))["global_B"]
P, T0, T14, FL, DEP = TR["P"], TR["T0"], TR["T14_d"], TR["flat_frac"], TR["depth"]
EPH2 = json.load(open("/home/user/Space/tess-search/vetting_2026-10-16/results/214311906/ephemeris_v2.json"))


def linfit(k, tc, err):
    W = 1 / err; A = np.column_stack([np.ones(len(k)), k])
    c, *_ = np.linalg.lstsq(A * W[:, None], tc * W, rcond=None)
    cov = np.linalg.pinv((A * W[:, None]).T @ (A * W[:, None])); r = tc - A @ c
    return c, cov, r, float(np.sum((r / err) ** 2))


def main():
    PT = pd.read_csv(f"{ROOT}/out/per_transit_M3V.csv")
    out = {}
    M = pd.read_csv(f"{ROOT}/lc_M3V_all.csv")
    # ---- 2. injection calibration of timing errors (per sector)
    calib = {}
    for s in SECTORS:
        d = M[(M.sector == s) & (M.quality_ok == 1) & (M.flare_flag == 0)]
        t = d.time_btjd_tdb.values; f = d.flux_norm.values; ex = d.cadence_s.values[0] / 86400
        sd = robust_std(np.diff(f)) / np.sqrt(2)
        realw = np.abs((((t - T0) / P + 0.5) % 1 - 0.5) * P) < 0.35
        pulls, dpulls = [], []
        tries = 0
        while len(pulls) < 60 and tries < 600:
            tries += 1
            tc = RNG.uniform(t.min() + 0.3, t.max() - 0.3)
            if np.min(np.abs((((tc - T0) / P + 0.5) % 1 - 0.5) * P)) < 0.6: continue
            w = np.abs(t - tc) < 0.3
            if np.sum(np.abs(t - tc) < T14 / 2) < 0.8 * T14 / ex or np.any(realw[w]): continue
            f2 = f * (1 - DEP * trap_frac(t, tc, T14, FL, ex))
            tc0 = tc + RNG.uniform(-0.01, 0.01)      # start offset (as for real events, predicted time differs)
            r = fit_single(t[w], f2[w], np.full(w.sum(), sd), ex, tc0, T14, FL, win=0.3, dt_max=0.08)
            if r is None: continue
            e = max(r["tc_err"], r["tc_err_pb"])
            pulls.append((r["tc"] - tc) / e); dpulls.append((r["depth"] - DEP) / max(r["depth_err"], r["depth_err_pb"]))
        pulls, dpulls = np.array(pulls), np.array(dpulls)
        calib[s] = dict(n=len(pulls), pull_std=float(np.std(pulls)), pull_rstd=float(robust_std(pulls)), pull_mean=float(np.mean(pulls)),
                        frac_gt3=float(np.mean(np.abs(pulls) > 3)), depth_pull_std=float(np.std(dpulls)), depth_pull_mean=float(np.mean(dpulls)))
        print(f"S{s}: timing pulls n={len(pulls)} std {np.std(pulls):.2f} robust {robust_std(pulls):.2f} mean {np.mean(pulls):+.2f} |>3| {np.mean(np.abs(pulls)>3):.2f}; "
              f"depth pulls std {np.std(dpulls):.2f} mean {np.mean(dpulls):+.2f}", flush=True)
    out["injection_calibration"] = calib
    # ---- 1. linear ephemeris on our times
    sel = PT[(PT.model == "B") & (~PT.relaxed) & (PT.cov_frac >= 0.5)].copy()
    sel["err"] = sel[["tc_err", "tc_err_pb"]].max(axis=1)
    sel["err_cal"] = sel.err * sel.sector.map(lambda s: max(1.0, calib[s]["pull_std"]))
    for lab, d in (("all", sel), ("snr>=3", sel[sel.snr >= 3])):
        for ecol in ("err", "err_cal"):
            c, cov, r, chi = linfit(d.epoch.values, d.tc.values, d[ecol].values)
            dof = len(d) - 2
            key = f"{lab}_{ecol}"
            out[key] = dict(n=len(d), T0=c[0], P=c[1], T0_err=float(np.sqrt(cov[0, 0])), P_err=float(np.sqrt(cov[1, 1])), chi2=chi, dof=dof,
                            p=float(__import__("scipy").stats.chi2.sf(chi, dof)), oc_rms_min=float(np.std(r) * 1440),
                            worst=[(int(e), float(x * 1440), float(x / er)) for e, x, er in sorted(zip(d.epoch, r, d[ecol]), key=lambda z: -abs(z[1] / z[2]))[:3]])
            print(f"ours [{key}]: n {len(d)}, P {c[1]:.7f}+-{np.sqrt(cov[1,1]):.7f}, T0 {c[0]:.5f}+-{np.sqrt(cov[0,0]):.5f}, chi2 {chi:.1f}/{dof} "
                  f"(p {out[key]['p']:.3g}), O-C rms {np.std(r)*1440:.1f} min; worst (epoch, min, sigma) {[(w[0], round(w[1],1), round(w[2],1)) for w in out[key]['worst']]}")
    # ---- 3. QLP comparison
    q = pd.DataFrame(EPH2["events"])
    q = q.rename(columns={"k": "epoch_qlp", "tc": "tc_qlp", "err": "err_qlp"})
    # QLP epochs are counted from the QLP T0 (same as ours to within 0.002 d)
    q["epoch"] = np.round((q.tc_qlp - T0) / P).astype(int)
    mm = sel.merge(q[["epoch", "tc_qlp", "err_qlp", "depth", "snr"]].rename(columns={"depth": "depth_qlp", "snr": "snr_qlp"}), on="epoch", how="outer")
    mm["diff_min"] = (mm.tc - mm.tc_qlp) * 1440
    mm["diff_sig"] = (mm.tc - mm.tc_qlp) / np.hypot(mm.err_cal, mm.err_qlp)
    mm["oc_ours_min"] = (mm.tc - (T0 + mm.epoch * P)) * 1440
    mm["oc_qlp_min"] = (mm.tc_qlp - (T0 + mm.epoch * P)) * 1440
    mm = mm.sort_values("epoch")
    cols = ["sector", "epoch", "tc", "err_cal", "oc_ours_min", "tc_qlp", "err_qlp", "oc_qlp_min", "diff_min", "diff_sig", "depth", "depth_qlp", "snr", "snr_qlp"]
    print(mm[cols].round(4).to_string())
    mm.to_csv(f"{ROOT}/out/timing_table.csv", index=False, float_format="%.6f")
    both = mm.dropna(subset=["tc", "tc_qlp"])
    out["qlp_comparison"] = dict(n_common=len(both), mean_diff_min=float(both.diff_min.mean()), rms_diff_min=float(both.diff_min.std()),
                                 chi2_diff=float(np.sum(both.diff_sig ** 2)), max_abs_sig=float(np.max(np.abs(both.diff_sig))))
    # QLP-only linear fit on its own times, errors as published, and with our calibrated error scale
    cq, covq, rq, chiq = linfit(q.epoch.values, q.tc_qlp.values, q.err_qlp.values)
    out["qlp_own_fit"] = dict(chi2=chiq, dof=len(q) - 2, oc_rms_min=float(np.std(rq) * 1440), median_err_min=float(np.median(q.err_qlp) * 1440))
    print(f"QLP times alone: chi2 {chiq:.1f}/{len(q)-2}, O-C rms {np.std(rq)*1440:.1f} min, median quoted error {np.median(q.err_qlp)*1440:.1f} min; "
          f"common events: mean diff (ours-QLP) {both.diff_min.mean():+.1f} min, rms {both.diff_min.std():.1f} min, chi2 {np.sum(both.diff_sig**2):.1f}/{len(both)}")
    save_json(out, f"{ROOT}/out/timing.json")
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.errorbar(mm.epoch, mm.oc_ours_min, mm.err_cal * 1440, fmt="o", label="deblended M3V (this work)")
    ax.errorbar(mm.epoch + 0.8, mm.oc_qlp_min, mm.err_qlp * 1440, fmt="s", label="QLP (corrected times)", alpha=0.7)
    ax.axhline(0, color="k", lw=0.5); ax.set_xlabel("epoch"); ax.set_ylabel("O - C (min)"); ax.legend()
    ax.set_title(f"O-C vs global ephemeris P={P:.6f}, T0={T0:.4f}")
    fig.tight_layout(); fig.savefig(f"{ROOT}/figs/oc.png", dpi=90)


if __name__ == "__main__":
    main()
