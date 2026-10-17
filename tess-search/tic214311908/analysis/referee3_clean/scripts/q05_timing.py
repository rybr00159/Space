"""q05: independent per-transit mid-times, depths, timing-error calibration, TTV chi2 and linear ephemeris from the
analyst's cleaned light curves (lc_M3V_INV_clean.csv, lc_M3V_REF_clean.csv).

Fit per window (+-0.30 d INV; REF strict windows as delivered): fixed batman shape (k 0.117, b 0.52, a/R* 25.3,
u 0.125/0.533, exposure-integrated), free mid-time (grid 0.5 min + parabola, +-58 min), free depth, linear baseline.
Timing-error calibration (LC level): the same fit on fake windows with the transit model INJECTED at a known time
(15.2 ppt): INV fake windows at phases 0.12-0.88 (step 0.04) of the real LC; REF the analyst's 'fake' REF windows
(clean/data/win_ref_fake.npz, the referee-method pipeline run at fake epochs, no transit) with the model added.
Per sector: sigma_t = std(recovered - injected) (and robust std); per-transit error = that sector value.
Ephemeris: WLS on epochs; combined = WLS on the per-event mean of INV and REF (equal-weight when both exist, error
from the per-sector fake-epoch covariance of the two pipelines is not available for REF-vs-INV pairs, so the
combination is done (i) treating them as fully correlated (conservative, error = mean of the two) and (ii) as
independent; both are reported).
Output: out/q05_timing.json, out/q05_per_transit.csv
"""
import json, numpy as np, pandas as pd, sys
sys.path.insert(0, "tic214311908/analysis/referee3_clean/scripts")
from q02_bigap import shape
R3 = "tic214311908/analysis/referee3_clean"
CL = "tic214311908/analysis/clean"
T0A, PA = 3840.91488, 5.6915020          # analyst final (used only to place windows)
DEP = 15.2e-3
rng = np.random.default_rng(31)


def fit_tc(t, f, tc0, ex, dmax=0.04, step=0.5 / 1440):
    tau0 = t - tc0
    grid = np.arange(-dmax, dmax + 1e-12, step)
    chi = np.full(len(grid), np.inf); cs = []
    for j, d in enumerate(grid):
        s = shape(t, tc0 + d, ex)
        if np.sum(s > 0.5) < 2: cs.append(None); continue
        M = np.column_stack([-s, np.ones_like(t), tau0])
        c, *_ = np.linalg.lstsq(M, f, rcond=None)
        r = f - M @ c; chi[j] = r @ r; cs.append(c)
    j = int(np.argmin(chi))
    if not np.isfinite(chi[j]): return np.nan, np.nan
    tc = tc0 + grid[j]
    if 0 < j < len(grid) - 1 and np.all(np.isfinite(chi[j - 1:j + 2])):
        y0, y1, y2 = chi[j - 1:j + 2]; den = y0 - 2 * y1 + y2
        if den > 0: tc += 0.5 * (y0 - y2) / den * step
    return tc, cs[j][0]


def accept(t, tc, half, ex):
    tau = t - tc; oot = np.abs(tau) > 0.07; nexp = 2 * (half - 0.07) / ex
    return np.sum(oot & (tau < 0)) >= 0.25 * nexp / 2 and np.sum(oot & (tau > 0)) >= 0.25 * nexp / 2 and np.sum(np.abs(tau) < 0.028) >= 2


def main():
    inv = pd.read_csv(f"{CL}/lc_M3V_INV_clean.csv", usecols=["time_btjd_tdb", "sector", "cadence_s", "quality_ok", "flux_raw_norm", "flare_flag"])
    inv = inv[(inv.quality_ok == 1) & np.isfinite(inv.flux_raw_norm) & (inv.flare_flag == 0)]
    ref = pd.read_csv(f"{CL}/lc_M3V_REF_clean.csv")
    rows = []; cal = {}
    # ---------------- INV real + injected fakes
    for sec, g in inv.groupby("sector"):
        t = g.time_btjd_tdb.values; y = g.flux_raw_norm.values; ex = g.cadence_s.values[0] / 86400
        e0 = int(np.floor((t.min() - T0A) / PA)); e1 = int(np.ceil((t.max() - T0A) / PA))
        dts = []
        for e in range(e0, e1 + 1):
            for fph in np.r_[0.0, np.arange(0.12, 0.881, 0.04)]:
                tc = T0A + (e + fph) * PA
                m = np.abs(t - tc) < 0.30
                if m.sum() < 8 or not accept(t[m], tc, 0.30, ex): continue
                if fph == 0.0:
                    tcf, d = fit_tc(t[m], y[m], tc, ex)
                    rows.append(dict(pl="inv", sector=sec, epoch=e, tc=tcf, depth=d * 1e3, n=int(m.sum())))
                else:
                    tin = tc + rng.uniform(-0.01, 0.01)
                    yi = y[m] - DEP * shape(t[m], tin, ex) * 1.0
                    tcf, d = fit_tc(t[m], yi, tc, ex)
                    if np.isfinite(tcf): dts.append((tcf - tin) * 1440)
        dts = np.array(dts)
        cal[f"inv_{sec}"] = dict(n=len(dts), std_min=float(np.std(dts)), rstd_min=float(1.4826 * np.median(np.abs(dts - np.median(dts)))),
                                 frac_gt_15min=float(np.mean(np.abs(dts) > 15)))
    # ---------------- REF real
    for (sec, ep), g in ref.groupby(["sector", "epoch"]):
        t = g.time_btjd_tdb.values; y = g.flux_norm.values; ex = g.cadence_s.values[0] / 86400
        tc = T0A + ep * PA
        tcf, d = fit_tc(t, y, tc, ex)
        rows.append(dict(pl="ref", sector=sec, epoch=int(ep), tc=tcf, depth=d * 1e3, n=len(t)))
    # ---------------- REF fakes (analyst's fake windows) with injection
    z = np.load(f"{CL}/data/win_ref_fake.npz")
    for sec in np.unique(z["sector"]):
        dts = []
        for wi in np.where(z["sector"] == sec)[0]:
            m = z["w"] == wi
            t = z["t"][m]; f = z["f"][m]; ex = z["ex"][m][0]; tc = z["tc0"][wi]
            tin = tc + rng.uniform(-0.01, 0.01)
            fi = f - DEP * shape(t, tin, ex)
            if np.sum(np.abs(t - tin) < 0.028) < 2: continue
            tcf, d = fit_tc(t, fi, tc, ex)
            if np.isfinite(tcf): dts.append((tcf - tin) * 1440)
        dts = np.array(dts)
        cal[f"ref_{sec}"] = dict(n=len(dts), std_min=float(np.std(dts)), rstd_min=float(1.4826 * np.median(np.abs(dts - np.median(dts)))),
                                 frac_gt_15min=float(np.mean(np.abs(dts) > 15)))
    df = pd.DataFrame(rows)
    df["err_min"] = [cal[f"{p_}_{s_}"]["std_min"] for p_, s_ in zip(df["pl"], df["sector"])]
    out = dict(calibration=cal)

    def wls(E, T, s):
        A = np.column_stack([np.ones(len(E)), E]) / s[:, None]
        c, *_ = np.linalg.lstsq(A, T / s, rcond=None); cov = np.linalg.inv(A.T @ A)
        r = (T - (c[0] + c[1] * E)) / s
        return c, cov, float(r @ r), r

    def eph(sub, name):
        E = sub.epoch.values.astype(float); T = sub.tc.values; s = sub.err_min.values / 1440
        c, cov, chi2, r = wls(E, T, s)
        dof = len(E) - 2
        # Monte Carlo p of chi2 under Gaussian errors
        pred = []
        for Ef in (100, 130, 160):
            v = cov[0, 0] + Ef ** 2 * cov[1, 1] + 2 * Ef * cov[0, 1]
            pred.append(dict(E=Ef, T=float(c[0] + c[1] * Ef), sig_min=float(np.sqrt(v) * 1440)))
        res = dict(T0_E0=float(c[0]), T0_err=float(np.sqrt(cov[0, 0])), P=float(c[1]), P_err=float(np.sqrt(cov[1, 1])),
                   corr=float(cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])), chi2=chi2, dof=dof, n=len(E), oc_rms_min=float(np.std(r * s) * 1440),
                   max_abs_pull=float(np.max(np.abs(r))), predictions=pred)
        print(f"{name:14s} T0 {c[0]:.5f}+-{np.sqrt(cov[0,0]):.5f} P {c[1]:.7f}+-{np.sqrt(cov[1,1]):.7f} chi2 {chi2:.1f}/{dof} "
              f"O-C rms {res['oc_rms_min']:.1f} min, max pull {res['max_abs_pull']:.1f}; pred E130 +-{pred[1]['sig_min']:.1f} min")
        return res
    df = df[np.isfinite(df.tc)]
    out["inv"] = eph(df[df.pl == "inv"], "INV")
    out["ref"] = eph(df[df.pl == "ref"], "REF")
    out["ref_noE58"] = eph(df[(df.pl == "ref") & ~((df.sector == 103) & (df.epoch == 58))], "REF no E58")
    # combined: per event mean
    pv = df.pivot_table(index="epoch", columns="pl", values=["tc", "err_min"])
    comb = []
    for ep, row in pv.iterrows():
        ti, tr = row[("tc", "inv")], row[("tc", "ref")]
        ei, er = row[("err_min", "inv")], row[("err_min", "ref")]
        if np.isfinite(ti) and np.isfinite(tr):
            w = np.array([1 / ei ** 2, 1 / er ** 2]); w /= w.sum()
            tcm = w[0] * ti + w[1] * tr
            e_ind = 1 / np.sqrt(1 / ei ** 2 + 1 / er ** 2)
            e_cor = w[0] * ei + w[1] * er          # fully correlated
            comb.append(dict(epoch=ep, tc=tcm, err_ind=e_ind, err_cor=e_cor, dt_inv_ref_min=(ti - tr) * 1440))
        elif np.isfinite(ti):
            comb.append(dict(epoch=ep, tc=ti, err_ind=ei, err_cor=ei, dt_inv_ref_min=np.nan))
        else:
            comb.append(dict(epoch=ep, tc=tr, err_ind=er, err_cor=er, dt_inv_ref_min=np.nan))
    cb = pd.DataFrame(comb)
    for k in ("ind", "cor"):
        cb["err_min"] = cb[f"err_{k}"]
        out[f"comb_{k}"] = eph(cb, f"comb ({k})")
    d = cb.dt_inv_ref_min.dropna()
    out["inv_minus_ref_time_min"] = dict(n=int(len(d)), mean=float(d.mean()), rms=float(np.std(d)))
    df.to_csv(f"{R3}/out/q05_per_transit.csv", index=False)
    json.dump(out, open(f"{R3}/out/q05_timing.json", "w"), indent=1)
    print(json.dumps(cal, indent=0))


if __name__ == "__main__":
    main()
