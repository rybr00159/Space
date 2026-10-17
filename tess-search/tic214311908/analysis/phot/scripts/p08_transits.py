"""Transit measurements on the deblended light curves.

Inputs: lc_M3V_all.csv (model B, primary, plus model A columns), lc_pair_all.csv, lc_Mira_all.csv, data/lc_s*.npz.
1. Global trapezoid fit to the M3V (model B): depth, T14, flat fraction, P, T0 (least squares, exposure-integrated).
2. Per-sector depth at the global shape/ephemeris for every light curve (M3V B/A, pair Bp/A, Mira B/A); errors:
   white (formal) and empirical 'phase-scramble' (std of depths measured at 150 random phases of the same sector,
   real transit windows excluded), which includes red noise; adopted = max of the two.
3. Per-transit fits (free mid-time, depth, linear baseline in +-0.3 d; fixed T14, flat fraction) on the M3V (B and A);
   formal, prayer-bead and chi2-profile errors; events in rejected (scattered-light) data measured on a relaxed
   light curve and flagged.
4. Odd/even depths, depth at phase 0.5 and a phase scan (all phases, step 0.002) with its empirical false-alarm level.
5. Out-of-transit variability: Lomb-Scargle of the undetrended (raw) M3V flux per sector after removing the part
   correlated with the Mira's flux; flare candidates (>= 2 consecutive cadences > 3 sigma).
Outputs: out/transits.json, out/per_transit_M3V.csv, out/per_sector_depths.csv, figs/*.png
"""
import sys; sys.path.insert(0, "tic214311908/analysis/phot/scripts")
from photlib import *
from tlib import trap_frac, ephem_frac, lin_depth, fit_single
from scipy.optimize import least_squares
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

RNG = np.random.default_rng(42)


def load():
    M = pd.read_csv(f"{ROOT}/lc_M3V_all.csv")
    Pq = pd.read_csv(f"{ROOT}/lc_pair_all.csv"); Mi = pd.read_csv(f"{ROOT}/lc_Mira_all.csv")
    st = json.load(open(f"{ROOT}/out/lc_stats.json"))
    return M, Pq, Mi, st


def series(df, col, errcol, flarecol=None):
    ok = (df.quality_ok == 1) & np.isfinite(df[col]) & np.isfinite(df[errcol])
    if flarecol is not None and flarecol in df: ok &= df[flarecol] == 0
    d = df[ok]
    return d.time_btjd_tdb.values, d[col].values, d[errcol].values, d.sector.values, d.cadence_s.values / 86400.0


def global_fit(t, f, e, exp, P0, T00):
    """Parameters relative to the reference ephemeris (dT0, dP) so that scipy's diff_step is an absolute step."""
    def res(p):
        dT0, dP, T14, fl, dep = p
        y = np.zeros(len(t))
        for ex in np.unique(exp):
            m = exp == ex; y[m] = ephem_frac(t[m], P0 + dP, T00 + dT0, T14, fl, ex)
        return (f - (1 - dep * y)) / e
    best = None
    for T14 in (0.06, 0.067, 0.075):
        for fl in (0.3, 0.6, 0.8):
            p0 = [0.0, 0.0, T14, fl, 0.015]
            q = least_squares(res, p0, bounds=([-0.02, -0.0005, 0.03, 0.0, 0.0], [0.02, 0.0005, 0.12, 0.98, 0.1]),
                              x_scale=[0.001, 1e-5, 0.003, 0.05, 0.002], diff_step=[2e-5, 2e-8, 2e-4, 2e-3, 2e-5])
            if best is None or q.cost < best.cost: best = q
    q = best
    J = q.jac; cov = np.linalg.pinv(J.T @ J) * max(1, 2 * q.cost / (len(t) - 5))
    x = q.x.copy(); x[0] += T00; x[1] += P0
    return x, np.sqrt(np.diag(cov)), cov, 2 * q.cost


def frac_all(t, exp, P, T0, T14, fl):
    y = np.zeros(len(t))
    for ex in np.unique(exp):
        m = exp == ex; y[m] = ephem_frac(t[m], P, T0, T14, fl, ex)
    return y


def depth_with_scramble(t, f, e, exp, P, T0, T14, fl, nscr=150, exclude_real=True):
    y = frac_all(t, exp, P, T0, T14, fl)
    d, de, chi, _ = lin_depth(f, e, y)
    real = np.abs((((t - T0_REF) / P_REF + 0.5) % 1 - 0.5) * P_REF) < 0.1 if exclude_real else np.zeros(len(t), bool)
    nulls = []
    for k in range(nscr):
        ph = RNG.uniform(0.08, 0.92)
        y2 = frac_all(t, exp, P, T0 + ph * P, T14, fl)
        m = ~real
        if np.sum(y2[m] > 0.5) < 3: continue
        nulls.append(lin_depth(f[m], e[m], y2[m])[0])
    nulls = np.array(nulls)
    return d, de, (float(np.std(nulls)) if len(nulls) > 10 else np.nan), (float(np.mean(nulls)) if len(nulls) > 10 else np.nan)


def main():
    global P_REF, T0_REF
    M, Pq, Mi, st = load()
    out = {}
    t, f, e, sec, exp = series(M, "flux_norm", "flux_err", "flare_flag")
    # ---- 1. global fit (model B M3V)
    x, xe, cov, chi2 = global_fit(t, f, e, exp, P_EPH, T0_EPH)
    T0g, Pg, T14g, flg, depg = x
    P_REF, T0_REF = Pg, T0g
    out["global_B"] = dict(T0=T0g, T0_err=xe[0], P=Pg, P_err=xe[1], T14_d=T14g, T14_err=xe[2], flat_frac=flg, flat_err=xe[3],
                           depth=depg, depth_err=xe[4], chi2=chi2, n=len(t))
    print(f"global (B): T0 {T0g:.5f}+-{xe[0]:.5f}, P {Pg:.7f}+-{xe[1]:.7f}, T14 {T14g*24:.3f}+-{xe[2]*24:.3f} h, flat {flg:.2f}+-{xe[3]:.2f}, "
          f"depth {depg*1e3:.2f}+-{xe[4]*1e3:.2f} ppt, chi2 {chi2:.0f}/{len(t)}", flush=True)
    # V-shape vs flat: fixed fl = 0 fit
    def res_v(p):
        y = frac_all(t, exp, Pg + p[1], T0g + p[0], p[2], 0.0); return (f - (1 - p[3] * y)) / e
    qv = least_squares(res_v, [0.0, 0.0, T14g, depg], x_scale=[0.001, 1e-5, 0.003, 0.002], diff_step=[2e-5, 2e-8, 2e-4, 2e-5])
    out["global_B"]["chi2_V"] = 2 * qv.cost
    print(f"  V-shaped (flat=0) chi2 {2*qv.cost:.1f} vs {chi2:.1f}: delta {2*qv.cost-chi2:.1f}")
    # model A global (depth only, same shape/ephemeris)
    tA, fA, eA, secA, expA = series(M, "flux_norm_modelA", "flux_err_modelA")
    # ---- 2. per-sector depths for all light curves
    rows = []
    sets = [("M3V_B", M, "flux_norm", "flux_err", "flare_flag"), ("M3V_A", M, "flux_norm_modelA", "flux_err_modelA", None),
            ("pair_Bp", Pq, "flux_norm", "flux_err", None), ("pair_A", Pq, "flux_norm_modelA", "flux_err", None),
            ("Mira_B", Mi, "flux_norm", "flux_err", None), ("Mira_A", Mi, "flux_norm_modelA", "flux_err", None)]
    for key, df, col, ecol, fcol in sets:
        star = key.split("_")[0]
        tt, ff, ee, ss, xx = series(df, col, ecol, fcol)
        if key.endswith("_A") and star != "M3V":
            # model-A errors: scale from the scatter of the model-A light curve itself
            pass
        for s in SECTORS + ("all",):
            m = (ss == s) if s != "all" else np.ones(len(tt), bool)
            if m.sum() < 50: continue
            # use the empirical per-sector scatter for weights
            w_e = np.zeros(m.sum())
            for s2 in np.unique(ss[m]):
                k = ss[m] == s2; w_e[k] = robust_std(ff[m][k]) * np.sqrt(np.clip(ee[m][k] / np.median(ee[m][k]), 1, None))
            d, de, dscr, dmean = depth_with_scramble(tt[m], ff[m], w_e, xx[m], Pg, T0g, T14g, flg)
            fref = np.mean([st[str(s2)]["fref"][star] for s2 in np.unique(ss[m])])
            rows.append(dict(lc=key, sector=s, depth_ppt=d * 1e3, err_white_ppt=de * 1e3, err_scramble_ppt=dscr * 1e3, null_mean_ppt=dmean * 1e3,
                             err_adopted_ppt=max(de, dscr if np.isfinite(dscr) else 0) * 1e3, fref=fref,
                             dip_e_s=d * fref, err_e_s=max(de, dscr if np.isfinite(dscr) else 0) * fref))
    D = pd.DataFrame(rows); D.to_csv(f"{ROOT}/out/per_sector_depths.csv", index=False, float_format="%.4f")
    print(D.pivot(index="sector", columns="lc", values="dip_e_s").round(1).to_string())
    print(D.pivot(index="sector", columns="lc", values="err_e_s").round(1).to_string())
    # chi2 of sector depths around the weighted mean (M3V B, A)
    for key in ("M3V_B", "M3V_A"):
        d = D[(D.lc == key) & (D.sector != "all")]; w = 1 / d.err_adopted_ppt ** 2
        mu = np.sum(w * d.depth_ppt) / w.sum(); c2 = np.sum(w * (d.depth_ppt - mu) ** 2)
        out[f"sector_consistency_{key}"] = dict(mean_ppt=mu, err=1 / np.sqrt(w.sum()), chi2=c2, dof=len(d) - 1)
        print(f"{key}: sector depths weighted mean {mu:.2f}+-{1/np.sqrt(w.sum()):.2f} ppt, chi2 {c2:.1f}/{len(d)-1}")
    # ---- 3. per-transit fits
    ptr = []
    allM = M.copy()
    for key, col, ecol in (("B", "flux_norm", "flux_err"), ("A", "flux_norm_modelA", "flux_err_modelA")):
        for relaxed in (False, True):
            if relaxed:
                # all finite cadences except extreme scattered light (bkg excess < 800 e/s/px), undetrended flux with a
                # local linear baseline (the global trend is not defined in rejected stretches)
                rcol = "flux_raw_norm" if key == "B" else "flux_raw_norm_modelA"
                bmed = allM.groupby("sector").bkg_e_s_px.transform("median")
                d = allM[np.isfinite(allM[rcol]) & ((allM.bkg_e_s_px - bmed) < 800)]
                ff_col = rcol
            else:
                d = allM[(allM.quality_ok == 1) & np.isfinite(allM[col]) & (allM.flare_flag == 0)]
                ff_col = col
            tt, ff, ee, ss = d.time_btjd_tdb.values, d[ff_col].values, d[ecol].values, d.sector.values
            ex = d.cadence_s.values / 86400
            for s in SECTORS:
                m = ss == s
                if m.sum() < 30: continue
                sd = robust_std(np.diff(ff[m])) / np.sqrt(2)
                ks = np.arange(np.ceil((tt[m].min() - T0g) / Pg - 0.02), np.floor((tt[m].max() - T0g) / Pg + 0.02) + 1)
                for k in ks:
                    tc0 = T0g + k * Pg
                    w = m & (np.abs(tt - tc0) < 0.3)
                    n_in = np.sum(m & (np.abs(tt - tc0) < T14g / 2))
                    exp_in = T14g / ex[m][0]
                    if n_in < max(2, 0.5 * exp_in) or w.sum() < 10: continue
                    r = fit_single(tt[w], ff[w], np.full(w.sum(), sd), ex[w][0], tc0, T14g, flg, win=0.3, dt_max=0.08)
                    if r is None: continue
                    r.update(model=key, relaxed=relaxed, sector=int(s), epoch=int(k), tc_pred=tc0, cov_frac=float(n_in / exp_in), snr=r["depth"] / r["depth_err"])
                    ptr.append(r)
    PT = pd.DataFrame(ptr)
    PT["tc_err_adopt"] = PT[["tc_err", "tc_err_pb"]].max(axis=1)
    PT["depth_err_adopt"] = PT[["depth_err", "depth_err_pb"]].max(axis=1)
    PT["oc_min"] = (PT.tc - PT.tc_pred) * 1440
    PT.to_csv(f"{ROOT}/out/per_transit_M3V.csv", index=False, float_format="%.6f")
    strict = PT[(PT.model == "B") & (~PT.relaxed)]
    print(strict[["sector", "epoch", "tc", "oc_min", "tc_err", "tc_err_pb", "depth", "depth_err", "depth_err_pb", "snr", "n_in", "cov_frac", "redchi"]].round(4).to_string())
    # per-transit odd/even (strict data), with jackknife over transits
    for mdl in ("B", "A"):
        d = PT[(PT.model == mdl) & (~PT.relaxed) & (PT.cov_frac > 0.5)].copy()
        d["par"] = d.epoch % 2
        def oe_stat(dd):
            r = {}
            for par in (0, 1):
                x = dd[dd.par == par]; w = 1 / x.depth_err_adopt ** 2
                r[par] = (np.sum(w * x.depth) / w.sum(), 1 / np.sqrt(w.sum()))
            return r[1][0] - r[0][0], np.hypot(r[1][1], r[0][1]), r
        dd_, de_, r_ = oe_stat(d)
        w = 1 / d.depth_err_adopt ** 2; mu = np.sum(w * d.depth) / w.sum(); chi_all = float(np.sum(w * (d.depth - mu) ** 2))
        jk = []
        for i in d.index:
            x_, e_, _ = oe_stat(d.drop(index=i)); jk.append((int(d.epoch[i]), x_ / e_))
        jk = np.array(jk)
        out[f"odd_even_pertransit_{mdl}"] = dict(odd=r_[1], even=r_[0], diff=dd_, diff_err=de_, nsigma=dd_ / de_, n_odd=int((d.par == 1).sum()),
                                                 n_even=int((d.par == 0).sum()), chi2_all_depths=chi_all, dof=len(d) - 1,
                                                 jackknife_nsigma_range=[float(jk[:, 1].min()), float(jk[:, 1].max())],
                                                 jackknife_min_epoch=int(jk[np.argmin(jk[:, 1]), 0]), jackknife_max_epoch=int(jk[np.argmax(jk[:, 1]), 0]))
        print(f"per-transit odd/even ({mdl}): odd {r_[1][0]*1e3:.2f}+-{r_[1][1]*1e3:.2f}, even {r_[0][0]*1e3:.2f}+-{r_[0][1]*1e3:.2f} ppt; "
              f"diff {dd_*1e3:.2f}+-{de_*1e3:.2f} ({dd_/de_:.1f} sigma); all-depth chi2 {chi_all:.1f}/{len(d)-1}; jackknife sigma range "
              f"{jk[:,1].min():.1f}..{jk[:,1].max():.1f}")
    # ---- 4. odd/even, phase 0.5, scan
    y = frac_all(t, exp, Pg, T0g, T14g, flg); ep = np.round((t - T0g) / Pg)
    sdv = np.zeros(len(t))
    for s in np.unique(sec): k = sec == s; sdv[k] = robust_std(f[k])
    oe = {}
    for par, lab in ((0, "even"), (1, "odd")):
        m = (ep % 2 == par) | (y == 0)
        d, de, dscr, _ = depth_with_scramble(t[m], f[m], sdv[m], exp[m], Pg, T0g, T14g, flg, nscr=100)
        oe[lab] = (d * 1e3, max(de, dscr) * 1e3)
    diff = oe["odd"][0] - oe["even"][0]; derr = np.hypot(oe["odd"][1], oe["even"][1])
    out["odd_even_ppt"] = dict(odd=oe["odd"], even=oe["even"], diff=diff, diff_err=derr, nsigma=abs(diff) / derr)
    print(f"odd {oe['odd'][0]:.2f}+-{oe['odd'][1]:.2f}, even {oe['even'][0]:.2f}+-{oe['even'][1]:.2f} ppt: diff {diff:.2f}+-{derr:.2f} ({abs(diff)/derr:.1f} sigma)")
    y5 = frac_all(t, exp, Pg, T0g + 0.5 * Pg, T14g, flg)
    m = y == 0
    d5, de5, _, _ = lin_depth(f[m], sdv[m], y5[m])
    scan = []
    for ph in np.arange(0.03, 0.97, 0.002):
        y2 = frac_all(t, exp, Pg, T0g + ph * Pg, T14g, flg)
        if np.sum(y2[m] > 0.5) < 10: continue
        dd, ee_, _, _ = lin_depth(f[m], sdv[m], y2[m]); scan.append((ph, dd, ee_))
    scan = np.array(scan); snr = scan[:, 1] / scan[:, 2]
    rs = robust_std(snr)
    i = int(np.argmax(snr))
    out["secondary"] = dict(phase05_ppt=d5 * 1e3, phase05_err_white=de5 * 1e3, phase05_err_scan=rs * de5 * 1e3,
                            scan_max_snr=float(snr[i]), scan_max_phase=float(scan[i, 0]), scan_max_depth_ppt=float(scan[i, 1] * 1e3),
                            scan_snr_rstd=float(rs), scan_min_snr=float(snr.min()))
    print(f"phase 0.5: {d5*1e3:.2f}+-{de5*1e3:.2f} (white) / {rs*de5*1e3:.2f} (scan-calibrated) ppt; scan max S/N {snr[i]:.2f} at {scan[i,0]:.3f} "
          f"({scan[i,1]*1e3:.2f} ppt); scan S/N robust std {rs:.2f}; min S/N {snr.min():.2f}")
    np.save(f"{ROOT}/out/phase_scan.npy", scan)
    # ---- 5. variability: raw M3V flux vs Mira; Lomb-Scargle; flares
    from astropy.timeseries import LombScargle
    var = {}
    Mi_ = Mi.set_index(["sector", "time_btjd_tdb"])
    for s in SECTORS:
        d = M[(M.sector == s) & (M.quality_ok == 1)]
        dm = Mi[(Mi.sector == s) & (Mi.quality_ok == 1)]
        if len(d) < 100: continue
        tt = d.time_btjd_tdb.values; raw = d.flux_raw_norm.values
        mira = np.interp(tt, dm.time_btjd_tdb.values, dm.dflux_e_s.values)
        A = np.column_stack([np.ones(len(tt)), mira, tt - tt.mean()]); c, *_ = np.linalg.lstsq(A, raw, rcond=None)
        resid = raw - A[:, :2] @ c[:2]
        good = np.isfinite(resid) & ~((np.abs((((tt - T0g) / Pg + 0.5) % 1 - 0.5) * Pg) < 0.1))
        fr, pw = LombScargle(tt[good], resid[good]).autopower(minimum_frequency=1 / 13, maximum_frequency=10, samples_per_peak=10)
        j = int(np.argmax(pw)); fap = LombScargle(tt[good], resid[good]).false_alarm_probability(pw[j])
        amp = np.std(resid[good]) * np.sqrt(2)
        var[s] = dict(best_period_d=float(1 / fr[j]), power=float(pw[j]), fap=float(fap), coeff_mira_per_eps=float(c[1]),
                      raw_rms_ppt=float(np.std(raw[good]) * 1e3), resid_rms_ppt=float(np.std(resid[good]) * 1e3))
        print(f"S{s} variability: raw rms {np.std(raw[good])*1e3:.1f} ppt; Mira coupling {c[1]*1e3:.3f} ppt/(e/s); LS best P {1/fr[j]:.2f} d power {pw[j]:.3f} FAP {fap:.2g}")
    out["variability"] = var
    # flares: >= 2 consecutive cadences > 3 sigma in the detrended (B) flux
    fl = []
    for s in SECTORS:
        d = M[(M.sector == s) & (M.quality_ok == 1)].reset_index(drop=True)
        if len(d) < 100: continue
        sd = robust_std(d.flux_norm.values); hi = (d.flux_norm.values - 1) > 3 * sd
        run = 0
        for i_, h in enumerate(hi):
            run = run + 1 if h else 0
            if run == 2: fl.append(dict(sector=s, t=float(d.time_btjd_tdb[i_ - 1]), peak_ppt=float((d.flux_norm[i_ - 1: i_ + 3].max() - 1) * 1e3), sigma=float(sd * 1e3)))
    out["flare_candidates"] = fl
    print(f"flare candidates (>=2 consecutive >3 sigma): {len(fl)}", [(f_['sector'], round(f_['t'], 3), round(f_['peak_ppt'], 1)) for f_ in fl][:20])
    save_json(out, f"{ROOT}/out/transits.json")
    # ---- figures
    fig, axs = plt.subplots(2, 4, figsize=(18, 8))
    for ax, (key, df, col, ecol, fcol) in zip(axs.ravel(), sets + [("M3V_B relaxed", M, "flux_norm", "flux_err", None)]):
        tt, ff, ee, ss, xx = series(df, col, ecol, fcol) if "relaxed" not in key else (df.time_btjd_tdb.values, df[col].values, df[ecol].values, df.sector.values, df.cadence_s.values)
        ph = (((tt - T0g) / Pg + 0.5) % 1 - 0.5) * Pg * 24
        m = np.abs(ph) < 6
        b = np.arange(-6, 6.01, 0.25); idx = np.digitize(ph[m], b)
        yb = np.array([np.nanmedian(ff[m][idx == i]) if np.sum(idx == i) > 5 else np.nan for i in range(1, len(b))])
        eb = np.array([1.2533 * robust_std(ff[m][idx == i]) / np.sqrt(np.sum(idx == i)) if np.sum(idx == i) > 5 else np.nan for i in range(1, len(b))])
        ax.errorbar(b[:-1] + 0.125, (yb - 1) * 1e3, eb * 1e3, fmt="o", ms=3)
        xm = np.linspace(-6, 6, 500) / 24
        ax.plot(xm * 24, -depg * 1e3 * trap_frac(xm, 0.0, T14g, flg, 0.0014), "r-", lw=1)
        ax.set_title(key); ax.axhline(0, color="k", lw=0.5); ax.set_xlabel("hours from mid-transit"); ax.set_ylabel("ppt")
    axs.ravel()[-1].plot(scan[:, 0], scan[:, 1] / scan[:, 2], "k-", lw=0.7); axs.ravel()[-1].set_title("phase scan S/N (M3V B)")
    fig.tight_layout(); fig.savefig(f"{ROOT}/figs/fold_all.png", dpi=80); plt.close(fig)


if __name__ == "__main__":
    main()
