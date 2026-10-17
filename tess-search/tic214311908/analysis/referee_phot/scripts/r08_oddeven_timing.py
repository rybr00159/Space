"""Referee r08: odd/even, timing and secondary eclipse on the referee local-window photometry.
Odd/even: (a) weighted per-transit depths (r07, calibrated errors; strict selection, coverage > 0.5, excluding the
non-detection S12 E-388 only in a variant), with leave-one-out jackknife; (b) folded fits of odd and even windows
at the global shape (r06) with errors from the fake-window scatter.
Timing: linear ephemeris to per-transit times with calibrated errors; chi2; O-C; also with formal errors; robust
check using the chi2-profile half widths.
Secondary: windows centred at phase 0.5 processed by r05 (strict), depth at the global shape, error from fake windows.
Phase scan: depth at 400 trial phases (windows centred at T0 + (k + phase) P), S/N distribution.
Output: out/oddeven_timing.json, out/phase_scan.csv"""
import sys; sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
from rlib import *
from r05_local import Sector, trap_shape
from scipy.stats import chi2 as C2
G = json.load(open(f"{R}/out/transit_fit.json"))["global"]
T0, P, T14, FL = G["T0"], G["P"], G["T14_h"] / 24, G["flat"]
PT = pd.read_csv(f"{R}/out/per_transit.csv")


def wmean(v, e):
    w = 1 / e ** 2; return float(np.sum(w * v) / w.sum()), float(1 / np.sqrt(w.sum()))


def oe_stats(d, ecol):
    od = d[d.epoch % 2 != 0]; ev = d[d.epoch % 2 == 0]
    mo, eo = wmean(od.depth.values, od[ecol].values); me, ee = wmean(ev.depth.values, ev[ecol].values)
    return mo, eo, me, ee, mo - me, float(np.hypot(eo, ee))


def depth_fixed(S, tcs, half_shift=0.0, bmax=60.0):
    """depth (e/s) of each window centred at tcs, fixed shape centred at window centre; returns list of (depth, err_white)"""
    out = []
    for tc in tcs:
        r = S.window(tc, bmax=bmax)
        if r is None: continue
        y = trap_shape(r["t"], tc, T14, FL, S.exptime)
        if np.sum(y > 0.5) < 2: continue
        sg = rstd(r["aM"][r["oot"]]); A = np.column_stack([np.ones(len(y)), -y])
        c, *_ = np.linalg.lstsq(A / sg, r["aM"] / sg, rcond=None)
        cov = np.linalg.pinv((A / sg).T @ (A / sg))
        out.append((tc, c[1], np.sqrt(cov[1, 1])))
    return out


def main():
    out = {}
    base = PT[(PT.sel == "strict") & (PT["cov"] > 0.5)].copy()
    for lab, d in (("all_strict", base), ("strict_no_S12E-388", base[~((base.sector == 12) & (base.epoch == -388))]),
                   ("strict_no_S12", base[base.sector != 12]),
                   ("relaxed_no_S12E-388", PT[(PT.sel == "relaxed") & (PT["cov"] > 0.5) & ~((PT.sector == 12) & (PT.epoch == -388))])):
        for ecol in ("depth_err", "depth_err_cal"):
            mo, eo, me, ee, df, de = oe_stats(d, ecol)
            jk = []
            for i in d.index:
                x = oe_stats(d.drop(index=i), ecol); jk.append((int(d.epoch[i]), x[4] / x[5]))
            jk = np.array(jk)
            mu, er = wmean(d.depth.values, d[ecol].values); c2 = float(np.sum(((d.depth - mu) / d[ecol]) ** 2))
            out[f"oddeven_{lab}_{ecol}"] = dict(odd=(mo, eo), even=(me, ee), diff=df, diff_err=de, nsig=df / de, n_odd=int((d.epoch % 2 != 0).sum()), n_even=int((d.epoch % 2 == 0).sum()),
                                                jack_range=(float(jk[:, 1].min()), float(jk[:, 1].max())), depth_chi2=c2, dof=len(d) - 1, p_depth=float(C2.sf(c2, len(d) - 1)))
            print(f"odd/even [{lab}, {ecol}]: odd {mo:.2f}+-{eo:.2f}, even {me:.2f}+-{ee:.2f} e/s, diff {df:+.2f}+-{de:.2f} ({df/de:+.2f} sigma); jackknife {jk[:,1].min():.2f}..{jk[:,1].max():.2f}; "
                  f"all-depth chi2 {c2:.1f}/{len(d)-1} (p {C2.sf(c2, len(d)-1):.3f})")
    # ---- timing
    tm = base[~((base.sector == 12) & (base.epoch == -388))].copy()
    for ecol in ("tc_err", "tc_err_cal", "tc_prof_half"):
        e = np.maximum(tm[ecol].values, 1e-4); A = np.column_stack([np.ones(len(tm)), tm.epoch.values]) / e[:, None]
        c, *_ = np.linalg.lstsq(A, tm.tc.values / e, rcond=None); cov = np.linalg.pinv(A.T @ A)
        r = tm.tc.values - (c[0] + c[1] * tm.epoch.values); c2 = float(np.sum((r / e) ** 2)); dof = len(tm) - 2
        worst = sorted(zip(tm.sector, tm.epoch, r * 1440, r / e), key=lambda z: -abs(z[3]))[:3]
        out[f"timing_{ecol}"] = dict(T0=c[0], P=c[1], T0_err=float(np.sqrt(cov[0, 0])), P_err=float(np.sqrt(cov[1, 1])), chi2=c2, dof=dof, p=float(C2.sf(c2, dof)),
                                     oc_rms_min=float(np.std(r) * 1440), oc_mad_min=float(rstd(r) * 1440), worst=[(int(a), int(b), float(x), float(y)) for a, b, x, y in worst])
        print(f"timing [{ecol}]: n {len(tm)} T0 {c[0]:.5f}+-{np.sqrt(cov[0,0]):.5f} P {c[1]:.7f}+-{np.sqrt(cov[1,1]):.7f} chi2 {c2:.1f}/{dof} (p {C2.sf(c2,dof):.3f}); O-C rms {np.std(r)*1440:.1f} min "
              f"(robust {rstd(r)*1440:.1f}); worst {[(a, b, round(x, 1), round(y, 1)) for a, b, x, y in worst]}")
    # ---- secondary eclipse and phase scan
    scan_rows = []
    sec_rows = []
    phases = np.round(np.arange(0.06, 0.9401, 0.0025), 4)
    for s in SECTORS:
        S = Sector(s)
        t = S.ttdb[S.qok]; k = np.arange(np.floor((t.min() - T0) / P) - 1, np.ceil((t.max() - T0) / P) + 1)
        for ph in phases:
            res = depth_fixed(S, T0 + (k + ph) * P)
            for tc, dd, de in res: scan_rows.append(dict(sector=s, phase=ph, tc=tc, depth=dd, err=de))
    SC = pd.DataFrame(scan_rows); SC.to_csv(f"{R}/out/phase_scan_windows.csv", index=False)
    # empirical single-window error per sector: robust std of single-window depths over the scan (excluding |phase-0.5|<0.02)
    sig_s = {s: rstd(SC[(SC.sector == s)].depth) for s in SECTORS}
    SC["err_emp"] = SC.sector.map(sig_s)
    agg = SC.groupby("phase").apply(lambda d: pd.Series(dict(depth=np.sum(d.depth / d.err_emp ** 2) / np.sum(1 / d.err_emp ** 2), err=1 / np.sqrt(np.sum(1 / d.err_emp ** 2)), n=len(d)))).reset_index()
    agg["snr"] = agg.depth / agg.err
    agg.to_csv(f"{R}/out/phase_scan.csv", index=False)
    p5 = agg.iloc[np.argmin(np.abs(agg.phase - 0.5))]
    i = int(np.argmax(agg.snr))
    out["secondary"] = dict(phase05_depth_e_s=float(p5.depth), phase05_err_e_s=float(p5.err), scan_snr_rstd=float(rstd(agg.snr)), scan_snr_std=float(np.std(agg.snr)),
                            scan_max_snr=float(agg.snr[i]), scan_max_phase=float(agg.phase[i]), scan_min_snr=float(agg.snr.min()),
                            single_window_err_e_s=sig_s)
    print(f"phase 0.5: {p5.depth:.2f}+-{p5.err:.2f} e/s (S/N {p5.snr:.2f}); scan S/N std {np.std(agg.snr):.2f} (robust {rstd(agg.snr):.2f}), max {agg.snr[i]:.2f} at {agg.phase[i]:.4f}, min {agg.snr.min():.2f}; single-window err per sector {np.round(list(sig_s.values()),2).tolist()}")
    # folded odd/even with the same empirical single-window errors (real windows)
    D = []
    for s in SECTORS:
        S = Sector(s)
        tcs = PT[(PT.sector == s) & (PT.sel == "strict") & (PT["cov"] > 0.5)]
        for r in tcs.itertuples():
            if s == 12 and r.epoch == -388: continue
            res = depth_fixed(S, [T0 + r.epoch * P])
            for tc, dd, de in res: D.append(dict(sector=s, epoch=r.epoch, depth=dd, err=sig_s[s]))
    D = pd.DataFrame(D)
    mo, eo = wmean(D[D.epoch % 2 != 0].depth.values, D[D.epoch % 2 != 0].err.values); me, ee = wmean(D[D.epoch % 2 == 0].depth.values, D[D.epoch % 2 == 0].err.values)
    out["oddeven_fixed_ephemeris_emp_err"] = dict(odd=(mo, eo), even=(me, ee), diff=mo - me, diff_err=float(np.hypot(eo, ee)), nsig=(mo - me) / np.hypot(eo, ee),
                                                  per_window=D.round(3).values.tolist())
    print(f"odd/even at fixed ephemeris with empirical single-window errors: odd {mo:.2f}+-{eo:.2f} even {me:.2f}+-{ee:.2f} diff {(mo-me):+.2f}+-{np.hypot(eo,ee):.2f} ({(mo-me)/np.hypot(eo,ee):+.2f} sigma)")
    save_json(out, f"{R}/out/oddeven_timing.json")


if __name__ == "__main__":
    main()
