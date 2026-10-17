"""Referee r20: where does the investigator's odd/even difference come from?
(a) Recompute their quoted per-transit odd-even (model B, strict, from phot/out/per_transit_M3V.csv) with weighted
    means, to check the arithmetic.
(b) Measure per-transit depths on their deliverable lc_M3V_all.csv (flux_norm, quality_ok & no flare) with the
    referee code: fixed shape at their global ephemeris (T0 3840.91414, P 5.6915087, T14 0.0764 d, flat 0.47), local
    constant baseline within +-0.3 d (+-0.6 d for S12), error = OOT robust scatter / sqrt(n_in-weight); and an
    empirical error from 100 fake epochs per sector on the same LC. Odd-even with both errors, with and without S12.
(c) Same with the referee's local-window LC values (r08 per_window, empirical errors) on the SAME event list, to
    separate the effect of the event selection from that of the photometry.
Output: out/inv_oddeven.json"""
import sys; sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
from rlib import *
from r05_local import trap_shape
T0i, Pi, T14i, FLi = 3840.91414, 5.6915087, 0.0764, 0.47
out = {}


def wm(v, e):
    w = 1 / np.asarray(e) ** 2; return float(np.sum(w * v) / w.sum()), float(1 / np.sqrt(w.sum()))


# (a)
th = pd.read_csv(f"{R}/../phot/out/per_transit_M3V.csv"); tb = th[(th.model == "B") & (~th.relaxed)].copy()
tb["e"] = np.maximum(tb.depth_err, tb.depth_err_pb)
od = tb[tb.epoch % 2 != 0]; ev = tb[tb.epoch % 2 == 0]
mo, eo = wm(od.depth, od.e); me, ee = wm(ev.depth, ev.e)
out["a_quoted_table"] = dict(odd=(mo * 1e3, eo * 1e3), even=(me * 1e3, ee * 1e3), diff_ppt=(mo - me) * 1e3, err_ppt=float(np.hypot(eo, ee) * 1e3),
                            nsig=float((mo - me) / np.hypot(eo, ee)), events=tb[["sector", "epoch", "depth", "e"]].round(5).values.tolist())
print(f"(a) investigator table, model B strict: odd {mo*1e3:.2f}+-{eo*1e3:.2f}, even {me*1e3:.2f}+-{ee*1e3:.2f} ppt, diff {(mo-me)*1e3:+.2f}+-{np.hypot(eo,ee)*1e3:.2f} ({(mo-me)/np.hypot(eo,ee):.2f} sigma)")
# (b)
M = pd.read_csv(f"{R}/../phot/lc_M3V_all.csv")
ok = (M.quality_ok == 1) & (M.flare_flag == 0) & np.isfinite(M.flux_norm)
rows = []; rng = np.random.default_rng(11); emp = {}
for s in SECTORS:
    d = M[ok & (M.sector == s)]; t = d.time_btjd_tdb.values; f = d.flux_norm.values; ex = d.cadence_s.values[0] / 86400
    half = 0.6 if s == 12 else 0.3

    def one(tc):
        w = np.abs(t - tc) < half
        y = trap_shape(t[w], tc, T14i, FLi, ex)
        if np.sum(y > 0.5) < 2 or np.sum((y == 0) & (t[w] < tc)) < 5 or np.sum((y == 0) & (t[w] > tc)) < 5: return None
        A = np.column_stack([np.ones(w.sum()), -y]); c, *_ = np.linalg.lstsq(A, f[w], rcond=None)
        sg = rstd(f[w][y == 0]); cv = np.linalg.pinv(A.T @ A) * sg ** 2
        return c[1], np.sqrt(cv[1, 1])
    k = np.arange(np.ceil((t.min() - T0i) / Pi), np.floor((t.max() - T0i) / Pi) + 1).astype(int)
    for kk in k:
        r = one(T0i + kk * Pi)
        if r is not None: rows.append(dict(sector=s, epoch=int(kk), depth=r[0], err_white=r[1]))
    nul = []
    for i in range(400):
        tc = rng.uniform(t.min() + half, t.max() - half)
        if abs(((tc - T0i) / Pi + 0.5) % 1 - 0.5) * Pi < half + 0.05: continue
        r = one(tc)
        if r is not None: nul.append(r[0])
        if len(nul) >= 100: break
    emp[s] = float(np.std(nul)); emp[f"{s}_robust"] = float(rstd(nul)); emp[f"{s}_mean"] = float(np.mean(nul))
D = pd.DataFrame(rows); D["err_emp"] = D.sector.map(lambda s: emp[s]); D["err"] = np.maximum(D.err_white, D.err_emp)
out["b_fake_single_err"] = emp
for lab, d in (("all", D), ("no_S12", D[D.sector != 12]), ("table_events", D.merge(tb[["sector", "epoch"]], on=["sector", "epoch"]))):
    for ecol in ("err_white", "err"):
        od = d[d.epoch % 2 != 0]; ev = d[d.epoch % 2 == 0]
        mo, eo = wm(od.depth, od[ecol]); me, ee = wm(ev.depth, ev[ecol])
        out[f"b_{lab}_{ecol}"] = dict(odd=(mo * 1e3, eo * 1e3), even=(me * 1e3, ee * 1e3), diff_ppt=(mo - me) * 1e3, err_ppt=float(np.hypot(eo, ee) * 1e3), nsig=float((mo - me) / np.hypot(eo, ee)), n=len(d))
        print(f"(b) investigator LC, {lab:12s} [{ecol:9s}]: n {len(d)} odd {mo*1e3:.2f}+-{eo*1e3:.2f} even {me*1e3:.2f}+-{ee*1e3:.2f} ppt, diff {(mo-me)*1e3:+.2f}+-{np.hypot(eo,ee)*1e3:.2f} ({(mo-me)/np.hypot(eo,ee):+.2f} sigma)")
out["b_events"] = D.round(5).values.tolist()
print("empirical single-transit errors (ppt):", {k: round(v * 1e3, 2) for k, v in emp.items() if isinstance(k, int)})
# (c) referee LC on the investigator's strict event list (r08 per-window depths at fixed ephemeris, empirical errors)
OE = json.load(open(f"{R}/out/oddeven_timing.json"))["oddeven_fixed_ephemeris_emp_err"]["per_window"]
Rw = pd.DataFrame(OE, columns=["sector", "epoch", "depth", "err"]); Rw[["sector", "epoch"]] = Rw[["sector", "epoch"]].astype(int)
for lab, d in (("referee_all", Rw), ("referee_on_table_events", Rw.merge(tb[["sector", "epoch"]], on=["sector", "epoch"]))):
    od = d[d.epoch % 2 != 0]; ev = d[d.epoch % 2 == 0]
    mo, eo = wm(od.depth, od.err); me, ee = wm(ev.depth, ev.err)
    out[f"c_{lab}"] = dict(odd=(mo, eo), even=(me, ee), diff_e_s=mo - me, err=float(np.hypot(eo, ee)), nsig=float((mo - me) / np.hypot(eo, ee)), n=len(d))
    print(f"(c) referee LC, {lab:24s}: n {len(d)} odd {mo:.2f}+-{eo:.2f} even {me:.2f}+-{ee:.2f} e/s, diff {mo-me:+.2f}+-{np.hypot(eo,ee):.2f} ({(mo-me)/np.hypot(eo,ee):+.2f} sigma)")
# event-by-event: investigator LC depth (x 979 e/s scale-free: use x935 to compare with their e/s) vs referee
J = D.merge(Rw, on=["sector", "epoch"], suffixes=("_inv", "_ref"))
J["d_inv_e_s"] = J.depth_inv * 935.0; J["diff"] = J.d_inv_e_s - J.depth_ref
out["event_compare"] = J.round(4).values.tolist(); out["event_compare_cols"] = list(J.columns)
out["event_diff_rms_e_s"] = float(J["diff"].std()); out["event_diff_mean_e_s"] = float(J["diff"].mean())
out["event_corr"] = float(np.corrcoef(J.d_inv_e_s, J.depth_ref)[0, 1])
print(J[["sector", "epoch", "d_inv_e_s", "depth_ref", "diff"]].round(1).to_string())
print(f"event depth inv(LC x 935) - referee: mean {J['diff'].mean():+.2f} rms {J['diff'].std():.2f} e/s; correlation {out['event_corr']:.2f}")
save_json(out, f"{R}/out/inv_oddeven.json")
