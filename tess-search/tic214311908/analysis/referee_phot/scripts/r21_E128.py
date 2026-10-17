"""Referee r21: the investigator's timing outlier S66 E-128 (their O-C -23 min; referee r07 +1.7 min; QLP 19 min
later than the investigator). On identical cadences, compare the investigator's deliverable (lc_M3V_all.csv, model B
flux_norm and model A) with the referee local-window M3V flux (r05 'M' fit) in 20-min bins, list background, quality
and pointing around the event, and refit the mid-time on each light curve with the same code (fixed referee shape,
free tc and depth, local constant + slope baseline within +-0.3 d). Also: S104 E60/E63, S103 E59, which carry most of
the investigator's odd/even difference. Output: out/E128.json"""
import sys; sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
from rlib import *
from r05_local import Sector, trap_shape
G = json.load(open(f"{R}/out/transit_fit.json"))["global"]; T0, P, T14, FL = G["T0"], G["P"], G["T14_h"] / 24, G["flat"]
M = pd.read_csv(f"{R}/../phot/lc_M3V_all.csv")


def fit_tc(t, f, e, ex, tc0, slope=True):
    grid = np.arange(-0.05, 0.0501, 0.0005); best = None; prof = []
    for dt in grid:
        y = trap_shape(t, tc0 + dt, T14, FL, ex)
        if np.sum(y > 0.5) < 2: prof.append(np.inf); continue
        cols = [np.ones(len(t)), -y] + ([t - tc0] if slope else [])
        A = np.column_stack(cols) / e[:, None]; c, *_ = np.linalg.lstsq(A, f / e, rcond=None)
        c2 = float(np.sum((A @ c - f / e) ** 2)); prof.append(c2)
        if best is None or c2 < best[0]: best = (c2, dt, c[1])
    prof = np.array(prof); ok = prof <= prof.min() + 1
    return dict(dt_min=best[1] * 1440, depth=best[2], lo_min=grid[ok].min() * 1440, hi_min=grid[ok].max() * 1440, chi2=best[0], n=len(t))


out = {}
for s, ep in ((66, -128), (66, -130), (104, 60), (104, 63), (103, 59), (103, 56)):
    tc0 = T0 + ep * P
    S = Sector(s); r = S.window(tc0)
    d = M[(M.sector == s) & (np.abs(M.time_btjd_tdb - tc0) < 0.35)]
    okB = (d.quality_ok == 1) & (d.flare_flag == 0)
    ex = S.exptime
    res = {}
    if r is not None:
        sg = rstd(r["aM"][r["oot"]])
        res["referee_aM"] = fit_tc(r["t"], -r["aM"] / 979.0, np.full(len(r["t"]), sg / 979.0), ex, tc0)
        res["referee_n"] = len(r["t"])
    dB = d[okB]
    for col, lab in (("flux_norm", "inv_B"), ("flux_norm_modelA", "inv_A")):
        x = dB[np.isfinite(dB[col])]
        f = 1 - x[col].values; e = np.full(len(x), rstd(f[np.abs(x.time_btjd_tdb.values - tc0) > 0.07]))
        res[lab] = fit_tc(x.time_btjd_tdb.values, -f * -1 if False else f, e, ex, tc0)
    res["inv_n"] = int(okB.sum()); res["inv_n_all"] = len(d)
    # cadences in the referee window but not in the investigator's strict set and vice versa
    if r is not None:
        tr = np.round(r["t"], 5); ti = np.round(dB.time_btjd_tdb.values, 5)
        res["only_referee"] = int(np.sum(~np.isin(tr, ti))); res["only_investigator"] = int(np.sum(~np.isin(ti, tr)))
        # in-transit only
        intr = np.abs(r["t"] - tc0) < T14 / 2; inti = np.abs(dB.time_btjd_tdb.values - tc0) < T14 / 2
        res["n_in_transit_referee"] = int(intr.sum()); res["n_in_transit_investigator"] = int(inti.sum())
    # 20-min binned comparison table
    edges = np.arange(-0.3, 0.3001, 1 / 72)
    tab = []
    for a, b in zip(edges[:-1], edges[1:]):
        mi = (dB.time_btjd_tdb - tc0 >= a) & (dB.time_btjd_tdb - tc0 < b)
        row = dict(t_min=round((a + b) / 2 * 1440), inv_B_ppt=float((dB.flux_norm[mi].mean() - 1) * 1e3) if mi.sum() else np.nan,
                   inv_A_ppt=float((dB.flux_norm_modelA[mi].mean() - 1) * 1e3) if mi.sum() else np.nan,
                   bkg=float(d.bkg_e_s_px[(d.time_btjd_tdb - tc0 >= a) & (d.time_btjd_tdb - tc0 < b)].mean()) if mi.sum() else np.nan,
                   n_inv=int(mi.sum()))
        if r is not None:
            mr = (r["t"] - tc0 >= a) & (r["t"] - tc0 < b)
            row["ref_ppt"] = float(r["aM"][mr].mean() / 979.0 * 1e3) if mr.sum() else np.nan; row["n_ref"] = int(mr.sum())
        tab.append(row)
    res["binned"] = tab
    out[f"S{s}_E{ep}"] = res
    print(f"S{s} E{ep}: " + "; ".join(f"{k} dt {v['dt_min']:+.1f} min [{v['lo_min']:+.1f},{v['hi_min']:+.1f}] depth {v['depth']*1e3:.1f} ppt" for k, v in res.items() if isinstance(v, dict) and "dt_min" in v)
          + f"; n inv {res['inv_n']} ref {res.get('referee_n')}, only-ref {res.get('only_referee')}, only-inv {res.get('only_investigator')}, in-transit inv/ref {res.get('n_in_transit_investigator')}/{res.get('n_in_transit_referee')}", flush=True)
    if (s, ep) == (66, -128):
        print(pd.DataFrame(tab).round(1).to_string())
save_json(out, f"{R}/out/E128.json")
