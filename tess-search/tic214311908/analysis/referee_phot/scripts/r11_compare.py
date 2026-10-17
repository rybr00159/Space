"""Referee r11: event-by-event comparison of the referee (r07) and investigator (phot/out/per_transit_M3V.csv, model B
strict) per-transit times and depths, and QLP event times (phot/out/timing_table.csv), to measure the method-induced
scatter on identical pixel data. Output: out/compare_events.csv, out/compare_events.json"""
import sys; sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
from rlib import *
me = pd.read_csv(f"{R}/out/per_transit.csv"); me = me[me.sel == "strict"]
th = pd.read_csv(f"{R}/../phot/out/per_transit_M3V.csv"); th = th[(th.model == "B") & (~th.relaxed)]
tt = pd.read_csv(f"{R}/../phot/out/timing_table.csv")
m = me.merge(th[["sector", "epoch", "tc", "tc_err", "tc_err_pb", "depth", "depth_err"]], on=["sector", "epoch"], suffixes=("_ref", "_inv"))
m = m.merge(tt[["epoch", "tc_qlp", "err_qlp", "err_cal"]], on="epoch", how="left")
m["dt_ref_inv_min"] = (m.tc_ref - m.tc_inv) * 1440
m["dt_qlp_ref_min"] = (m.tc_qlp - m.tc_ref) * 1440
m["dt_qlp_inv_min"] = (m.tc_qlp - m.tc_inv) * 1440
m["depth_inv_e_s"] = m.depth_inv * 935.0
m["ddepth_e_s"] = m.depth_ref - m.depth_inv_e_s
m.to_csv(f"{R}/out/compare_events.csv", index=False, float_format="%.5f")
k = m[m.sector != 12]; q = m.dropna(subset=["tc_qlp"])
out = dict(n=len(k), dt_ref_inv_rms_min=float(k.dt_ref_inv_min.std()), dt_ref_inv_mean_min=float(k.dt_ref_inv_min.mean()), dt_ref_inv_robust_min=float(rstd(k.dt_ref_inv_min)),
           ddepth_rms_e_s=float(k.ddepth_e_s.std()), ddepth_mean_e_s=float(k.ddepth_e_s.mean()),
           ref_formal_tc_err_med_min=float(np.median(k.tc_err_ref) * 1440), ref_cal_tc_err_med_min=float(np.median(k.tc_err_cal) * 1440),
           inv_formal_tc_err_med_min=float(np.median(np.maximum(k.tc_err_inv, k.tc_err_pb)) * 1440), inv_cal_tc_err_med_min=float(np.nanmedian(k.err_cal) * 1440),
           n_qlp=len(q), dt_qlp_ref_rms_min=float(q.dt_qlp_ref_min.std()), dt_qlp_ref_mean_min=float(q.dt_qlp_ref_min.mean()),
           dt_qlp_inv_rms_min=float(q.dt_qlp_inv_min.std()), dt_qlp_inv_mean_min=float(q.dt_qlp_inv_min.mean()),
           ref_depth_err_med=float(np.median(k.depth_err_cal)), inv_depth_err_med_e_s=float(np.median(k.depth_err_inv) * 935))
save_json(out, f"{R}/out/compare_events.json")
print(json.dumps(out, indent=0))
