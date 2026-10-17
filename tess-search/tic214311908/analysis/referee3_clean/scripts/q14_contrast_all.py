"""q14: late-minus-early depth contrast (S103+S104 minus S39+S66) and sector chi2 for every photometry, with this
referee's per-transit depths (fixed shape, linear baseline, fake-epoch per-sector errors):
INV delivered (q04 inv_l), REF delivered (fixed-tc fit of lc_M3V_REF_clean.csv; errors from win_ref_fake), and the raw
apertures (q04 box3/m25). Output: out/q14_contrast_all.json"""
import json, sys, numpy as np, pandas as pd
sys.path.insert(0, "tic214311908/analysis/referee3_clean/scripts")
from q02_bigap import shape
R3 = "tic214311908/analysis/referee3_clean"
CL = "tic214311908/analysis/clean"
T0A, PA = 3840.91488, 5.6915020
r = pd.read_csv(f"{R3}/out/q04_paired_per_transit.csv"); oe = json.load(open(f"{R3}/out/q10_oddeven.json"))["per_sector_err_ppt"]
def lin(t, f, tc, ex):
    s = shape(t, tc, ex); M = np.column_stack([-s, np.ones_like(t), t - tc]); c, *_ = np.linalg.lstsq(M, f, rcond=None); return c[0] * 1e3
ref = pd.read_csv(f"{CL}/lc_M3V_REF_clean.csv")
rr = pd.DataFrame([dict(sector=s, epoch=int(e), d=lin(g.time_btjd_tdb.values, g.flux_norm.values, T0A + e * PA, g.cadence_s.values[0] / 86400)) for (s, e), g in ref.groupby(["sector", "epoch"])])
sets = {"INV_delivered": r[["sector", "epoch"]].assign(d=r.inv_l, e=[oe["INV"][str(s)] for s in r.sector]),
        "REF_delivered": rr.assign(e=[oe["REF"][str(s)] for s in rr.sector]),
        "REF_delivered_noE58": rr[~((rr.sector == 103) & (rr.epoch == 58))].assign(e=lambda x: [oe["REF"][str(s)] for s in x.sector]),
        "AP_box3": r[["sector", "epoch"]].assign(d=r.L_box3_pt_ppt, e=[oe["AP"]["L_box3_pt"][str(s)] for s in r.sector]),
        "AP_m25": r[["sector", "epoch"]].assign(d=r.L_m25_pt_ppt, e=[oe["AP"]["L_m25_pt"][str(s)] for s in r.sector])}
out = {}
for k, d in sets.items():
    d = d[np.isfinite(d.d)]
    def grp(ss):
        g = d[d.sector.isin(ss)]; w = 1 / g.e ** 2; return float((g.d * w).sum() / w.sum()), float(1 / np.sqrt(w.sum()))
    sec = {int(s): grp([s]) for s in sorted(d.sector.unique())}
    a, ae = grp([39, 66]); b, be = grp([103, 104]); m, me = grp(list(sec))
    chi2 = float(sum(((v[0] - m) / v[1]) ** 2 for v in sec.values()))
    out[k] = dict(per_sector=sec, global_=[m, me], contrast=[b - a, float(np.hypot(ae, be)), (b - a) / float(np.hypot(ae, be))], chi2=[chi2, len(sec) - 1])
    print(f"{k:22s} global {m:5.2f}+-{me:4.2f}  contrast {b-a:+5.2f}+-{np.hypot(ae,be):4.2f} ({(b-a)/np.hypot(ae,be):+.2f} sig)  chi2 {chi2:4.1f}/{len(sec)-1}  " +
          " ".join(f"{s}:{v[0]:5.1f}+-{v[1]:3.1f}" for s, v in sec.items()))
json.dump(out, open(f"{R3}/out/q14_contrast_all.json", "w"), indent=1)
