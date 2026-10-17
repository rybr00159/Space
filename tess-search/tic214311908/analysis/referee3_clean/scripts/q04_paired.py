"""q04: paired comparison of the analyst's cleaned INV PRF light curve (lc_M3V_INV_clean.csv, flux_raw_norm,
quality_ok == 1) with the referee3 raw-pixel apertures (q02) at IDENTICAL real and fake transit centres.

Same fit for both (fixed batman shape normalised to 1 at mid-transit, quadratic baseline over +-0.30 d; also a
linear baseline for INV as in the analyst's fits). For each sector:
  * mean real depth (ppt) of INV and of each aperture (EE-corrected), and their difference;
  * error of the difference from the paired fake-epoch differences (same tc in both curves), so the shared noise
    of the two measurements is taken into account;
  * the late-minus-early contrast (S103+S104 minus S39+S66) per method and the difference of the contrasts with
    paired errors; same for a 'fake contrast' null built from fake epochs.
Output: out/q04_paired.json, out/q04_paired_per_transit.csv
"""
import json, numpy as np, pandas as pd, sys
sys.path.insert(0, "tic214311908/analysis/referee3_clean/scripts")
from q02_bigap import shape, accept, HALF
R3 = "tic214311908/analysis/referee3_clean"
CL = "tic214311908/analysis/clean"
tag = sys.argv[1] if len(sys.argv) > 1 else ""
inv = pd.read_csv(f"{CL}/lc_M3V_INV_clean.csv", usecols=["time_btjd_tdb", "sector", "cadence_s", "quality_ok", "flux_raw_norm", "flux_norm"])
inv = inv[(inv.quality_ok == 1) & np.isfinite(inv.flux_raw_norm)]
r = pd.read_csv(f"{R3}/out/q02_bigap_per_transit{tag}.csv"); f = pd.read_csv(f"{R3}/out/q02_bigap_fake{tag}.csv")
S = json.load(open(f"{R3}/out/q02_bigap{tag}.json"))
APS = [c for c in r.columns if c.startswith("L_") and c.endswith("_pt")] + ["L_box3_noreg", "L_m25_noreg"]


def fit(t, y, tc, ex, order):
    tau = t - tc
    s = shape(t, tc, ex)
    M = np.column_stack([-s] + [tau ** o for o in range(order + 1)])
    c, *_ = np.linalg.lstsq(M, y, rcond=None)
    return c[0]


def inv_depths(df):
    out = []
    for _, row in df.iterrows():
        g = inv[inv.sector == row.sector]
        t = g.time_btjd_tdb.values
        m = np.abs(t - row.tc) < HALF
        ex = g.cadence_s.values[0] / 86400.0
        if m.sum() < 8 or not accept(t[m], row.tc, ex):
            out.append((np.nan, np.nan, np.nan)); continue
        y = g.flux_raw_norm.values[m]
        yd = g.flux_norm.values[m]
        out.append((fit(t[m], y, row.tc, ex, 2) * 1e3, fit(t[m], y, row.tc, ex, 1) * 1e3,
                    fit(t[m][np.isfinite(yd)], yd[np.isfinite(yd)], row.tc, ex, 1) * 1e3 if np.isfinite(yd).sum() > 8 else np.nan))
    a = np.array(out)
    df = df.copy(); df["inv_q"], df["inv_l"], df["invdet_l"] = a[:, 0], a[:, 1], a[:, 2]
    return df


r = inv_depths(r); f = inv_depths(f)
for col in APS:
    ap = col.split("_")[1]
    for df in (r, f):
        df[col + "_ppt"] = [row[col] / (979 * S[str(int(row.sector))]["EE"][ap]) * 1e3 for _, row in df.iterrows()]
r.to_csv(f"{R3}/out/q04_paired_per_transit{tag}.csv", index=False)
res = {}
EARLY, LATE = (39, 66), (103, 104)
for col in APS:
    c2 = col + "_ppt"
    for invc in ("inv_q", "inv_l"):
        key = f"{invc}_vs_{col}"
        per = {}
        for sec in sorted(r.sector.unique()):
            rr = r[(r.sector == sec) & np.isfinite(r[invc]) & np.isfinite(r[c2])]
            ff = f[(f.sector == sec) & np.isfinite(f[invc]) & np.isfinite(f[c2])]
            n = len(rr)
            if n == 0: continue
            sd_inv = np.std(ff[invc]); sd_ap = np.std(ff[c2]); sd_diff = np.std(ff[invc] - ff[c2])
            rho = np.corrcoef(ff[invc], ff[c2])[0, 1]
            per[int(sec)] = dict(n=n, inv=float(rr[invc].mean()), inv_err=float(sd_inv / np.sqrt(n)), ap=float(rr[c2].mean()),
                                 ap_err=float(sd_ap / np.sqrt(n)), diff=float((rr[invc] - rr[c2]).mean()), diff_err=float(sd_diff / np.sqrt(n)),
                                 rho_fake=float(rho), nfake=int(len(ff)))
        def grp(ss, k):
            v = np.array([per[s][k] for s in ss if s in per]); e = np.array([per[s][k + "_err"] for s in ss if s in per])
            return float(np.sum(v / e ** 2) / np.sum(1 / e ** 2)), float(1 / np.sqrt(np.sum(1 / e ** 2)))
        con = {}
        for k in ("inv", "ap"):
            a, ae = grp(EARLY, k); b, be = grp(LATE, k)
            con[k] = [b - a, float(np.hypot(ae, be))]
        # contrast difference with paired errors (simple unweighted average of sector means to keep pairing exact)
        def cd(ss):
            return np.mean([per[s]["diff"] for s in ss]), np.sqrt(np.sum([per[s]["diff_err"] ** 2 for s in ss])) / len(ss)
        e1, e1e = cd(EARLY); l1, l1e = cd(LATE)
        con["contrast_inv_minus_ap_unweighted"] = [float(l1 - e1), float(np.hypot(e1e, l1e)), float((l1 - e1) / np.hypot(e1e, l1e))]
        con["global_diff_inv_minus_ap"] = grp(per.keys(), "diff")
        res[key] = dict(per_sector=per, contrast=con)
        print(f"{key:28s} " + " ".join(f"{s}:{per[s]['inv']:5.1f}/{per[s]['ap']:5.1f} d{per[s]['diff']:+5.1f}+-{per[s]['diff_err']:3.1f}" for s in per) +
              f" | contrast INV {con['inv'][0]:+.1f}+-{con['inv'][1]:.1f} AP {con['ap'][0]:+.1f}+-{con['ap'][1]:.1f}; (INV-AP) late-early {con['contrast_inv_minus_ap_unweighted'][0]:+.1f}+-{con['contrast_inv_minus_ap_unweighted'][1]:.1f}"
              f" ({con['contrast_inv_minus_ap_unweighted'][2]:+.2f} sig); glob diff {con['global_diff_inv_minus_ap'][0]:+.2f}+-{con['global_diff_inv_minus_ap'][1]:.2f}")
json.dump(res, open(f"{R3}/out/q04_paired{tag}.json", "w"), indent=1)
