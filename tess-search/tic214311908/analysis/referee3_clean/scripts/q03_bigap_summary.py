"""q03: summarise q02 (large-aperture flux losses): per-aperture fake-epoch scatter, per-sector depths in ppt
(EE-corrected), global depth, S39+S66 vs S103+S104 contrast, correlation with the analyst's per-transit depths.
Output: out/q03_bigap_summary.json"""
import json, numpy as np, pandas as pd
R3 = "tic214311908/analysis/referee3_clean"
CL = "tic214311908/analysis/clean"
import sys
tag = sys.argv[1] if len(sys.argv) > 1 else ""
r = pd.read_csv(f"{R3}/out/q02_bigap_per_transit{tag}.csv"); f = pd.read_csv(f"{R3}/out/q02_bigap_fake{tag}.csv")
S = json.load(open(f"{R3}/out/q02_bigap{tag}.json"))
an = pd.read_csv(f"{CL}/out/per_transit_clean.csv")
out = {}
for col in [c for c in r.columns if c.startswith("L_")]:
    ap = col.split("_")[1]
    res = {}
    persec = {}
    for sec, g in r.groupby("sector"):
        if str(int(sec)) not in S: continue
        fk = f[f.sector == sec][col]
        sd = float(np.std(fk)); rsd = float(1.4826 * np.median(np.abs(fk - np.median(fk))))
        ee = S[str(sec)]["EE"][ap]
        # per-sector mean flux loss, error from fake std / sqrt(n); also an empirical 'fake sector' check of the mean
        L = g[col].values
        persec[int(sec)] = dict(n=int(len(L)), L_mean=float(L.mean()), L_err=sd / np.sqrt(len(L)), fake_std=sd, fake_rstd=rsd,
                                 fake_mean=float(fk.mean()), nfake=int(len(fk)), EE=ee,
                                 depth_ppt=float(L.mean() / (979 * ee) * 1e3), depth_err=float(sd / np.sqrt(len(L)) / (979 * ee) * 1e3),
                                 per_transit_L=[float(x) for x in L], epochs=[int(x) for x in g.epoch])
    secs = sorted(persec)
    d = np.array([persec[s]["depth_ppt"] for s in secs]); e = np.array([persec[s]["depth_err"] for s in secs])
    wm = np.sum(d / e ** 2) / np.sum(1 / e ** 2); we = 1 / np.sqrt(np.sum(1 / e ** 2))
    chi2 = float(np.sum(((d - wm) / e) ** 2))
    def grp(ss):
        dd = np.array([persec[s]["depth_ppt"] for s in ss]); ee_ = np.array([persec[s]["depth_err"] for s in ss])
        return np.sum(dd / ee_ ** 2) / np.sum(1 / ee_ ** 2), 1 / np.sqrt(np.sum(1 / ee_ ** 2))
    a, ae = grp([39, 66]); b, be = grp([103, 104])
    ea, eae = grp([12, 39, 66]); la, lae = grp([93, 103, 104])
    res = dict(per_sector=persec, global_depth=[float(wm), float(we)], chi2_sectors=[chi2, len(secs) - 1],
               late_minus_early_39_66_vs_103_104=[float(b - a), float(np.hypot(ae, be)), float((b - a) / np.hypot(ae, be))],
               late_minus_early_12_39_66_vs_93_103_104=[float(la - ea), float(np.hypot(eae, lae)), float((la - ea) / np.hypot(eae, lae))])
    out[col] = res
    print(f"{col:14s} glob {wm:6.2f}+-{we:4.2f} chi2 {chi2:5.1f}/5  S39+66 {a:5.1f}+-{ae:3.1f} S103+104 {b:5.1f}+-{be:3.1f} diff {b-a:+5.1f}+-{np.hypot(ae,be):3.1f} | " +
          " ".join(f"{s}:{persec[s]['depth_ppt']:5.1f}+-{persec[s]['depth_err']:3.1f}" for s in secs) +
          " | fake std e/s " + " ".join(f"{persec[s]['fake_std']:4.1f}" for s in secs) + " fake mean " + " ".join(f"{persec[s]['fake_mean']:+4.1f}" for s in secs))
# per-transit comparison with analyst's per-transit depths
best = "L_g40_pt" if "L_g40_pt" in r.columns else None
cmp_ = {}
for col in [c for c in r.columns if c.startswith("L_")]:
    ap = col.split("_")[1]
    rr = r.copy(); rr["dep"] = [row[col] / (979 * S[str(int(row.sector))]["EE"][ap]) * 1e3 for _, row in rr.iterrows()]
    for pl in ("inv", "ref"):
        a2 = an[an.pipeline == pl].merge(rr[["sector", "epoch", "dep"]], on=["sector", "epoch"])
        if len(a2) > 3:
            cmp_[f"{col}_{pl}"] = dict(n=int(len(a2)), corr=float(np.corrcoef(a2.depth_ppt, a2.dep)[0, 1]),
                                       mean_diff=float(np.mean(a2.dep - a2.depth_ppt)), rms_diff=float(np.std(a2.dep - a2.depth_ppt)))
out["compare_with_analyst_per_transit"] = cmp_
for k, v in cmp_.items():
    if "pt" in k: print(k, v)
json.dump(out, open(f"{R3}/out/q03_bigap_summary{tag}.json", "w"), indent=1)
