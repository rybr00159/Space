"""q15: split-independent (omnibus) test of sector-dependent INV-minus-aperture depth differences, from the paired
per-sector differences of q04 (all variants): chi2 of the six per-sector differences about their weighted mean,
p from chi2(5). Output: out/q15_omnibus.json"""
import json, glob, numpy as np
from scipy.stats import chi2 as C2
R3 = "tic214311908/analysis/referee3_clean"
out = {}
for fn in sorted(glob.glob(f"{R3}/out/q04_paired*.json")):
    tag = fn.split("q04_paired")[1].replace(".json", "") or "_default"
    d = json.load(open(fn))
    for key in ("inv_l_vs_L_box3_pt", "inv_l_vs_L_m25_pt"):
        ps = d[key]["per_sector"]
        x = np.array([v["diff"] for v in ps.values()]); e = np.array([v["diff_err"] for v in ps.values()])
        w = 1 / e ** 2; m = np.sum(x * w) / np.sum(w); c = float(np.sum(((x - m) / e) ** 2))
        out[f"{tag}:{key}"] = dict(chi2=c, dof=len(x) - 1, p=float(C2.sf(c, len(x) - 1)), mean=float(m), contrast=d[key]["contrast"]["contrast_inv_minus_ap_unweighted"])
        print(f"{tag:10s} {key:22s} chi2 {c:5.1f}/{len(x)-1} p {C2.sf(c, len(x)-1):.2f}  mean INV-AP {m:+.2f}  split contrast diff {d[key]['contrast']['contrast_inv_minus_ap_unweighted'][0]:+.1f}+-{d[key]['contrast']['contrast_inv_minus_ap_unweighted'][1]:.1f}")
json.dump(out, open(f"{R3}/out/q15_omnibus.json", "w"), indent=1)
