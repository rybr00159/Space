"""q12: secondary eclipse at phase 0.50 (and 0.48/0.52) in the referee3 raw-pixel apertures (q02 fake table):
weighted mean over all events of the fitted 'depth' at that phase (transit-shaped template of the same duration),
errors from the per-sector scatter over all fake phases.  Also the INV delivered curve with the same windows (q04).
Output: out/q12_secondary_ap.json"""
import json, numpy as np, pandas as pd
R3 = "tic214311908/analysis/referee3_clean"
f = pd.read_csv(f"{R3}/out/q02_bigap_fake.csv"); S = json.load(open(f"{R3}/out/q02_bigap.json"))
out = {}
for col in ("L_box3_pt", "L_m25_pt"):
    ap = col.split("_")[1]
    f["d"] = [row[col] / (979 * S[str(int(row.sector))]["EE"][ap]) * 1e3 for _, row in f.iterrows()]
    sd = f.groupby("sector").d.std()
    for ph in (0.48, 0.50, 0.52):
        g = f[np.isclose(f.phase, ph)]
        w = 1 / sd.loc[g.sector].values ** 2
        m = float(np.sum(g.d.values * w) / np.sum(w)); e = float(1 / np.sqrt(np.sum(w)))
        out[f"{col}_ph{ph}"] = dict(n=int(len(g)), depth_ppt=m, err=e, ul3=m + 3 * e if m > 0 else 3 * e)
        print(col, ph, round(m, 2), "+-", round(e, 2), "n", len(g))
json.dump(out, open(f"{R3}/out/q12_secondary_ap.json", "w"), indent=1)
