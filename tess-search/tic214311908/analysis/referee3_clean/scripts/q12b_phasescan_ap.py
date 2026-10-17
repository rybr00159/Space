"""q12b: occultation-phase scan 0.40-0.60 (step 0.0025) in the referee3 raw-pixel apertures (q02 run with
Q2_PHRANGE, tag _ph49). Depth (transit-shaped template, same duration) per phase = weighted mean over events;
error per event = per-sector std of all fake phases of the main q02 run (0.12-0.88). S/N calibrated by the scan's
robust std. Output: out/q12b_phasescan_ap.json"""
import json, numpy as np, pandas as pd
R3 = "tic214311908/analysis/referee3_clean"
f = pd.read_csv(f"{R3}/out/q02_bigap_fake_ph49.csv"); f0 = pd.read_csv(f"{R3}/out/q02_bigap_fake.csv"); S = json.load(open(f"{R3}/out/q02_bigap.json"))
out = {}
for col in ("L_box3_pt", "L_m25_pt"):
    ap = col.split("_")[1]
    for d_ in (f, f0):
        d_["d"] = [row[col] / (979 * S[str(int(row.sector))]["EE"][ap]) * 1e3 for _, row in d_.iterrows()]
    sd = f0.groupby("sector").d.std()
    rows = []
    for ph, g in f.groupby("phase"):
        w = 1 / sd.loc[g.sector].values ** 2
        rows.append((ph, np.sum(g.d.values * w) / np.sum(w), 1 / np.sqrt(np.sum(w))))
    a = np.array(rows); sn = a[:, 1] / a[:, 2]; rs = 1.4826 * np.median(np.abs(sn - np.median(sn)))
    sel = (a[:, 0] >= 0.4899) & (a[:, 0] <= 0.4951)
    out[col] = dict(scan_robust_std=float(rs), at_0490_0495=[[float(x) for x in r] for r in a[sel]], at_05=[float(x) for x in a[np.isclose(a[:, 0], 0.5)][0]],
                    max_sn=float(np.max(sn / rs)), phase_max=float(a[np.argmax(sn), 0]))
    print(col, "robust std of S/N", round(rs, 2), "0.490-0.495:", np.round(a[sel], 2).tolist(), "0.5:", np.round(a[np.isclose(a[:, 0], 0.5)][0], 2), "max S/N cal", round(out[col]["max_sn"], 2), "at", out[col]["phase_max"])
json.dump(out, open(f"{R3}/out/q12b_phasescan_ap.json", "w"), indent=1)
