import sys; sys.path.insert(0, ".")
from photlib import *
from p06_lc import detrend, in_transit_windows
edges = [0, 5, 10, 20, 40, 80, 160, 320, 640, 1e9]
for s in SECTORS:
    z = np.load(f"{ROOT}/data/phot_s{s}.npz"); J = json.load(open(f"{ROOT}/out/phot_s{s}.json"))
    t = z["t"]; bm = np.median(z["bkg"][z["q0"]]); ex = z["bkg"] - bm
    ok0 = z["qok"] & z["jit_ok"] & np.isfinite(z["F_B"][:, 1]) & (ex < 1000)
    raw = 1 + z["F_B"][:, 1] / 940.0
    tr = detrend(t, raw, ok0 & (ex < 40), in_transit_windows(t), window=0.75)
    det = raw / tr
    chi = z["chi2_B"] / np.nanmedian(z["chi2_B"][z["q0"]])
    row = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = ok0 & (ex >= a) & (ex < b) & np.isfinite(det) & ~in_transit_windows(t)
        if m.sum() > 20: row.append(f"[{a:.0f},{b:.0f}) n={m.sum()} rstd={robust_std(det[m])*1e3:.1f} med={np.median(det[m]-1)*1e3:+.1f} chi={np.median(chi[m]):.1f}")
    print(f"S{s}: bkg median {bm:.0f}; " + " | ".join(row))
