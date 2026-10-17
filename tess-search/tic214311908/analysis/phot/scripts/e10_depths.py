import sys; sys.path.insert(0, ".")
from photlib import *
def tmpl(t, P=P_EPH, T0=T0_EPH, T14=0.067, fl=0.78):
    ph = ((t - T0) / P + 0.5) % 1 - 0.5; a = np.abs(ph * P); T23 = fl * T14
    return np.where(a <= T23 / 2, 1.0, np.where(a >= T14 / 2, 0.0, (T14 / 2 - a) / ((T14 - T23) / 2)))
tot = {}
for s in SECTORS:
    L = np.load(f"{ROOT}/data/lc_s{s}.npz"); t = L["t_tdb"]; ok = L["quality_ok"]
    y = tmpl(t); row = []
    ep = np.round((t - T0_EPH) / P_EPH); inn = ok & (y > 0)
    cov = {int(e): int(np.sum(inn & (ep == e))) for e in np.unique(ep[inn])}
    for key in ("M3V_B", "M3V_A", "pair_Bp", "pair_A", "Mira_B", "Mira_A"):
        d = L[f"det_{key}"]; m = ok & np.isfinite(d) & ~L[f"flare_{key}"]
        fref = json.load(open(f"{ROOT}/out/lc_stats.json"))[str(s)]["fref"][key.split("_")[0]]
        A = np.column_stack([np.ones(m.sum()), -y[m]]); c, *_ = np.linalg.lstsq(A, d[m], rcond=None)
        sd = robust_std(d[m & (y == 0)]); e = sd / np.sqrt(np.sum((y[m] - y[m].mean()) ** 2))
        row.append(f"{key} {c[1]*1e3:5.1f}+-{e*1e3:4.1f} ppt ({c[1]*fref:5.1f}+-{e*fref:4.1f} e/s)")
        tot.setdefault(key, []).append((c[1] * fref, e * fref))
    print(f"S{s}: " + " | ".join(row) + f" | in-transit cadences per epoch {cov}")
for k, v in tot.items():
    v = np.array(v); w = 1 / v[:, 1] ** 2; m = np.sum(w * v[:, 0]) / w.sum(); e = 1 / np.sqrt(w.sum()); chi = np.sum(w * (v[:, 0] - m) ** 2)
    print(f"{k}: weighted mean {m:.1f} +- {e:.1f} e/s; chi2 {chi:.1f}/5")
