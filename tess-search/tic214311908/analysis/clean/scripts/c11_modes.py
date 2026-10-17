"""c11: what the cleaning does to the transit. The same per-transit measurement (c05.fit_window: fixed batman shape,
depth at the predicted mid-time, linear baseline, +-0.30 d window, identical acceptance rules) applied to every
version of the light curves:
  INV modes none (= phase-1 phot, reproduced exactly), rr3, all (adopted), allreg; flux versions raw and det.
  REF modes orig (= the referee's own photometry) and all (adopted).
Errors: per-window formal error (window OOT scatter) x beta_s, beta_s = std(depth/formal) of the same fit at 60 fake
epochs per sector (>= 1). Reported: per-transit depth changes,
weighted global depth, per-sector depths and the S39+S66 vs S103+S104 and early/late contrasts, odd-even.
Also: the per-transit depth change (all - none) against the RRab leak bias predicted by referee2_search r01 (slope).
Outputs: out/modes_compare.json, out/modes_per_transit.csv
"""
import sys
sys.path.insert(0, "tic214311908/analysis/clean/scripts")
from cleanlib import *
from c05_inject import fit_window, accept
import pandas as pd

rng = np.random.default_rng(1111)
EPH = json.load(open(f"{CL}/out/ephemeris_final.json"))["combined"]
P, T0 = EPH["P_d"], EPH["T0_epoch0_BTJD_TDB"]
FAKE = {}


def centres(t, s):
    tf = t[np.isfinite(t)]
    E = np.arange(np.ceil((tf.min() - T0) / P - 0.1), np.floor((tf.max() - T0) / P + 0.1) + 1).astype(int)
    if s not in FAKE:
        fk = []
        while len(fk) < 60:
            tc = rng.uniform(tf.min() + 0.3, tf.max() - 0.3)
            if abs((((tc - T0) / P + 0.5) % 1 - 0.5) * P) < 0.45: continue
            if fk and min(abs(f - tc) for f in fk) < 0.15: continue
            fk.append(tc)
        FAKE[s] = fk
    return [(T0 + e * P, int(e)) for e in E], [(tc, -9999) for tc in FAKE[s]]


def inv_rows(mode, ver):
    rows = []
    for s in SECTORS:
        z = np.load(f"{CL}/data/lcinv_{mode}_s{s}.npz")
        t = z["t_tdb"]; f = z[f"{ver}_M3V_B"]; ok = z["quality_ok"].astype(bool) & np.isfinite(f)
        det = z["det_M3V_B"]; sd = robust_std(det[ok & ~z["intransit_window"].astype(bool)])
        hi = ok & (det > 1 + 4 * sd); ok &= ~(hi | np.roll(hi, 1) | np.roll(hi, -1))
        real, fake = centres(t[ok], s)
        for kind, cen in (("real", real), ("fake", fake)):
            for tc, ep in cen:
                m = ok & (np.abs(t - tc) < 0.3)
                if m.sum() < 8 or not accept(t[m], tc, 0.3, float(z["exptime"])): continue
                r = fit_window(t[m], f[m], tc, float(z["exptime"]), free=False)
                if r: rows.append(dict(version=f"inv_{mode}_{ver}", sector=s, kind=kind, epoch=ep, tc=tc, depth=r["depth_fix"], formal=r["depth_fix_err"]))
    return rows


def ref_rows(mode):
    rows = []
    for s in SECTORS:
        z = np.load(f"{CL}/data/refwin_{mode}_s{s}.npz")
        ex = float(z["exptime"])
        for kind in ("strict", "fake"):
            for j, tc in enumerate(z[f"{kind}_tc"]):
                m = z[f"{kind}_win"] == j
                tt = z[f"{kind}_t"][m]; ff = 1 + z[f"{kind}_aM"][m] / F_M3V
                if not accept(tt, tc, {12: 0.6, 39: 0.4}.get(s, 0.35), ex): continue
                r = fit_window(tt, ff, tc, ex, free=False)
                ep = int(z["strict_epoch"][j]) if kind == "strict" else -9999
                if r: rows.append(dict(version=f"ref_{mode}", sector=s, kind="real" if kind == "strict" else "fake", epoch=ep, tc=tc, depth=r["depth_fix"], formal=r["depth_fix_err"]))
    return rows


def summarise(d):
    real = d[d.kind == "real"].copy(); fake = d[d.kind == "fake"]
    # per-window error = formal (window's own OOT scatter) x beta_s, beta_s = std of fake-epoch depth / formal (sector)
    beta = {s: float(np.std(fake[fake.sector == s].depth / fake[fake.sector == s].formal)) for s in SECTORS}
    err = {s: float(np.std(fake[fake.sector == s].depth)) for s in SECTORS}
    real["err"] = real.formal * real.sector.map(beta).clip(lower=1.0)
    w = 1 / real.err ** 2
    mu = float(np.sum(w * real.depth) / np.sum(w)); emu = float(1 / np.sqrt(np.sum(w)))
    chi2 = float(np.sum(w * (real.depth - mu) ** 2))
    def wm(sel):
        x = real[sel]; ww = 1 / x.err ** 2
        return float(np.sum(ww * x.depth) / np.sum(ww)), float(1 / np.sqrt(np.sum(ww)))
    ps = {int(s): wm(real.sector == s) for s in SECTORS if np.any(real.sector == s)}
    a = wm(real.sector.isin([39, 66])); b = wm(real.sector.isin([103, 104]))
    e1 = wm(real.sector.isin([12, 39, 66])); l1 = wm(real.sector.isin([93, 103, 104]))
    o = wm(real.epoch % 2 != 0); ev = wm(real.epoch % 2 == 0)
    return dict(n=len(real), epochs=real.epoch.tolist(), depth=mu, depth_err=emu, chi2=chi2, dof=len(real) - 1, per_sector=ps,
                single_err=err, beta=beta, S39S66=a, S103S104=b, late_minus_early_pair=[b[0] - a[0], float(np.hypot(a[1], b[1]))],
                late_minus_early_all=[l1[0] - e1[0], float(np.hypot(e1[1], l1[1]))], odd_minus_even=[o[0] - ev[0], float(np.hypot(o[1], ev[1]))])


def main():
    rows = []
    for mode in ("none", "rr3", "all", "allreg"):
        for ver in ("raw", "det"):
            rows += inv_rows(mode, ver)
    for mode in ("orig", "all"):
        rows += ref_rows(mode)
    D = pd.DataFrame(rows)
    D.to_csv(f"{CL}/out/modes_per_transit.csv", index=False, float_format="%.4f")
    out = {}
    for v, g in D.groupby("version"):
        out[v] = summarise(g)
        r = out[v]
        print(f"{v:18s}: n {r['n']:2d} depth {r['depth']:.2f}+-{r['depth_err']:.2f} chi2 {r['chi2']:.1f}/{r['dof']}; S39+66 {r['S39S66'][0]:.1f}+-{r['S39S66'][1]:.1f} "
              f"S103+104 {r['S103S104'][0]:.1f}+-{r['S103S104'][1]:.1f} (diff {r['late_minus_early_pair'][0]:+.1f}+-{r['late_minus_early_pair'][1]:.1f}); "
              f"late-early(all) {r['late_minus_early_all'][0]:+.1f}+-{r['late_minus_early_all'][1]:.1f}; odd-even {r['odd_minus_even'][0]:+.1f}+-{r['odd_minus_even'][1]:.1f}", flush=True)
    # per-transit change from cleaning vs the predicted RRab leak bias (referee2_search r01)
    tab = pd.DataFrame(json.load(open(f"{BASE}/referee2_search/out/r01_rrab_leak.json"))["D"]["table"])[["epoch", "bias_ppt"]]
    res = {}
    for a, b in (("inv_none_raw", "inv_all_raw"), ("inv_none_det", "inv_all_det"), ("ref_orig", "ref_all")):
        x = D[(D.version == a) & (D.kind == "real")].set_index("epoch").depth
        y = D[(D.version == b) & (D.kind == "real")].set_index("epoch").depth
        j = pd.concat([x.rename("before"), y.rename("after")], axis=1, join="inner").join(tab.set_index("epoch"), how="inner")
        ch = j.after - j.before
        A = np.column_stack([np.ones(len(j)), j.bias_ppt.values]); c, *_ = np.linalg.lstsq(A, ch.values, rcond=None)
        r_ = ch.values - A @ c; s2 = np.sum(r_ ** 2) / max(len(j) - 2, 1); cov = s2 * np.linalg.inv(A.T @ A)
        # regression of 'after' depths on the predicted bias (should be ~0 if the leak is gone)
        A2 = np.column_stack([np.ones(len(j)), j.bias_ppt.values]); c2, *_ = np.linalg.lstsq(A2, j.after.values, rcond=None)
        r2 = j.after.values - A2 @ c2; cov2 = np.sum(r2 ** 2) / max(len(j) - 2, 1) * np.linalg.inv(A2.T @ A2)
        A3 = np.column_stack([np.ones(len(j)), j.bias_ppt.values]); c3, *_ = np.linalg.lstsq(A3, j.before.values, rcond=None)
        r3 = j.before.values - A3 @ c3; cov3 = np.sum(r3 ** 2) / max(len(j) - 2, 1) * np.linalg.inv(A3.T @ A3)
        res[f"{a}->{b}"] = dict(n=len(j), change_rms=float(np.std(ch)), change_max=float(np.max(np.abs(ch))),
                                slope_change_vs_bias=[float(c[1]), float(np.sqrt(cov[1, 1]))],
                                slope_before_vs_bias=[float(c3[1]), float(np.sqrt(cov3[1, 1]))],
                                slope_after_vs_bias=[float(c2[1]), float(np.sqrt(cov2[1, 1]))],
                                per_epoch={int(k): [round(float(v1), 2), round(float(v2), 2), round(float(bb), 2)] for k, v1, v2, bb in zip(j.index, j.before, j.after, j.bias_ppt)})
        print(f"{a}->{b}: n {len(j)} change rms {np.std(ch):.2f} max {np.max(np.abs(ch)):.2f} ppt; slope(change vs predicted leak) {c[1]:+.2f}+-{np.sqrt(cov[1,1]):.2f} "
              f"(expected -1 for full removal of an INV-like leak); slope before {c3[1]:+.2f}+-{np.sqrt(cov3[1,1]):.2f}, after {c2[1]:+.2f}+-{np.sqrt(cov2[1,1]):.2f}", flush=True)
    out["leak_regression"] = res
    jd(out, f"{CL}/out/modes_compare.json")


if __name__ == "__main__":
    main()
