"""q10: odd/even depth difference with fake-epoch per-sector errors, for the INV delivered light curve, the REF
delivered windows and the referee3 raw-pixel apertures (box3, r <= 2.5 px), at the analyst's final ephemeris
(epoch parity of the phase-1 numbering, E = 0 at BTJD 3840.915).
Per-transit depth: fixed shape (q02 shape), linear (INV/REF) or as in q02/q04; per-transit error = per-sector std of
the same fit at fake epochs (INV, apertures: q04 fake windows; REF: clean/data/win_ref_fake.npz).
Also a permutation null: the odd-even statistic for 2000 random relabellings of transits into two groups of the
same sizes (gives an empirical p that does not rely on the error model).
Output: out/q10_oddeven.json"""
import json, sys, numpy as np, pandas as pd
sys.path.insert(0, "tic214311908/analysis/referee3_clean/scripts")
from q02_bigap import shape
R3 = "tic214311908/analysis/referee3_clean"
CL = "tic214311908/analysis/clean"
T0A, PA = 3840.91488, 5.6915020
rng = np.random.default_rng(10)
tag = sys.argv[1] if len(sys.argv) > 1 else ""
r = pd.read_csv(f"{R3}/out/q04_paired_per_transit{tag}.csv")
import glob
fk = pd.read_csv(f"{R3}/out/q02_bigap_fake{tag}.csv")


def lin(t, f, tc, ex):
    s = shape(t, tc, ex); tau = t - tc
    M = np.column_stack([-s, np.ones_like(t), tau]); c, *_ = np.linalg.lstsq(M, f, rcond=None); return c[0] * 1e3


# REF per-transit (fixed tc) and fake errors
ref = pd.read_csv(f"{CL}/lc_M3V_REF_clean.csv")
rr = []
for (sec, ep), g in ref.groupby(["sector", "epoch"]):
    rr.append(dict(sector=sec, epoch=int(ep), d=lin(g.time_btjd_tdb.values, g.flux_norm.values, T0A + ep * PA, g.cadence_s.values[0] / 86400)))
rr = pd.DataFrame(rr)
z = np.load(f"{CL}/data/win_ref_fake.npz"); refsd = {}
for sec in np.unique(z["sector"]):
    v = []
    for wi in np.where(z["sector"] == sec)[0]:
        m = z["w"] == wi; t = z["t"][m]
        if np.sum(np.abs(t - z["tc0"][wi]) < 0.028) < 2: continue
        v.append(lin(t, z["f"][m], z["tc0"][wi], z["ex"][m][0]))
    refsd[int(sec)] = float(np.std(v))
rr["err"] = [refsd[s] for s in rr.sector]

# INV per-transit from q04 (linear baseline) needs fake errors: recompute INV fake depths from q04 is not stored -> use
# the q04 fake table recomputation: q04 stored only real; recompute fake INV std from the delivered LC here
inv = pd.read_csv(f"{CL}/lc_M3V_INV_clean.csv", usecols=["time_btjd_tdb", "sector", "cadence_s", "quality_ok", "flux_raw_norm"])
inv = inv[(inv.quality_ok == 1) & np.isfinite(inv.flux_raw_norm)]
invsd = {}
for sec, g in inv.groupby("sector"):
    t = g.time_btjd_tdb.values; f = g.flux_raw_norm.values; ex = g.cadence_s.values[0] / 86400; v = []
    for tc in fk[fk.sector == sec].tc.values:
        m = np.abs(t - tc) < 0.30
        if m.sum() < 8 or np.sum(np.abs(t[m] - tc) < 0.028) < 2: continue
        v.append(lin(t[m], f[m], tc, ex))
    invsd[int(sec)] = float(np.std(v))
S = json.load(open(f"{R3}/out/q02_bigap{tag}.json"))
apsd = {}
for col in ("L_box3_pt", "L_m25_pt"):
    ap = col.split("_")[1]
    apsd[col] = {int(s): float(np.std(fk[fk.sector == s][col]) / (979 * S[str(int(s))]["EE"][ap]) * 1e3) for s in np.unique(fk.sector)}

sets = {"INV": (r[["sector", "epoch"]].assign(d=r.inv_l, err=[invsd[s] for s in r.sector])),
        "REF": rr, "REF_noE58": rr[~((rr.sector == 103) & (rr.epoch == 58))],
        "AP_box3": r[["sector", "epoch"]].assign(d=r.L_box3_pt_ppt, err=[apsd["L_box3_pt"][s] for s in r.sector]),
        "AP_m25": r[["sector", "epoch"]].assign(d=r.L_m25_pt_ppt, err=[apsd["L_m25_pt"][s] for s in r.sector])}
out = {}
for k, d in sets.items():
    d = d[np.isfinite(d.d)]
    odd = (d.epoch % 2 != 0).values
    w = 1 / d.err.values ** 2
    def stat(o):
        mo = np.sum(d.d.values[o] * w[o]) / np.sum(w[o]); me = np.sum(d.d.values[~o] * w[~o]) / np.sum(w[~o])
        return mo - me, np.sqrt(1 / np.sum(w[o]) + 1 / np.sum(w[~o]))
    oe, oee = stat(odd)
    perm = np.array([stat(rng.permutation(odd))[0] for _ in range(4000)])
    out[k] = dict(n=int(len(d)), n_odd=int(odd.sum()), odd_minus_even=float(oe), err=float(oee), z=float(oe / oee),
                  p_perm_two_sided=float(np.mean(np.abs(perm) >= abs(oe))), perm_std=float(np.std(perm)))
    print(k, {kk: (round(v, 3) if isinstance(v, float) else v) for kk, v in out[k].items()})
out["per_sector_err_ppt"] = dict(INV=invsd, REF=refsd, AP=apsd)
json.dump(out, open(f"{R3}/out/q10_oddeven{tag}.json", "w"), indent=1)
