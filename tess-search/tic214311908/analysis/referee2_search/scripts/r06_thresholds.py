"""r06: false-alarm thresholds and completeness, re-evaluated.

1. Analyst's null maxima (search/out/null_*_clean_*.json, search/logs/s02_inv_roll_1_killed.log):
   - the recurring 119.5 maximum of the inverted+rolled nulls is diagnosed in r06b (it is NOT a wrap junction: it is
     a 4-bin, 45-ppt positive excursion in S12 at BTJD 1634.41 that straddles an 8-h data gap in the real data;
     rolling re-assigns those four values to contiguous time stamps, so the inverted copy exposes it in full);
   - Gumbel fits (method of moments and ML) to: all 60 maxima as used by the analyst; the 30 non-inverted rolled
     maxima only; inverted maxima with the gap-joined 119.5 value removed.
   - per-band (0.3-1, 1-3, 3-10, 10-30 d) Gumbel thresholds at FAP 0.25 % per band (4 bands -> 1 % global).
2. Recovery fractions of the analyst's 800 injections (search/out/inject_*.csv, column dchi2_match) with the global
   threshold 147.5 and with the per-band thresholds -> R50 per period bin (logistic fit in log R).
3. The same for referee-2's own injections (out/r05_inject_*.csv) against referee-2's own nulls (out/r03_null_*.json),
   if present.
Output: out/r06_thresholds.json
"""
import sys, glob, re
sys.path.insert(0, "tic214311908/analysis/referee2_search/scripts")
from rlib2 import *
from scipy import stats, optimize

S = f"{BASE}/search"
BANDS = ["0.3-1.0", "1.0-3.0", "3.0-10.0", "10.0-30.0"]
BEDGE = [(0.3, 1.0), (1.0, 3.0), (3.0, 10.0), (10.0, 30.0)]


def gumbel_thr(x, fap):
    x = np.asarray(x, float)
    mu, beta = stats.gumbel_r.fit(x)
    return float(stats.gumbel_r.ppf(1 - fap, mu, beta)), float(mu), float(beta)


def r50(df, thr_of_P, Pcol, Rcol, dcol, bins):
    out = {}
    for a, b in bins:
        m = (df[Pcol] >= a) & (df[Pcol] < b)
        x = np.log(df[Rcol][m].values); det = (df[dcol][m].values >= np.array([thr_of_P(p) for p in df[Pcol][m].values])).astype(float)
        if m.sum() < 10 or det.min() == det.max():
            out[f"{a}-{b}"] = dict(n=int(m.sum()), frac=float(det.mean())); continue
        def nll(p):
            z = (x - p[0]) / max(p[1], 1e-3)
            pr = np.clip(1 / (1 + np.exp(-z)), 1e-9, 1 - 1e-9)
            return -np.sum(det * np.log(pr) + (1 - det) * np.log(1 - pr))
        r = optimize.minimize(nll, [np.log(3.0), 0.1], method="Nelder-Mead")
        R50 = float(np.exp(r.x[0])); R90 = float(np.exp(r.x[0] + abs(r.x[1]) * np.log(9)))
        out[f"{a}-{b}"] = dict(n=int(m.sum()), R50=R50, R90=R90, frac=float(det.mean()))
    return out


def main():
    out = {}
    # ---------------- 1 analyst nulls
    fam = {}
    for fn in glob.glob(f"{S}/out/null_*_clean_*.json"):
        d = json.load(open(fn)); k = d["kind"]
        if k in ("roll", "inv_roll"):
            fam.setdefault(k, []).extend(d["res"])
    killed = [float(l.split()[1]) for l in open(f"{S}/logs/s02_inv_roll_1_killed.log") if re.match(r"^\d+ ", l)]
    roll = np.array([r["max"] for r in fam["roll"]]); inv = np.array([r["max"] for r in fam["inv_roll"]])
    allmax = np.concatenate([roll, inv, killed])
    out["analyst_null_counts"] = dict(roll=len(roll), inv_roll=len(inv), killed_inv=len(killed),
                                      n_equal_119p5=int(np.sum(np.abs(np.concatenate([inv, killed]) - 119.5) < 0.05)))
    th = {}
    th["all60_as_analyst"] = gumbel_thr(allmax, 0.01)
    th["roll30_only"] = gumbel_thr(roll, 0.01)
    inv_nj = np.array([v for v in np.concatenate([inv, killed]) if abs(v - 119.5) > 0.05])
    th["inverted_without_gapjoined_value"] = gumbel_thr(inv_nj, 0.01) if len(inv_nj) >= 5 else None
    th["roll_plus_inv_without_gapjoined"] = gumbel_thr(np.concatenate([roll, inv_nj]), 0.01)
    out["analyst_global_thresholds_fap1"] = th
    # per band (only the 46 realisations with band maxima)
    bands = {}
    for b in BANDS:
        v = np.array([r["bands"][b] for r in fam["roll"]] + [r["bands"][b] for r in fam["inv_roll"]])
        vr = np.array([r["bands"][b] for r in fam["roll"]])
        bands[b] = dict(n=len(v), median=float(np.median(v)), max=float(v.max()), thr_fap0p25=gumbel_thr(v, 0.0025)[0],
                        thr_fap0p25_roll_only=gumbel_thr(vr, 0.0025)[0])
    out["analyst_band_thresholds"] = bands
    real = np.load(f"{S}/data/bls_real_all_clean.npz"); Pr = 1 / real["freq"]
    out["analyst_real_band_max"] = {f"{a}-{b}": float(real["dchi2"][(Pr >= a) & (Pr < b)].max()) for a, b in BEDGE}
    print(json.dumps(out, indent=0), flush=True)
    # ---------------- 2 analyst injections
    inj = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{S}/out/inject_*.csv"))], ignore_index=True)
    gthr = 147.4747783814821

    def band_thr_fn(bt):
        def f(P_):
            for (a, b), k in zip(BEDGE, BANDS):
                if a <= P_ < b: return bt[k]
            return bt[BANDS[-1]]
        return f
    bt = {k: bands[k]["thr_fap0p25"] for k in BANDS}
    pbins = [(0.5, 1), (1, 2), (2, 4), (4, 8), (8, 15), (15, 30)]
    out["analyst_inj_global"] = r50(inj, lambda p: gthr, "P", "rp", "dchi2_match", pbins)
    out["analyst_inj_band"] = r50(inj, band_thr_fn(bt), "P", "rp", "dchi2_match", pbins)
    out["analyst_inj_band_thresholds_used"] = bt
    print("analyst injections, global:", out["analyst_inj_global"])
    print("analyst injections, band  :", out["analyst_inj_band"])
    # ---------------- 3 referee-2 nulls and injections
    mine = []
    for fn in glob.glob(f"{ME}/out/r03_null_*.json"):
        d = json.load(open(fn)); mine += [dict(kind=d["kind"], **r) for r in d["res"]]
    if len(mine) >= 8:
        mx = np.array([r["max"] for r in mine])
        g = gumbel_thr(mx, 0.01)
        mb = {}
        for b in BANDS:
            v = np.array([r["bands"][b] for r in mine])
            mb[b] = dict(n=len(v), median=float(np.median(v)), max=float(v.max()), thr_fap0p25=gumbel_thr(v, 0.0025)[0])
        rr = json.load(open(f"{ME}/out/r03_real.json"))
        p_emp = float(np.mean(mx >= rr["real"]["peaks"][0]["dchi2"]))
        out["mine_nulls"] = dict(n=len(mine), kinds={k: int(sum(1 for r in mine if r["kind"] == k)) for k in set(r["kind"] for r in mine)},
                                 maxima=mx.tolist(), median=float(np.median(mx)), thr_fap1=g[0], gumbel=g[1:], bands=mb,
                                 real_top=rr["real"]["peaks"][0]["dchi2"], p_emp_real_top=p_emp,
                                 real_band_max=rr["real"]["bands"],
                                 p_emp_real_bands={b: float(np.mean(np.array([r["bands"][b] for r in mine]) >= rr["real"]["bands"][b])) for b in BANDS},
                                 by_kind={k: dict(median=float(np.median([r["max"] for r in mine if r["kind"] == k])),
                                                  max=float(np.max([r["max"] for r in mine if r["kind"] == k]))) for k in set(r["kind"] for r in mine)})
        print("mine nulls:", json.dumps(out["mine_nulls"], indent=0))
        files = sorted(glob.glob(f"{ME}/out/r05_inject_*.csv"))
        if files:
            I = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
            mbt = {k: mb[k]["thr_fap0p25"] for k in BANDS}
            tab = {}
            for (Pc, Rc), gdf in I.groupby(["Pcell", "Rcell"]):
                f = band_thr_fn(mbt)
                tab[f"P{Pc}_R{Rc}"] = dict(n=len(gdf), frac_global=float(np.mean(gdf.d2_match >= g[0])),
                                           frac_band=float(np.mean(gdf.d2_match >= np.array([f(p) for p in gdf.P]))),
                                           median_d2=float(gdf.d2_match.median()))
            out["mine_injections"] = dict(global_thr=g[0], band_thr=mbt, cells=tab)
            print(pd.DataFrame(tab).T.to_string())
    jdump(out, f"{ME}/out/r06_thresholds.json")


if __name__ == "__main__":
    main()
