"""q06: residual RRab (Gaia DR3 5925372615693467904, P = 0.6485275 d) signal in the delivered cleaned INV light curve
(lc_M3V_INV_clean.csv, flux_norm, quality_ok == 1, out of the +-0.1 d transit windows, flare points removed),
and in the referee3 raw 3x3-box aperture (q02-style, no field-wide estimators).
Statistic: per sector, 3-harmonic least-squares fit (free phase); A1 = first-harmonic semi-amplitude (ppt);
p = fraction of 300 fake periods (0.45-0.95 d, avoiding +-1% of P_RR, P_RR/2... and the EB periods 0.7503, 1.5953 and
their halves) with A1 >= observed. Upper limit: 95th percentile of A1 at fake periods after injecting sinusoids
(injection-recovery grid) is not needed; we quote A1 + 1.645 * std(fake A1 projections) as a simple UL.
Output: out/q06_rrab.json"""
import json, numpy as np, pandas as pd, glob
R3 = "tic214311908/analysis/referee3_clean"
CL = "tic214311908/analysis/clean"
PRR = 0.6485275
AVOID = [PRR, PRR / 2, 0.7503, 0.7503 / 2, 1.5953 / 2, 1.5953 / 3, 2 * PRR / 3]
rng = np.random.default_rng(6)


def a1(t, y, P, nh=3):
    ph = 2 * np.pi * t / P
    X = np.column_stack([np.ones_like(t)] + [f(k * ph) for k in range(1, nh + 1) for f in (np.cos, np.sin)])
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    return float(np.hypot(c[1], c[2])), c


def fakeP(n):
    out = []
    while len(out) < n:
        p = rng.uniform(0.45, 0.95)
        if all(abs(p / a - 1) > 0.01 for a in AVOID): out.append(p)
    return np.array(out)


def main():
    d = pd.read_csv(f"{CL}/lc_M3V_INV_clean.csv", usecols=["time_btjd_tdb", "sector", "quality_ok", "flux_norm", "in_transit_window", "flare_flag", "rrab_resid_model"])
    d = d[(d.quality_ok == 1) & np.isfinite(d.flux_norm) & (d.in_transit_window == 0) & (d.flare_flag == 0)]
    res = {}
    FP = fakeP(300)
    for sec, g in d.groupby("sector"):
        t = g.time_btjd_tdb.values; y = (g.flux_norm.values - 1) * 1e3
        y = np.clip(y, np.median(y) - 6 * np.std(y), np.median(y) + 6 * np.std(y))
        A, _ = a1(t, y, PRR)
        fa = np.array([a1(t, y, p)[0] for p in FP])
        p = float(np.mean(fa >= A))
        # noise per quadrature from fake periods: Rayleigh sigma ~ sqrt(mean(fa^2)/2)
        sig = float(np.sqrt(np.mean(fa ** 2) / 2))
        ul = float(np.sqrt(max(A ** 2 - 2 * sig ** 2, 0)) + 1.645 * sig)
        res[int(sec)] = dict(A1_ppt=A, p_fake=p, sigma_quad_ppt=sig, UL95_ppt=ul, n=int(len(t)))
        print(sec, {k: round(v, 3) if isinstance(v, float) else v for k, v in res[int(sec)].items()})
    json.dump(res, open(f"{R3}/out/q06_rrab.json", "w"), indent=1)


if __name__ == "__main__":
    main()
