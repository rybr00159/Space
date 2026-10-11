"""Injection-recovery test: add limb-darkened transits to real light curves and rerun the search."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("NUMBA_NUM_THREADS", "1")
import sys, warnings
from multiprocessing import Pool
import numpy as np, pandas as pd, batman

warnings.filterwarnings("ignore")
from bls import search_robust, detrend
from run_search import stars


def transit_model(t, P, t0, rp, a, inc):
    p = batman.TransitParams()
    p.t0, p.per, p.rp, p.a, p.inc, p.ecc, p.w = t0, P, rp, a, inc, 0.0, 90.0
    p.u, p.limb_dark = [0.4, 0.25], "quadratic"
    return batman.TransitModel(p, t, supersample_factor=5, exp_time=600 / 86400).light_curve(p)


def one(args):
    tic, t, f, seed = args
    rng = np.random.default_rng(seed)
    P = 10 ** rng.uniform(np.log10(0.7), np.log10(12.0))
    depth = 10 ** rng.uniform(np.log10(150e-6), np.log10(5000e-6))
    rp = np.sqrt(depth); a = 4.2 * P ** (2 / 3)  # Sun-like star: a/R* ~ 4.2 P^(2/3)
    b = rng.uniform(0, 0.8); inc = np.degrees(np.arccos(b / a))
    t0 = t.min() + rng.uniform(0, P)
    f = f * transit_model(t.astype(float), P, t0, rp, a, inc)
    dur = P / np.pi / a  # approximate full duration
    ep = np.round((t - t0) / P); intr = np.abs(t - t0 - ep * P) < dur / 2
    cnt = np.unique(ep[intr], return_counts=True)[1]
    n_inj = int((cnt >= 0.5 * dur / (600 / 86400)).sum())
    t, f = detrend(t.astype(float), f.astype(float))
    try:
        r, _ = search_robust(t, f)
    except Exception:
        return None
    ratio = r["P"] / P
    hit = any(abs(ratio - k) < 0.01 * k for k in (1, 0.5, 2)) and abs(((r["t0"] - t0) / P + 0.5) % 1 - 0.5) < 0.05 * max(1, ratio)
    exact = abs(ratio - 1) < 0.01
    return dict(tic=int(tic), P=P, depth=depth, b=b, P_found=r["P"], snr_red=r["snr_red"], n_tr=r["n_tr"],
                passed=r["passed"], n_inj=n_inj, hit=hit, exact=exact, sigma=r["sigma"])


if __name__ == "__main__":
    npz, n, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    jobs = [(tic, t, f, i) for i, (tic, t, f) in enumerate(stars(npz)) if i < n]
    import multiprocessing as mp
    with mp.get_context("spawn").Pool(4) as pool:
        res = [r for r in pool.imap(one, jobs, chunksize=10) if r]
    df = pd.DataFrame(res); df.to_parquet(out)
    df["snr_exp"] = df.depth / df.sigma * np.sqrt(df.n_tr.clip(1) * 6)  # rough expectation
    df["dbin"] = pd.cut(df.depth * 1e6, [150, 300, 500, 800, 1200, 2000, 5000])
    df = df[df.n_inj >= 2]
    print("injections with >=2 transits in the data:", len(df))
    print(df.groupby("dbin", observed=True).agg(n=("hit", "size"), recovered=("hit", "mean"), exact=("exact", "mean")).round(2))
