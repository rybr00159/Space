"""Search stars observed in two sectors (S104 + S105) on their combined light curves."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("NUMBA_NUM_THREADS", "1")
import sys, glob, json, time, warnings
import numpy as np, pandas as pd

warnings.filterwarnings("ignore")
from bls import search_robust, detrend
from run_search import KEYS


def index(datadir):
    idx = {}
    for f in sorted(glob.glob(f"{datadir}/chunk_*.npz")):
        z = np.load(f)
        for i, t in enumerate(z["tic"]): idx[int(t)] = (f, i)
    return idx


def get(z, i, wins):
    off = z["off"]; t = z["time"][off[i]:off[i + 1]].astype(float); f = z["flux"][off[i]:off[i + 1]].astype(float)
    m = np.ones(len(t), bool)
    for a, b in wins.get(int(z["camccd"][i][0]), []): m &= ~((t > a) & (t < b))
    return detrend(t[m], f[m])


def one(args):
    tic, t, f = args
    try:
        b, _ = search_robust(t, f, max_mask=2, pmax=min(25.0, (t.max() - t.min()) / 2))
        r = {k: b.get(k, np.nan) for k in KEYS}; r["tic"] = tic
        return r
    except Exception as e:
        return {"tic": tic, "error": str(e)[:200]}


if __name__ == "__main__":
    dA, dB, wA, wB, out = sys.argv[1:6]
    WA = {int(k): v["windows"] for k, v in json.load(open(wA)).items()}
    WB = {int(k): v["windows"] for k, v in json.load(open(wB)).items()}
    iA, iB = index(dA), index(dB)
    both = sorted(set(iA) & set(iB))
    print(len(both), "stars in both sectors", flush=True)
    os.makedirs(out, exist_ok=True)
    import multiprocessing as mp
    cache = {}
    def z_of(path):
        if path not in cache:
            if len(cache) > 6: cache.pop(next(iter(cache)))
            cache[path] = {k: v for k, v in np.load(path).items()}
        return cache[path]
    # process in groups ordered by sector-A chunk to keep file loads low
    both.sort(key=lambda t: (iA[t][0], iB[t][0]))
    with mp.get_context("spawn").Pool(4) as pool:
        for g0 in range(0, len(both), 4000):
            fn = f"{out}/group_{g0 // 4000:04d}.parquet"
            if os.path.exists(fn): continue
            t0 = time.time(); jobs = []
            for tic in both[g0:g0 + 4000]:
                ta, fa = get(z_of(iA[tic][0]), iA[tic][1], WA); tb, fb = get(z_of(iB[tic][0]), iB[tic][1], WB)
                jobs.append((tic, np.concatenate([ta, tb]), np.concatenate([fa, fb])))
            res = [r for r in pool.imap(one, jobs, chunksize=10) if r]
            pd.DataFrame(res).to_parquet(fn)
            print(f"{fn}: {len(res)} in {time.time() - t0:.0f}s", flush=True)
