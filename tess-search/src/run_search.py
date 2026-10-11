"""Run the robust BLS search on every star in downloaded chunk files (multiprocessing)."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("NUMBA_NUM_THREADS", "1")
import sys, os, glob, time, warnings
from multiprocessing import Pool
import numpy as np, pandas as pd

warnings.filterwarnings("ignore")
from bls import search_robust, detrend

KEYS = ["P", "t0", "dur", "snr_w", "sde", "n_tr", "n_tr_any", "max_event_frac", "depth", "depth_err", "sigma_dur",
        "snr_red", "red_ratio", "oddeven_sig", "sec_depth", "sec_sig", "shape_ratio", "sigma", "npts", "baseline",
        "n_masked", "passed", "rank", "alias_of"]


def one(args):
    tic, t, f, wins = args
    try:
        t = t.astype(float); f = f.astype(float)
        if wins:
            m = np.ones(len(t), bool)
            for a, b in wins: m &= ~((t > a) & (t < b))
            t, f = t[m], f[m]
        t, f = detrend(t, f)
        if len(t) < 500: return None
        b, _ = search_robust(t, f, max_mask=2)
        r = {k: b.get(k, np.nan) for k in KEYS}; r["tic"] = int(tic)
        return r
    except Exception as e:
        return {"tic": int(tic), "error": str(e)[:200]}


def stars(npz):
    z = np.load(npz)
    off, tics, T, F = z["off"], z["tic"], z["time"], z["flux"]  # read each array once
    for i, tic in enumerate(tics):
        yield tic, T[off[i]:off[i + 1]], F[off[i]:off[i + 1]]


def main(datadir, outdir, procs=4, badwin=None):
    import json
    os.makedirs(outdir, exist_ok=True)
    bw = {int(k): v["windows"] for k, v in json.load(open(badwin)).items()} if badwin else {}
    done = set()
    import multiprocessing as mp
    with mp.get_context("spawn").Pool(procs) as pool:
        while True:
            files = sorted(glob.glob(f"{datadir}/chunk_*.npz"))
            todo = [f for f in files if f not in done and not os.path.exists(f"{outdir}/{os.path.basename(f)[:-4]}.parquet")]
            if not todo:
                # keep polling while the downloader is still running
                break
            f = todo[0]; t0 = time.time()
            z = np.load(f); camccd = dict(zip(z["tic"], map(tuple, z["camccd"])))
            excl_file = f"{datadir.rstrip('/')}_exclude.txt"  # stars handled by the multi-sector search instead
            excl = set(np.loadtxt(excl_file, dtype=np.int64).tolist()) if os.path.exists(excl_file) else set()
            jobs = ((tic, t, fl, bw.get(int(camccd[tic][0]), [])) for tic, t, fl in stars(f) if int(tic) not in excl)
            res = [r for r in pool.imap(one, jobs, chunksize=20) if r is not None]
            df = pd.DataFrame(res)
            df["cam"] = df.tic.map(lambda x: camccd[x][0]); df["ccd"] = df.tic.map(lambda x: camccd[x][1])
            df.to_parquet(f"{outdir}/{os.path.basename(f)[:-4]}.parquet")
            done.add(f)
            print(f"{f}: {len(df)} stars in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4] if len(sys.argv) > 4 else None)
