"""Find sector-wide systematic time windows from the population of light curves.

For each camera, every star's detrended flux is normalised by its own scatter and binned in time; bins
where the median across stars is strongly negative, or where many stars dip at once, are flagged.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys, glob, json, warnings
import numpy as np

warnings.filterwarnings("ignore")
from bls import detrend, robust_std


def windows(datadir, nstar=2500, binw=0.02, thr_med=-0.6, thr_frac=0.12):
    per_cam = {}
    for f in sorted(glob.glob(f"{datadir}/chunk_*.npz")):
        z = np.load(f); off, T, F, cc = z["off"], z["time"], z["flux"], z["camccd"]
        for i in range(len(z["tic"])):
            cam = int(cc[i][0])
            if len(per_cam.setdefault(cam, [])) >= nstar: continue
            t, fl = detrend(T[off[i]:off[i + 1]].astype(float), F[off[i]:off[i + 1]].astype(float))
            s = robust_std(fl - 1)
            if s <= 0 or s > 0.01: continue
            per_cam[cam].append((t, (fl - 1) / s))
        if all(len(v) >= nstar for v in per_cam.values()) and len(per_cam) >= 4: break
    out = {}
    for cam, lst in per_cam.items():
        tmin = min(t.min() for t, _ in lst); tmax = max(t.max() for t, _ in lst)
        edges = np.arange(tmin, tmax + binw, binw); nb = len(edges) - 1
        sums = [[] for _ in range(nb)]
        M = np.full((len(lst), nb), np.nan)
        for j, (t, z) in enumerate(lst):
            k = np.clip(((t - tmin) / binw).astype(int), 0, nb - 1)
            c = np.bincount(k, minlength=nb); s = np.bincount(k, z, minlength=nb)
            m = c > 0; M[j, m] = s[m] / c[m]
        med = np.nanmedian(M, axis=0)
        frac = np.nanmean(M < -3.0, axis=0)  # fraction of stars with a >3 sigma dip in this bin
        nvalid = np.sum(np.isfinite(M), axis=0)
        bad = (nvalid > 50) & ((med < thr_med) | (frac > thr_frac))
        wins = []
        for k in np.where(bad)[0]:
            a, b = edges[k] - binw, edges[k + 1] + binw
            if wins and a <= wins[-1][1]: wins[-1][1] = b
            else: wins.append([a, b])
        out[cam] = dict(windows=wins, nstars=len(lst), median_frac=float(np.nanmedian(frac)))
    return out


if __name__ == "__main__":
    res = windows(sys.argv[1])
    json.dump({str(k): v for k, v in res.items()}, open(sys.argv[2], "w"), indent=1)
    for cam, v in sorted(res.items()):
        tot = sum(b - a for a, b in v["windows"])
        print(f"camera {cam}: {v['nstars']} stars, {len(v['windows'])} windows, {tot:.2f} d masked:",
              [(round(a, 2), round(b, 2)) for a, b in v["windows"]][:20])
