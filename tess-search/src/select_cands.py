"""Stage-1 candidate selection for one sector.

1. Find sector-wide systematic times: moments where far more stars than usual have a 'transit'.
2. Re-score each candidate with those windows masked.
3. Drop signals already in the TOI / CTOI / confirmed-planet / SPOC-TCE lists.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys, glob, warnings
import numpy as np, pandas as pd

warnings.filterwarnings("ignore")
from bls import detrend, clean, diagnostics
import known


def bad_windows(res, snr_min=8.0, binw=0.05, nsig=5.0):
    x = res[(res.snr_red >= snr_min) & res.passed]
    t = []
    for r in x.itertuples():
        k = np.arange(-60, 61); tt = r.t0 + k * r.P
        t.append(tt)
    t = np.concatenate(t)
    lo, hi = np.nanpercentile(res.t0, 0.1) - 30, np.nanpercentile(res.t0, 99.9) + 30
    t = t[(t > lo) & (t < hi)]
    edges = np.arange(t.min(), t.max() + binw, binw)
    h, e = np.histogram(t, edges)
    med = np.median(h[h > 0]); mad = 1.4826 * np.median(np.abs(h[h > 0] - med))
    bad = h > med + nsig * max(mad, 1.0)
    wins = [(e[i] - binw, e[i + 1] + binw) for i in np.where(bad)[0]]
    return wins, h, e


def rescore(cands, datadir, wins, camwins=None):
    out = []
    by_chunk = {}
    for f in sorted(glob.glob(f"{datadir}/chunk_*.npz")):
        z = np.load(f); by_chunk[f] = (set(z["tic"].tolist()), z)
    for r in cands.itertuples():
        for f, (ids, z) in by_chunk.items():
            if r.tic in ids: break
        else:
            continue
        i = int(np.where(z["tic"] == r.tic)[0][0]); off = z["off"]
        t = z["time"][off[i]:off[i + 1]].astype(float); fl = z["flux"][off[i]:off[i + 1]].astype(float)
        m = np.ones(len(t), bool)
        extra = (camwins or {}).get(int(z["camccd"][i][0]), [])
        for a, b in list(wins) + list(extra): m &= ~((t > a) & (t < b))
        t, fl = detrend(t[m], fl[m]); t, fl, s = clean(t, fl)
        d = diagnostics(t, fl, r.P, r.t0, r.dur, s)
        out.append(dict(tic=r.tic, **{f"m_{k}": v for k, v in d.items()}))
    return pd.DataFrame(out)


if __name__ == "__main__":
    sector, resdir, datadir, out = int(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4]
    import json
    camwins = {int(k): v["windows"] for k, v in json.load(open(sys.argv[5])).items()} if len(sys.argv) > 5 else {}
    res = pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(f"{resdir}/chunk_*.parquet"))], ignore_index=True)
    res = res[res.get("error", pd.Series(index=res.index, dtype=object)).isna()] if "error" in res else res
    wins, h, e = bad_windows(res)
    print(f"{len(res)} stars searched; {len(wins)} systematic windows: {[(round(a, 2), round(b, 2)) for a, b in wins][:30]}")
    c = res[(res.snr_red >= 8.5) & res.passed].copy()
    print("candidates before masking:", len(c))
    m = rescore(c, datadir, wins, camwins)
    c = c.merge(m, on="tic")
    keep = (c.m_snr_red >= 9) & (c.m_n_tr >= 3) & (c.m_max_event_frac <= 0.5) & (c.m_depth > 0)
    c = c[keep].copy()
    print("after masking (snr_red>=9, >=3 transits, no dominant event):", len(c))
    ks = c.apply(lambda r: known.known_signals(int(r.tic), r.P), axis=1)
    c["known_on_star"] = ks.map(lambda x: ";".join(x[0][:6])); c["known_match"] = ks.map(lambda x: ";".join(x[1][:6]))
    c["is_known"] = c.known_match.str.len() > 0
    c["sector"] = sector
    print("already known (period match):", int(c.is_known.sum()), " new:", int((~c.is_known).sum()))
    c.sort_values("m_snr_red", ascending=False).to_csv(out, index=False)
