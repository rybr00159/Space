"""Red-noise factor beta of the PSF pair light curve (out-of-transit and out-of-occultation, biweight-detrended 'f' as
used by the transit/twin fits): rms of bin means / (white rms / sqrt(n)) for 0.5, 1, 2 and 4.5-h bins, per sector and
pooled (n-weighted mean of beta^2). Output: results/beta_pair.json"""
import os, json
import numpy as np
ROOT = "/home/user/work/tic87/neighbour"
lc = np.load(os.path.join(ROOT, "data/lc_pair.npy"), allow_pickle=True).item()
P, T0 = 4.1361919, 3271.53062
out = {}
for s, v in lc.items():
    t, f = v["t"], v["f"]
    ph = ((t - T0) / P) % 1
    m = (np.minimum(ph, 1 - ph) * P > 0.2) & (np.abs(ph - 0.5) * P > 0.2)
    t, r = t[m], f[m] - np.median(f[m])
    sw = np.std(r); dt = np.median(np.diff(t))
    out[s] = {}
    for hb in (0.5, 1.0, 2.0, 4.5):
        nb = max(1, int(round(hb / 24 / dt)))
        k = np.floor((t - t[0]) / (hb / 24)).astype(int)
        u, inv, cnt = np.unique(k, return_inverse=True, return_counts=True)
        means = np.bincount(inv, r) / cnt
        ok = cnt >= 0.7 * nb
        exp = sw / np.sqrt(cnt[ok])
        beta = float(np.sqrt(np.mean((means[ok] / exp) ** 2)))
        out[s][f"{hb}h"] = beta
    print(s, out[s])
pool = {hb: float(np.sqrt(np.mean([out[s][hb] ** 2 for s in out]))) for hb in out[13]}
out["pooled_rms_over_sectors"] = pool
print(pool)
json.dump(out, open(os.path.join(ROOT, "results/beta_pair.json"), "w"), indent=1, default=str)
