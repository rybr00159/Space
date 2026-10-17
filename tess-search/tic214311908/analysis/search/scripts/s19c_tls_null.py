"""s19c: empirical null for the S103+S104 TLS search (s19): 12 realisations with the cleaned flux cyclically rolled
within each sector (non-inverted, keeps dips but destroys any 16.79-d phasing across the two sectors / within) and
12 inverted + rolled. Same TLS settings as s19. Output: out/tls_null_S103_104.json (saved after each run)."""
import sys; sys.path.insert(0, "tic214311908/analysis/search/scripts")
from slib import *
from transitleastsquares import transitleastsquares
rng = np.random.default_rng(193)
z = np.load(f"{ROOT}/data/m3v_clean.npz")
t, f, s, cad, intr = z["t"], z["det_clean"], z["sector"], z["cad"].astype(float), z["intr"]
m = ~intr & np.isfinite(f)
tb, fb, nb_, sb, cb = bin_lc(t[m], f[m], s[m], cad[m], 10.0)
k = np.isin(sb, [103, 104]); tt, ff, ss = tb[k], fb[k] / np.median(fb[k]), sb[k]
err = np.full(k.sum(), robust_std(ff)); span = tt.max() - tt.min()
out = []
for i in range(24):
    y = ff.copy()
    for sec in (103, 104):
        j = np.where(ss == sec)[0]; y[j] = np.roll(ff[j], rng.integers(len(j) // 10, len(j) - len(j) // 10))
    sign = 1 if i % 2 == 0 else -1
    y = 1 + sign * (y - 1)
    r = transitleastsquares(tt, y, err).power(R_star=RSTAR, R_star_min=0.3, R_star_max=0.5, M_star=MSTAR, M_star_min=0.3, M_star_max=0.45,
        u=[0.25, 0.35], period_min=0.3, period_max=span / 2, oversampling_factor=3, use_threads=1, show_progress_bar=False)
    out.append(dict(sign=sign, SDE=float(r.SDE), P=float(r.period), n_transits=int(r.distinct_transit_count)))
    print(i, out[-1], flush=True)
    save_json(out, f"{ROOT}/out/tls_null_S103_104.json")
