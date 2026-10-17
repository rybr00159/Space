"""s19b: TLS sanity check - the same per-sector TLS set-up as s19 but WITHOUT masking the 5.69-d transits
(S39, S93, S103+S104): the known signal must be recovered. Output: out/tls_sanity.json"""
import sys; sys.path.insert(0, "tic214311908/analysis/search/scripts")
from slib import *
from transitleastsquares import transitleastsquares
z = np.load(f"{ROOT}/data/m3v_clean.npz")
t, f, s, cad = z["t"], z["det_clean"], z["sector"], z["cad"].astype(float)
m = np.isfinite(f)
tb, fb, nb_, sb, cb = bin_lc(t[m], f[m], s[m], cad[m], 10.0)
out = {}
for g, ss in {"S39": [39], "S93": [93], "S103_104": [103, 104]}.items():
    k = np.isin(sb, ss); tt, ff = tb[k], fb[k] / np.median(fb[k]); err = np.full(k.sum(), robust_std(ff))
    r = transitleastsquares(tt, ff, err).power(R_star=RSTAR, R_star_min=0.3, R_star_max=0.5, M_star=MSTAR, M_star_min=0.3, M_star_max=0.45,
        u=[0.25, 0.35], period_min=0.3, period_max=(tt.max() - tt.min()) / 2, oversampling_factor=3, use_threads=1, show_progress_bar=False)
    out[g] = dict(SDE=float(r.SDE), P=float(r.period), depth_ppt=float((1 - r.depth) * 1e3), duration_h=float(r.duration * 24), snr=float(r.snr))
    print(g, out[g], flush=True)
save_json(out, f"{ROOT}/out/tls_sanity.json")
