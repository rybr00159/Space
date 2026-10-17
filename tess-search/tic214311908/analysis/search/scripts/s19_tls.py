"""s19: Transit Least Squares (Hippke & Heller 2019) per sector on the cleaned M3V light curve (5.69-d transits
masked, 10-min bins; S12 30 min), real and inverted (empirical null), plus S103+S104 together.
Stellar parameters R = 0.381, M = 0.370 (single-star case), quadratic limb darkening u = (0.25, 0.35);
period 0.3 d to half the span; oversampling 3; one thread.
TLS' own SDE-FAP mapping is for white noise; the inverted runs give the empirical comparison (a real dip signal
should exceed the inverted SDE of the same sector).
Output: out/tls.json
"""
import sys; sys.path.insert(0, "tic214311908/analysis/search/scripts")
from slib import *
from transitleastsquares import transitleastsquares


def main():
    z = np.load(f"{ROOT}/data/m3v_clean.npz")
    t, f, s, cad, intr = z["t"], z["det_clean"], z["sector"], z["cad"].astype(float), z["intr"]
    m = ~intr & np.isfinite(f)
    tb, fb, nb_, sb, cb = bin_lc(t[m], f[m], s[m], cad[m], 10.0)
    out = {}
    groups = {"S12": [12], "S39": [39], "S66": [66], "S93": [93], "S103": [103], "S104": [104], "S103_104": [103, 104]}
    for g, ss in groups.items():
        k = np.isin(sb, ss)
        tt, ff = tb[k], fb[k] / np.median(fb[k])
        err = np.full(k.sum(), robust_std(ff))
        span = tt.max() - tt.min()
        for sign in (1, -1):
            yy = 1 + sign * (ff - 1)
            model = transitleastsquares(tt, yy, err)
            r = model.power(R_star=RSTAR, R_star_min=0.3, R_star_max=0.5, M_star=MSTAR, M_star_min=0.3, M_star_max=0.45,
                            u=[0.25, 0.35], period_min=0.3, period_max=span / 2, oversampling_factor=3, use_threads=1,
                            show_progress_bar=False)
            key = g + ("" if sign == 1 else "_inv")
            out[key] = dict(SDE=float(r.SDE), P=float(r.period), T0=float(r.T0), depth_ppt=float((1 - r.depth) * 1e3),
                            duration_h=float(r.duration * 24), snr=float(r.snr), n_transits=int(r.distinct_transit_count),
                            FAP_tls=float(r.FAP) if r.FAP is not None else None, n_periods=len(r.periods),
                            odd_even_mismatch=float(r.odd_even_mismatch) if np.isfinite(r.odd_even_mismatch) else None)
            print(key, json.dumps(out[key]), flush=True)
            save_json(out, f"{ROOT}/out/tls.json")


if __name__ == "__main__":
    main()
