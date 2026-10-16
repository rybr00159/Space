"""Run the exact stage-1 search used on Sectors 104/105 on TIC 4206066's own QLP light curves."""
import sys, json, warnings
import os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tess-search", "src"))  # the search pipeline, unchanged
warnings.filterwarnings("ignore")
import numpy as np
import fetch_qlp
from bls import detrend, search_robust, clean, diagnostics, _bls, freq_grid, DUR_H, _local_peaks
from classify import tmax_hours

TIC = 4206066
EPH = json.load(open(os.path.join(HERE, "..", "sources", "ephemeris_final.json")))
A, B = EPH["A"], EPH["B"]

def stage1(t, f, label):
    t, f = detrend(t, f)
    best, (tt, ff) = search_robust(t, f, max_mask=2)
    pre = best["snr_red"] >= 8.5 and best["passed"]
    post = best["snr_red"] >= 9 and best["n_tr"] >= 3 and best["max_event_frac"] <= 0.5 and best["depth"] > 0
    print(f"\n== {label}: {len(t)} points, baseline {best['baseline']:.1f} d, sigma {best['sigma']*1e6:.0f} ppm per bin")
    print(f"   BLS winner: P={best['P']:.5f} d  depth={best['depth']*1e6:.0f} ppm  dur={best['dur']*24:.2f} h  "
          f"n_tr={best['n_tr']}  snr_w={best['snr_w']:.1f}  snr_red={best['snr_red']:.2f}  SDE={best['sde']:.1f}  "
          f"masked={best['n_masked']}  alias_of={best.get('alias_of')}")
    print(f"   stage-1 pass (snr_red>=8.5 & passed): {pre};  final cut (snr_red>=9, >=3 transits, no dominant event): {pre and post}")
    # what our metric gives at the published ephemerides
    tc, fc, s = clean(t, f)
    for name, e in (("A 3.18 d", A), ("B 11.13 d", B)):
        d = diagnostics(tc, fc, e["P"], e["T0"], e["T14"], s)
        print(f"   at published {name}: depth={d['depth']*1e6:.0f}±{d['depth_err']*1e6:.0f} ppm  n_tr={d['n_tr']}  "
              f"snr_red={d['snr_red']:.2f}  red_ratio={d['red_ratio']:.2f}  oddeven={d['oddeven_sig']:.1f}  "
              f"max_event_frac={d['max_event_frac']:.2f}")
    # rank of the 3.18 d peak in the coarse BLS spectrum
    y = fc - 1.0; w = np.full(len(tc), 1 / s ** 2); y -= np.sum(w * y) / np.sum(w)
    freqs = freq_grid(tc.max() - tc.min())
    pw, dp, t0, du = _bls(tc - tc.min(), y, w, freqs, DUR_H / 24, 10 / 1440)
    peaks = _local_peaks(pw, freqs, 15)
    Ps = [1 / freqs[i] for i in peaks]
    hit = [k for k, P in enumerate(Ps) if any(abs(P / (A["P"] * m) - 1) < 0.01 for m in (0.5, 1, 2))]
    print("   top BLS peaks (d):", ", ".join(f"{P:.3f}" for P in Ps[:10]))
    print("   rank of 3.18 d (or 2x/0.5x) among top 15 peaks:", [h + 1 for h in hit] or "not in top 15")
    return best

for sector, binsize in ((98, 3),):  # S6/S32 QLP files have no DET_FLUX; see would_detect2.py
    r = fetch_qlp.load(TIC, sector, binsize=binsize)
    if r is None:
        print(f"\nS{sector}: no QLP light curve"); continue
    t, f, n, meta = r
    stage1(t, f, f"Sector {sector} QLP DET_FLUX (camera/CCD {meta})")

print("\nduration check: predicted max duration for a central transit =",
      round(float(tmax_hours(np.array([A['P']]), np.array([0.648]), np.array([0.63]))[0]), 2), "h; published 1.87 h")
