"""s18: semi-amplitude of the TARS ~8.7-d modulation and of the 5.69-d transit in the same TARS light curves (fixed
trapezoid, T14 1.6 h, flat 70 %, local linear baselines), and the implied M3V semi-amplitude if the modulation were on
the M3V: A_M3V = (A_TARS / depth_TARS) x 15.4 e-/s (both diluted identically by the 1-px aperture, assuming TARS'
processing preserves an 8.7-d signal as it preserves the transit; a lower bound if TARS' detrending suppresses it).
Errors: sinusoid amplitudes at 200 fake periods 3-14 d; transit depth from 200 fake epochs.
Output: out/tars_amp.json
"""
import sys; sys.path.insert(0, "tic214311908/analysis/search/scripts")
from slib import *
from astropy.io import fits
from s12_alias import trap
from s14_rotation import sin_vec
RNG = np.random.default_rng(18)
pp = json.load(open(f"{ROOT}/out/public_periodograms.json"))
out = {}
for s in (12, 39, 66, 93):
    f = f"{BASE}/ctx/lc/hlsp_tars_tess_ffi_s{s:04d}-0000000214311908_tess_v01_lc.fits"
    d = fits.open(f)[1].data; t = np.array(d["TIME"], float); y = np.array(d["FLUX"], float)
    m = np.isfinite(t) & np.isfinite(y); t, y = t[m], y[m] / np.nanmedian(y[m]) - 1
    cad = np.full(len(t), np.median(np.diff(t)) * 86400)
    intr = transit_mask(t)
    Pt = pp[f"TARS_M3V_S{s}"]["top"][0][0]
    a = float(np.hypot(*sin_vec(t[~intr], y[~intr], Pt)))
    fk = [np.hypot(*sin_vec(t[~intr], y[~intr], pf)) for pf in RNG.uniform(3, 14, 300) if abs(pf - Pt) / Pt > 0.1][:200]
    # transit depth with local linear baselines (+-0.3 d), stacked
    def depth(T0):
        num = den = 0.0
        for e in range(-500, 100):
            tc = T0 + e * P_B; k = np.abs(t - tc) < 0.3
            if k.sum() < 10: continue
            tm = trap(t[k], P_B, tc, cad[k])
            if (tm > 0.5).sum() < 2: continue
            X = np.column_stack([np.ones(k.sum()), t[k] - tc, -tm]); c, *_ = np.linalg.lstsq(X, y[k], rcond=None)
            w = (tm > 0.5).sum(); num += c[2] * w; den += w
        return num / den if den else np.nan
    dep = depth(T0_B); fd = [depth(T0_B + RNG.uniform(0.1, 0.9) * P_B) for _ in range(150)]
    out[s] = dict(P_tars=Pt, amp_ppt=a * 1e3, amp_fake_median_ppt=float(np.median(fk) * 1e3), amp_p_emp=float(np.mean(np.array(fk) >= a)),
                  depth_ppt=dep * 1e3, depth_err_emp_ppt=float(np.nanstd(fd) * 1e3),
                  implied_M3V_amp_e_s=float(a / dep * 15.4), implied_M3V_amp_pct=float(a / dep * 15.4 / 979 * 100))
    print(s, {k: round(v, 3) for k, v in out[s].items()})
save_json({str(k): v for k, v in out.items()}, f"{ROOT}/out/tars_amp.json")
