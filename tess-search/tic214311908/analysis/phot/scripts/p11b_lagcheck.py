"""Sensitivity check of the QLP-TESScut cross-correlation lag (p11): repeat for S93, S103 and S104 WITHOUT the QLP
time correction and with a lag grid of +-5 cadences; expect about -5.2 (S103) and -7.5 min (S104), ~0 for S93."""
import sys; sys.path.insert(0, "tic214311908/analysis/phot/scripts")
from photlib import *
from astropy.io import fits
from scipy.ndimage import median_filter
TC = json.load(open(f"{ROOT}/out/timecorr.json")); res = {}
for s in (93, 103, 104):
    C = Cutout(s); W = json.load(open(f"{ROOT}/out/wcs_s{s}.json"))
    g = C.gaia(15); r = g[g.source_id == IDS["M3V"]].iloc[0]
    yy, xx = np.mgrid[0:C.ny, 0:C.nx]
    ap = np.hypot(xx - r.x - W["global_shift"][0], yy - r.y - W["global_shift"][1]) <= 1.5
    ok = C.fin & (C.q == 0)
    med = np.nanmedian(C.flux[ok], 0); faint = med < np.percentile(med, 15)
    bk = np.nanmedian(C.flux[:, faint], 1); fa = C.flux[:, ap].sum(1) - bk * ap.sum()
    ta = C.t     # TESScut times as delivered (FFI-direction barycentred)
    with fits.open(f"{LCDIR}/hlsp_qlp_tess_ffi_s{s:04d}-0000000214311908_tess_v01_llc.fits") as h:
        d = h[1].data; tq = np.array(d["TIME"], float); fq = np.array(d["SAP_FLUX"], float); qq = np.array(d["QUALITY"])
    mq = (qq == 0) & np.isfinite(tq) & np.isfinite(fq); m2 = ok & np.isfinite(fa)
    def hp(t, f):
        n = max(3, int(round(0.3 / np.median(np.diff(t)))) | 1); return f - median_filter(f, size=n, mode="nearest")
    ha = hp(ta[m2], fa[m2] / np.nanmedian(fa[m2])); hq = hp(tq[mq], fq[mq] / np.nanmedian(fq[mq]))
    ha = np.clip(ha, -5 * robust_std(ha), 5 * robust_std(ha)); hq = np.clip(hq, -5 * robust_std(hq), 5 * robust_std(hq))
    lags = np.arange(-5, 5.0001, 0.05) * C.exptime; cc = []
    for L in lags:
        v = np.interp(ta[m2], tq[mq] - L, hq, left=np.nan, right=np.nan); k = np.isfinite(v)
        cc.append(np.corrcoef(ha[k], v[k])[0, 1])
    cc = np.array(cc); i = int(np.nanargmax(cc)); a_, b_, c_ = np.polyfit(lags[i - 1:i + 2], cc[i - 1:i + 2], 2); best = -b_ / (2 * a_)
    res[s] = dict(lag_uncorrected_min=float(best * 1440), timecorr_median_min=float(np.nanmedian(C.timecorr) * 1440))
    print(f"S{s}: uncorrected QLP - TESScut lag {best*1440:+.2f} min (TESScut TIMECORR median {np.nanmedian(C.timecorr)*1440:.2f} min)")
save_json(res, f"{ROOT}/out/public_timing_lag_uncorrected.json")
