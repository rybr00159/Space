"""Referee r02b: QLP time stamps vs TESScut cadence times (nearest cadence within 0.4 cadence), for three references:
TESScut TIME (FFI barycentric), TESScut TIME - TIMECORR (spacecraft TDB), and TESScut TIME + target correction (r02).
Output: out/qlp_stamps.json"""
import sys; sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
from rlib import *
out = {}
for s in SECTORS:
    tq = np.array(fits.open(f"/home/user/work/lc_cache/hlsp_qlp_tess_ffi_s{s:04d}-0000000214311908_tess_v01_llc.fits")[1].data["TIME"], float)
    tq = tq[np.isfinite(tq)]
    C = Cut(s); m = np.isfinite(C.t) & np.isfinite(C.tcorr)
    z = np.load(f"{R}/data/tcorr_s{s}.npy")
    refs = dict(ffi_bary=C.t[m], spacecraft=C.t[m] - C.tcorr[m], target_bary=C.t[m] + np.interp(C.t[m], z[0], z[1]) / 86400)
    row = {}
    for k, tr in refs.items():
        j = np.clip(np.searchsorted(tr, tq), 1, len(tr) - 1)
        dd = np.where(np.abs(tq - tr[j - 1]) < np.abs(tq - tr[j]), tq - tr[j - 1], tq - tr[j]) * 86400
        ok = np.abs(dd) < 0.4 * C.dt * 86400
        row[k] = dict(frac_matched=float(ok.mean()), median_s=float(np.median(dd[ok])) if ok.sum() else np.nan,
                      rstd_s=float(rstd(dd[ok])) if ok.sum() > 3 else np.nan)
    out[s] = row
    print(f"S{s}: " + "; ".join(f"QLP - {k}: matched {v['frac_matched']:.2f}, median {v['median_s']:+.1f} s" for k, v in row.items()))
save_json({str(k): v for k, v in out.items()}, f"{R}/out/qlp_stamps.json")
