"""Referee r02: independent check of the target-specific barycentric correction of TESScut times, using astropy's
Earth barycentric ephemeris (TESS is within ~1.3 light-s of Earth, so the Earth-based Roemer delay is good to ~1.5 s)
instead of the JPL Horizons TESS vectors used by the investigator. Also compares TESScut and QLP time stamps
cadence-by-cadence (QLP TIME vs TESScut TIME at matching cadences).
Output: out/timecheck.json, data/tcorr_s{sector}.npy (TESScut TIME, seconds to add)."""
import sys; sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
from rlib import *
from astropy.time import Time
from astropy.coordinates import get_body_barycentric, SkyCoord
import astropy.units as u
RA, DE = 260.970485, -51.324202
n = SkyCoord(RA, DE, unit="deg").cartesian.xyz.value
LC = "/home/user/work/lc_cache"
inv = json.load(open(f"{R}/../phot/out/timecorr.json"))
out = {}
for s in SECTORS:
    C = Cut(s); m = np.isfinite(C.t) & np.isfinite(C.tcorr)
    tsc = C.t[m] - C.tcorr[m]
    tt = Time(tsc + 2457000, format="jd", scale="tdb")
    pe = get_body_barycentric("earth", tt).xyz.to(u.m).value.T
    roem = pe @ n / 299792458.0
    add = roem - C.tcorr[m] * 86400
    np.save(f"{R}/data/tcorr_s{s}.npy", np.vstack([C.t[m], add]))
    # QLP comparison at matching cadences (QLP files: TIME column; S103/S104 known to lack the barycentric correction)
    f = f"{LC}/hlsp_qlp_tess_ffi_s{s:04d}-0000000214311908_tess_v01_llc.fits"
    with fits.open(f) as h:
        tq = np.array(h[1].data["TIME"], float); cq = np.array(h[1].data["CADENCENO"], int) if "CADENCENO" in h[1].data.names else None
        hq = h[1].header
    tq = tq[np.isfinite(tq)]
    tcut = C.t[m] + add / 86400          # referee target-barycentred TESScut times
    j = np.clip(np.searchsorted(tcut, tq), 1, len(tcut) - 1)
    near = np.where(np.abs(tq - tcut[j - 1]) < np.abs(tq - tcut[j]), j - 1, j)
    dd = (tq - tcut[near]) * 86400
    ok = np.abs(dd) < 0.4 * C.dt * 86400
    out[s] = dict(add_mean_s=float(add.mean()), add_min_s=float(add.min()), add_max_s=float(add.max()),
                  investigator_add_mean_s=inv[str(s)]["add_s_mean"], diff_vs_investigator_s=float(add.mean() - inv[str(s)]["add_s_mean"]),
                  timecorr_mean_min=float(C.tcorr[m].mean() * 1440),
                  qlp_minus_referee_median_s=float(np.median(dd[ok])) if ok.sum() else np.nan, qlp_minus_referee_rstd_s=float(rstd(dd[ok])) if ok.sum() > 3 else np.nan,
                  n_match=int(ok.sum()), n_qlp=len(tq))
    print(f"S{s}: add to TESScut TIME {add.mean():+.1f} s (range {add.min():+.1f}..{add.max():+.1f}); investigator {inv[str(s)]['add_s_mean']:+.1f} s; "
          f"QLP - (TESScut+corr) at matching cadences: median {out[s]['qlp_minus_referee_median_s']:+.1f} s (rstd {out[s]['qlp_minus_referee_rstd_s']:.2f} s, n {ok.sum()}/{len(tq)})", flush=True)
save_json({str(k): v for k, v in out.items()}, f"{R}/out/timecheck.json")
