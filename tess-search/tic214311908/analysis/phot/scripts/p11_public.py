"""Cross-check with independent public light curves of TIC 214311908 (and the Mira / pair where available).

Products on MAST (queried 2026-10-16, data/mast/mast_timeseries.csv): QLP (S12, 39, 66, 93, 103, 104), TARS
(S12, 39, 66, 93; flux only), GSFC-ELEANOR-LITE (S12; also for the Mira and the pair). No TGLC, TESS-SPOC, SPOC
2-min, CDIPS or PATHOS light curves exist for these stars.
For every product: same detrending (biweight 0.75 d, transit windows masked), per-sector fractional depth at the
global ephemeris/shape of p08 (phase-scramble errors), and the depth converted to e-/s with the product's
normalisation where meaningful. QLP times in S103/S104 are corrected by the FFI TIMECORR median (as in
ephemeris_v2) plus our target-specific barycentric term.
Timing-offset check: cross-correlation of the high-pass filtered QLP SAP flux of TIC 214311908 with a TESScut
aperture sum (r = 1.5 px around the M3V, faint-pixel background removed) of the same pixels' systematics
(scattered light, jitter) gives the QLP - TESScut time offset per sector.
Outputs: out/public_depths.csv, out/public_timing_lag.json
"""
import sys; sys.path.insert(0, "tic214311908/analysis/phot/scripts")
from photlib import *
from tlib import ephem_frac, lin_depth
from p06_lc import detrend
from astropy.io import fits

DL = f"{ROOT}/data/mast/dl"
TR = json.load(open(f"{ROOT}/out/transits.json"))["global_B"]
P, T0, T14, FL = TR["P"], TR["T0"], TR["T14_d"], TR["flat_frac"]
EPH2 = json.load(open("/home/user/Space/tess-search/vetting_2026-10-16/results/214311906/ephemeris_v2.json"))
TC = json.load(open(f"{ROOT}/out/timecorr.json"))
RNG = np.random.default_rng(5)


def depth_scr(t, f, ex, n=150):
    y = ephem_frac(t, P, T0, T14, FL, ex)
    real = np.abs((((t - T0) / P + 0.5) % 1 - 0.5) * P) < 0.1
    sd = robust_std(f[~real]); e = np.full(len(t), sd)
    d, de, _, _ = lin_depth(f, e, y)
    nulls = []
    for k in range(n):
        y2 = ephem_frac(t, P, T0 + RNG.uniform(0.08, 0.92) * P, T14, FL, ex)
        if np.sum(y2[~real] > 0.5) < 3: continue
        nulls.append(lin_depth(f[~real], e[~real], y2[~real])[0])
    return d, de, float(np.std(nulls)), sd


def prep(t, f, ok):
    m = ok & np.isfinite(t) & np.isfinite(f) & (f > 0)
    t, f = t[m], f[m] / np.nanmedian(f[m])
    tr = detrend(t, f, np.ones(len(t), bool), np.abs((((t - T0) / P + 0.5) % 1 - 0.5) * P) < 0.1)
    d = f / tr; k = np.isfinite(d)
    return t[k], d[k], f[k]


rows = []
# ---- QLP for all three stars
for tic, star in ((214311908, "M3V"), (214311910, "pair"), (214311906, "Mira")):
    for s in SECTORS:
        fn = f"{LCDIR}/hlsp_qlp_tess_ffi_s{s:04d}-{tic:016d}_tess_v01_llc.fits"
        with fits.open(fn) as h:
            d = h[1].data; t = np.array(d["TIME"], float); q = np.array(d["QUALITY"])
            sap = np.array(d["SAP_FLUX"], float)
            det = np.array(d["DET_FLUX"] if "DET_FLUX" in d.columns.names else d["KSPSAP_FLUX"], float)
        shift = EPH2["qlp_shift_days"].get(str(s), 0.0) + TC[str(s)]["add_s_mean"] / 86400.0
        t = t + shift
        ex = float(np.nanmedian(np.diff(t)))
        for col, fl in (("SAP", sap), ("DET", det)):
            tt, ff, raw = prep(t, fl, q == 0)
            dep, de, dscr, sd = depth_scr(tt, ff, ex)
            rows.append(dict(product="QLP", col=col, star=star, sector=s, depth_ppt=dep * 1e3, err_white=de * 1e3, err_scramble=dscr * 1e3,
                             rstd_ppt=sd * 1e3, n=len(tt), time_shift_s=shift * 86400))
# ---- TARS (M3V)
for s in (12, 39, 66, 93):
    fn = f"{DL}/hlsp_tars_tess_ffi_s{s:04d}-0000000214311908_tess_v01_lc.fits"
    with fits.open(fn) as h:
        d = h[1].data; t = np.array(d["TIME"], float); f = np.array(d["FLUX"], float)
        pers = [h[1].header.get(f"PER{i}") for i in range(1, 6)]
    tt, ff, raw = prep(t + TC[str(s)]["add_s_mean"] / 86400.0, f, np.ones(len(t), bool))
    ex = float(np.nanmedian(np.diff(tt)))
    dep, de, dscr, sd = depth_scr(tt, ff, ex)
    rows.append(dict(product="TARS", col="FLUX", star="M3V", sector=s, depth_ppt=dep * 1e3, err_white=de * 1e3, err_scramble=dscr * 1e3,
                     rstd_ppt=sd * 1e3, n=len(tt), tars_periods=str(pers)))
# ---- eleanor-lite S12 (three stars)
for tic, star in ((214311908, "M3V"), (214311910, "pair"), (214311906, "Mira")):
    fn = f"{DL}/hlsp_gsfc-eleanor-lite_tess_ffi_s0012-{tic:016d}_tess_v1.0_lc.fits"
    with fits.open(fn) as h:
        d = h[1].data; t = np.array(d["TIME"], float); q = np.array(d["QUALITY"])
        cols = {c: np.array(d[c], float) for c in ("RAW_FLUX", "CORR_FLUX", "PCA_FLUX")}
    for c, fl in cols.items():
        tt, ff, raw = prep(t + TC["12"]["add_s_mean"] / 86400.0, fl, q == 0)
        if len(tt) < 100: continue
        dep, de, dscr, sd = depth_scr(tt, ff, 1800 / 86400)
        rows.append(dict(product="eleanor-lite", col=c, star=star, sector=12, depth_ppt=dep * 1e3, err_white=de * 1e3, err_scramble=dscr * 1e3,
                         rstd_ppt=sd * 1e3, n=len(tt)))
D = pd.DataFrame(rows)
D["err_adopted"] = D[["err_white", "err_scramble"]].max(axis=1); D["snr"] = D.depth_ppt / D.err_adopted
D.to_csv(f"{ROOT}/out/public_depths.csv", index=False, float_format="%.3f")
pd.set_option("display.width", 220)
print(D[["product", "col", "star", "sector", "depth_ppt", "err_adopted", "snr", "rstd_ppt", "n"]].round(2).to_string())
for (pr, col, star), g in D.groupby(["product", "col", "star"]):
    w = 1 / g.err_adopted ** 2; mu = np.sum(w * g.depth_ppt) / w.sum()
    print(f"{pr} {col} {star}: weighted mean {mu:.2f} +- {1/np.sqrt(w.sum()):.2f} ppt over {len(g)} sectors; chi2 {np.sum(w*(g.depth_ppt-mu)**2):.1f}/{len(g)-1}")

# ---- timing offset QLP vs TESScut (cross-correlation of high-pass filtered systematics)
lag = {}
for s in SECTORS:
    C = Cutout(s); W = json.load(open(f"{ROOT}/out/wcs_s{s}.json"))
    g = C.gaia(15); r = g[g.source_id == IDS["M3V"]].iloc[0]
    yy, xx = np.mgrid[0:C.ny, 0:C.nx]
    ap = np.hypot(xx - r.x - W["global_shift"][0], yy - r.y - W["global_shift"][1]) <= 1.5
    ok = C.fin & (C.q == 0)
    med = np.nanmedian(C.flux[ok], 0); faint = med < np.percentile(med, 15)
    bk = np.nanmedian(C.flux[:, faint], 1)
    fa = C.flux[:, ap].sum(1) - bk * ap.sum()
    ta = C.t + TC[str(s)]["add_s_mean"] / 86400.0
    fn = f"{LCDIR}/hlsp_qlp_tess_ffi_s{s:04d}-0000000214311908_tess_v01_llc.fits"
    with fits.open(fn) as h:
        d = h[1].data; tq = np.array(d["TIME"], float); fq = np.array(d["SAP_FLUX"], float); qq = np.array(d["QUALITY"])
    tq = tq + EPH2["qlp_shift_days"].get(str(s), 0.0) + TC[str(s)]["add_s_mean"] / 86400.0
    mq = (qq == 0) & np.isfinite(tq) & np.isfinite(fq)
    ex = C.exptime
    # high-pass: subtract running median over 0.3 d, on each series' own grid
    def hp(t, f):
        from scipy.ndimage import median_filter
        n = max(3, int(round(0.3 / np.median(np.diff(t)))) | 1)
        return f - median_filter(f, size=n, mode="nearest")
    m2 = ok & np.isfinite(fa)
    ha = hp(ta[m2], fa[m2] / np.nanmedian(fa[m2])); hq = hp(tq[mq], fq[mq] / np.nanmedian(fq[mq]))
    # clip outliers
    ha = np.clip(ha, -5 * robust_std(ha), 5 * robust_std(ha)); hq = np.clip(hq, -5 * robust_std(hq), 5 * robust_std(hq))
    lags = np.arange(-2.0, 2.0001, 0.05) * ex
    cc = []
    for L in lags:
        v = np.interp(ta[m2], tq[mq] - L, hq, left=np.nan, right=np.nan)
        k = np.isfinite(v) & np.isfinite(ha)
        # do not interpolate across QLP gaps
        cc.append(np.corrcoef(ha[k], v[k])[0, 1])
    cc = np.array(cc); i = int(np.nanargmax(cc))
    # parabola refinement
    if 0 < i < len(cc) - 1:
        a_, b_, c_ = np.polyfit(lags[i - 1:i + 2], cc[i - 1:i + 2], 2); best = -b_ / (2 * a_)
    else: best = lags[i]
    lag[s] = dict(lag_min=float(best * 1440), peak_corr=float(cc[i]), cadence_min=float(ex * 1440),
                  note="positive lag: QLP time stamps late relative to TESScut (after the corrections applied)")
    print(f"S{s}: QLP-TESScut lag {best*1440:+.2f} min (peak corr {cc[i]:.3f}; cadence {ex*1440:.1f} min)", flush=True)
save_json(lag, f"{ROOT}/out/public_timing_lag.json")
