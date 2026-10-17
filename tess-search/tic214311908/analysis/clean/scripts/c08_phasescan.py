"""c08: phase scans of the cleaned light curves (both pipelines) to re-check (i) the 1.3-1.9 ppt bump at phase
0.490-0.495 and (ii) the 3.3-3.6 sigma excursions near phase 0.57 reported on the uncleaned light curves; plus the
secondary-eclipse depth at phase 0.5.

Statistic: at each scan phase phi (0.15-0.85, step 0.0025; ephemeris = ephemeris_final 'combined'), events at
T0 + P (E + phi) for every epoch E with data; per event a +-0.30 d window (INV: undetrended flux of the cleaned INV
curve, quality_ok, flares removed; REF: a dedicated referee-method window centred on the event, built from the cleaned
REF photometry, i.e. per-pixel baselines with the centre masked, exactly as for transits); a common 'dip' depth with a
fixed trapezoid shape (T14 = 1.70 h, T23/T14 = 0.72 = the planet's occultation shape) and per-window linear baselines.
Positive depth = dimming. Errors: formal (per-window OOT robust sigma), then the scan S/N is divided by its robust std
over all scan phases (empirical calibration, as in phase 1).
Variant without S12 (30-min data). Outputs: out/phase_scan_{inv,ref}[_noS12].csv, out/phase_scan.json
"""
import sys
sys.path.insert(0, "tic214311908/analysis/clean/scripts")
from cleanlib import *
from refclean import SectorC
import pandas as pd

EPH = json.load(open(f"{CL}/out/ephemeris_final.json"))["combined"]
P, T0 = EPH["P_d"], EPH["T0_epoch0_BTJD_TDB"]
T14, FL = 1.70 / 24, 0.72
PHASES = np.round(np.arange(0.15, 0.8501, 0.0025), 4)
HALF = 0.30


def accept(tt, tc, half, ex):
    oot = np.abs(tt - tc) > 0.07
    nexp = 2 * (half - 0.07) / ex
    return np.sum(oot & (tt < tc)) >= 0.25 * nexp / 2 and np.sum(oot & (tt > tc)) >= 0.25 * nexp / 2 and np.sum(np.abs(tt - tc) < 0.028) >= 2


def common_depth(wins):
    """wins: list of (t, f, tc, ex). Common depth with per-window linear baselines."""
    if not wins: return np.nan, np.nan, 0
    nw = len(wins); ncol = 1 + 2 * nw
    rows = []; ys = []
    for i, (t, f, tc, ex) in enumerate(wins):
        oot = np.abs(t - tc) > 0.07
        sg = robust_std(f[oot] - np.polyval(np.polyfit(t[oot] - tc, f[oot], 1), t[oot] - tc))
        X = np.zeros((len(t), ncol)); X[:, 0] = -trap_frac(t, tc, T14, FL, ex); X[:, 1 + 2 * i] = 1; X[:, 2 + 2 * i] = t - tc
        rows.append(X / sg); ys.append(f / sg)
    X = np.vstack(rows); y = np.concatenate(ys)
    c, *_ = np.linalg.lstsq(X, y, rcond=None); cov = np.linalg.pinv(X.T @ X)
    return c[0] * 1e3, np.sqrt(cov[0, 0]) * 1e3, nw


def inv_scan(secs=SECTORS):
    data = {}
    for s in secs:
        z = np.load(f"{CL}/data/lcinv_all_s{s}.npz")
        t = z["t_tdb"]; f = z["raw_M3V_B"]; ok = z["quality_ok"].astype(bool) & np.isfinite(f) & np.isfinite(t)
        det = z["det_M3V_B"]; sd = robust_std(det[ok & ~z["intransit_window"].astype(bool)])
        hi = ok & (det > 1 + 4 * sd); ok &= ~(hi | np.roll(hi, 1) | np.roll(hi, -1))
        data[s] = (t[ok], f[ok], float(z["exptime"]))
    rows = []
    for ph in PHASES:
        wins = []
        for s, (t, f, ex) in data.items():
            E = np.arange(np.floor((t.min() - T0) / P) - 1, np.ceil((t.max() - T0) / P) + 1)
            for e in E:
                tc = T0 + P * (e + ph); m = np.abs(t - tc) < HALF
                if m.sum() >= 8 and accept(t[m], tc, HALF, ex): wins.append((t[m], f[m], tc, ex))
        d, de, n = common_depth(wins)
        rows.append(dict(phase=ph, depth_ppt=d, err_ppt=de, n_events=n))
    return pd.DataFrame(rows)


def ref_scan(secs=SECTORS):
    rows = {ph: [] for ph in PHASES}
    for s in secs:
        S = SectorC(s, "all")
        t = S.ttdb[S.qok]
        E = np.arange(np.floor((t.min() - T0) / P) - 1, np.ceil((t.max() - T0) / P) + 1)
        for ph in PHASES:
            for e in E:
                tc = T0 + P * (e + ph)
                if tc < t.min() or tc > t.max(): continue
                w = S.window(tc, half=HALF)
                if w is None: continue
                tt = w["t"]; ff = 1 + w["aM"] / F_M3V
                if accept(tt, tc, HALF, S.exptime): rows[ph].append((tt, ff, tc, S.exptime))
        print(f"REF scan S{s} done", flush=True)
    out = []
    for ph in PHASES:
        d, de, n = common_depth(rows[ph]); out.append(dict(phase=ph, depth_ppt=d, err_ppt=de, n_events=n))
    return pd.DataFrame(out)


def summarise(df, lab):
    df["snr"] = df.depth_ppt / df.err_ppt
    far = np.abs(((df.phase + 0.5) % 1) - 0.5) > 0.03
    rs = robust_std(df.snr[far]); df["snr_cal"] = df.snr / rs
    df["err_cal_ppt"] = df.err_ppt * rs
    b = df[(df.phase >= 0.4899) & (df.phase <= 0.4951)]
    n57 = df[(df.phase >= 0.555) & (df.phase <= 0.595)]
    i_max = int(np.nanargmax(np.abs(df.snr_cal.values)))
    ntr = (0.85 - 0.15) / (T14 / P)
    zmax = float(np.abs(df.snr_cal.values[i_max]))
    from scipy import stats
    p_glob = float(1 - (1 - 2 * stats.norm.sf(zmax)) ** ntr)
    s05 = df.iloc[int(np.argmin(np.abs(df.phase - 0.5)))]
    res = dict(scan_snr_robust_std=float(rs), n_independent_trials=float(ntr),
               bump_0490_0495=dict(depth_ppt=b.depth_ppt.round(3).tolist(), snr_cal=b.snr_cal.round(2).tolist(), phases=b.phase.tolist()),
               near_057=dict(min_snr_cal=float(n57.snr_cal.min()), phase_min=float(n57.phase.values[np.argmin(n57.snr_cal.values)]),
                             max_snr_cal=float(n57.snr_cal.max()), phase_max=float(n57.phase.values[np.argmax(n57.snr_cal.values)]),
                             depth_at_0576=float(df.depth_ppt.values[np.argmin(np.abs(df.phase.values - 0.576))]),
                             snr_at_0576=float(df.snr_cal.values[np.argmin(np.abs(df.phase.values - 0.576))]),
                             depth_at_05725=float(df.depth_ppt.values[np.argmin(np.abs(df.phase.values - 0.5725))]),
                             snr_at_05725=float(df.snr_cal.values[np.argmin(np.abs(df.phase.values - 0.5725))])),
               global_max=dict(phase=float(df.phase.values[i_max]), snr_cal=float(df.snr_cal.values[i_max]),
                               depth_ppt=float(df.depth_ppt.values[i_max]), p_global_gauss=p_glob),
               secondary_phase05=dict(depth_ppt=float(s05.depth_ppt), err_cal_ppt=float(s05.err_cal_ppt),
                                      ul3sigma_ppt=float(max(s05.depth_ppt, 0) + 3 * s05.err_cal_ppt)))
    print(lab, json.dumps(res)[:1500], flush=True)
    return res


if __name__ == "__main__":
    out = {}
    di = inv_scan(); out["inv"] = summarise(di, "INV"); di.to_csv(f"{CL}/out/phase_scan_inv.csv", index=False, float_format="%.5f")
    dr = ref_scan(); out["ref"] = summarise(dr, "REF"); dr.to_csv(f"{CL}/out/phase_scan_ref.csv", index=False, float_format="%.5f")
    m = di.merge(dr, on="phase", suffixes=("_inv", "_ref"))
    out["inv_ref_scan_correlation"] = float(np.corrcoef(m.depth_ppt_inv, m.depth_ppt_ref)[0, 1])
    # variant without the 30-min S12 data (its few-point windows at gap edges give the REF single-step spike at 0.5375)
    s5 = tuple(s for s in SECTORS if s != 12)
    di5 = inv_scan(s5); out["inv_noS12"] = summarise(di5, "INV noS12"); di5.to_csv(f"{CL}/out/phase_scan_inv_noS12.csv", index=False, float_format="%.5f")
    dr5 = ref_scan(s5); out["ref_noS12"] = summarise(dr5, "REF noS12"); dr5.to_csv(f"{CL}/out/phase_scan_ref_noS12.csv", index=False, float_format="%.5f")
    jd(out, f"{CL}/out/phase_scan.json")
