"""c05: pixel-level injection of a CONSTANT 15.7-ppt transit (15.37 e-/s at mid-transit; batman shape of the adopted
solution: k 0.1195, b 0.52, a/R* 25.4, u 0.125/0.533, exposure-integrated) at the M3V position into the raw TESScut
cubes, at 10 fake ephemerides per sector (P = 5.6915087 d, phase offsets 0.13-0.87, so >= 0.74 d from any real
transit), and re-run of BOTH cleaned pipelines:
  INV : cleanlib.photometry (mode 'all'; fake windows +-0.1 d masked in the pixel regression, as for real transits)
        + make_lc (fake windows masked in the detrending), then per-window fits as for the real data.
  REF : refclean.SectorC (mode 'all') windows centred on each fake transit.
Injected images ('variants'; the fitting templates are never changed):
  invprf : the INV fitting template (photlib PRF at the INV M3V position)          -> matched for INV
  refprf : the REF fitting template (PRF.locate at the REF M3V position)           -> matched for REF
  sx+, sx-, sy+, sy- : invprf displaced by +-0.07 px in column / row               -> template error (WCS/PRF)
  blur   : invprf convolved with a 0.15-px Gaussian                                -> template error (focus/jitter)
For each fake transit the same template is fitted (fixed batman shape; depth at the true mid-time, and depth + free
mid-time); a no-injection run with identical masks is subtracted (paired difference = pipeline response).
Outputs: out/inject_pix_s{sector}.csv (one row per fake transit and variant), logs.
"""
import sys
sys.path.insert(0, "tic214311908/analysis/clean/scripts")
from cleanlib import *
from refclean import SectorC
from tfit import Model, mid_depth, single_fit
import pandas as pd

K, B_, AR, U1, U2 = 0.1195, 0.52, 25.4, 0.125, 0.533
D0 = mid_depth(K, B_, AR, U1, U2, P1)
DEP_PPT = 15.7
AMP = DEP_PPT * 1e-3 * F_M3V           # e-/s at mid-transit
FAKE_F = (0.13, 0.21, 0.29, 0.37, 0.45, 0.55, 0.63, 0.71, 0.79, 0.87)
HALF = 0.30


def shape_frac(t, T0, ex):
    """Fraction of mid-transit depth (exposure-integrated) for an ephemeris (T0, P1)."""
    M = Model(t, ex)
    return (1 - M.flux(T0, P1, K, B_, AR, U1, U2)) / D0


def fit_window(tt, ff, tc, ex, order=1, free=True):
    """depth at fixed tc (mid-transit units) and (if free) a free-tc fit."""
    oot = np.abs(tt - tc) > 0.07
    if oot.sum() < 6: return None
    sg = robust_std(ff[oot] - np.polyval(np.polyfit(tt[oot] - tc, ff[oot], 1), tt[oot] - tc))
    s = np.full(len(tt), max(sg, 1e-6))
    y = shape_frac(tt, tc, ex)
    X = np.column_stack([-y] + [(tt - tc) ** o for o in range(order + 1)]) / s[:, None]
    c, *_ = np.linalg.lstsq(X, ff / s, rcond=None)
    cov = np.linalg.pinv(X.T @ X)
    sf = single_fit(tt, ff, s, ex, tc, lambda t_, tc_, ex_: shape_frac(t_, tc_, ex_), dt_max=0.05, order=order, step=0.001) if free else None
    return dict(depth_fix=c[0] * 1e3, depth_fix_err=np.sqrt(cov[0, 0]) * 1e3,
                depth_free=(sf["depth"] * 1e3 if sf else np.nan), tc_free=(sf["tc"] if sf else np.nan),
                tc_err=(sf["tc_err"] if sf else np.nan), sig=sg)


def accept(tt, tc, half, exptime):
    oot = np.abs(tt - tc) > 0.07
    nexp = 2 * (half - 0.07) / exptime
    return (np.sum(oot & (tt < tc)) >= 0.25 * nexp / 2) and (np.sum(oot & (tt > tc)) >= 0.25 * nexp / 2) and np.sum(np.abs(tt - tc) < 0.028) >= 2


def run(sector):
    S = setup(sector, "all"); C = S["C"]
    t_tdb = t_tdb_of(sector, C.t, C.timecorr)
    m3 = S["g"][S["g"].name == "M3V"].iloc[0]; prf = S["prf"]
    R_ = SectorC(sector, "all")
    imgs = {"invprf": prf.image(m3.x, m3.y), "refprf": R_.P_M3V,
            "sx+": prf.image(m3.x + 0.07, m3.y), "sx-": prf.image(m3.x - 0.07, m3.y),
            "sy+": prf.image(m3.x, m3.y + 0.07), "sy-": prf.image(m3.x, m3.y - 0.07),
            "blur": prf.image(m3.x, m3.y, blur=np.hypot(prf.blur, 0.15))}
    flat0 = R_.flat.copy()
    rows = []
    tfin = np.isfinite(t_tdb)
    for f in FAKE_F:
        T0f = T01 + f * P1
        k = np.arange(np.floor((np.nanmin(t_tdb) - T0f) / P1) - 1, np.ceil((np.nanmax(t_tdb) - T0f) / P1) + 2)
        fake_ep = T0f + k * P1
        fake_ep = fake_ep[(fake_ep > np.nanmin(t_tdb) - 0.2) & (fake_ep < np.nanmax(t_tdb) + 0.2)]
        amp = np.zeros(len(t_tdb)); amp[tfin] = AMP * shape_frac(t_tdb[tfin], T0f, C.exptime)
        # REF times (referee's own barycentric correction; agrees within 0.6 s)
        ampR = np.zeros(len(R_.ttdb)); fr = np.isfinite(R_.ttdb); ampR[fr] = AMP * shape_frac(R_.ttdb[fr], T0f, R_.exptime)
        res = {}
        for var in ("none",) + tuple(imgs):
            # ---- INV
            inj = None if var == "none" else [dict(img=imgs[var], amp=amp)]
            R = photometry(sector, S, inject=inj, extra_mask_times=fake_ep, fit_models=("B",))
            extra = in_windows(t_tdb, fake_ep, 0.1)
            L = make_lc(sector, R, extra_mask=extra, models=(("B", "M3V"),))
            ok = L["quality_ok"] & np.isfinite(L["raw_M3V_B"])
            # ---- REF
            R_.flat = flat0 if var == "none" else flat0 - ampR[:, None] * imgs[var].ravel()[None]
            for j, tc in enumerate(fake_ep):
                key = (j, var)
                out = {}
                m = ok & (np.abs(L["t_tdb"] - tc) < HALF)
                if m.sum() >= 8 and accept(L["t_tdb"][m], tc, HALF, C.exptime):
                    for ver in ("raw", "det"):
                        r = fit_window(L["t_tdb"][m], L[f"{ver}_M3V_B"][m], tc, C.exptime)
                        if r: out.update({f"inv_{ver}_{kk}": vv for kk, vv in r.items()})
                w = R_.window(tc)
                if w is not None:
                    tt = w["t"]; ff = 1 + w["aM"] / F_M3V
                    if accept(tt, tc, {12: 0.6, 39: 0.4}.get(sector, 0.35), R_.exptime):
                        r = fit_window(tt, ff, tc, R_.exptime)
                        if r: out.update({f"ref_{kk}": vv for kk, vv in r.items()})
                res[key] = out
        R_.flat = flat0
        for j, tc in enumerate(fake_ep):
            base = res[(j, "none")]
            for var in imgs:
                o = res[(j, var)]
                row = dict(sector=sector, fake_phase=f, tc=tc, variant=var)
                for kk, vv in o.items():
                    row[kk] = vv
                    if kk.startswith(("inv_raw_depth", "inv_det_depth", "ref_depth")) and kk in base:
                        row[kk + "_resp"] = vv - base[kk]
                    if kk.endswith("tc_free"):
                        row[kk.replace("tc_free", "dt_min")] = (vv - tc) * 1440
                for kk, vv in base.items():
                    if "depth_fix" in kk and not kk.endswith("err"): row[kk + "_noinj"] = vv
                rows.append(row)
        d = pd.DataFrame([r for r in rows if r["fake_phase"] == f and r["variant"] in ("invprf", "refprf")])
        if len(d):
            print(f"S{sector} f{f}: n {len(d)//2}; INV raw resp (invprf) {np.nanmean(d[d.variant=='invprf'].get('inv_raw_depth_fix_resp', np.nan)):.2f}; "
                  f"REF resp (refprf) {np.nanmean(d[d.variant=='refprf'].get('ref_depth_fix_resp', np.nan)):.2f} ppt", flush=True)
    pd.DataFrame(rows).to_csv(f"{CL}/out/inject_pix_s{sector}.csv", index=False)


if __name__ == "__main__":
    for s in [int(a) for a in sys.argv[1:]] or SECTORS:
        run(s)
