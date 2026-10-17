"""Referee r07: per-transit mid-times and depths on the referee local-window M3V light curve, with errors calibrated by
pixel-level injections.
Per transit: least squares (tc free within +-0.06 d of the window centre, depth free; T14 and flat fraction fixed at
the r06 global values; the r05 baseline already removed), formal errors scaled by sqrt(chi2_red >= 1), and a
chi2-profile interval for tc.
Calibration: in each sector, 80 fake windows (r05 fake centres) receive a pixel-level injected transit (PRF of the
M3V, 15 e-/s, global shape, true time = window centre + U(-0.015, 0.015) d), are re-processed by r05 and fitted
exactly as the real events. Pull distributions -> per-sector error scale factors (std and robust).
A second injection set uses a deliberately mismatched template (PRF displaced by (+0.10,-0.10) px and blurred by
0.15 px) to measure template-mismatch depth bias.
Outputs: out/per_transit.csv, out/inject_timing.json"""
import sys; sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
from rlib import *
from r05_local import Sector, trap_shape, real_centres, fake_centres
from scipy.optimize import least_squares
G = json.load(open(f"{R}/out/transit_fit.json"))["global"]
T0, P, T14, FL = G["T0"], G["P"], G["T14_h"] / 24, G["flat"]


def fit_one(t, v, ex, tc0, sg=None):
    oot = np.abs(t - tc0) > 0.07
    if sg is None: sg = rstd(v[oot])
    grid = np.arange(-0.06, 0.0601, 0.001)
    best = None; prof = []
    for dt in grid:
        y = trap_shape(t, tc0 + dt, T14, FL, ex)
        if np.sum(y > 0.5) < 2: prof.append(np.inf); continue
        A = np.column_stack([np.ones(len(t)), -y]); c, *_ = np.linalg.lstsq(A / sg, v / sg, rcond=None)
        c2 = np.sum(((v - A @ c) / sg) ** 2); prof.append(c2)
        if best is None or c2 < best[0]: best = (c2, dt, c)
    if best is None: return None
    prof = np.array(prof)
    def rr(p): return (v - (p[2] - p[1] * trap_shape(t, tc0 + p[0], T14, FL, ex))) / sg
    q = least_squares(rr, [best[1] + 1e-5, best[2][1], best[2][0]], x_scale=[2e-3, 2, 1], diff_step=[1e-2, 1e-3, 1e-3])
    cov = np.linalg.pinv(q.jac.T @ q.jac); red = max(1.0, 2 * q.cost / max(1, len(t) - 3))
    ok = prof <= prof.min() + 1.0 * red
    lo, hi = grid[ok].min(), grid[ok].max()
    return dict(tc=tc0 + q.x[0], depth=q.x[1], tc_err=float(np.sqrt(cov[0, 0] * red)), depth_err=float(np.sqrt(cov[1, 1] * red)),
                tc_prof_half=float((hi - lo) / 2), redchi=float(2 * q.cost / max(1, len(t) - 3)), n=len(t),
                n_in=int(np.sum(np.abs(t - tc0 - q.x[0]) < T14 / 2)))


def main():
    rows = []; calib = {}
    rng = np.random.default_rng(99)
    for s in SECTORS:
        S = Sector(s)
        cen, ks = real_centres(S, P, T0)
        for lab, bmax in (("strict", 60.0), ("relaxed", 800.0)):
            for tc, k in zip(cen, ks):
                r = S.window(tc, bmax=bmax)
                if r is None: continue
                f = fit_one(r["t"], r["aM"], S.exptime, tc)
                if f is None: continue
                # three-star amplitudes for this transit at the fitted shape (fixed tc = fitted)
                y = trap_shape(r["t"], f["tc"], T14, FL, S.exptime)
                amps = []
                for j in range(3):
                    v = r["aT"][:, j]; sg = rstd(v[r["oot"]]); A = np.column_stack([np.ones(len(v)), -y])
                    c, *_ = np.linalg.lstsq(A / sg, v / sg, rcond=None); amps.append(c[1])
                rows.append(dict(sector=s, epoch=int(k), sel=lab, B_max=float(np.max(r["B"]) - S.Bmed), cov=float(np.sum(np.abs(r["t"] - f["tc"]) < T14 / 2) / (T14 / S.exptime)),
                                 oc_min=(f["tc"] - tc) * 1440, **f, T_M3V=amps[0], T_pair=amps[1], T_Mira=amps[2]))
        # ---- injection calibration
        fk = fake_centres(S, 80, rng, P, T0)
        mis = S.C.psf(S.pos["M3V"][0] + 0.10, S.pos["M3V"][1] - 0.10)
        # blur by averaging shifted PRFs
        n_ = np.array([-2, -1, 0, 1, 2]) * 0.15 * 0.8; ww = np.exp(-0.5 * (n_ / 0.15) ** 2); ww /= ww.sum()
        mis = sum(ww[i] * ww[j] * S.C.psf(S.pos["M3V"][0] + 0.10 + n_[i], S.pos["M3V"][1] - 0.10 + n_[j]) for i in range(5) for j in range(5))
        pulls = []; dp = []; dmis = []; dpx = []
        for tcw in fk:
            tt = tcw + rng.uniform(-0.015, 0.015)
            r = S.window(tcw, inject=[(S.P_M3V, 15.0, tt, T14, FL)])
            if r is None: continue
            f = fit_one(r["t"], r["aM"], S.exptime, tcw)
            r0 = S.window(tcw)
            f0 = fit_one(r0["t"], r0["aM"], S.exptime, tcw) if r0 is not None else None
            if f is None: continue
            pulls.append((f["tc"] - tt) / f["tc_err"]); dp.append((f["depth"] - 15.0) / f["depth_err"]); dpx.append(f["depth"])
            rm = S.window(tcw, inject=[(mis, 15.0, tt, T14, FL)])
            if rm is not None:
                y = trap_shape(rm["t"], tt, T14, FL, S.exptime); y0 = trap_shape(r0["t"], tt, T14, FL, S.exptime)
                A = np.column_stack([np.ones(len(y)), -y]); c, *_ = np.linalg.lstsq(A, rm["aM"], rcond=None)
                A0 = np.column_stack([np.ones(len(y0)), -y0]); c0, *_ = np.linalg.lstsq(A0, r0["aM"], rcond=None)
                dmis.append((c[1] - c0[1]) / 15.0)
        pulls = np.array(pulls); dp = np.array(dp)
        calib[s] = dict(n=len(pulls), pull_std=float(np.std(pulls)), pull_rstd=float(rstd(pulls)), pull_mean=float(np.mean(pulls)), frac_gt3=float(np.mean(np.abs(pulls) > 3)),
                        depth_pull_std=float(np.std(dp)), depth_pull_mean=float(np.mean(dp)), depth_rec_mean=float(np.mean(dpx)),
                        mismatch_ratio_mean=float(np.mean(dmis)), mismatch_ratio_std=float(np.std(dmis)))
        print(f"S{s}: timing pulls n {len(pulls)} std {np.std(pulls):.2f} robust {rstd(pulls):.2f} mean {np.mean(pulls):+.2f} |>3| {np.mean(np.abs(pulls)>3):.2f}; depth pulls std {np.std(dp):.2f} mean {np.mean(dp):+.2f} "
              f"(recovered {np.mean(dpx):.2f} of 15); mismatched-template recovery {np.mean(dmis):.3f}+-{np.std(dmis):.3f}", flush=True)
    PT = pd.DataFrame(rows)
    PT["tc_err_cal"] = PT.tc_err * PT.sector.map(lambda s: max(1.0, calib[s]["pull_std"]))
    PT["depth_err_cal"] = PT.depth_err * PT.sector.map(lambda s: max(1.0, calib[s]["depth_pull_std"]))
    PT.to_csv(f"{R}/out/per_transit.csv", index=False, float_format="%.6f")
    save_json({str(k): v for k, v in calib.items()}, f"{R}/out/inject_timing.json")
    pd.set_option("display.width", 250)
    print(PT[["sector", "epoch", "sel", "B_max", "cov", "oc_min", "tc_err", "tc_err_cal", "tc_prof_half", "depth", "depth_err", "depth_err_cal", "redchi", "T_M3V", "T_pair", "T_Mira"]].round(4).to_string())


if __name__ == "__main__":
    main()
