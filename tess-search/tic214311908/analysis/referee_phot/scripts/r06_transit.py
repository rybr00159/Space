"""Referee r06: transit measurements on the referee local-window M3V light curve (r05, 'M' fit, flux change in e-/s).
1. global trapezoid fit (T0, P, T14, flat fraction, flux loss) on all strict windows of the six sectors;
   V-shape (flat=0) comparison.
2. per-sector flux loss at the global shape; errors: white (per-window OOT scatter) and empirical (the same fit at
   the fake-epoch windows of the same sector: std of single-window depths / sqrt(n_real_windows)).
3. per-transit fits: free tc (+-0.06 d) and depth, fixed global T14/flat; error on tc and depth calibrated with
   pixel-level injections into fake windows (r07).
4. odd/even depths (per transit and folded), depth at phase 0.5 (separate windows), 'three-star' fit amplitudes.
Output: out/transit_fit.json, out/per_transit.csv, out/per_sector.csv"""
import sys; sys.path.insert(0, "tic214311908/analysis/referee_phot/scripts")
from rlib import *
from r05_local import trap_shape, P_G, T0_G
from scipy.optimize import least_squares


def load(s, lab="strict"):
    z = np.load(f"{R}/data/win_s{s}.npz", allow_pickle=True)
    d = {k[len(lab) + 1:]: z[k] for k in z.files if k.startswith(lab + "_")}
    d["exptime"] = float(z["exptime"]); d["sector"] = s
    return d


def window_sigma(d, key="aM"):
    """per-cadence error: robust std of OOT values of each window."""
    v = d[key] if key == "aM" else d["aT"][:, key]
    sg = np.zeros(len(v))
    for i in np.unique(d["win"]):
        m = d["win"] == i; sg[m] = rstd(v[m & d["oot"]])
    return sg


def model_all(p, t, exps):
    T0, P, T14, fl, dep = p
    y = np.zeros(len(t))
    for ex in np.unique(exps):
        m = exps == ex
        tc = T0 + np.round((t[m] - T0) / P) * P
        y[m] = trap_shape(t[m], tc, T14, fl, ex) if np.ndim(tc) == 0 else np.array([0.0])
    return -dep * y


def trapv(t, tc, T14, fl, ex):
    # vectorised over per-cadence tc
    nsub = 9; y = np.zeros(len(t)); T23 = fl * T14
    for o in (np.arange(nsub) + 0.5) / nsub - 0.5:
        a = np.abs(t + o * ex - tc)
        y += np.where(a <= T23 / 2, 1.0, np.where(a >= T14 / 2, 0.0, (T14 / 2 - a) / max((T14 - T23) / 2, 1e-9)))
    return y / nsub


def shape_eph(t, ex, T0, P, T14, fl):
    tc = T0 + np.round((t - T0) / P) * P
    return trapv(t, tc, T14, fl, ex)


def main():
    out = {}
    D = {s: load(s) for s in SECTORS}
    F = {s: load(s, "fake") for s in SECTORS}
    t = np.concatenate([D[s]["t"] for s in SECTORS]); f = np.concatenate([D[s]["aM"] for s in SECTORS])
    e = np.concatenate([window_sigma(D[s]) for s in SECTORS]); ex = np.concatenate([np.full(len(D[s]["t"]), D[s]["exptime"]) for s in SECTORS])
    sec = np.concatenate([np.full(len(D[s]["t"]), s) for s in SECTORS])
    def res(p):
        dT0, dP, T14, fl, dep = p
        y = np.zeros(len(t))
        for x_ in np.unique(ex):
            m = ex == x_; y[m] = shape_eph(t[m], x_, T0_G + dT0, P_G + dP, T14, fl)
        return (f + dep * y) / e
    best = None
    for T14 in (0.065, 0.072, 0.08):
        for fl in (0.3, 0.6, 0.85):
            q = least_squares(res, [0, 0, T14, fl, 15], bounds=([-0.02, -5e-4, 0.03, 0, 0], [0.02, 5e-4, 0.12, 0.98, 60]),
                              x_scale=[1e-3, 1e-5, 3e-3, 0.05, 1], diff_step=[1e-4, 1e-7, 1e-3, 1e-2, 1e-3])
            if best is None or q.cost < best.cost: best = q
    q = best; cov = np.linalg.pinv(q.jac.T @ q.jac); chi2 = 2 * q.cost; n = len(t)
    pe = np.sqrt(np.diag(cov) * max(1, chi2 / (n - 5)))
    T0, P, T14, FL, DEP = T0_G + q.x[0], P_G + q.x[1], q.x[2], q.x[3], q.x[4]
    qv = least_squares(lambda p: res([p[0], p[1], p[2], 0.0, p[3]]), [q.x[0], q.x[1], q.x[2], q.x[4]], x_scale=[1e-3, 1e-5, 3e-3, 1], diff_step=[1e-4, 1e-7, 1e-3, 1e-3])
    out["global"] = dict(T0=T0, T0_err=pe[0], P=P, P_err=pe[1], T14_h=T14 * 24, T14_err_h=pe[2] * 24, flat=FL, flat_err=pe[3], dF_e_s=DEP, dF_err=pe[4],
                         chi2=chi2, n=n, chi2_V=2 * qv.cost, dchi2_V=2 * qv.cost - chi2, T14_V_h=qv.x[2] * 24)
    print(f"global: T0 {T0:.5f}+-{pe[0]:.5f} P {P:.7f}+-{pe[1]:.7f} T14 {T14*24:.2f}+-{pe[2]*24:.2f} h flat {FL:.2f}+-{pe[3]:.2f} dF {DEP:.2f}+-{pe[4]:.2f} e/s; chi2 {chi2:.0f}/{n}; V-shape dchi2 {2*qv.cost-chi2:.1f}", flush=True)
    # profile over flat fraction
    prof = []
    for flv in (0.0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9):
        qq = least_squares(lambda p: res([p[0], p[1], p[2], flv, p[3]]), [q.x[0], q.x[1], q.x[2], q.x[4]], x_scale=[1e-3, 1e-5, 3e-3, 1], diff_step=[1e-4, 1e-7, 1e-3, 1e-3])
        prof.append((flv, 2 * qq.cost - chi2, qq.x[2] * 24, qq.x[3]))
    out["flat_profile"] = prof
    print("flat profile (flat, dchi2, T14 h, dF):", [tuple(np.round(p_, 2)) for p_ in prof])
    # ---- 2. per-sector depths with empirical errors from fake windows
    rows = []
    for s in SECTORS:
        d = D[s]; y = shape_eph(d["t"], d["exptime"], T0, P, T14, FL); sg = window_sigma(d)
        def lin(v, yy, ss):
            A = np.column_stack([-yy]) / ss[:, None]; c, *_ = np.linalg.lstsq(A, v / ss, rcond=None)
            return c[0], 1 / np.sqrt(np.sum((yy / ss) ** 2))
        dep, dwh = lin(d["aM"], y, sg)
        # empirical: single-window depth at fake centres (shape centred on each fake window centre)
        fk = F[s]; sgf = window_sigma(fk); singles = []
        for i in np.unique(fk["win"]):
            m = fk["win"] == i; yy = trapv(fk["t"][m], fk["tc"][i], T14, FL, fk["exptime"])
            if np.sum(yy > 0.5) < 2: continue
            singles.append(lin(fk["aM"][m], yy, sgf[m])[0])
        singles = np.array(singles)
        reals = []
        for i in np.unique(d["win"]):
            m = d["win"] == i; reals.append(lin(d["aM"][m], y[m], sg[m])[0])
        nw = len(np.unique(d["win"]))
        demp = float(np.std(singles) / np.sqrt(nw)); demp_r = float(1.4826 * np.median(np.abs(singles - np.median(singles))) / np.sqrt(nw))
        # three-star amplitudes
        amps = {}
        for j, nm in enumerate(("M3V", "pair", "Mira")):
            sgj = window_sigma(d, j); a_, ae_ = lin(d["aT"][:, j], y, sgj)
            nul = []
            for i in np.unique(fk["win"]):
                m = fk["win"] == i; yy = trapv(fk["t"][m], fk["tc"][i], T14, FL, fk["exptime"])
                if np.sum(yy > 0.5) < 2: continue
                sgf_j = rstd(fk["aT"][m & fk["oot"], j]); nul.append(lin(fk["aT"][m, j], yy, np.full(m.sum(), sgf_j))[0])
            amps[nm] = (float(a_), float(max(ae_, np.std(nul) / np.sqrt(nw))))
        rows.append(dict(sector=s, n_windows=nw, dF=dep, err_white=dwh, err_fake=demp, err_fake_robust=demp_r, null_mean=float(np.mean(singles)),
                         n_fake=len(singles), err_adopt=max(dwh, demp), single_depths=np.round(reals, 2).tolist(),
                         T_M3V=amps["M3V"][0], T_M3V_err=amps["M3V"][1], T_pair=amps["pair"][0], T_pair_err=amps["pair"][1], T_Mira=amps["Mira"][0], T_Mira_err=amps["Mira"][1]))
    PS = pd.DataFrame(rows); PS.to_csv(f"{R}/out/per_sector.csv", index=False)
    w = 1 / PS.err_adopt ** 2; mu = float(np.sum(w * PS.dF) / w.sum()); c2 = float(np.sum(w * (PS.dF - mu) ** 2))
    out["per_sector_mean"] = dict(dF=mu, err=float(1 / np.sqrt(w.sum())), chi2=c2, dof=len(PS) - 1)
    for nm in ("M3V", "pair", "Mira"):
        ww = 1 / PS[f"T_{nm}_err"] ** 2; out[f"three_{nm}"] = (float(np.sum(ww * PS[f"T_{nm}"]) / ww.sum()), float(1 / np.sqrt(ww.sum())))
    pd.set_option("display.width", 250)
    print(PS.drop(columns=["single_depths"]).round(2).to_string())
    print(f"sector mean dF {mu:.2f}+-{1/np.sqrt(w.sum()):.2f} e/s, chi2 {c2:.1f}/{len(PS)-1}; three-star fit: " + ", ".join(f"{nm} {out['three_'+nm][0]:.1f}+-{out['three_'+nm][1]:.1f}" for nm in ("M3V", "pair", "Mira")))
    save_json(out, f"{R}/out/transit_fit.json")
    return out


if __name__ == "__main__":
    main()
