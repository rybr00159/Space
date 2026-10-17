"""r01: independent check of the 0.6485-d 'RRab leak' in the phot model-B M3V light curve and its effect on the
5.69-d transit depths and odd/even.

A. Lomb-Scargle (astropy) per sector of the transit-masked, 10-min-binned detrended flux (0.1-15 d): top peaks and
   an empirical 1% level from 200 1-d-block-shuffled copies.
B. Coherent common-shape 3-harmonic scan (S39-S104, per-sector offsets) on a fine grid around 0.6485 d; best period vs
   the Gaia DR3 RRab period 0.6485275 d.
C. Per-sector 3-harmonic fits at the Gaia period with free per-sector coefficients; cleaned light curve.
D. Predicted leak bias of each transit depth (mean model in the flat part minus the mean over a +-0.3 d local
   baseline, for a fixed trapezoid T14 1.6 h, flat 0.7), and regressions of
     (i)  phot model-B strict per-transit depths (phot/out/per_transit_M3V.csv),
     (ii) the referee's independent local-window depths (referee_phot/out/per_transit.csv, strict, e-/s),
   on that bias: slope 1 means the photometry carries the leak, slope 0 means it does not.
E. Own local per-transit depths (trapezoid + linear baseline, +-0.3 d) on det and on the cleaned flux, empirical
   per-sector errors (100 fake epochs), odd-even.
Outputs: out/r01_rrab_leak.json, data/m3v_clean_r2.npz
"""
import sys, os
sys.path.insert(0, "tic214311908/analysis/referee2_search/scripts")
from rlib2 import *
from astropy.timeseries import LombScargle
rng = np.random.default_rng(101)


def main():
    out = {}
    d = load_m3v("B")
    oot = ~d.intr.values
    t, y, s, cad = d.t.values, d.det.values - 1, d.s.values, d.cad.values
    # ---------------- A
    A = {}
    fr = np.linspace(1 / 15, 10.0, 60000)
    for sec in SECTORS:
        m = (s == sec) & oot
        tb, yb, nb, sb = binlc(t[m], y[m], s[m], cad[m], 10.0)
        p = LombScargle(tb, yb).power(fr)
        order = np.argsort(p)[::-1]; top = []
        for i in order:
            if len(top) >= 4: break
            if all(abs(fr[i] - fr[j]) > 0.02 * fr[j] for j in top): top.append(i)
        # null: shuffle 1-d blocks
        kb = np.floor((tb - tb.min()) / 1.0).astype(int); ub = np.unique(kb)
        mx = []
        for _ in range(200):
            perm = rng.permutation(len(ub))
            y2 = np.concatenate([yb[kb == ub[i]] for i in perm])
            mx.append(LombScargle(tb, y2).power(fr[::4]).max())
        # power at the Gaia RR frequency
        pr = LombScargle(tb, yb).power(np.array([1 / P_RR_GAIA]))[0]
        A[sec] = dict(top=[(float(1 / fr[i]), float(p[i])) for i in top], null_p99=float(np.percentile(mx, 99)),
                      null_p50=float(np.median(mx)), power_at_PRR=float(pr))
        print(sec, A[sec], flush=True)
    out["A_lombscargle"] = A
    # ---------------- B coherent scan
    sel = np.isin(s, [39, 66, 93, 103, 104]) & oot
    ts, ys, ss = t[sel], y[sel], s[sel]
    wsec = {sec: 1 / rstd(ys[ss == sec]) ** 2 for sec in np.unique(ss)}
    w = np.array([wsec[q] for q in ss])
    off = np.column_stack([(ss == q).astype(float) for q in np.unique(ss)])

    def chi(Pp, coh=True):
        if coh:
            X = np.column_stack([off, harmX(ts, Pp)[:, 1:]])
        else:
            X = np.column_stack([off] + [harmX(ts, Pp)[:, 1:] * (ss == q)[:, None] for q in np.unique(ss)])
        sw = np.sqrt(w)
        c, *_ = np.linalg.lstsq(X * sw[:, None], ys * sw, rcond=None)
        return float(np.sum(((X @ c) - ys) ** 2 * w))
    Pg = np.arange(0.6460, 0.6510, 2.5e-7)
    cc = np.array([chi(p_) for p_ in Pg])
    j = int(np.argmin(cc))
    chi0 = float(np.sum(w * (ys - (off @ np.linalg.lstsq(off * np.sqrt(w)[:, None], ys * np.sqrt(w), rcond=None)[0])) ** 2))
    mins = [i for i in range(1, len(cc) - 1) if cc[i] < cc[i - 1] and cc[i] <= cc[i + 1]]
    mins = sorted(mins, key=lambda i: cc[i])[:5]
    # fake-period distribution of the coherent dchi2 (signal vs none), periods 0.5-0.9 d away from P_RR and harmonics
    fk = []
    for p_ in rng.uniform(0.5, 0.9, 120):
        if min(abs(p_ - P_RR_GAIA), abs(p_ - 2 * P_RR_GAIA), abs(p_ - P_RR_GAIA / 2 * 1.0)) < 0.01: continue
        fk.append(chi0 - chi(p_))
    B = dict(P_best=float(Pg[j]), dchi2_signal=float(chi0 - cc[j]), fake_dchi2_median=float(np.median(fk)),
             fake_dchi2_max=float(np.max(fk)), n_fake=len(fk), P_gaia=P_RR_GAIA,
             dchi2_at_gaia=float(chi(P_RR_GAIA) - cc[j]),
             teeth=[(float(Pg[i]), float(cc[i] - cc[j])) for i in mins],
             dchi2_incoherent_minus_coherent_at_best=float(chi(Pg[j], False) - cc[j]))
    print("B", B, flush=True)
    out["B_coherent"] = B
    # ---------------- C per-sector fits at the Gaia period
    C = {}
    model = np.zeros(len(t))
    for sec in SECTORS:
        m = (s == sec) & oot
        X = harmX(t[m], P_RR_GAIA)
        c, *_ = np.linalg.lstsq(X, y[m], rcond=None)
        res = y[m] - X @ c
        # block-bootstrap amplitude error (0.5-d blocks)
        kb = np.floor((t[m] - t[m].min()) / 0.5).astype(int); ub = np.unique(kb)
        idxs = [np.where(kb == u)[0] for u in ub]
        a1 = []
        for _ in range(150):
            ii = np.concatenate([idxs[i] for i in rng.integers(0, len(ub), len(ub))])
            cb, *_ = np.linalg.lstsq(X[ii], y[m][ii], rcond=None); a1.append(np.hypot(cb[1], cb[2]))
        amp = [float(np.hypot(c[1 + 2 * h], c[2 + 2 * h])) for h in range(3)]
        mf = s == sec
        model[mf] = harmX(t[mf], P_RR_GAIA) @ c - c[0]
        ptp = float(np.ptp(harmX(np.linspace(0, P_RR_GAIA, 500), P_RR_GAIA) @ c))
        C[sec] = dict(amp_ppt=[a * 1e3 for a in amp], amp1_err_ppt=float(np.std(a1) * 1e3), ptp_ppt=ptp * 1e3,
                      phase1=float(np.arctan2(c[2], c[1])))
        print(sec, C[sec], flush=True)
    out["C_per_sector"] = C
    clean = y - model
    np.savez_compressed(f"{ME}/data/m3v_clean_r2.npz", t=t, s=s, cad=cad, det=d.det.values, clean=1 + clean,
                        model=model, intr=d.intr.values)
    # ---------------- D predicted bias per transit and regressions
    rows = []
    for e in range(-400, 80):
        tc = T0 + e * P
        k = np.abs(t - tc) < 0.3
        if k.sum() < 10 or not np.any(np.abs(t[k] - tc) < 0.03): continue
        tm = trapz_model(t[k], tc, cad_s=cad[k])
        inn = tm > 0.99; base = tm == 0
        if inn.sum() < 1: continue
        # bias in the depth fitted with a linear baseline: regress model on [1, t, -tm]
        X = np.column_stack([np.ones(k.sum()), t[k] - tc, -tm])
        bias = np.linalg.lstsq(X, model[k], rcond=None)[0][2]
        rows.append(dict(epoch=e, sector=int(s[k][0]), bias_ppt=float(bias * 1e3)))
    Dm = pd.DataFrame(rows)
    PB = pd.read_csv(f"{BASE}/phot/out/per_transit_M3V.csv")
    print(PB.columns.tolist())
    out["D"] = {}
    PBs = PB[(PB.model == "B") & (PB.relaxed == False)].copy()
    dcol = "depth" if "depth" in PBs.columns else None
    PBs = PBs.merge(Dm, on="epoch", how="inner")
    R = pd.read_csv(f"{BASE}/referee_phot/out/per_transit.csv")
    R = R[(R.sel == "strict") & (R.epoch != -388)].merge(Dm, on="epoch", how="inner")

    def regress(dep, err, bias):
        # dep = a + b * bias ; weighted
        wv = 1 / err ** 2
        X = np.column_stack([np.ones(len(dep)), bias])
        cov = np.linalg.inv((X * wv[:, None]).T @ X)
        c = cov @ ((X * wv[:, None]).T @ dep)
        chi2 = float(np.sum(wv * (dep - X @ c) ** 2))
        return dict(intercept=float(c[0]), slope=float(c[1]), slope_err=float(np.sqrt(cov[1, 1])), chi2=chi2, n=len(dep),
                    slope_err_scaled=float(np.sqrt(cov[1, 1] * max(1, chi2 / (len(dep) - 2)))))
    # phot depths are fractional (935 normalisation): convert to ppt at 979
    dep_ph = PBs[dcol].values * 1e3 * SCALE
    err_ph = PBs["depth_err"].values * 1e3 * SCALE * 1.3
    out["D"]["phot_B_strict"] = regress(dep_ph, err_ph, PBs.bias_ppt.values)
    dep_rf = R.depth.values / FREF * 1e3; err_rf = R.depth_err_cal.values / FREF * 1e3
    out["D"]["referee_strict"] = regress(dep_rf, err_rf, R.bias_ppt.values)
    out["D"]["table"] = Dm.merge(R[["epoch", "depth", "depth_err_cal"]].rename(columns=dict(depth="ref_depth_e_s", depth_err_cal="ref_err_e_s")), on="epoch", how="left") \
        .merge(PBs[["epoch", dcol, "depth_err"]].rename(columns={dcol: "phot_depth_frac935", "depth_err": "phot_err_frac935"}), on="epoch", how="left").to_dict(orient="records")
    print("D", {k: v for k, v in out["D"].items() if k != "table"}, flush=True)
    # ---------------- E own local per-transit depths (det and clean) + odd/even
    def local_depth(yy, tc):
        k = np.abs(t - tc) < 0.3
        tm = trapz_model(t[k], tc, cad_s=cad[k])
        if (tm > 0.5).sum() < 2 or k.sum() < 10: return np.nan
        X = np.column_stack([np.ones(k.sum()), t[k] - tc, -tm])
        return float(np.linalg.lstsq(X, yy[k], rcond=None)[0][2])
    E = []
    for e in Dm.epoch.values:
        tc = T0 + e * P
        E.append(dict(epoch=int(e), sector=int(s[np.argmin(np.abs(t - tc))]), det=local_depth(y, tc) * 1e3, clean=local_depth(clean, tc) * 1e3))
    E = pd.DataFrame(E).dropna()
    errs = {}
    for sec in SECTORS:
        m = s == sec; fk = []
        while len(fk) < 100:
            tc = rng.uniform(t[m].min() + 0.3, t[m].max() - 0.3)
            if abs(phase_dist(tc)) < 0.35: continue
            v = local_depth(clean, tc)
            if np.isfinite(v): fk.append(v)
        errs[sec] = float(np.std(fk) * 1e3)
    E["err"] = E.sector.map(errs)

    def oe(col, df):
        wv = 1 / df.err ** 2; odd = (df.epoch % 2) != 0
        mo = np.sum(wv[odd] * df[col][odd]) / wv[odd].sum(); me = np.sum(wv[~odd] * df[col][~odd]) / wv[~odd].sum()
        mu = np.sum(wv * df[col]) / wv.sum()
        return dict(odd_minus_even=float(mo - me), err=float(np.sqrt(1 / wv[odd].sum() + 1 / wv[~odd].sum())),
                    mean=float(mu), mean_err=float(1 / np.sqrt(wv.sum())), chi2=float(np.sum(wv * (df[col] - mu) ** 2)), n=len(df))
    out["E_own_depths"] = dict(err_ppt=errs, det=oe("det", E), clean=oe("clean", E),
                               rms_change=float(np.std(E.clean - E.det)), max_change=float(np.max(np.abs(E.clean - E.det))),
                               table=E.to_dict(orient="records"))
    # referee depths: odd/even raw and after subtracting the predicted bias (should not be needed if slope ~0)
    Rr = R.copy(); Rr["err"] = Rr.depth_err_cal / FREF * 1e3; Rr["dppt"] = Rr.depth / FREF * 1e3
    Rr["dppt_corr"] = Rr.dppt - Rr.bias_ppt
    out["E_referee_oddeven"] = dict(raw=oe("dppt", Rr), bias_subtracted=oe("dppt_corr", Rr))
    # phot strict depths, raw and bias-subtracted
    Pp = PBs.copy(); Pp["err"] = err_ph; Pp["dppt"] = dep_ph; Pp["dppt_corr"] = dep_ph - Pp.bias_ppt.values
    out["E_phot_oddeven"] = dict(raw=oe("dppt", Pp), bias_subtracted=oe("dppt_corr", Pp))
    print("E", {k: (v if k != "table" else None) for k, v in out["E_own_depths"].items()}, flush=True)
    print("E ref", out["E_referee_oddeven"]); print("E phot", out["E_phot_oddeven"])
    # per-sector mean depth on the cleaned light curve (S103+S104 vs S39+S66)
    grp = {}
    for nm, ss_ in (("S39_S66", [39, 66]), ("S93", [93]), ("S103_S104", [103, 104]), ("S12", [12])):
        k = E.sector.isin(ss_)
        for col in ("det", "clean"):
            wv = 1 / E.err[k] ** 2
            grp[f"{nm}_{col}"] = (float(np.sum(wv * E[col][k]) / wv.sum()), float(1 / np.sqrt(wv.sum())))
    out["E_groups"] = grp
    print(grp)
    jdump(out, f"{ME}/out/r01_rrab_leak.json")


if __name__ == "__main__":
    main()
