"""c07: depth versus epoch on the cleaned light curves, and its tests.

1. Per-sector depths (fixed shape of the c04 'single' INV/REF posterior medians, per-window linear baselines, beta-scaled
   errors = fake-epoch calibrated) for INV and REF; group comparisons; omnibus chi2; linear trend in time.
   Look-elsewhere: 200000 null draws of the six sector depths (common mean, the measured per-sector errors); the
   observed maximum |z| is compared with the null maximum over (A) the 5 contiguous early/late splits, (B) all 31
   bipartitions, (C) all pairs of disjoint sector subsets with >= 2 sectors each (includes the post-hoc S39+S66 vs
   S103+S104).
2. Pixel-level injections (c05): recovered/injected depth per sector, pipeline and injected-image variant.
3. Per-transit and per-sector correlates of the depth:
   (a) Mira flux in the window (aperture flux Fap of the INV pipeline, e-/s), its slope, and the pulsation phase from the
       VSX ephemeris (ASASSN-V J172354.64-511925.3: P = 294 d, epoch of maximum JD 2457654.0);
   (b) Mira PRF-model residual near the M3V (static mean image fitted with a free-position SPOC PRF for the Mira, fixed
       M3V/pair PRFs and a plane; residual rms within 1.5 px of the M3V and its projection on the M3V PRF);
       and the correlation between the M3V PRF and the Mira's empirical image over the fit pixels (template degeneracy);
   (c) M3V pixel phase, camera/CCD, M3V-Mira separation and position angle in pixels;
   (d) background (faint-pixel level) in the window;
   (e) expected per-sector M3V flux in PRF units (field-star PRF/TIC scale per sector, referee_phot r17) -> depth
       renormalised per sector.
   Weighted regressions of per-transit depth on each continuous variable (slope / error and permutation p-value), and
   of per-sector depths on sector-level variables.
Outputs: out/depth_epoch.json, out/per_transit_correlates.csv, out/inject_summary.json
"""
import sys, itertools
sys.path.insert(0, "tic214311908/analysis/clean/scripts")
from cleanlib import *
from c03_windows import load_windows
from c04_global import load
from tfit import Model, mid_depth, fit_depth_fixed_shape
import pandas as pd
from scipy import stats

rng = np.random.default_rng(707)
YEAR = 365.25
SEC_T = {}


def shape_from(pipe):
    g = json.load(open(f"{CL}/out/global_{pipe}_single.json"))
    m = g["map"]
    k, b, aR, u1, u2 = (g[x]["med"] for x in ("k", "b", "aR", "u1", "u2"))
    return dict(T0=g["T0"]["med"], P=g["P"]["med"], k=k, b=b, aR=aR, u1=u1, u2=u2)


def sector_depths(pipe, drop_epochs=()):
    W = load(pipe)
    if drop_epochs:
        keep = [i for i in range(W.nw) if W.epoch[i] not in drop_epochs]; W = W.subset(keep)     # sigma already beta-scaled
    sh = shape_from(pipe if pipe in ("inv", "ref") else "ref")
    D0 = mid_depth(sh["k"], sh["b"], sh["aR"], sh["u1"], sh["u2"], sh["P"])
    M = Model(W.t, W.ex)
    frac = (1 - M.flux(sh["T0"], sh["P"], sh["k"], sh["b"], sh["aR"], sh["u1"], sh["u2"])) / D0
    dep, err, chi2, ug = fit_depth_fixed_shape(W, frac, groups=W.sector)
    dall, eall, _, _ = fit_depth_fixed_shape(W, frac)
    tmid = {s: float(np.mean(W.t[W.sector[W.w] == s])) for s in ug}
    return {int(s): (float(d * 1e3), float(e * 1e3)) for s, d, e in zip(ug, dep, err)}, (float(dall[0] * 1e3), float(eall[0] * 1e3)), tmid, W, frac


def group_z(d, A, B):
    def wm(g):
        v = np.array([d[s][0] for s in g]); e = np.array([d[s][1] for s in g]); w = 1 / e ** 2
        return np.sum(w * v) / np.sum(w), 1 / np.sqrt(np.sum(w))
    a, ea = wm(A); b, eb = wm(B)
    return (b - a) / np.hypot(ea, eb), b - a, np.hypot(ea, eb)


def lee(d, observed_pairs):
    secs = sorted(d)
    e = np.array([d[s][1] for s in secs])
    fams = {}
    fams["A_contiguous"] = [(secs[:i], secs[i:]) for i in range(1, len(secs))]
    bip = []
    for r in range(1, len(secs)):
        for A in itertools.combinations(secs, r):
            B = tuple(s for s in secs if s not in A)
            if (B, A) not in bip: bip.append((A, B))
    fams["B_bipartitions"] = bip
    pairs = []
    for ra in range(2, len(secs) - 1):
        for A in itertools.combinations(secs, ra):
            rest = [s for s in secs if s not in A]
            for rb in range(2, len(rest) + 1):
                for B in itertools.combinations(rest, rb):
                    if (B, A) not in pairs: pairs.append((A, B))
    fams["C_disjoint_pairs_ge2"] = pairs
    N = 200000
    sims = rng.normal(0, 1, (N, len(secs))) * e[None]
    out = {}
    w = 1 / e ** 2
    for name, fam in fams.items():
        obs = max(abs(group_z(d, A, B)[0]) for A, B in fam)
        # vectorised null maxima
        mx = np.zeros(N)
        for A, B in fam:
            ia = [secs.index(s) for s in A]; ib = [secs.index(s) for s in B]
            ma = (sims[:, ia] * w[ia]).sum(1) / w[ia].sum(); mb = (sims[:, ib] * w[ib]).sum(1) / w[ib].sum()
            z = np.abs(mb - ma) / np.sqrt(1 / w[ia].sum() + 1 / w[ib].sum())
            mx = np.maximum(mx, z)
        out[name] = dict(n_tests=len(fam), observed_max_z=float(obs), p_global=float(np.mean(mx >= obs)))
        for lab, (A, B) in observed_pairs.items():
            z = abs(group_z(d, A, B)[0])
            out[name][f"p_global_of_{lab}"] = float(np.mean(mx >= z))
    return out


def inject_summary():
    rows = []
    for s in SECTORS:
        d = pd.read_csv(f"{CL}/out/inject_pix_s{s}.csv")
        for var, g in d.groupby("variant"):
            r = dict(sector=s, variant=var, n=len(g))
            for pipe, col in (("inv", "inv_raw_depth_fix_resp"), ("inv_det", "inv_det_depth_fix_resp"), ("ref", "ref_depth_fix_resp")):
                v = g[col].dropna() / 15.7
                r[f"{pipe}_resp_mean"] = float(v.mean()); r[f"{pipe}_resp_std"] = float(v.std()); r[f"{pipe}_n"] = int(len(v))
            for pipe, col in (("inv", "inv_raw_depth_fix"), ("ref", "ref_depth_fix")):
                v = g[col].dropna()
                r[f"{pipe}_abs_mean_ppt"] = float(v.mean()); r[f"{pipe}_abs_err_ppt"] = float(v.std() / np.sqrt(max(len(v), 1)))
            rows.append(r)
    T = pd.DataFrame(rows)
    summ = {}
    for var in T.variant.unique():
        x = T[T.variant == var]
        summ[var] = {p: dict(per_sector={int(s): round(float(v), 4) for s, v in zip(x.sector, x[f"{p}_resp_mean"])},
                             spread_max_minus_min=float(x[f"{p}_resp_mean"].max() - x[f"{p}_resp_mean"].min()),
                             late_over_early=float(x[x.sector.isin([103, 104])][f"{p}_resp_mean"].mean() / x[x.sector.isin([39, 66])][f"{p}_resp_mean"].mean()))
                     for p in ("inv", "inv_det", "ref")}
    # most adverse combination: each sector takes any variant (template error differs by sector)
    adv = {}
    for p in ("inv", "ref"):
        P_ = T.pivot(index="sector", columns="variant", values=f"{p}_resp_mean")
        late = P_.loc[[103, 104]]; early = P_.loc[[39, 66]]
        adv[p] = dict(max_late_over_early=float(late.max(1).mean() / early.min(1).mean()),
                      min_late_over_early=float(late.min(1).mean() / early.max(1).mean()))
    # absolute recovered depth from the noisy injected runs (matched template), per sector
    absd = {}
    for p, var in (("inv", "invprf"), ("ref", "refprf")):
        x = T[T.variant == var]
        absd[p] = {int(s): (round(m, 2), round(e, 2)) for s, m, e in zip(x.sector, x[f"{p}_abs_mean_ppt"], x[f"{p}_abs_err_ppt"])}
    T.to_csv(f"{CL}/out/inject_summary.csv", index=False, float_format="%.4f")
    return dict(by_variant=summ, adverse_combination=adv, absolute_recovered_matched=absd)


def static_residual(sector):
    """Static mean image around the Mira: Mira = SPOC PRF (free x, y, flux; + optional blur), M3V/pair PRFs at fixed
    fluxes (979 x PRF-scale, pair from TIC x 1.107), plane; residual near the M3V. Also corr(M3V PRF, Mira b) over the
    model-B fit pixels."""
    from scipy.optimize import least_squares
    S = setup(sector, "all"); C = S["C"]; prf = S["prf"]; g = S["g"]
    R = photometry(sector, S, fit_models=("B",))
    a = R["a_img"]
    mira = g[g.name == "Mira"].iloc[0]; m3 = g[g.name == "M3V"].iloc[0]; pr = g[g.name == "pair"].iloc[0]
    yy, xx = np.mgrid[0:C.ny, 0:C.nx]
    reg = np.hypot(xx - mira.x, yy - mira.y) <= 3.5
    T3 = prf.image(m3.x, m3.y); TP = prf.image(pr.x, pr.y)
    others = np.zeros_like(a)
    ff = S["W"]["free_fluxes"]
    for r_ in g.itertuples():
        if r_.name in ("Mira", "M3V", "pair"): continue
        if np.hypot(r_.x - mira.x, r_.y - mira.y) > 9: continue
        others += ff.get(str(int(r_.source_id)), tmag_flux(r_.T) * S["W"]["faint_scale"]) * prf.image(r_.x, r_.y)
    fixed = 979.0 * T3 + 1.107 * (tmag_flux(13.3749) + tmag_flux(14.56)) * TP + others
    u = (xx - mira.x) / 3; v = (yy - mira.y) / 3
    wgt = 1 / np.sqrt(np.abs(a) + 25.0)
    def res(p):
        dx, dy, bl = p
        Tm = prf.image(mira.x + dx, mira.y + dy, blur=abs(bl))
        A = np.column_stack([Tm[reg], np.ones(reg.sum()), u[reg], v[reg]]) * wgt[reg][:, None]
        y = (a - fixed)[reg] * wgt[reg]
        c, *_ = np.linalg.lstsq(A, y, rcond=None)
        return y - A @ c
    best = None
    for b0 in (0.0, 0.3):
        q = least_squares(res, [0.0, 0.0, b0], diff_step=[0.01, 0.01, 0.01], bounds=([-1.5, -1.5, 0], [1.5, 1.5, 2.0]))
        if best is None or q.cost < best.cost: best = q
    dx, dy, bl = best.x
    Tm = prf.image(mira.x + dx, mira.y + dy, blur=abs(bl))
    A = np.column_stack([Tm[reg], np.ones(reg.sum()), u[reg], v[reg]]) * wgt[reg][:, None]
    c, *_ = np.linalg.lstsq(A, (a - fixed)[reg] * wgt[reg], rcond=None)
    model = fixed + c[0] * Tm + c[1] + c[2] * u + c[3] * v
    resid = a - model
    near = np.hypot(xx - m3.x, yy - m3.y) <= 1.5
    proj = float(np.sum(resid[reg] * T3[reg]) / np.sum(T3[reg] ** 2))     # e-/s of 'M3V-shaped' residual
    fm = R["fitmask"]; b = R["b_img"]
    corr = float(np.corrcoef(T3[fm], b[fm])[0, 1])
    sepx, sepy = m3.x - mira.x, m3.y - mira.y
    return dict(sector=sector, camera=C.camera, ccd=C.ccd, m3v_x=float(m3.x), m3v_y=float(m3.y),
                pix_phase_x=float((m3.x + 0.5) % 1), pix_phase_y=float((m3.y + 0.5) % 1),
                sep_px=float(np.hypot(sepx, sepy)), pa_deg=float(np.degrees(np.arctan2(sepy, sepx))),
                mira_fit_offset_px=[float(dx), float(dy)], mira_blur_px=float(abs(bl)), mira_flux_prf=float(c[0]),
                resid_rms_near_m3v=float(np.sqrt(np.mean(resid[near] ** 2))), resid_proj_m3v_e_s=proj,
                resid_proj_over_m3v=proj / 979.0, corr_m3vprf_mirab=corr,
                mira_mean_aperture_flux=float(R["mu"]["F"]))


def correlates(pipe, D, shape):
    rows = []
    for s in SECTORS:
        z = np.load(f"{CL}/data/lcinv_all_s{s}.npz")
        t = z["t_tdb"]; Fap = z["Fap"]; bkg = z["bkg"]; dx = z["dx"]; dy = z["dy"]; q = z["quality_ok"].astype(bool)
        for _, r in D[(D.pipeline == pipe) & (D.sector == s)].iterrows():
            m = q & (np.abs(t - r.tc_pred) < 0.3) & np.isfinite(Fap)
            mm = q & (np.abs(t - r.tc_pred) < 0.05)
            if m.sum() < 5: continue
            sl = np.polyfit(t[m] - r.tc_pred, Fap[m], 1)[0]
            phase = (((r.tc_pred + 2457000.0) - 2457654.0) / 294.0) % 1
            rows.append(dict(pipeline=pipe, sector=s, epoch=int(r.epoch), depth=r.depth_ppt, err=r.depth_err_cal, tc=r.tc_pred,
                             mira_flux=float(np.median(Fap[m])), mira_slope_per_d=float(sl), mira_phase=float(phase),
                             bkg=float(np.median(bkg[mm]) if mm.sum() else np.median(bkg[m])), bkg_max=float(np.max(bkg[m])),
                             jitter_rms=float(np.std(np.hypot(dx[m] - np.median(dx[m]), dy[m] - np.median(dy[m]))))))
    return pd.DataFrame(rows)


def wreg(x, y, e, nperm=20000):
    w = 1 / e ** 2
    X = np.column_stack([np.ones(len(x)), x - np.average(x, weights=w)])
    A = X * np.sqrt(w)[:, None]; c, *_ = np.linalg.lstsq(A, y * np.sqrt(w), rcond=None)
    cov = np.linalg.inv(A.T @ A); r = (y - X @ c) * np.sqrt(w); chi2 = float(r @ r)
    sc = max(1.0, np.sqrt(chi2 / max(len(x) - 2, 1)))
    slope, se = c[1], np.sqrt(cov[1, 1]) * sc
    perm = []
    for _ in range(nperm):
        xp = rng.permutation(x); Xp = np.column_stack([np.ones(len(x)), xp - np.average(xp, weights=w)]) * np.sqrt(w)[:, None]
        cp, *_ = np.linalg.lstsq(Xp, y * np.sqrt(w), rcond=None); perm.append(cp[1])
    perm = np.array(perm)
    return dict(slope=float(slope), slope_err=float(se), z=float(slope / se), p_perm=float(np.mean(np.abs(perm) >= abs(slope))),
                spearman_r=float(stats.spearmanr(x, y).correlation), spearman_p=float(stats.spearmanr(x, y).pvalue), n=int(len(x)))


def main():
    out = {}
    D = pd.read_csv(f"{CL}/out/per_transit_clean.csv")
    pairs = {"S39S66_vs_S103S104": ((39, 66), (103, 104)), "early_vs_late": ((12, 39, 66), (93, 103, 104)),
             "S39S66_vs_S93S103S104": ((39, 66), (93, 103, 104))}
    for pipe, drop in (("inv", ()), ("ref", ()), ("ref", (58,)), ("reforig", ())):
        lab = pipe + ("_noE58" if drop else "")
        d, dall, tmid, W, frac = sector_depths(pipe, drop)
        SEC_T.update(tmid)
        secs = sorted(d)
        v = np.array([d[s][0] for s in secs]); e = np.array([d[s][1] for s in secs])
        mu = np.sum(v / e ** 2) / np.sum(1 / e ** 2); chi2 = float(np.sum(((v - mu) / e) ** 2))
        tt = np.array([tmid[s] for s in secs]) / YEAR
        A = np.column_stack([np.ones(len(tt)), tt - tt.mean()]) / e[:, None]
        c, *_ = np.linalg.lstsq(A, v / e, rcond=None); cov = np.linalg.inv(A.T @ A)
        grp = {k: dict(zip(("z", "diff_ppt", "err_ppt"), map(float, group_z(d, *p)))) for k, p in pairs.items()}
        out[lab] = dict(per_sector=d, all=dall, chi2=chi2, dof=len(secs) - 1, p_chi2=float(stats.chi2.sf(chi2, len(secs) - 1)),
                        trend_ppt_per_yr=[float(c[1]), float(np.sqrt(cov[1, 1]))], trend_z=float(c[1] / np.sqrt(cov[1, 1])),
                        groups=grp, look_elsewhere=lee(d, pairs))
        print(f"{lab}: sectors " + ", ".join(f"S{s} {d[s][0]:.1f}+-{d[s][1]:.1f}" for s in secs) + f"; all {dall[0]:.2f}+-{dall[1]:.2f}; chi2 {chi2:.1f}/5 "
              f"(p {stats.chi2.sf(chi2, 5):.3f}); trend {c[1]:+.2f}+-{np.sqrt(cov[1,1]):.2f} ppt/yr; " +
              "; ".join(f"{k} {g['diff_ppt']:+.1f}+-{g['err_ppt']:.1f} (z {g['z']:.2f})" for k, g in grp.items()), flush=True)
        print("   LEE:", {k: (v_["n_tests"], round(v_["observed_max_z"], 2), round(v_["p_global"], 3), round(v_["p_global_of_S39S66_vs_S103S104"], 3)) for k, v_ in out[lab]["look_elsewhere"].items()}, flush=True)
    out["injection"] = inject_summary()
    print("injection matched INV per sector:", out["injection"]["by_variant"]["invprf"]["inv"]["per_sector"], flush=True)
    print("injection matched REF per sector:", out["injection"]["by_variant"]["refprf"]["ref"]["per_sector"], flush=True)
    print("adverse:", out["injection"]["adverse_combination"], flush=True)
    # sector-level geometry / residuals
    geo = [static_residual(s) for s in SECTORS]
    out["sector_geometry"] = geo
    for g_ in geo:
        print({k: (round(v_, 3) if isinstance(v_, float) else v_) for k, v_ in g_.items()}, flush=True)
    # (e) per-sector PRF/TIC scale
    inv_scale = json.load(open(f"{RF}/out/invscale.json"))["per_sector"]
    sc = {int(s): v_["inv_over_tic_median"] for s, v_ in inv_scale.items()}
    msc = np.mean(list(sc.values()))
    out["renormalised_per_sector"] = {}
    for pipe in ("inv", "ref"):
        d = out[pipe]["per_sector"]
        dn = {s: (d[s][0] * msc / sc[s], d[s][1] * msc / sc[s]) for s in d}
        out["renormalised_per_sector"][pipe] = dict(depths=dn, groups={k: dict(zip(("z", "diff_ppt", "err_ppt"), map(float, group_z(dn, *p)))) for k, p in pairs.items()})
    print("renormalised:", {p: {k: round(v_["z"], 2) for k, v_ in x["groups"].items()} for p, x in out["renormalised_per_sector"].items()}, flush=True)
    # sector-level regressions (per-sector depths vs geometry variables)
    out["sector_regressions"] = {}
    for pipe in ("inv", "ref"):
        d = out[pipe]["per_sector"]; secs = sorted(d)
        y = np.array([d[s][0] for s in secs]); e = np.array([d[s][1] for s in secs])
        res = {}
        G = {g_["sector"]: g_ for g_ in geo}
        for var in ("pix_phase_x", "pix_phase_y", "sep_px", "pa_deg", "resid_rms_near_m3v", "resid_proj_over_m3v", "corr_m3vprf_mirab",
                    "mira_mean_aperture_flux", "mira_blur_px"):
            x = np.array([G[s][var] for s in secs], float)
            res[var] = wreg(x, y, e, nperm=5000)
        res["prf_scale"] = wreg(np.array([sc[s] for s in secs]), y, e, nperm=5000)
        out["sector_regressions"][pipe] = res
        print(pipe, "sector-level:", {k: (round(v_["z"], 2), round(v_["p_perm"], 3)) for k, v_ in res.items()}, flush=True)
    # per-transit correlates
    allc = []
    out["transit_regressions"] = {}
    for pipe in ("inv", "ref"):
        C = correlates(pipe, D, None)
        allc.append(C)
        res = {}
        for var in ("mira_flux", "mira_slope_per_d", "bkg", "bkg_max", "jitter_rms", "tc"):
            res[var] = wreg(C[var].values, C.depth.values, C.err.values)
        ph = C.mira_phase.values
        res["mira_phase_sincos"] = dict(sin=wreg(np.sin(2 * np.pi * ph), C.depth.values, C.err.values),
                                        cos=wreg(np.cos(2 * np.pi * ph), C.depth.values, C.err.values))
        out["transit_regressions"][pipe] = res
        print(pipe, "per-transit:", {k: (round(v_["z"], 2), round(v_["p_perm"], 3)) for k, v_ in res.items() if "z" in v_}, flush=True)
    pd.concat(allc).to_csv(f"{CL}/out/per_transit_correlates.csv", index=False, float_format="%.5f")
    jd(out, f"{CL}/out/depth_epoch.json")


if __name__ == "__main__":
    main()
