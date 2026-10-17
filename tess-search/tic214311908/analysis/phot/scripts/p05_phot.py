"""Deblended difference-image PSF photometry of the Mira, the M3V (TIC 214311908) and the M-dwarf pair (TIC 214311910 A+B).

Per sector:
1. Pixel regression (quality==0 cadences outside the transit windows of the real ephemeris and of any injected
   signal, extreme background/jitter removed, 1-h bins): pixel_j(t) = a_j + b_j (F_ap(t) - <F_ap>) + c_j (B(t) - <B>)
   + d_j (dx(t) - <dx>) + e_j (dy(t) - <dy>), with F_ap the background-subtracted r = 2.5 px sum around the Mira
   (Mira-flux proxy), B the faint-pixel background and dx, dy the field-star pointing offsets (p04_jitter).
   b = the Mira's empirical image (its true, colour-affected PSF); d, e = the scene's empirical jitter response.
2. Difference image of every cadence: D(t) = pixel(t) - a - c (B - <B>) - d (dx - <dx>) - e (dy - <dy>).
3. Per cadence weighted linear fit of D over pixels within RFIT px of the Mira or the M3V:
   D = dF_Mira * b + dF_M3V * PRF(M3V) + dF_pair * PRF(pair) + plane, PRF positions from Gaia + field-star WCS.
   The M dwarfs keep the SPOC PRF model (p02 blur); variants allow a perturbed PRF.
4. Static normalisation: fit of the mean image a with b-image (Mira), PRFs of M3V, pair and every other star within
   7 px (fluxes of others fixed to p02 values), plane background -> mean fluxes F_M3V, F_pair, F_Mira.
Outputs (function photometry returns a dict; main writes data/phot_s{sector}.npz and out/phot_s{sector}.json).
"""
import sys; sys.path.insert(0, "tic214311908/analysis/phot/scripts")
from photlib import *
from p03_scene import build_stars

RFIT = 3.0
QUAL_LOOSE = 2048           # straylight flag only: accepted as 'quality_ok' if the background is not extreme
REG_BKG_MAX = 150.0         # max background excess (e/s/px) of cadences used in the pixel regression


def in_windows(t, epochs_t, half):
    if len(epochs_t) == 0: return np.zeros(len(t), bool)
    e = np.sort(np.asarray(epochs_t))
    i = np.clip(np.searchsorted(e, t), 1, len(e) - 1)
    return np.minimum(np.abs(t - e[i - 1]), np.abs(t - e[i])) < half


def transit_times(t, P, T0):
    k = np.arange(np.floor((np.nanmin(t) - T0) / P) - 1, np.ceil((np.nanmax(t) - T0) / P) + 2)
    return T0 + k * P


def binned(t, Y, width):
    k = np.floor((t - t[0]) / width).astype(int)
    u, inv = np.unique(k, return_inverse=True)
    n = np.bincount(inv)
    out = np.zeros((len(u),) + Y.shape[1:])
    np.add.at(out, inv, Y)
    return out / n.reshape((-1,) + (1,) * (Y.ndim - 1)), n


def setup(sector):
    C = Cutout(sector); W = json.load(open(f"{ROOT}/out/wcs_s{sector}.json"))
    J = np.load(f"{ROOT}/data/jitter_s{sector}.npz")
    g, rBA = build_stars(C, W)
    prf = PRFModel(C); prf.blur = W["blur"]
    return C, W, J, g, prf, rBA


def photometry(sector, S=None, inject=None, extra_mask_times=(), prf_perturb=None, rfit=RFIT, mira_template="empirical",
               free_others=False, verbose=True, controls=None, return_D=False, jitter_order=1):
    """inject: list of dicts {xy:(x,y) or 'b' (Mira empirical image), amp: array (e/s removed per cadence)} applied to the cube.
    prf_perturb: dict(blur=, dx=, dy=) applied to the M-dwarf fitting templates (to test template mismatch).
    mira_template: 'empirical' (b image) or 'prf' (SPOC PRF at the b-image fitted position)."""
    C, W, J, g, prf, rBA = S if S is not None else setup(sector)
    cube = C.flux.copy()
    mira = g[g.name == "Mira"].iloc[0]; m3v = g[g.name == "M3V"].iloc[0]; pair = g[g.name == "pair"].iloc[0]
    if inject:
        for inj in inject:
            if isinstance(inj["img"], str): continue
            cube -= inj["amp"][:, None, None] * inj["img"][None]
    yy, xx = np.mgrid[0:C.ny, 0:C.nx]
    # ---------------- cadence selection
    fin = C.fin & np.isfinite(J["dx"]) & np.isfinite(J["dy"]) & np.isfinite(J["bkg"])
    q0 = fin & (C.q == 0)
    qok = fin & ((C.q == 0) | (C.q == QUAL_LOOSE))
    B = J["bkg"]
    bl = np.median(B[q0]) + 5 * robust_std(B[q0])
    jit_ok = np.ones(len(C.t), bool)
    for k in ("dx", "dy"):
        v = J[k]; jit_ok &= np.abs(v - np.nanmedian(v[q0])) < 6 * robust_std(v[q0]) + 0.02
    win = in_windows(C.t, transit_times(C.t, P_EPH, T0_EPH), 0.1)
    if len(extra_mask_times): win |= in_windows(C.t, extra_mask_times, 0.1)
    t = C.t
    flat = cube.reshape(len(t), -1)
    # background plane from the faintest 15% pixels (by sector median) of every cadence: level and gradients;
    # subtracted from every pixel with unit response (uniform additive scattered light), residual response regressed
    med = np.nanmedian(flat[q0], 0); faint = med < np.percentile(med, 15)
    uf = ((xx.ravel()[faint]) - 10) / 10; vf = ((yy.ravel()[faint]) - 10) / 10
    Af = np.column_stack([np.ones(faint.sum()), uf, vf])
    BG = np.full((len(t), 3), np.nan)
    cf, *_ = np.linalg.lstsq(Af, flat[fin][:, faint].T, rcond=None); BG[fin] = cf.T
    uall = (xx.ravel() - 10) / 10; vall = (yy.ravel() - 10) / 10
    flat = flat - (BG[:, 0:1] + BG[:, 1:2] * uall[None] + BG[:, 2:3] * vall[None])
    bexc = B - np.median(B[q0])
    reg = qok & (bexc < REG_BKG_MAX) & jit_ok & ~win
    # ---------------- step 1: pixel regression (1 h bins)
    ap = (np.hypot(xx - mira.x, yy - mira.y) <= 2.5).ravel()
    Fap = np.full(len(t), np.nan)
    Fap[fin] = flat[fin][:, ap].sum(1)
    mu = dict(F=np.mean(Fap[reg]), B=np.mean(BG[reg, 0]), dx=np.mean(J["dx"][reg]), dy=np.mean(J["dy"][reg]))
    mu.update(dx2=np.mean((J["dx"][reg] - mu["dx"]) ** 2), dy2=np.mean((J["dy"][reg] - mu["dy"]) ** 2),
              dxy=np.mean((J["dx"][reg] - mu["dx"]) * (J["dy"][reg] - mu["dy"])))
    def regressors(sel):
        ddx, ddy = J["dx"][sel] - mu["dx"], J["dy"][sel] - mu["dy"]
        cols = [np.ones(sel.sum()), Fap[sel] - mu["F"], BG[sel, 0] - mu["B"], ddx, ddy]
        if jitter_order >= 2:       # second-order pointing response (curvature of the bright Mira's PSF)
            cols += [ddx ** 2 - mu["dx2"], ddy ** 2 - mu["dy2"], ddx * ddy - mu["dxy"]]
        return np.column_stack(cols)
    Xb, n = binned(t[reg], regressors(reg), 1 / 24); Yb, _ = binned(t[reg], flat[reg], 1 / 24)
    kb = n >= max(2, 0.5 * np.median(n)); Xb, Yb = Xb[kb], Yb[kb]
    coef, *_ = np.linalg.lstsq(Xb, Yb, rcond=None)
    rb = Yb - Xb @ coef
    covu = np.linalg.pinv(Xb.T @ Xb)
    berr = np.sqrt((rb ** 2).sum(0) / (len(Yb) - Xb.shape[1]) * covu[1, 1]).reshape(C.ny, C.nx)
    a_img, b_img = coef[0].reshape(C.ny, C.nx), coef[1].reshape(C.ny, C.nx)
    # ---------------- step 2: difference images for all finite cadences
    D = np.full(flat.shape, np.nan)
    cz = coef.copy(); cz[1] = 0.0
    D[fin] = flat[fin] - regressors(fin) @ cz
    # note: regressors(fin) column 0 is 1 -> subtracts a; Fap column multiplied by zero (Mira variation kept in D)
    # ---------------- step 3: per-cadence fit
    fm = (np.hypot(xx - mira.x, yy - mira.y) <= rfit) | (np.hypot(xx - m3v.x, yy - m3v.y) <= rfit)
    fmr = fm.ravel()
    pp = prf_perturb or {}
    def tmpl(x, y):
        return prf.image(x + pp.get("dx", 0.0), y + pp.get("dy", 0.0), blur=np.hypot(prf.blur, pp.get("blur", 0.0)))
    if mira_template == "empirical":
        TM = b_img.copy()
    else:
        P_ = json.load(open(f"{ROOT}/out/mirapsf_s{sector}.json"))["fits"]["2"]["p"]
        TM = prf.image(mira.x + P_[0], mira.y + P_[1], blur=np.hypot(prf.blur, P_[2]))
        TM = TM / TM.ravel()[ap].sum()
    T3 = tmpl(m3v.x, m3v.y); TP = tmpl(pair.x, pair.y)
    u = (xx - mira.x) / 3; v = (yy - mira.y) / 3
    base = dict(Mira=TM, M3V=T3, pair=TP, bkg0=np.ones_like(u), bkgx=u, bkgy=v)
    # model definitions: B = primary for M3V and Mira (pair not fitted); Bp = pair light curve (M3V not fitted);
    # A = all three free (noisier, strongly correlated templates)
    MODELS = dict(B=["Mira", "M3V", "bkg0", "bkgx", "bkgy"], Bp=["Mira", "pair", "bkg0", "bkgx", "bkgy"],
                  A=["Mira", "M3V", "pair", "bkg0", "bkgx", "bkgy"])
    if free_others:
        for r_ in g.itertuples():
            if r_.name in ("Mira", "M3V", "pair"): continue
            if r_.T < 15.5 and np.hypot(r_.x - mira.x, r_.y - mira.y) < rfit + 1.5:
                base[str(r_.source_id)] = tmpl(r_.x, r_.y)
                for k in MODELS: MODELS[k].append(str(r_.source_id))
    Dg = D[q0][:, fmr]
    sig = 1.4826 * np.nanmedian(np.abs(np.diff(Dg, axis=0)), axis=0) / np.sqrt(2) + 0.5
    Y = (D[fin][:, fmr] / sig).T
    FIT = {}
    for mk, nms in MODELS.items():
        A = np.column_stack([base[n_].ravel()[fmr] for n_ in nms])
        Aw = A / sig[:, None]
        cov = np.linalg.pinv(Aw.T @ Aw)
        F = np.full((len(t), A.shape[1]), np.nan); chi2 = np.full(len(t), np.nan)
        c_ = cov @ Aw.T @ Y
        F[fin] = c_.T
        chi2[fin] = ((Aw @ c_ - Y) ** 2).sum(0) / (fmr.sum() - A.shape[1])
        FIT[mk] = dict(F=F, names=nms, ferr=np.sqrt(np.diag(cov)), chi2=chi2, corr=(cov / np.outer(np.sqrt(np.diag(cov)), np.sqrt(np.diag(cov)))))
    F, names, chi2 = FIT["A"]["F"], FIT["A"]["names"], FIT["A"]["chi2"]
    ferr = FIT["A"]["ferr"]
    # ---------------- optional control stars: per-cadence PRF + plane fit to D within 2.5 px of each control star
    CTRL = {}
    if controls:
        for nm_, (cx_, cy_) in controls.items():
            cm = (np.hypot(xx - cx_, yy - cy_) <= 2.5).ravel()
            uc = ((xx.ravel() - cx_) / 3)[cm]; vc = ((yy.ravel() - cy_) / 3)[cm]
            Ac = np.column_stack([tmpl(cx_, cy_).ravel()[cm], np.ones(cm.sum()), uc, vc])
            sc = 1.4826 * np.nanmedian(np.abs(np.diff(D[q0][:, cm], axis=0)), axis=0) / np.sqrt(2) + 0.5
            Acw = Ac / sc[:, None]; covc = np.linalg.pinv(Acw.T @ Acw)
            Fc = np.full(len(t), np.nan); Fc[fin] = (covc @ Acw.T @ (D[fin][:, cm] / sc).T)[0]
            CTRL[nm_] = Fc
    # ---------------- step 4: static normalisation from the mean image a
    others = np.zeros((C.ny, C.nx)); ff = W["free_fluxes"]
    for r_ in g.itertuples():
        if r_.name in ("Mira", "M3V", "pair"): continue
        if np.hypot(r_.x - mira.x, r_.y - mira.y) > rfit + 7: continue
        others += ff.get(str(int(r_.source_id)), tmag_flux(r_.T) * W["faint_scale"]) * prf.image(r_.x, r_.y)
    sm = np.hypot(xx - mira.x, yy - mira.y) <= rfit + 0.5
    As = np.column_stack([TM[sm], T3[sm], TP[sm], np.ones(sm.sum()), u[sm], v[sm]])
    ws = 1 / (0.02 * np.abs(a_img[sm]) + 5)
    cs, *_ = np.linalg.lstsq(As * ws[:, None], (a_img - others)[sm] * ws, rcond=None)
    rs = (As @ cs - (a_img - others)[sm]) * ws
    static = dict(Mira=float(cs[0]), M3V=float(cs[1]), pair=float(cs[2]), bkg=float(cs[3]), chi2=float(np.sum(rs ** 2)), npix=int(sm.sum()))
    # alternative static fit with the pair fixed at its Gaia-predicted flux
    exp_pair = tmag_flux(pair["T"]); exp_m3v = tmag_flux(m3v["T"])
    As2 = np.column_stack([TM[sm], T3[sm], np.ones(sm.sum()), u[sm], v[sm]])
    cs2, *_ = np.linalg.lstsq(As2 * ws[:, None], (a_img - others - exp_pair * TP)[sm] * ws, rcond=None)
    static["M3V_pairfixed"] = float(cs2[1]); static["Mira_pairfixed"] = float(cs2[0])
    static["expected"] = dict(M3V=float(exp_m3v), pair=float(exp_pair))
    out = dict(t=t, q=C.q, cadn=C.cadn, timecorr=C.timecorr, fin=fin, q0=q0, qok=qok, reg=reg, bkg=B, bkg_lim=bl, jit_ok=jit_ok,
               dx=J["dx"], dy=J["dy"], Fap=Fap, F=F, names=names, ferr=ferr, chi2=chi2, static=static, b_img=b_img, berr=berr,
               a_img=a_img, fitmask=fm, sig=sig, mu=mu, exptime=C.exptime, FIT=FIT, CTRL=CTRL,
               templates=dict(Mira=TM, M3V=T3, pair=TP), pos={k: (float(r_.x), float(r_.y)) for k, r_ in (("Mira", mira), ("M3V", m3v), ("pair", pair))})
    if return_D: out["D"] = D.reshape(len(t), C.ny, C.nx)
    if verbose:
        print(f"S{sector}: reg cadences {reg.sum()} / fin {fin.sum()} (q0 {q0.sum()}, qok {qok.sum()}); fit pixels {fmr.sum()}; "
              f"formal per-cadence err [A] Mira {ferr[0]:.1f} M3V {ferr[1]:.1f} pair {ferr[2]:.1f}; [B] Mira {FIT['B']['ferr'][0]:.1f} M3V {FIT['B']['ferr'][1]:.1f}; [Bp] pair {FIT['Bp']['ferr'][1]:.1f} e/s; static fluxes Mira {cs[0]:.0f} M3V {cs[1]:.0f} pair {cs[2]:.0f} "
              f"(pair fixed {exp_pair:.0f}: M3V {cs2[1]:.0f}); expected M3V {exp_m3v:.0f}; static chi2 {static['chi2']:.0f}/{sm.sum()}; median chi2r {np.nanmedian(chi2[q0]):.2f}", flush=True)
    return out


def save(sector, R, tag=""):
    keep = {k: v for k, v in R.items() if isinstance(v, np.ndarray)}
    for mk, d in R["FIT"].items():
        keep[f"F_{mk}"] = d["F"]; keep[f"chi2_{mk}"] = d["chi2"]; keep[f"ferr_{mk}"] = d["ferr"]; keep[f"names_{mk}"] = np.array(d["names"])
    np.savez_compressed(f"{ROOT}/data/phot_s{sector}{tag}.npz", names=np.array(R["names"]), **keep)
    save_json(dict(sector=sector, names=R["names"], ferr=R["ferr"], static=R["static"], mu=R["mu"], bkg_lim=R["bkg_lim"],
                   models={mk: dict(names=d["names"], ferr=d["ferr"], corr=d["corr"]) for mk, d in R["FIT"].items()},
                   exptime=R["exptime"], pos=R["pos"]), f"{ROOT}/out/phot_s{sector}{tag}.json")


if __name__ == "__main__":
    for s in ([int(a) for a in sys.argv[1:]] or SECTORS):
        R = photometry(s); save(s, R)
