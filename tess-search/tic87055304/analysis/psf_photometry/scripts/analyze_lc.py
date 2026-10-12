"""Transit depths in the PSF light curves of the target and the neighbour (per sector and combined) for one photometry tag.

Series per sector (masked: high-background cadences > median + 6 * 1.48 MAD, chi2 outliers > 5x median):
  t_abs, n_abs : fitted fluxes (e-/s) -> A = flux decrement Delta F (e-/s)
  t_rel, n_rel : fluxes / median -> fractional depth in the PSF-normalised light curve
  sum_rel      : (Ft + Fn) / median -> depth of the blended pair
Errors: max(white, per-transit scatter, control-phase joint, control-phase single-epoch / sqrt(N)) (depthfit.py).
Combination: per sector the vector (dFt, dFn) with control-phase covariance C_s;  model dFt = f S_s, dFn = (1-f) S_s;
  GLS profile likelihood for f (S_s free per sector); chi2 at f = 1 (target hosts) and f = 0 (neighbour hosts).
Optional regressors (--regs): field registration dx, dy and background in the per-window baseline.
"""
import os, sys, json, argparse
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import *
from depthfit import summarize, detrended, depth_fit, WIN


BKG_NMAD = float(os.environ.get("BKG_NMAD", "6"))


def masks(d):
    bk = d["bkg"][:, 0]; med = np.median(bk); mad = 1.4826 * np.median(np.abs(bk - med))
    c2 = d["chi2"]
    return (bk < med + BKG_NMAD * mad) & (c2 < 5 * np.median(c2))


def sector_series(d, js):
    ft = d["flux"][:, int(d["jt"])]; fn = d["flux"][:, int(d["jn"])]
    return ft, fn


def fit_f(per):
    """per: list of (dFt, dFn, cov2x2). Profile likelihood in f with free S per sector."""
    fs = np.linspace(-1.0, 2.0, 3001); chi = []
    for f in fs:
        c = 0
        for a, b, C in per:
            v = np.array([f, 1 - f]); Ci = np.linalg.inv(C); y = np.array([a, b])
            S = (v @ Ci @ y) / (v @ Ci @ v)
            r = y - S * v; c += r @ Ci @ r
        chi.append(c)
    chi = np.array(chi); k = np.argmin(chi)
    inside = fs[chi < chi[k] + 1]
    return dict(f_best=float(fs[k]), f_lo=float(inside.min()), f_hi=float(inside.max()), chi2_min=float(chi[k]),
                chi2_f1=float(np.interp(1.0, fs, chi)), chi2_f0=float(np.interp(0.0, fs, chi)), fs=fs, chi=chi)


def run(tag, sectors=SECTORS, use_regs=False, plot=True, P_=P, t0=T0, label=None, quiet=False):
    res = {}
    per = []
    od = os.path.join(ROOT, "lc", tag)
    if plot:
        fig, axs = plt.subplots(2, len(sectors), figsize=(4 * len(sectors), 7.5), sharey="row")
    for j, s in enumerate(sectors):
        d = np.load(os.path.join(od, f"psf_s{s}.npz"))
        js = json.load(open(os.path.join(ROOT, "static", f"static_s{s}.json")))
        t = d["t"]; ex = float(d["exptime"])
        ft, fn = sector_series(d, js)
        m = masks(d)
        regs = None
        if use_regs:
            regs = [d["reg_field"][0][m], d["reg_field"][1][m], d["bkg"][:, 0][m]]
        tt = t[m]; Ft = ft[m]; Fn = fn[m]
        mt, mn = np.median(Ft), np.median(Fn); msum = np.median(Ft + Fn)
        names = ["t_abs", "n_abs", "t_rel", "n_rel", "sum_rel"]
        r = summarize(tt, [Ft, Fn, Ft / mt, Fn / mn, (Ft + Fn) / msum], ex, names, P_=P_, t0=t0, regs=regs)
        rg = js["ratio_gaia"]
        Fn_cat = msum * rg / (1 + rg); Ft_cat = msum / (1 + rg)
        C = np.array(r["_cov_control"])[:2, :2]
        # inflate covariance to the adopted errors (keep correlation)
        rho = C[0, 1] / np.sqrt(C[0, 0] * C[1, 1])
        et, en = r["t_abs"]["err"], r["n_abs"]["err"]
        Cad = np.array([[et ** 2, rho * et * en], [rho * et * en, en ** 2]])
        per.append((r["t_abs"]["A"], r["n_abs"]["A"], Cad))
        S = r["t_abs"]["A"] + r["n_abs"]["A"]
        # f_target per sector with delta-method error
        g = np.array([r["n_abs"]["A"], -r["t_abs"]["A"]]) / S ** 2
        f_err = float(np.sqrt(g @ Cad @ g))
        res[s] = dict(n_tr=r["t_abs"]["n_tr"], N=int(m.sum()), masked=int((~m).sum()), F_t_med=float(mt), F_n_med=float(mn),
                      ratio_psf=float(mn / mt), ratio_gaia=rg,
                      dF_t=r["t_abs"]["A"], dF_t_err=et, dF_n=r["n_abs"]["A"], dF_n_err=en, rho_tn=float(rho),
                      depth_t_psfnorm=r["t_rel"]["A"], depth_t_psfnorm_err=r["t_rel"]["err"],
                      depth_n_psfnorm=r["n_rel"]["A"], depth_n_psfnorm_err=r["n_rel"]["err"],
                      depth_t_catnorm=float(r["t_abs"]["A"] / Ft_cat), depth_t_catnorm_err=float(et / Ft_cat),
                      depth_n_catnorm=float(r["n_abs"]["A"] / Fn_cat), depth_n_catnorm_err=float(en / Fn_cat),
                      depth_pair=r["sum_rel"]["A"], depth_pair_err=r["sum_rel"]["err"],
                      f_target=float(r["t_abs"]["A"] / S), f_target_err=f_err,
                      err_detail={k: {kk: vv for kk, vv in v.items() if kk.startswith("err") or kk == "mean_control_joint"}
                                  for k, v in r.items() if not k.startswith("_")},
                      A_epoch_t=r["t_abs"].get("A_epoch"), A_epoch_n=r["n_abs"].get("A_epoch"))
        if not quiet:
            q = res[s]
            print(f"S{s} [{tag}{' +regs' if use_regs else ''}] ntr={q['n_tr']} dFt={q['dF_t']:.2f}+-{et:.2f} dFn={q['dF_n']:.2f}+-{en:.2f} "
                  f"(rho {rho:.2f}) f_t={q['f_target']:.3f}+-{f_err:.3f} | depth t(psf) {q['depth_t_psfnorm']*1e3:.2f}+-{q['depth_t_psfnorm_err']*1e3:.2f} ppt "
                  f"n(psf) {q['depth_n_psfnorm']*1e3:.1f}+-{q['depth_n_psfnorm_err']*1e3:.1f} ppt; n(cat) {q['depth_n_catnorm']*1e3:.1f}+-{q['depth_n_catnorm_err']*1e3:.1f}; "
                  f"pair {q['depth_pair']*1e3:.2f}+-{q['depth_pair_err']*1e3:.2f} ppt; Fn/Ft psf {q['ratio_psf']:.3f}", flush=True)
        if plot:
            for row, (y, nm, ref) in enumerate([(Ft / mt, "target", mt), (Fn / mn, "neighbour", mn)]):
                dd = detrended(tt, y, ex, P_, t0, regs=regs)
                ax = axs[row, j]
                if dd is None: continue
                dtk, yd, A, keep = dd
                hrs = dtk * 24
                ax.plot(hrs, yd * 1e3, ".", color="0.75", ms=1.5, zorder=1)
                bins = np.arange(-WIN * 24, WIN * 24 + 0.01, 0.5)
                idx = np.digitize(hrs, bins)
                bx = [hrs[idx == i].mean() for i in range(1, len(bins)) if (idx == i).sum() > 2]
                by = [yd[idx == i].mean() for i in range(1, len(bins)) if (idx == i).sum() > 2]
                be = [yd[idx == i].std() / np.sqrt((idx == i).sum()) for i in range(1, len(bins)) if (idx == i).sum() > 2]
                ax.errorbar(bx, np.array(by) * 1e3, np.array(be) * 1e3, fmt="o", color="k", ms=3, zorder=3)
                tg = np.linspace(-WIN, WIN, 400) + t0
                tm = template(tg, ex, P_, t0)
                dp = r["sum_rel"]["A"]
                if nm == "target":
                    pred_t = dp * msum / mt; pred_n = 0.0
                    ax.plot((tg - t0) * 24, -pred_t * tm * 1e3, "C0-", lw=1.5, label="if target hosts")
                    ax.plot((tg - t0) * 24, 0 * tm, "C3--", lw=1.5, label="if neighbour hosts")
                else:
                    pred_n = dp * msum / mn
                    ax.plot((tg - t0) * 24, 0 * tm, "C0-", lw=1.5, label="if target hosts")
                    ax.plot((tg - t0) * 24, -pred_n * tm * 1e3, "C3--", lw=1.5, label="if neighbour hosts")
                ax.set_title(f"S{s} {nm} PSF flux")
                if j == 0: ax.set_ylabel("relative flux - 1 (ppt)")
                if row == 1: ax.set_xlabel("hours from mid-transit")
                if j == 0 and row == 0: ax.legend(fontsize=7, loc="lower left")
        axs[0, 0].set_ylim(-12, 6); axs[1, 0].set_ylim(-60, 30)
    fit = fit_f(per)
    if plot:
        plt.suptitle(f"PSF light curves folded on P={P_:.6f} d [{label or tag}{' +regs' if use_regs else ''}]; "
                     f"combined f_target = {fit['f_best']:.2f} [{fit['f_lo']:.2f},{fit['f_hi']:.2f}]", fontsize=10)
        plt.tight_layout()
        plt.savefig(os.path.join(ROOT, f"fold_psf_{tag}{'_regs' if use_regs else ''}{f'_bkg{BKG_NMAD:g}' if BKG_NMAD != 6 else ''}.png"), dpi=90)
        plt.close()
    # combined relative depths (inverse-variance weighted)
    def wmean(keyv, keye):
        v = np.array([res[s][keyv] for s in sectors]); e = np.array([res[s][keye] for s in sectors])
        w = 1 / e ** 2; return float(np.sum(w * v) / np.sum(w)), float(1 / np.sqrt(np.sum(w))), float(np.sum(w * (v - np.sum(w * v) / np.sum(w)) ** 2))
    comb = {}
    for k in ["depth_t_psfnorm", "depth_n_psfnorm", "depth_t_catnorm", "depth_n_catnorm", "depth_pair", "dF_t", "dF_n"]:
        m_, e_, c2 = wmean(k, k + "_err"); comb[k] = dict(value=m_, err=e_, chi2_between_sectors=c2, dof=len(sectors) - 1)
    comb["f_target"] = {k: v for k, v in fit.items() if k not in ("fs", "chi")}
    comb["delta_chi2_target_minus_neighbour"] = fit["chi2_f1"] - fit["chi2_f0"]
    if not quiet:
        print(f"COMBINED [{tag}{' +regs' if use_regs else ''}]: f_target = {fit['f_best']:.3f} [{fit['f_lo']:.3f}, {fit['f_hi']:.3f}] (1-sigma); "
              f"chi2(f=1) {fit['chi2_f1']:.1f}, chi2(f=0) {fit['chi2_f0']:.1f}, chi2min {fit['chi2_min']:.1f} (dof {len(sectors)-1}); "
              f"depth_t(psf) {comb['depth_t_psfnorm']['value']*1e3:.2f}+-{comb['depth_t_psfnorm']['err']*1e3:.2f} ppt; "
              f"depth_n(psf) {comb['depth_n_psfnorm']['value']*1e3:.1f}+-{comb['depth_n_psfnorm']['err']*1e3:.1f} ppt; "
              f"depth_n(cat) {comb['depth_n_catnorm']['value']*1e3:.1f}+-{comb['depth_n_catnorm']['err']*1e3:.1f} ppt; "
              f"pair {comb['depth_pair']['value']*1e3:.2f}+-{comb['depth_pair']['err']*1e3:.2f} ppt", flush=True)
    out = dict(tag=tag, regs=use_regs, sectors={str(k): v for k, v in res.items()}, combined=comb)
    sfx = ("_regs" if use_regs else "") + (f"_bkg{BKG_NMAD:g}" if BKG_NMAD != 6 else "")
    out["bkg_nmad"] = BKG_NMAD
    json.dump(out, open(os.path.join(od, f"depths{sfx}.json"), "w"), indent=1, default=float)
    return out, fit


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="nominal")
    ap.add_argument("--regs", action="store_true")
    ap.add_argument("--sectors", default="13,39,66,93,104")
    ap.add_argument("--noplot", action="store_true")
    a = ap.parse_args()
    run(a.tag, [int(v) for v in a.sectors.split(",")], a.regs, plot=not a.noplot)
