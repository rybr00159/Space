"""Transit-depth estimator with red-noise-aware uncertainties.

depth_fit(): within +-WIN days of each transit (P, t0), model  y = c0_e + c1_e dt + c2_e dt^2 - A * tmpl(t)  (per-epoch
quadratic baseline, shared amplitude A; tmpl = unit-depth batman template from common.template).  y may be a normalised
flux (A = fractional depth) or an absolute flux in e-/s (A = flux decrement).  Epochs need >= 50% of the in-transit
cadences and >= 25% of the out-of-transit window present.
control(): the same estimator at fake epochs t0 + phi*P (phi grid avoiding the real transit +- window and phase 0.5),
both the joint estimator (std over phi = sector-level red-noise error) and single-epoch amplitudes (std/sqrt(N_tr)).
"""
import numpy as np
from common import template, phase_dt, P, T0, DUR

WIN = 0.45   # days (half window), ~2.4 transit durations


def design(t, y, exptime, P_=P, t0=T0, win=WIN, poly=2, min_in=0.5, min_out=0.25, mask=None, regs=None):
    dt, ep = phase_dt(t, P_, t0)
    sel = (np.abs(dt) < win) & np.isfinite(y)
    if mask is not None: sel &= mask
    cad = exptime / 86400
    n_in_exp = DUR / cad; n_out_exp = (2 * win - DUR) / cad
    eps = []
    for e in np.unique(ep[sel]):
        m = sel & (ep == e)
        nin = np.sum(np.abs(dt[m]) < DUR / 2); nout = np.sum(np.abs(dt[m]) > DUR / 2 + 0.02)
        if nin >= min_in * n_in_exp and nout >= min_out * n_out_exp: eps.append(e)
    if not eps: return None
    keep = sel & np.isin(ep, eps)
    # transit shape always from the true orbital period (P_ may be 2P for odd/even)
    tm = template(t0 + dt[keep], exptime, P, t0)
    cols = [-tm]
    epk = ep[keep]; dtk = dt[keep]
    for e in eps:
        m = (epk == e).astype(float)
        for p in range(poly + 1): cols.append(m * (dtk / win) ** p)
    if regs is not None:
        for rg in regs:
            v = rg[keep]; cols.append((v - np.mean(v)) / (np.std(v) + 1e-12))
    A = np.vstack(cols).T
    return A, y[keep], keep, np.array(eps), epk, dtk, tm


def depth_fit(t, y, exptime, P_=P, t0=T0, win=WIN, poly=2, sig=None, mask=None, per_epoch=False, regs=None):
    r = design(t, y, exptime, P_, t0, win, poly, mask=mask, regs=regs)
    if r is None: return dict(A=np.nan, err=np.nan, n=0)
    A, yy, keep, eps, epk, dtk, tm = r
    c, *_ = np.linalg.lstsq(A, yy, rcond=None)
    res = yy - A @ c
    s = np.std(res) if sig is None else sig
    cov = np.linalg.pinv(A.T @ A) * s ** 2
    out = dict(A=float(c[0]), err=float(np.sqrt(cov[0, 0])), n=len(eps), rms=float(np.std(res)), epochs=eps.tolist())
    if per_epoch:
        Ae, Ee = [], []
        for k, e in enumerate(eps):
            m = epk == e
            cols = [-tm[m]] + [(dtk[m] / win) ** p for p in range(poly + 1)]
            B = np.vstack(cols).T
            ce, *_ = np.linalg.lstsq(B, yy[m], rcond=None)
            re = yy[m] - B @ ce
            Ae.append(ce[0]); Ee.append(np.sqrt(np.linalg.pinv(B.T @ B)[0, 0]) * np.std(re))
        out["A_epoch"] = np.array(Ae); out["err_epoch"] = np.array(Ee)
    return out


def control(t, ys, exptime, P_=P, t0=T0, win=WIN, poly=2, phis=None, mask=None, exclude_secondary=True, regs=None):
    """ys: list of series (fitted jointly at the same fake epochs so their covariance can be estimated).
    Returns array (nphi, nseries) of joint amplitudes and list of single-epoch amplitude arrays."""
    if phis is None:
        lo = (win + DUR / 2 + 0.05) / P_
        phis = np.arange(lo, 1 - lo + 1e-9, 0.01)
        if exclude_secondary:
            phis = phis[np.abs(phis - 0.5) > (DUR / 2 + 0.05) / P_ + 0.0]
    joint, single = [], []
    for ph in phis:
        row, srow = [], []
        for y in ys:
            r = depth_fit(t, y, exptime, P_, t0 + ph * P_, win, poly, mask=mask, per_epoch=True, regs=regs)
            row.append(r["A"]); srow.append(r.get("A_epoch", np.array([])))
        joint.append(row); single.append(srow)
    return np.array(phis), np.array(joint), single


def summarize(t, ys, exptime, names, P_=P, t0=T0, win=WIN, poly=2, mask=None, do_control=True, regs=None):
    out = {}
    real = [depth_fit(t, y, exptime, P_, t0, win, poly, mask=mask, per_epoch=True, regs=regs) for y in ys]
    if do_control:
        phis, J, S = control(t, ys, exptime, P_, t0, win, poly, mask=mask, regs=regs)
        ok = np.all(np.isfinite(J), axis=1)
        covJ = np.cov(J[ok].T) if ok.sum() > 3 else None
    for k, (nm, r) in enumerate(zip(names, real)):
        d = dict(A=r["A"], err_white=r["err"], n_tr=r["n"])
        if "A_epoch" in r and r["n"] > 1:
            d["err_epoch_scatter"] = float(np.std(r["A_epoch"], ddof=1) / np.sqrt(r["n"]))
            d["A_epoch"] = r["A_epoch"].tolist()
        if do_control and ok.sum() > 3:
            d["err_control_joint"] = float(np.std(J[ok, k], ddof=1))
            d["mean_control_joint"] = float(np.mean(J[ok, k]))
            sing = np.concatenate([S[i][k] for i in range(len(S)) if ok[i]])
            d["err_control_single"] = float(np.std(sing, ddof=1) / np.sqrt(max(r["n"], 1)))
            d["n_control"] = int(ok.sum())
            # adopted: max of white, epoch scatter, control estimates
            d["err"] = float(np.nanmax([d["err_white"], d.get("err_epoch_scatter", 0), d["err_control_joint"], d["err_control_single"]]))
        else:
            d["err"] = float(np.nanmax([d["err_white"], d.get("err_epoch_scatter", 0)]))
        out[nm] = d
    if do_control and ok.sum() > 3:
        out["_cov_control"] = covJ.tolist(); out["_names"] = names
    return out


def detrended(t, y, exptime, P_=P, t0=T0, win=WIN, poly=2, mask=None, regs=None):
    """Joint fit, then return dt and y minus the baseline (poly + regressors) divided by the mean baseline level,
    i.e. a 'relative' series with the transit kept (for plotting)."""
    r = design(t, y, exptime, P_, t0, win, poly, mask=mask, regs=regs)
    if r is None: return None
    A, yy, keep, eps, epk, dtk, tm = r
    c, *_ = np.linalg.lstsq(A, yy, rcond=None)
    base = A[:, 1:] @ c[1:]
    return dtk, yy - base, c[0], keep
