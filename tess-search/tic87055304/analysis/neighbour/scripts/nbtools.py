"""Shared helpers for the neighbour-host analysis (neighbour/).

template(): unit-depth batman template (shape defaults to the neighbour-host shape: k 0.141, b 0.64, a/R* 6.25,
            u = (0.35, 0.25); exposure-integrated).
depth_fit(): joint fit within +-WIN d of each event: per-epoch quadratic baseline + shared amplitude A of the template.
control_err(): same estimator at fake phases (avoiding the transit and phase 0.5) -> red-noise-aware error.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

ROOT = "/home/user/work/tic87/neighbour"
P, T0, DUR = 4.136195, 1658.41325, 0.1875
SECTORS = [13, 39, 66, 93, 104]
WIN = 0.45
SHAPE_NB = dict(k=0.141, b=0.64, a=6.25, u1=0.35, u2=0.25)
SHAPE_TG = dict(k=0.0668, b=0.85, a=4.30, u1=0.30, u2=0.25)


def template(t, exptime, P_=P, t0=T0, shape=None):
    import batman
    sh = SHAPE_NB if shape is None else shape
    p = batman.TransitParams()
    p.t0, p.per, p.rp, p.a = t0, P_, sh["k"], sh["a"]
    p.inc = np.degrees(np.arccos(sh["b"] / sh["a"])); p.ecc, p.w = 0.0, 90.0
    p.u, p.limb_dark = [sh["u1"], sh["u2"]], "quadratic"
    ss = max(1, int(round(exptime / 60)))
    lc = batman.TransitModel(p, np.asarray(t, float), supersample_factor=ss, exp_time=exptime / 86400).light_curve(p)
    lc0 = batman.TransitModel(p, np.array([t0])).light_curve(p)[0]
    return (1 - lc) / (1 - lc0)


def phase_dt(t, P_=P, t0=T0):
    ep = np.round((t - t0) / P_)
    return t - t0 - ep * P_, ep


def design(t, y, exptime, P_=P, t0=T0, win=WIN, poly=2, min_in=0.5, min_out=0.25, shape=None, Ptmpl=None, epochs=None):
    dt, ep = phase_dt(t, P_, t0)
    sel = (np.abs(dt) < win) & np.isfinite(y)
    cad = exptime / 86400
    n_in_exp = DUR / cad; n_out_exp = (2 * win - DUR) / cad
    eps = []
    for e in np.unique(ep[sel]):
        if epochs is not None and e not in epochs: continue
        m = sel & (ep == e)
        nin = np.sum(np.abs(dt[m]) < DUR / 2); nout = np.sum(np.abs(dt[m]) > DUR / 2 + 0.02)
        if nin >= min_in * n_in_exp and nout >= min_out * n_out_exp: eps.append(e)
    if not eps: return None
    keep = sel & np.isin(ep, eps)
    tm = template(t0 + dt[keep], exptime, P if Ptmpl is None else Ptmpl, t0, shape)
    cols = [-tm]
    epk = ep[keep]; dtk = dt[keep]
    for e in eps:
        m = (epk == e).astype(float)
        for p in range(poly + 1): cols.append(m * (dtk / win) ** p)
    return np.vstack(cols).T, y[keep], keep, np.array(eps), epk, dtk, tm


def depth_fit(t, y, exptime, P_=P, t0=T0, **kw):
    r = design(t, y, exptime, P_, t0, **kw)
    if r is None: return dict(A=np.nan, err=np.nan, n=0, epochs=[])
    A, yy, keep, eps, epk, dtk, tm = r
    c, *_ = np.linalg.lstsq(A, yy, rcond=None)
    res = yy - A @ c
    s = np.std(res)
    cov = np.linalg.pinv(A.T @ A) * s ** 2
    # per-epoch amplitudes (for the scatter-based error)
    Ae = []
    for e in eps:
        m = epk == e
        B = np.vstack([-tm[m]] + [(dtk[m] / kw.get("win", WIN)) ** p for p in range(kw.get("poly", 2) + 1)]).T
        ce, *_ = np.linalg.lstsq(B, yy[m], rcond=None); Ae.append(ce[0])
    Ae = np.array(Ae)
    err_sc = np.std(Ae, ddof=1) / np.sqrt(len(Ae)) if len(Ae) > 2 else np.nan
    return dict(A=float(c[0]), err_white=float(np.sqrt(cov[0, 0])), err_scatter=float(err_sc), n=len(eps),
                epochs=eps.tolist(), per_epoch=Ae.tolist(), rms=float(s))


def control_phases(n=60, avoid=0.09):
    ph = np.linspace(0.0, 1.0, n, endpoint=False) + 0.5 / n
    ok = (np.abs(ph - 0.0) > avoid) & (np.abs(ph - 1.0) > avoid) & (np.abs(ph - 0.5) > avoid)
    return ph[ok]


def depth_with_errors(t, y, exptime, phase=0.0, ncontrol=60, **kw):
    """Amplitude at the given phase + adopted error = max(white, epoch scatter, control-phase scatter)."""
    r = depth_fit(t, y, exptime, P, T0 + phase * P, **kw)
    if r["n"] == 0: return r
    ctl = []
    for ph in control_phases(ncontrol):
        if abs(((ph - phase + 0.5) % 1) - 0.5) < 0.09: continue
        c = depth_fit(t, y, exptime, P, T0 + ph * P, Ptmpl=P, **kw)
        if c["n"] >= max(1, r["n"] // 2): ctl.append(c["A"])
    ctl = np.array(ctl)
    r["err_control"] = float(np.std(ctl)) if len(ctl) > 5 else np.nan
    r["mean_control"] = float(np.mean(ctl)) if len(ctl) > 5 else np.nan
    r["err"] = float(np.nanmax([r["err_white"], r["err_scatter"], r["err_control"]]))
    return r


def wmean(v, e):
    v, e = np.asarray(v, float), np.asarray(e, float)
    w = 1 / e ** 2
    m = np.sum(w * v) / np.sum(w)
    return float(m), float(1 / np.sqrt(np.sum(w))), float(np.sum(w * (v - m) ** 2))


def load_lc(name):
    return np.load(os.path.join(ROOT, f"data/lc_{name}.npy"), allow_pickle=True).item()
