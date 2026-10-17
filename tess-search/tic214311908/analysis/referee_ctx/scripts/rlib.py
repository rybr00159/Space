"""Referee library (independent of ctx/scripts): ephemeris, transit shape, local-baseline depth fits.

Shape: trapezoid T14 = 1.63 h, T23/T14 = 0.7 (same as the candidate's published shape), supersampled over the
exposure. Depths are fitted per transit window with an independent quadratic baseline in time (local, not a global
biweight), then combined by a joint linear least-squares fit (common depth, separate baselines).
"""
import json
import numpy as np

EPH = json.load(open("/home/user/Space/tess-search/vetting_2026-10-16/results/214311906/ephemeris_v2.json"))
P = float(EPH["P"]); T0 = float(EPH["T0"])
QLP_SHIFT = {int(k): float(v) for k, v in EPH["qlp_shift_days"].items()}
T14 = 1.63 / 24.0
FLAT = 0.7
EXPT = {12: 1800.0, 39: 600.0}


def exptime(sector):
    return EXPT.get(sector, 200.0) / 86400.0


def shape(t, tc, ex, t14=T14, flat=FLAT, nsub=11):
    off = ((np.arange(nsub) + 0.5) / nsub - 0.5) * ex
    tt = np.abs(t[:, None] + off[None, :] - tc)
    t23 = flat * t14
    y = np.clip((t14 / 2 - tt) / ((t14 - t23) / 2), 0, 1)
    return y.mean(axis=1)


def centres(tmin, tmax, phase=0.0):
    k0 = int(np.floor((tmin - T0) / P - phase)) - 1
    k1 = int(np.ceil((tmax - T0) / P - phase)) + 1
    return np.array([T0 + (k + phase) * P for k in range(k0, k1 + 1)]), np.arange(k0, k1 + 1)


def joint_fit(t, f, tcs, ex, half=0.30, deg=2, min_in=2, min_out=10):
    """Common depth (in units of f, positive = dimming) with a separate polynomial baseline per window.
    Returns depth, formal error (from residual scatter), n windows used, per-window list (tc, depth_i, err_i)."""
    blocks = []
    for tc in tcs:
        m = np.abs(t - tc) < half
        if m.sum() < min_out + min_in:
            continue
        s = shape(t[m], tc, ex)
        if np.sum(s > 0.5) < min_in or np.sum(s == 0) < min_out:
            continue
        blocks.append((m, s, tc))
    if not blocks:
        return np.nan, np.nan, 0, []
    nb = len(blocks); npar = 1 + nb * (deg + 1)
    rows = []; y = []
    for b, (m, s, tc) in enumerate(blocks):
        x = (t[m] - tc) / half
        A = np.zeros((m.sum(), npar)); A[:, 0] = -s
        for d in range(deg + 1):
            A[:, 1 + b * (deg + 1) + d] = x ** d
        rows.append(A); y.append(f[m])
    A = np.vstack(rows); y = np.concatenate(y)
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A @ coef; dof = max(len(y) - npar, 1)
    s2 = np.sum(r ** 2) / dof
    cov = s2 * np.linalg.pinv(A.T @ A)
    per = []
    for m, s, tc in blocks:
        x = (t[m] - tc) / half
        B = np.column_stack([-s] + [x ** d for d in range(deg + 1)])
        c, *_ = np.linalg.lstsq(B, f[m], rcond=None)
        rr = f[m] - B @ c; cv = np.sum(rr ** 2) / max(m.sum() - B.shape[1], 1) * np.linalg.pinv(B.T @ B)
        per.append((tc, c[0], np.sqrt(cv[0, 0])))
    return coef[0], np.sqrt(cov[0, 0]), nb, per


def fake_window_depths(t, f, ex, avoid_tcs, half=0.30, step=0.05, deg=2, guard=None):
    """Single-window depths at fake centres (every `step` d) whose window does not touch any real transit or
    phase-0.5 window. Returns centres and depths (one per fake window; neighbouring windows overlap, so effective
    sample size is ~ span/(2*half))."""
    guard = 2 * half if guard is None else guard
    cs = np.arange(t.min() + half, t.max() - half, step)
    out_c, out_d = [], []
    for c in cs:
        if np.min(np.abs(avoid_tcs - c)) < guard:
            continue
        d, e, n, _ = joint_fit(t, f, [c], ex, half, deg)
        if n == 1:
            out_c.append(c); out_d.append(d)
    return np.array(out_c), np.array(out_d)
