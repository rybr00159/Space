"""Transit-model helpers: exposure-integrated trapezoid, linear depth fits, per-transit fits with prayer-bead errors."""
import numpy as np
from scipy.optimize import least_squares


def trap_frac(t, tc, T14, fl, exptime, nsub=None):
    """Fraction of full depth (0..1) of a symmetric trapezoid (total duration T14, flat part fl*T14), averaged over
    the exposure. tc can be an array broadcastable with t (per-cadence nearest mid-time)."""
    if nsub is None: nsub = 1 if exptime < 0.003 else (5 if exptime < 0.008 else 11)
    y = np.zeros(np.shape(t))
    T23 = fl * T14
    for o in (np.arange(nsub) + 0.5) / nsub - 0.5:
        a = np.abs(t + o * exptime - tc)
        y = y + np.where(a <= T23 / 2, 1.0, np.where(a >= T14 / 2, 0.0, (T14 / 2 - a) / max((T14 - T23) / 2, 1e-9)))
    return y / nsub


def ephem_frac(t, P, T0, T14, fl, exptime):
    tc = T0 + np.round((t - T0) / P) * P
    return trap_frac(t, tc, T14, fl, exptime)


def lin_depth(f, e, y, extra=None):
    """Weighted LSQ of f = c0 - depth*y (+ extra columns); returns depth, err, chi2."""
    cols = [np.ones(len(f)), -y] + ([] if extra is None else list(extra))
    A = np.column_stack(cols) / e[:, None]
    c, *_ = np.linalg.lstsq(A, f / e, rcond=None)
    cov = np.linalg.pinv(A.T @ A)
    r = f / e - A @ c
    return c[1], np.sqrt(cov[1, 1]), float(np.sum(r ** 2)), c


def fit_single(t, f, e, exptime, tc0, T14, fl, win=0.3, dt_max=0.1, slope=True):
    """Fit one transit: free tc (within +-dt_max of tc0), depth, baseline (+slope). Grid then least squares.
    Returns dict with tc, depth, formal errors, chi2, n, residuals, and a prayer-bead error for tc and depth."""
    m = np.abs(t - tc0) < win
    tt, ff, ee = t[m], f[m], e[m]
    if m.sum() < 8: return None
    def model(p, tt_=tt):
        y = trap_frac(tt_, tc0 + p[0], T14, fl, exptime)
        base = p[2] + (p[3] * (tt_ - tc0) if slope else 0.0)
        return base * (1 - p[1] * y)
    grid = np.round(np.arange(-dt_max, dt_max + 1e-9, 0.002), 6)   # rounded: tiny float noise breaks scipy's finite differences
    best = None
    for dt in grid:
        y = trap_frac(tt, tc0 + dt, T14, fl, exptime)
        if np.sum(y > 0.5) < 2: continue
        cols = [np.ones(len(tt)), -y] + ([tt - tc0] if slope else [])
        A = np.column_stack(cols) / ee[:, None]; c, *_ = np.linalg.lstsq(A, ff / ee, rcond=None)
        c2 = np.sum((ff / ee - A @ c) ** 2)
        if c[1] > 0 and (best is None or c2 < best[0]): best = (c2, dt, c)
    if best is None: return None
    c2, dt, c = best
    # parameters: [dtc (d, relative to tc0), depth, baseline, slope]; all |x| < 1 so scipy's diff_step is absolute
    p0 = [dt if abs(dt) > 1e-9 else 1e-4, max(c[1] / max(c[0], 1e-6), 1e-4), c[0]] + ([c[2]] if slope else [0.0])
    def res(p): return (ff - model(p)) / ee
    lo = [-dt_max - 0.01, -0.5, 0.5, -1]; hi = [dt_max + 0.01, 0.5, 1.5, 1]
    if not slope: lo, hi = lo[:3] + [-1e-12], hi[:3] + [1e-12]
    DS = [2e-4, 1e-4, 1e-5, 1e-5]
    q = least_squares(res, p0, bounds=(lo, hi), x_scale=[0.002, 0.003, 0.003, 0.01], diff_step=DS)
    J = q.jac; dof = max(1, len(tt) - 4)
    cov = np.linalg.pinv(J.T @ J) * max(1.0, 2 * q.cost / dof)
    r = q.fun * ee
    # prayer-bead: cyclic shifts of the residuals added to the best model, refit
    mod = model(q.x)
    pb = []
    n = len(r)
    for k in np.unique(np.linspace(1, n - 1, min(n - 1, 60)).astype(int)):
        ff2 = mod + np.roll(r, k)
        q2 = least_squares(lambda p: (ff2 - model(p)) / ee, q.x, bounds=(lo, hi), x_scale=[0.002, 0.003, 0.003, 0.01], diff_step=DS)
        pb.append(q2.x[:2])
    pb = np.array(pb)
    # chi2 profile in tc for a robust (non-Gaussian) interval
    prof = []
    for dtc in np.arange(-dt_max, dt_max + 1e-9, 0.001):
        y = trap_frac(tt, tc0 + dtc, T14, fl, exptime)
        cols = [np.ones(len(tt)), -y] + ([tt - tc0] if slope else [])
        A = np.column_stack(cols) / ee[:, None]; c_, *_ = np.linalg.lstsq(A, ff / ee, rcond=None)
        prof.append((tc0 + dtc, np.sum((ff / ee - A @ c_) ** 2)))
    prof = np.array(prof)
    ok1 = prof[:, 1] <= prof[:, 1].min() + max(1.0, 2 * q.cost / dof)
    tci = (prof[ok1, 0].min(), prof[ok1, 0].max())
    return dict(tc=float(tc0 + q.x[0]), depth=float(q.x[1]), tc_err=float(np.sqrt(cov[0, 0])), depth_err=float(np.sqrt(cov[1, 1])),
                tc_err_pb=float(np.std(pb[:, 0])) if len(pb) > 3 else np.nan, depth_err_pb=float(np.std(pb[:, 1])) if len(pb) > 3 else np.nan,
                tc_prof_lo=float(tci[0]), tc_prof_hi=float(tci[1]),
                chi2=float(2 * q.cost), n=int(len(tt)), n_in=int(np.sum(trap_frac(tt, tc0 + q.x[0], T14, fl, exptime) > 0.5)),
                redchi=float(2 * q.cost / dof))
