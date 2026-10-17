"""tfit: transit-model and fitting helpers shared by the cleaned-light-curve analyses (c03-c08).

* batman quadratic-LD model, circular orbit, exposure-integrated with the midpoint rule (n_sub = 5 / 9 / 21 for the
  200 / 600 / 1800-s cadences; batman's own supersample_factor is not used, see referee2_shape).
* Data are handled as 'windows': one window per transit epoch (or per fake epoch), each with its own baseline
  (offset + slope by default, multiplicative: f = (c0 + c1 tau) * model), profiled analytically.
* Per-transit fits: fixed shape, free mid-time and depth.
"""
import numpy as np
import batman

G_CGS = 6.674e-8
LD0 = (0.125, 0.533)          # Claret (2017) PHOENIX TESS, 3350 K, logg 4.84 (shape/out/ld_claret17.json)
RHO_STAR = (9.3, 0.9)         # single-star density prior, g/cm3 (star/, referee_star/; shape used 9.41 +- 0.88)
R_STAR = (0.383, 0.012)       # Rsun
R_EARTH_RSUN = 0.0091577


def nsub_for(exptime_d):
    s = exptime_d * 86400
    return 5 if s < 300 else (9 if s < 900 else 21)


def a_over_r(rho, P):
    return (G_CGS * rho * (P * 86400.0) ** 2 / (3 * np.pi)) ** (1 / 3)


class Model:
    """Transit model evaluated at times t (BTJD_TDB) with exposure time(s) exptime (array, d)."""
    def __init__(self, t, exptime):
        t = np.asarray(t, float); ex = np.broadcast_to(np.asarray(exptime, float), t.shape)
        self.n = len(t)
        self.groups = []
        for e in np.unique(ex):
            m = ex == e; ns = nsub_for(e)
            off = ((np.arange(ns) + 0.5) / ns - 0.5) * e
            tt = (t[m][:, None] + off[None]).ravel()
            p = batman.TransitParams()
            p.t0, p.per, p.rp, p.a, p.inc, p.ecc, p.w, p.u, p.limb_dark = 0.0, 5.69, 0.12, 25.0, 89.0, 0.0, 90.0, list(LD0), "quadratic"
            self.groups.append((m, ns, tt, p, batman.TransitModel(p, tt)))

    def flux(self, T0, P, k, b, aR, u1, u2):
        out = np.ones(self.n)
        inc = np.degrees(np.arccos(np.clip(b / aR, 0, 1)))
        for m, ns, tt, p, bm in self.groups:
            p.t0, p.per, p.rp, p.a, p.inc, p.u = T0, P, k, aR, inc, [u1, u2]
            out[m] = bm.light_curve(p).reshape(-1, ns).mean(1)
        return out


def mid_depth(k, b, aR, u1, u2, P=5.6915):
    """Mid-transit depth of the instantaneous model."""
    p = batman.TransitParams()
    p.t0, p.per, p.rp, p.a, p.inc, p.ecc, p.w, p.u, p.limb_dark = 0.0, P, k, aR, np.degrees(np.arccos(min(b / aR, 1))), 0.0, 90.0, [u1, u2], "quadratic"
    return 1 - batman.TransitModel(p, np.array([0.0])).light_curve(p)[0]


def t14(P, k, b, aR):
    inc = np.arccos(b / aR)
    arg = np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) / (aR * np.sin(inc))
    return P / np.pi * np.arcsin(np.clip(arg, 0, 1))


def t23(P, k, b, aR):
    inc = np.arccos(b / aR)
    arg = np.sqrt(np.clip((1 - k) ** 2 - b ** 2, 0, None)) / (aR * np.sin(inc))
    return P / np.pi * np.arcsin(np.clip(arg, 0, 1))


# ------------------------------------------------------------------------------------------------ windows
class Windows:
    """Concatenated windows: t, f (normalised flux, ~1), s (per-point sigma), w (window index), tc0 (nominal centre per
    window), sector per window, epoch per window, exptime per point."""
    def __init__(self, t, f, s, w, tc0, sector, epoch, exptime, order=1):
        self.t, self.f, self.s, self.w = map(np.asarray, (t, f, s, w))
        self.tc0 = np.asarray(tc0, float); self.sector = np.asarray(sector); self.epoch = np.asarray(epoch)
        self.ex = np.broadcast_to(np.asarray(exptime, float), self.t.shape).copy()
        self.order = order
        self.nw = len(self.tc0)
        self.idx = [np.where(self.w == i)[0] for i in range(self.nw)]
        self.tau = self.t - self.tc0[self.w]
        self._M = None

    def subset(self, keep_windows):
        keep_windows = np.asarray(keep_windows)
        sel = np.isin(self.w, keep_windows)
        remap = -np.ones(self.nw, int); remap[keep_windows] = np.arange(len(keep_windows))
        return Windows(self.t[sel], self.f[sel], self.s[sel], remap[self.w[sel]], self.tc0[keep_windows], self.sector[keep_windows],
                       self.epoch[keep_windows], self.ex[sel], self.order)

    def model(self):
        if self._M is None: self._M = Model(self.t, self.ex)
        return self._M

    def profile(self, mod, scale=None):
        """chi2 with per-window multiplicative polynomial baseline profiled; returns chi2, baseline-corrected residuals."""
        s = self.s if scale is None else self.s * scale
        chi2 = 0.0; res = np.zeros(len(self.t)); base = np.zeros(len(self.t))
        for i, ii in enumerate(self.idx):
            X = np.column_stack([mod[ii] * self.tau[ii] ** j for j in range(self.order + 1)]) / s[ii, None]
            y = self.f[ii] / s[ii]
            c, *_ = np.linalg.lstsq(X, y, rcond=None)
            r = y - X @ c; chi2 += r @ r; res[ii] = r * s[ii]
            base[ii] = (X @ c) * s[ii] / mod[ii]
        return chi2, res, base


def fit_depth_fixed_shape(W, shape_frac, scale=None, groups=None):
    """Linear fit of a common depth (per group) with fixed shape 'shape_frac' (0 out, 1 at mid-transit-depth units),
    per-window polynomial baselines (additive). groups: array per window (group id), default one group.
    Returns depths, formal errors, chi2."""
    s = W.s if scale is None else W.s * scale
    g = np.zeros(W.nw, int) if groups is None else np.asarray(groups)
    ug = np.unique(g)
    ncol = len(ug) + W.nw * (W.order + 1)
    X = np.zeros((len(W.t), ncol))
    for j, gg in enumerate(ug):
        m = g[W.w] == gg; X[m, j] = -shape_frac[m]
    for i, ii in enumerate(W.idx):
        for o in range(W.order + 1):
            X[ii, len(ug) + i * (W.order + 1) + o] = W.tau[ii] ** o
    Xw = X / s[:, None]; yw = W.f / s
    c, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    cov = np.linalg.pinv(Xw.T @ Xw)
    r = yw - Xw @ c
    return c[:len(ug)], np.sqrt(np.diag(cov)[:len(ug)]), float(r @ r), ug


def single_fit(t, f, s, ex, tc0, shape_fn, dt_max=0.06, order=1, step=0.0005):
    """Fit one window: free tc (grid + parabolic refinement), depth (mid-transit units), polynomial baseline.
    shape_fn(t, tc, ex) -> fraction of full depth. Returns dict(tc, depth, depth_err, tc_err (curvature), chi2, n_in)."""
    tau = t - tc0
    base_cols = [tau ** o for o in range(order + 1)]
    grid = np.arange(-dt_max, dt_max + 1e-12, step)
    best = (np.inf, None, None)
    chis = np.full(len(grid), np.nan)
    for j, d in enumerate(grid):
        y = shape_fn(t, tc0 + d, ex)
        if np.sum(y > 0.5) < 2: continue
        X = np.column_stack([-y] + base_cols) / s[:, None]
        c, *_ = np.linalg.lstsq(X, f / s, rcond=None)
        r = f / s - X @ c; c2 = r @ r; chis[j] = c2
        if c2 < best[0]: best = (c2, d, c)
    if best[1] is None: return None
    c2, d, c = best
    j = int(np.argmin(np.where(np.isfinite(chis), chis, np.inf)))
    # parabolic refinement and curvature error
    tc = tc0 + d; tc_err = np.nan
    if 0 < j < len(grid) - 1 and np.all(np.isfinite(chis[j - 1:j + 2])):
        y0, y1, y2 = chis[j - 1:j + 2]; den = (y0 - 2 * y1 + y2)
        if den > 0:
            dd = 0.5 * (y0 - y2) / den * step; tc = tc0 + grid[j] + dd
            tc_err = step * np.sqrt(2 / den)
    # chi2 interval (delta chi2 = 1) from the grid, robust to non-parabolic profiles
    okc = np.isfinite(chis) & (chis <= np.nanmin(chis) + 1.0)
    lo, hi = grid[okc].min(), grid[okc].max()
    y = shape_fn(t, tc, ex)
    X = np.column_stack([-y] + base_cols) / s[:, None]
    cc, *_ = np.linalg.lstsq(X, f / s, rcond=None)
    cov = np.linalg.pinv(X.T @ X)
    r = f / s - X @ cc
    return dict(tc=float(tc), depth=float(cc[0]), depth_err=float(np.sqrt(cov[0, 0])), tc_err=float(tc_err),
                tc_lo=float(tc0 + lo), tc_hi=float(tc0 + hi), chi2=float(r @ r), n=int(len(t)), n_in=int(np.sum(y > 0.5)))


def wls_line(x, y, e):
    """Weighted linear fit y = a + b x; returns (a, b), cov, chi2."""
    A = np.column_stack([np.ones(len(x)), x]) / e[:, None]
    c, *_ = np.linalg.lstsq(A, y / e, rcond=None)
    cov = np.linalg.inv(A.T @ A)
    r = y / e - A @ c
    return c, cov, float(r @ r)
