"""Per-window image products for the two independent image methods.

DIFF ('difference image'): for a window centred at tc, in-transit cadences are those with exposure-integrated
  unit template >= 0.75 (flat bottom); out-of-transit (OOT) cadences lie 0.5*T14 + 0.5*exptime + 12 min to that
  + W_OOT on BOTH sides.  Per pixel, a straight line a + b (t - tc) is fitted to the OOT cadences (removes the
  Mira's slow variability and the background trend locally); D = (line(t_in) - flux_in).mean() / template_in.mean()
  (positive = the pixel dims; estimates the full transit depth in e-/s) with variance from the OOT scatter.
AMP ('amplitude image'): per pixel, global regression over each orbit segment on a cubic B-spline (knot spacing
  KNOT d) baseline plus the unit transit template.  By Frisch-Waugh-Lovell, alpha = tau^T M y / tau^T M tau with M the
  spline residual-maker, so per-window numerators/denominators can be cached and any set of windows combined.
Both are linear in the pixel data, so an injected signal s(t) * img adds exactly (response) * img.
"""
import numpy as np
from scipy.interpolate import BSpline
from l0lib import template, T14, P, T0, real_tcs

W_OOT = 0.25
KNOT = 0.5


def segments(t, gap=0.5):
    return np.concatenate([[0], np.cumsum(np.diff(t) > gap)])


def bspline_basis(t, knot, deg=3):
    lo, hi = t.min() - 1e-6, t.max() + 1e-6
    n = max(1, int(np.ceil((hi - lo) / knot)))
    inner = np.linspace(lo, hi, n + 1)
    k = np.concatenate([[lo] * deg, inner, [hi] * deg])
    nb = len(k) - deg - 1
    B = np.zeros((len(t), nb))
    for i in range(nb):
        c = np.zeros(nb); c[i] = 1
        B[:, i] = BSpline(k, c, deg, extrapolate=False)(t)
    return np.nan_to_num(B)


def background(cube, frac=0.3, clip=3.0):
    """Per-cadence background plane from the faintest pixels (by median image), minus its median level:
    returns (nt, ny, nx) background-variation cube and the pixel mask used."""
    nt, ny, nx = cube.shape
    med = np.median(cube, 0)
    sel = med <= np.percentile(med, 100 * frac)
    yy, xx = np.mgrid[0:ny, 0:nx]
    X = np.column_stack([np.ones(sel.sum()), xx[sel] - nx / 2, yy[sel] - ny / 2])
    Z = (cube - med[None])[:, sel]                          # (nt, nsel)
    w = np.ones_like(Z)
    for _ in range(4):
        coef = np.array([np.linalg.lstsq(X * w[i][:, None], Z[i] * w[i], rcond=None)[0] for i in range(nt)])
        r = Z - coef @ X.T
        s = 1.4826 * np.median(np.abs(r), 1, keepdims=True) + 1e-9
        w = (np.abs(r) < clip * s).astype(float)
    full = np.column_stack([np.ones(ny * nx), xx.ravel() - nx / 2, yy.ravel() - ny / 2])
    return (coef @ full.T).reshape(nt, ny, nx), sel


class Data:
    """Good cadences of one sector as (nt, npix); per-sector caches for both methods.  By default the per-cadence
    background variation (plane fitted to the faintest 30% pixels) is subtracted first."""
    def __init__(self, S, knot=KNOT, flux=None, bkg=True):
        self.S = S
        gd = S.good
        self.t = S.t[gd]
        cube = (S.flux[gd] if flux is None else flux[gd]).astype(float)
        if bkg:
            self.bkgcube, self.bkgsel = background(cube)
            cube = cube - self.bkgcube
        self.Y = cube.reshape(gd.sum(), -1)
        self.exptime = S.exptime
        self.seg = segments(self.t)
        self.knot = knot
        # spline residual-maker per segment, applied once to the cube
        self.Yres = np.zeros_like(self.Y)
        self.Q = []                       # orthonormal bases per segment (for M tau)
        for s in np.unique(self.seg):
            m = self.seg == s
            if m.sum() < 30:
                self.Yres[m] = 0.0; self.Q.append((m, None)); continue
            B = bspline_basis(self.t[m], knot)
            q, _ = np.linalg.qr(B)
            self.Q.append((m, q))
            self.Yres[m] = self.Y[m] - q @ (q.T @ self.Y[m])
        self.sig_amp = 1.4826 * np.median(np.abs(self.Yres - np.median(self.Yres, 0)), 0)

    def Mtau(self, tau):
        out = np.zeros_like(tau)
        for m, q in self.Q:
            if q is None: continue
            out[m] = tau[m] - q @ (q.T @ tau[m])
        return out

    # ---------------- AMP
    def amp_window(self, tc):
        tau = template(self.t, [tc], self.exptime)
        if tau.sum() < 0.5 * T14 / self.exptime * 0.85:          # need >= ~half of the transit covered
            return None
        mt = self.Mtau(tau)
        return dict(num=mt @ self.Y if False else tau @ self.Yres, den=float(tau @ mt), ncov=float(tau.sum()))

    # ---------------- DIFF
    def diff_window(self, tc, w_oot=W_OOT, min_in_frac=0.5, min_side_frac=0.3):
        t = self.t; e = self.exptime
        tau = template(t, [tc], e)
        a = np.abs(t - tc)
        inn = tau >= 0.75
        exp_in = max(1, int((T14 * 0.7 - e) / e))
        if inn.sum() < max(1, min_in_frac * exp_in):
            return None
        lo = 0.5 * T14 + 0.5 * e + 12 / 1440
        oot = (a >= lo) & (a <= lo + w_oot)
        left = oot & (t < tc); right = oot & (t > tc)
        nside = w_oot / e
        if left.sum() < min_side_frac * nside or right.sum() < min_side_frac * nside:
            return None
        dt = t[oot] - tc
        X = np.column_stack([np.ones(oot.sum()), dt])
        XtXi = np.linalg.inv(X.T @ X)
        coef = XtXi @ X.T @ self.Y[oot]
        res = self.Y[oot] - X @ coef
        s2 = (res ** 2).sum(0) / (oot.sum() - 2)
        xin = np.column_stack([np.ones(inn.sum()), t[inn] - tc])
        pred = (xin @ coef).mean(0)
        xbar = xin.mean(0)
        tin = tau[inn].mean()
        D = (pred - self.Y[inn].mean(0)) / tin
        V = s2 * (1.0 / inn.sum() + xbar @ XtXi @ xbar) / tin ** 2
        return dict(D=D, V=V, nin=int(inn.sum()), nout=int(oot.sum()))


def window_sets(D, margin=None, step=0.05):
    """Real transit centres with coverage, and fake centres (>= margin from every real transit; also avoiding the
    real transits' OOT windows)."""
    t = D.t
    tr = real_tcs(t)
    tr = tr[(tr > t.min() - 0.1) & (tr < t.max() + 0.1)]
    if margin is None:
        margin = T14 + 2 * (0.5 * T14 + W_OOT + 0.02) + 0.05     # fake window + its OOT never touch real ones
    fk = np.arange(t.min() + 0.3, t.max() - 0.3, step)
    ok = np.all(np.abs(fk[:, None] - tr[None, :]) > margin, axis=1)
    return tr, fk[ok]
