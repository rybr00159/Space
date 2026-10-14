"""Fast grid-based PRF fitting for difference images and direct images.

The TESS PRF (PRF package) is pre-evaluated on a 0.01-px grid of positions around the WCS target position.
For every grid position the linear parameters (amplitude, constant) are solved in closed form; the chi2
surface minimum is refined with a 2-D quadratic fit to the 3x3 neighbourhood.
Direct image model: F * [PRF(x,y) + sum_i r_i PRF(x+dx_i, y+dy_i)] + c   (neighbour fluxes tied to Gaia ratios)
Difference image model: A * PRF(x,y) + c
"""
import numpy as np
from ana import Ctx, circle_mask, RFIT

STEP = 0.01
HALF = 0.70


class Fast:
    def __init__(self, sec, psfkind='prf'):
        self.ctx = ctx = Ctx(sec, psfkind=psfkind)
        self.sec = sec
        self.mask = ctx.mask
        self.g = np.arange(-HALF, HALF + 1e-9, STEP)
        n = len(self.g)
        self.n = n
        m = self.mask
        self.M = np.empty((n, n, m.sum()))      # [ix, iy, pix]
        self.O = np.empty((n, n, m.sum()))
        for i, dx in enumerate(self.g):
            for j, dy in enumerate(self.g):
                x = ctx.xt0 + dx; y = ctx.yt0 + dy
                self.M[i, j] = ctx.psf(x, y)[m]
                self.O[i, j] = ctx.others_img(x, y)[m]
        self.Mf = self.M.reshape(n * n, -1)
        self.MOf = (self.M + self.O).reshape(n * n, -1)
        self.jac = ctx.jac

    def _solve(self, T, img, sig):
        """T: (P, npix) templates. Returns chi2 (P,), amp (P,), const (P,)."""
        m = self.mask
        w = 1.0 / sig[m]
        y = img[m] * w
        Tw = T * w[None, :]
        u = w
        TT = (Tw * Tw).sum(1); Tu = Tw @ u; Ty = Tw @ y
        uu = u @ u; uy = u @ y; yy = y @ y
        det = TT * uu - Tu ** 2
        a = (Ty * uu - Tu * uy) / det
        c = (TT * uy - Tu * Ty) / det
        chi2 = yy - a * Ty - c * uy
        return chi2, a, c

    def _refine(self, chi2):
        n = self.n
        C = chi2.reshape(n, n)
        k = int(np.argmin(chi2)); i, j = divmod(k, n)
        i = min(max(i, 1), n - 2); j = min(max(j, 1), n - 2)
        sub = C[i - 1:i + 2, j - 1:j + 2]
        # quadratic fit z = a + b x + c y + d x^2 + e xy + f y^2 on 3x3
        xx, yy = np.meshgrid([-1, 0, 1], [-1, 0, 1], indexing='ij')
        A = np.column_stack([np.ones(9), xx.ravel(), yy.ravel(), xx.ravel() ** 2, (xx * yy).ravel(), yy.ravel() ** 2])
        p, *_ = np.linalg.lstsq(A, sub.ravel(), rcond=None)
        H = np.array([[2 * p[3], p[4]], [p[4], 2 * p[5]]])
        try:
            d = -np.linalg.solve(H, p[1:3])
        except np.linalg.LinAlgError:
            d = np.zeros(2)
        if not np.all(np.abs(d) <= 1.0) or np.any(np.linalg.eigvalsh(H) <= 0):
            d = np.zeros(2)
        cmin = p[0] + p[1] * d[0] + p[2] * d[1] + p[3] * d[0] ** 2 + p[4] * d[0] * d[1] + p[5] * d[1] ** 2
        at_edge = (k // n in (0, n - 1)) or (k % n in (0, n - 1))
        # formal covariance (pixels) from curvature: chi2 = cmin + 0.5 d^T H' d, H' in px^-2
        Hp = H / STEP ** 2
        try:
            cov = 2.0 * np.linalg.inv(Hp)
        except np.linalg.LinAlgError:
            cov = np.full((2, 2), np.nan)
        return self.g[i] + d[0] * STEP, self.g[j] + d[1] * STEP, cmin, cov, at_edge

    def fit_diff(self, D, sig):
        chi2, a, c = self._solve(self.Mf, D, sig)
        dx, dy, cmin, cov, edge = self._refine(chi2)
        k = int(np.argmin(chi2))
        return dict(dx=dx, dy=dy, amp=float(a[k]), chi2=float(cmin), cov=cov, edge=edge, chi2_grid=chi2)

    def fit_direct(self, B, sig):
        chi2, a, c = self._solve(self.MOf, B, sig)
        dx, dy, cmin, cov, edge = self._refine(chi2)
        k = int(np.argmin(chi2))
        return dict(dx=dx, dy=dy, F=float(a[k]), c=float(c[k]), chi2=float(cmin), edge=edge)

    def chi2_at(self, D, sig, dx, dy):
        """chi2 of the difference image for a PRF at a given offset (px) from the WCS target position."""
        ctx = self.ctx; m = self.mask
        M = ctx.psf(ctx.xt0 + dx, ctx.yt0 + dy)[m][None, :]
        chi2, a, c = self._solve(M, D, sig)
        return float(chi2[0]), float(a[0])

    def offset_sky(self, rD, rB):
        v = np.array([rD['dx'] - rB['dx'], rD['dy'] - rB['dy']])
        return self.jac @ v, v

    def empirical_image(self, B, sig):
        """Unit-flux empirical target image from a direct image: (B - c - F*others)/F at the fitted position."""
        rB = self.fit_direct(B, sig)
        ctx = self.ctx
        x = ctx.xt0 + rB['dx']; y = ctx.yt0 + rB['dy']
        E = (B - rB['c'] - rB['F'] * ctx.others_img(x, y)) / rB['F']
        return E, rB

    def prf_img(self, dx, dy):
        ctx = self.ctx
        return ctx.psf(ctx.xt0 + dx, ctx.yt0 + dy)
