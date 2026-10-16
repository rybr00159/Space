"""Difference-image + PRF-fit localisation tools (referee, independent code)."""
import sys
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/five_referee/code')
import numpy as np
from scipy.optimize import least_squares
from reflib import *


def event_diff_images(cut, tcs, T14, fin=0.3, gap=0.6, wout=0.15, inject=None):
    """Per-event difference images (baseline - in-eclipse; positive = deficit).

    inject: optional (template_image, depth_fraction, shape_func) to inject a fake
    signal on the given pixel template at times tcs (for bias tests).
    Returns list of (tc, diff_img, var_img_estimate, n_in), and direct image.
    """
    t = cut['time']
    F = cut['flux'].copy()
    if inject is not None:
        tmpl, frac = inject
        # trapezoid-like V: depth weighting by fraction of time inside eclipse
        for tc in tcs:
            x = np.abs(t - tc)
            w = np.clip(1 - x / (T14 / 2.), 0, 1)  # V shape, peak 1 at centre
            F -= (w * frac)[:, None, None] * tmpl[None, :, :]
    out = []
    oot_all = np.ones(len(t), bool)
    for tc in tcs:
        dt = t - tc
        oot_all &= np.abs(dt) > gap * T14
        inn = np.abs(dt) < fin * T14
        oo = (np.abs(dt) > gap * T14) & (np.abs(dt) < gap * T14 + wout)
        if inn.sum() < 2 or (oo & (dt < 0)).sum() < 3 or (oo & (dt > 0)).sum() < 3:
            continue
        # per-pixel linear baseline fit vs time on oot samples
        A = np.vstack([np.ones(oo.sum()), dt[oo]]).T
        Y = F[oo].reshape(oo.sum(), -1)
        coef, *_ = np.linalg.lstsq(A, Y, rcond=None)
        res = Y - A @ coef
        sig2 = np.var(res, axis=0, ddof=2)
        Ai = np.vstack([np.ones(inn.sum()), dt[inn]]).T
        base_in = Ai @ coef
        din = base_in - F[inn].reshape(inn.sum(), -1)
        d = din.mean(axis=0).reshape(F.shape[1:])
        # variance of mean of in-eclipse residuals + baseline uncertainty (approx)
        v = (sig2 / inn.sum() + sig2 / oo.sum()).reshape(F.shape[1:])
        out.append((tc, d, v, inn.sum()))
    direct = np.nanmedian(F[oot_all], axis=0)
    return out, direct


def combine(evs, idx=None):
    if idx is None:
        idx = range(len(evs))
    D = np.array([evs[i][1] for i in idx])
    V = np.array([evs[i][2] for i in idx])
    d = D.mean(axis=0)
    # empirical error from event-to-event scatter, floor by formal
    if len(D) > 2:
        e = np.sqrt(np.maximum(D.var(axis=0, ddof=1) / len(D), V.mean(axis=0) / len(D)))
    else:
        e = np.sqrt(V.mean(axis=0) / len(D))
    return d, e


def fit_prf(img, err, prf, x0, y0, box=3, fix_xy=None, extra_sources=None):
    """Fit A*PRF(x,y)+c over a (2*box+1)^2 box centred near (x0,y0).
    extra_sources: list of (dx, dy, relflux) fixed relative to the main source (for direct image).
    Returns x, y, A, c, cov, chi2, npix."""
    ny, nx = img.shape
    x0, y0 = float(x0), float(y0)
    xc, yc = int(round(x0)), int(round(y0))
    ys, xs = slice(max(0, yc - box), min(ny, yc + box + 1)), slice(max(0, xc - box), min(nx, xc + box + 1))
    sub = img[ys, xs]
    se = err[ys, xs]

    def model(p):
        x, y, A, c = p
        m = A * prf(x, y)
        if extra_sources:
            for dx, dy, rf in extra_sources:
                m = m + A * rf * prf(x + dx, y + dy)
        return m[ys, xs] + c

    def resid(p):
        return ((sub - model(p)) / se).ravel()
    A0 = np.nansum(sub)
    r = least_squares(resid, [x0, y0, A0, 0.0], x_scale=[0.05, 0.05, abs(A0) * 0.05 + 1e-6, abs(A0) * 0.01 + 1e-6],
                      diff_step=[1e-3, 1e-3, 1e-4, 1e-4])
    J = r.jac
    chi2 = np.sum(r.fun ** 2)
    dof = r.fun.size - 4
    try:
        cov = np.linalg.inv(J.T @ J) * max(1.0, chi2 / dof)
    except np.linalg.LinAlgError:
        cov = np.full((4, 4), np.nan)
    return dict(x=r.x[0], y=r.x[1], A=r.x[2], c=r.x[3], cov=cov, chi2=chi2, dof=dof)


def pix_to_EN(cut, x, y, xr, yr):
    """Convert pixel offset (x-xr, y-yr) to (dE, dN) arcsec using the cutout WCS."""
    w = cut['wcs']
    ra1, de1 = w.pixel_to_world_values(xr, yr)
    ra2, de2 = w.pixel_to_world_values(x, y)
    dE = (ra2 - ra1) * np.cos(np.radians(de1)) * 3600.
    dN = (de2 - de1) * 3600.
    return dE, dN


def EN_to_pix(cut, ra0, dec0, dE, dN):
    w = cut['wcs']
    ra = ra0 + dE / 3600. / np.cos(np.radians(dec0))
    de = dec0 + dN / 3600.
    return w.world_to_pixel_values(ra, de)
