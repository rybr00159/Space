"""Two-sphere detached eclipsing-binary model with quadratic limb darkening and third light (v5, written fresh).

occ(z, p, u1, u2): fraction of the flux of a quadratically limb-darkened star of unit radius that is hidden by an
opaque disc of radius p whose centre is at projected distance z (both in units of the occulted star's radius).
Exact radial integral   blocked = int_0^1 I(r) 2 r alpha(r) dr / int_0^1 I(r) 2 r dr,
  alpha(r) = 1                                   for r < p - z,
           = arccos((r^2+z^2-p^2)/(2 r z)) / pi  for |z-p| < r < z+p,
           = 0                                   otherwise,
evaluated on [|z-p|, min(1,z+p)] by Gauss-Legendre after the substitution r = a + (b-a)(1-cos th)/2 (removes the
square-root endpoint singularities) plus the analytic fully covered core.  Valid for any p (p > 1 is needed for the
secondary eclipse, where the larger star occults the smaller one).

Light curve: star 1 is the star eclipsed at T0 (primary eclipse).  Parameters
  rsum = (R1+R2)/a, k = R2/R1, b0 = a cos i / (R1+R2)  [impact parameter in units of the radius sum, circular part],
  J = L2/(L1 k^2) = ratio of disc-averaged band surface brightness star2/star1,
  L3 = fraction of the total (EB + others) band flux that is not from the EB,
  ecosw, esinw, quadratic LD (u1,u2) of each star.
The model is computed on a fine grid of time since conjunction, integrated over the exposure with the cumulative
trapezoid (exact for the piecewise-linear grid model), and interpolated to the data times."""
import numpy as np
from numba import njit

NG = 48
XG, WG = np.polynomial.legendre.leggauss(NG)

@njit(cache=True)
def _F(x, u1, u2):
    # antiderivative in mu of I(mu)*2*mu, I = 1 - u1(1-mu) - u2(1-mu)^2
    return x * x - 2 * u1 * (x * x / 2 - x ** 3 / 3) - 2 * u2 * (x * x / 2 - 2 * x ** 3 / 3 + x ** 4 / 4)

@njit(cache=True)
def _cum(r, u1, u2):
    m = np.sqrt(max(0.0, 1.0 - r * r))
    return _F(1.0, u1, u2) - _F(m, u1, u2)

@njit(cache=True)
def occ(z, p, u1, u2, xg, wg):
    if p <= 0.0 or z >= 1.0 + p:
        return 0.0
    tot = _cum(1.0, u1, u2)
    if p - z >= 1.0:
        return 1.0
    blocked = 0.0
    if p - z > 0.0:
        blocked += _cum(p - z, u1, u2)
    ra = abs(z - p)
    rb = min(1.0, z + p)
    if rb > ra and z > 0.0:
        s = 0.0
        for k in range(xg.size):
            th = 0.5 * np.pi * (xg[k] + 1.0)
            r = ra + (rb - ra) * 0.5 * (1.0 - np.cos(th))
            dr = (rb - ra) * 0.5 * np.sin(th) * 0.5 * np.pi
            if r <= 0.0:
                continue
            c = (r * r + z * z - p * p) / (2.0 * r * z)
            if c >= 1.0:
                ang = 0.0
            elif c <= -1.0:
                ang = np.pi
            else:
                ang = np.arccos(c)
            mu = np.sqrt(max(0.0, 1.0 - r * r))
            I = 1.0 - u1 * (1.0 - mu) - u2 * (1.0 - mu) ** 2
            s += wg[k] * dr * I * 2.0 * r * ang / np.pi
        blocked += s
    return blocked / tot

@njit(cache=True)
def _kepE(M, e):
    E = M
    for _ in range(50):
        f = E - e * np.sin(E) - M
        E -= f / (1.0 - e * np.cos(E))
        if abs(f) < 1e-13:
            break
    return E

@njit(cache=True)
def flux_grid(tg, P, rsum, k, b0, J, ecosw, esinw, u1a, u2a, u1b, u2b, xg, wg):
    """EB flux (no third light; 1 = out of eclipse) at times tg measured from the primary conjunction."""
    n = tg.size
    out = np.empty(n)
    r1 = rsum / (1.0 + k); r2 = rsum * k / (1.0 + k)
    e = np.sqrt(ecosw * ecosw + esinw * esinw)
    w = np.arctan2(esinw, ecosw) if e > 0 else 0.5 * np.pi
    # b0 = a cos i/(R1+R2) -> cos i
    cosi = b0 * rsum
    sini = np.sqrt(max(0.0, 1.0 - cosi * cosi))
    nu0 = 0.5 * np.pi - w
    E0 = 2.0 * np.arctan(np.sqrt((1.0 - e) / (1.0 + e)) * np.tan(0.5 * nu0))
    M0 = E0 - e * np.sin(E0)
    L2 = k * k * J / (1.0 + k * k * J); L1 = 1.0 - L2
    for i in range(n):
        M = 2.0 * np.pi * tg[i] / P + M0
        if e > 1e-7:
            E = _kepE(M, e)
            nu = 2.0 * np.arctan2(np.sqrt(1.0 + e) * np.sin(0.5 * E), np.sqrt(1.0 - e) * np.cos(0.5 * E))
            rr = (1.0 - e * e) / (1.0 + e * np.cos(nu))
        else:
            nu = M; rr = 1.0
        th = nu + w
        st = np.sin(th)
        zs = rr * np.sqrt(max(0.0, 1.0 - sini * sini * st * st))
        f = 1.0
        if st > 0:   # star 2 in front: primary eclipse of star 1
            if zs < r1 + r2:
                f -= L1 * occ(zs / r1, k, u1a, u2a, xg, wg)
        else:
            if zs < r1 + r2:
                f -= L2 * occ(zs / r2, 1.0 / k, u1b, u2b, xg, wg)
        out[i] = f
    return out

def secondary_offset(P, ecosw, esinw):
    """time from primary to secondary conjunction (Sterne formula, exact for the conjunction geometry)."""
    e = np.hypot(ecosw, esinw); w = np.arctan2(esinw, ecosw) if e > 0 else 0.5 * np.pi
    def Mof(nu):
        E = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(nu / 2)); return E - e * np.sin(E)
    M1 = Mof(0.5 * np.pi - w); M2 = Mof(1.5 * np.pi - w)
    return ((M2 - M1) % (2 * np.pi)) / (2 * np.pi) * P

class GridModel:
    """Evaluate the EB model on a fine grid around both eclipses and integrate over each exposure time."""
    def __init__(self, P, hw=0.26, dtg=15.0 / 86400):
        self.P = P; self.hw = hw; self.dtg = dtg
        self.tg = np.arange(-hw, hw + dtg / 2, dtg)
    def eval(self, tphase, texp, P, rsum, k, b0, J, ecosw, esinw, u1a, u2a, u1b, u2b):
        """tphase: time since nearest primary conjunction modulo P, in [-P/4, 3P/4).  texp: exposure (d) per point
        (array, one value per point, few distinct values).  Returns EB flux without third light.
        The flux grid only spans the eclipse itself (half-width set by the radius sum, with margin); outside it the
        flux is exactly 1, so the cumulative integral is clamped correctly by np.interp."""
        dsec = secondary_offset(P, ecosw, esinw)
        e = np.hypot(ecosw, esinw)
        hwe = P / (2 * np.pi) * np.arcsin(min(1.0, 1.3 * rsum * (1 + e) / np.sqrt(1 - e * e))) + 0.003
        hwe = min(self.hw, hwe)
        tg = np.arange(-hwe, hwe + self.dtg / 2, self.dtg)
        out = np.ones(tphase.size)
        temax = np.max(texp)
        for centre in (0.0, dsec):
            tgc = tg + centre
            fg = flux_grid(tgc, P, rsum, k, b0, J, ecosw, esinw, u1a, u2a, u1b, u2b, XG, WG)
            cum = np.concatenate([[0.0], np.cumsum(0.5 * (fg[1:] + fg[:-1] - 2.0) * np.diff(tgc))])
            x = tphase - centre
            sel = np.abs(x) < hwe + 0.5 * temax + self.dtg
            if not np.any(sel):
                continue
            xs = x[sel]; te = texp[sel]
            hi = np.interp(xs + 0.5 * te, tg, cum); lo = np.interp(xs - 0.5 * te, tg, cum)
            with np.errstate(invalid='ignore', divide='ignore'):
                v = np.where(te > 0, (hi - lo) / np.where(te > 0, te, 1.0), np.interp(xs, tg, fg - 1.0))
            out[sel] += v
        return out
