"""Two-sphere eclipsing-binary light-curve model (no tidal distortion; spheres, quadratic LD).
occ_frac(z,p,u1,u2): fraction of the flux of a quadratically limb-darkened star of unit radius hidden
by an opaque disc of radius p at projected centre separation z (both in units of the occulted star's
radius).  Computed by radial integration: annulus of radius r is covered over an angle
2*arccos((r^2+z^2-p^2)/(2 r z)); Gauss-Legendre with cosine node clustering on the partial range
plus the analytic fully covered core.  Valid for any p (also p>1, needed for the secondary eclipse).
Validated against batman (p<1) and brute-force 2-D integration (p>1) in test_model.py."""
import numpy as np
from numba import njit, prange

_NG = 64
_xg, _wg = np.polynomial.legendre.leggauss(_NG)

@njit(cache=True)
def _Icum(r, u1, u2):
    # integral_0^r I(s) 2 s ds  for I(mu)=1-u1(1-mu)-u2(1-mu)^2, mu=sqrt(1-s^2)
    # use substitution mu: integral = int_{mu(r)}^{1} I(mu) 2 mu dmu
    m = np.sqrt(max(0.0, 1.0 - r * r))
    # I(mu) 2 mu = 2mu - 2u1 mu (1-mu) - 2u2 mu (1-mu)^2
    def F(x):
        return x * x - 2 * u1 * (x * x / 2 - x ** 3 / 3) - 2 * u2 * (x * x / 2 - 2 * x ** 3 / 3 + x ** 4 / 4)
    return F(1.0) - F(m)

@njit(cache=True)
def occ_frac(z, p, u1, u2, xg, wg):
    if z >= 1.0 + p or p <= 0.0:
        return 0.0
    tot = _Icum(1.0, u1, u2)
    blocked = 0.0
    # fully covered core r < p - z
    rin = p - z
    if rin >= 1.0:
        return 1.0
    if rin > 0:
        blocked += _Icum(rin, u1, u2)
    ra = abs(z - p)
    rb = min(1.0, z + p)
    if rb > ra:
        s = 0.0
        for k in range(xg.size):
            th = 0.5 * np.pi * (xg[k] + 1.0)  # 0..pi
            r = ra + (rb - ra) * 0.5 * (1.0 - np.cos(th))
            dr = (rb - ra) * 0.5 * np.sin(th) * 0.5 * np.pi
            if r <= 0:
                continue
            c = (r * r + z * z - p * p) / (2 * r * z) if z > 0 else -2.0
            if c >= 1.0:
                ang = 0.0
            elif c <= -1.0:
                ang = np.pi
            else:
                ang = np.arccos(c)
            mu = np.sqrt(max(0.0, 1.0 - r * r))
            I = 1.0 - u1 * (1 - mu) - u2 * (1 - mu) ** 2
            s += wg[k] * dr * I * 2 * r * ang / np.pi
        blocked += s
    return blocked / tot

@njit(cache=True)
def kepler_E(M, e):
    E = M + e * np.sin(M) if e < 0.8 else np.pi + 0.0 * M
    for _ in range(30):
        f = E - e * np.sin(E) - M
        E -= f / (1 - e * np.cos(E))
        if abs(f) < 1e-12:
            break
    return E

@njit(cache=True, parallel=True)
def eb_flux(t, T0, P, rsum, k, cosi, J, L3, ecosw, esinw, u1a, u2a, u1b, u2b, texp, nsub, xg, wg):
    """Normalised flux (1 = out of eclipse) of EB+third light.  Star 1 is eclipsed at T0 (primary).
    rsum=(R1+R2)/a, k=R2/R1, J = ratio of disc-averaged surface brightness (star2/star1) in band."""
    n = t.size
    out = np.empty(n)
    r1 = rsum / (1 + k); r2 = rsum * k / (1 + k)
    e = np.sqrt(ecosw ** 2 + esinw ** 2)
    w = np.arctan2(esinw, ecosw) if e > 0 else np.pi / 2
    sini = np.sqrt(1 - cosi * cosi)
    # time of periastron from primary-eclipse conjunction (true anomaly nu0 = pi/2 - w)
    nu0 = np.pi / 2 - w
    E0 = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(nu0 / 2))
    M0 = E0 - e * np.sin(E0)
    L2 = k * k * J / (1 + k * k * J); L1 = 1 - L2
    for i in prange(n):
        acc = 0.0
        for j in range(nsub):
            tt = t[i] + texp * ((j + 0.5) / nsub - 0.5)
            M = 2 * np.pi * (tt - T0) / P + M0
            M = M % (2 * np.pi)
            if e > 1e-6:
                E = kepler_E(M, e)
                nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
                rr = (1 - e * e) / (1 + e * np.cos(nu))
            else:
                nu = M; rr = 1.0
            th = nu + w
            # sky-projected separation in units of a
            zs = rr * np.sqrt(1 - sini * sini * np.sin(th) ** 2)
            front2 = np.sin(th) > 0   # star 2 in front -> star 1 eclipsed (primary)
            f = 1.0
            if front2:
                f -= L1 * occ_frac(zs / r1, r2 / r1, u1a, u2a, xg, wg)
            else:
                f -= L2 * occ_frac(zs / r2, r1 / r2, u1b, u2b, xg, wg)
            acc += f
        out[i] = (1 - L3) * acc / nsub + L3
    return out

def model(t, p, texp, nsub):
    return eb_flux(t, p['T0'], p['P'], p['rsum'], p['k'], p['cosi'], p['J'], p['L3'], p.get('ecosw', 0.0), p.get('esinw', 0.0),
                   p['u1a'], p['u2a'], p['u1b'], p['u2b'], texp, nsub, _xg, _wg)
