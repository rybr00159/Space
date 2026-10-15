"""Per-event eclipse timing with a fixed-shape template (global best-fit physical model, product-specific depth
scale) and a per-event multiplicative quadratic baseline profiled analytically.
measure(): chi2 on a grid of shifts (-25..+25 min, 0.2 min steps), parabolic refinement around the minimum,
formal sigma from the curvature (Delta chi2 = 1)."""
import numpy as np, sys
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6/lib')
import ebm

SHIFTS = np.arange(-25.0, 25.0001, 0.2) / 1440.0

class Template:
    def __init__(self, p, scale_of_key=None):
        self.p = dict(p); self.gm = ebm.GridModel(p['P'])
        self.scale = scale_of_key or {}
        self.dsec = ebm.secondary_offset(p['P'], p['ecosw'], p['esinw'])
        self._cache = {}
    def _profile(self, texp1, ph):
        k = (int(round(texp1 * 86400)), ph)
        if k not in self._cache:
            p = self.p
            xg = np.arange(-0.24, 0.24 + 1e-9, 4.0 / 86400)
            tph = xg + (self.dsec if ph == 0.5 else 0.0)
            f = self.gm.eval(tph, np.full(xg.size, texp1), p['P'], p['rsum'], p['k'], p['b0'], p['J'], p['ecosw'], p['esinw'],
                             p['u1a'], p['usa'] - p['u1a'], p['u1b'], p['usb'] - p['u1b'])
            self._cache[k] = (xg, (1 - p['L3']) * (1 - f))
        return self._cache[k]
    def depth(self, x, texp, ph, key):
        """eclipse depth profile (1 - normalised flux, diluted, product-scaled) at times x relative to the
        conjunction of the eclipse of type ph (0 primary, 0.5 secondary); texp is a scalar or a constant array"""
        te = float(np.atleast_1d(texp)[0])
        xg, dg = self._profile(te, ph)
        return self.scale.get(key, 1.0) * np.interp(x, xg, dg, left=0.0, right=0.0)

def _chi2_lin(y, m, dt, w):
    A = np.vstack([m, m * dt, m * dt * dt]).T
    Aw = A * w[:, None]
    c = np.linalg.solve(A.T @ Aw, Aw.T @ y)
    r = y - A @ c
    return float((w * r * r).sum())

def measure(t, y, sig, texp, tc, ph, key, tmpl, shifts=SHIFTS, depth_free=False):
    """time of the eclipse of type ph near predicted conjunction tc.  Returns (delta [d], sigma_formal [d], chi2min, n)"""
    w = 1.0 / sig ** 2
    dt = t - tc
    chi = np.empty(shifts.size)
    for i, s in enumerate(shifts):
        dep = tmpl.depth(dt - s, texp, ph, key)
        if depth_free:
            A = np.vstack([np.ones_like(dt), dt, dt * dt, -dep]).T
            Aw = A * w[:, None]; c = np.linalg.solve(A.T @ Aw, Aw.T @ y); r = y - A @ c; chi[i] = (w * r * r).sum()
        else:
            chi[i] = _chi2_lin(y, 1.0 - dep, dt, w)
    j = int(np.argmin(chi))
    if j < 3 or j > shifts.size - 4:
        return np.nan, np.nan, float(chi[j]), t.size
    sl = slice(j - 3, j + 4)
    c = np.polyfit(shifts[sl] - shifts[j], chi[sl], 2)
    if c[0] <= 0:
        return np.nan, np.nan, float(chi[j]), t.size
    d0 = -c[1] / (2 * c[0]) + shifts[j]
    return float(d0), float(1.0 / np.sqrt(c[0])), float(chi[j]), t.size
