"""Eclipse timing: per-event template fit (time shift + depth scale + quadratic baseline, linearised),
grid search in shift, parabolic refinement, formal error from delta-chi2=1 with per-sector sigma."""
import numpy as np, sys
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/final/lib')
import lclib as L, fit2
from ebmodel import eb_flux

def template(p, t, texp, nsub, scale=1.0):
    u1a, u2a = fit2.q2u(p['q1a'], p['q2a']); u1b, u2b = fit2.q2u(p['q1b'], p['q2b'])
    f = eb_flux(np.ascontiguousarray(t, dtype=float), p['T0'], p['P'], p['rsum'], p['k'], p['cosi'], p['J'], p['L3'], p['ecosw'], p['esinw'],
                u1a, u2a, u1b, u2b, texp, nsub, fit2.XG, fit2.WG)
    return 1 - scale * (1 - f)

def time_event(t, y, sig, tc, p, texp, nsub, scale=1.0, grid_min=12.0, step_min=0.1):
    """Return dict(dt_min, err_min, depth_scale, chi2min, curvature). tc = predicted centre (template eclipse time)."""
    dt = t - tc
    shifts = np.arange(-grid_min, grid_min + 1e-9, step_min) / 1440.0
    # template evaluated on a fine grid once, then interpolated for each shift
    tf = np.arange(dt.min() - grid_min / 1440 - 0.01, dt.max() + grid_min / 1440 + 0.01, 10 / 86400.0)
    mf = template(p, tf + tc, texp, nsub, scale)
    delta_f = 1 - mf
    w = 1 / sig
    X0 = np.vstack([np.ones_like(dt), dt, dt * dt]).T
    chi = np.empty(shifts.size); amp = np.empty(shifts.size)
    for i, s in enumerate(shifts):
        dlt = np.interp(dt - s, tf, delta_f)
        A = np.column_stack([X0, -dlt]) * w[:, None]
        c, *_ = np.linalg.lstsq(A, y * w, rcond=None)
        r = y * w - A @ c
        chi[i] = r @ r; amp[i] = c[3] / c[0]
    j = int(np.argmin(chi))
    if j == 0 or j == shifts.size - 1:
        return dict(dt_min=np.nan, err_min=np.nan, depth_scale=np.nan, chi2=np.nan, edge=True)
    # parabola through +-3 grid points
    sl = slice(max(0, j - 4), min(shifts.size, j + 5))
    cc = np.polyfit(shifts[sl] * 1440, chi[sl], 2)
    if cc[0] <= 0:
        return dict(dt_min=np.nan, err_min=np.nan, depth_scale=np.nan, chi2=np.nan, edge=True)
    x0 = -cc[1] / (2 * cc[0])
    err = 1 / np.sqrt(cc[0])
    return dict(dt_min=float(x0), err_min=float(err), depth_scale=float(amp[j]), chi2=float(chi[j]), npts=int(t.size), edge=False)
