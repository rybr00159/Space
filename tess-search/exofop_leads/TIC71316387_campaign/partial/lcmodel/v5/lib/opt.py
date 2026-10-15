"""Optimisation helpers: multi-start Nelder-Mead/Powell on the negative log posterior in a transformed space."""
import numpy as np
from scipy.optimize import minimize

def optimise(post, x0, steps, n_restart=3, maxiter=6000, rng=None, verbose=False, powell_fev=4000):
    """Repeated Nelder-Mead from x0 (with simplex scaled by steps), followed by Powell polish."""
    f = lambda x: -post(x)
    best_x = np.array(x0, float); best_f = f(best_x)
    if not np.isfinite(best_f):
        raise ValueError('start point outside prior')
    for r in range(n_restart):
        sim = [best_x] + [best_x + np.eye(len(x0))[i] * steps[i] * (1 if r % 2 == 0 else -1) for i in range(len(x0))]
        sim = np.array(sim)
        for i in range(1, len(sim)):
            if not np.isfinite(f(sim[i])): sim[i] = best_x - (sim[i] - best_x)
            if not np.isfinite(f(sim[i])): sim[i] = best_x + 0.1 * (sim[i] - best_x)
        o = minimize(f, best_x, method='Nelder-Mead', options=dict(initial_simplex=sim, maxiter=maxiter, maxfev=maxiter, xatol=1e-7, fatol=1e-4))
        if o.fun < best_f - 1e-6:
            best_x, best_f = o.x, o.fun
        if verbose: print(f'   restart {r}: {o.fun:.3f}', flush=True)
    o = minimize(f, best_x, method='Powell', options=dict(maxiter=50, maxfev=powell_fev, xtol=1e-6, ftol=1e-7))
    if o.fun < best_f: best_x, best_f = o.x, o.fun
    return best_x, best_f

def optimise_lbfgs(post, x0, steps, bounds=None, n_cycles=3, nm_polish=1500, eps=1e-3):
    """L-BFGS-B on the negative log posterior in a space scaled by `steps` (forward differences with step eps in
    scaled units), alternated with a short Nelder-Mead polish to escape kinks; returns (x, f)."""
    x0 = np.asarray(x0, float); s = np.asarray(steps, float)
    def g(u):
        v = -post(x0 + u * s)
        return v if np.isfinite(v) else 1e30
    if bounds is not None:
        bu = [((lo - a) / sc if lo is not None else None, (hi - a) / sc if hi is not None else None) for (lo, hi), a, sc in zip(bounds, x0, s)]
    else:
        bu = None
    u = np.zeros(x0.size); fbest = g(u)
    if fbest >= 1e30: raise ValueError('start point outside prior')
    for c in range(n_cycles):
        o = minimize(g, u, method='L-BFGS-B', bounds=bu, options=dict(eps=eps, maxiter=400, maxfun=20000, ftol=1e-10, gtol=1e-6))
        if o.fun < fbest: u, fbest = o.x, o.fun
        if nm_polish:
            sim = np.vstack([u] + [u + 0.3 * np.eye(u.size)[i] for i in range(u.size)])
            o2 = minimize(g, u, method='Nelder-Mead', options=dict(initial_simplex=sim, maxfev=nm_polish, xatol=1e-4, fatol=1e-4))
            if o2.fun < fbest - 1e-6: u, fbest = o2.x, o2.fun
            elif c > 0: break
    return x0 + u * s, fbest
