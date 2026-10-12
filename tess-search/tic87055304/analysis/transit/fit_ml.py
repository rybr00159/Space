"""Maximum-likelihood fit (global optimiser + Nelder-Mead polish) of the batman model; picks the reference epoch
that minimises the T0-P covariance; writes ml_start.json used to initialise the MCMC runs."""
import os, json, sys
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.optimize import minimize
import tmodel as tm

HERE = os.path.dirname(os.path.abspath(__file__))


def ml_fit(post, x0, n_restart=4, seed=1):
    rng = np.random.default_rng(seed)
    best = None
    f = lambda x: -post(x) if np.isfinite(post(x)) else 1e25
    for r in range(n_restart):
        x = np.array(x0, float)
        if r > 0:
            x[2] *= 1 + 0.05 * rng.standard_normal(); x[3] = np.clip(x[3] + 0.1 * rng.standard_normal(), 0.0, 1.0)
            x[4] = np.clip(x[4] * np.exp(0.3 * rng.standard_normal()), 0.01, 3)
        if not np.isfinite(post(x)):
            continue
        for it in range(3):
            res = minimize(f, x, method="Nelder-Mead", options=dict(maxiter=20000, maxfev=20000, xatol=1e-8, fatol=1e-4, adaptive=True))
            x = res.x
        res = minimize(f, x, method="Powell", options=dict(maxiter=20000, xtol=1e-6, ftol=1e-8))
        if best is None or res.fun < best.fun:
            best = res
    return best


def numeric_cov(post, x, idx, h):
    """Covariance of parameters idx from the numerical Hessian of -log post."""
    n = len(idx)
    H = np.zeros((n, n))
    f0 = post(x)
    for i in range(n):
        for j in range(i, n):
            xi = x.copy()
            def fx(di, dj):
                y = x.copy(); y[idx[i]] += di; y[idx[j]] += dj; return post(y)
            hi, hj = h[i], h[j]
            if i == j:
                H[i, i] = -(fx(hi, 0) - 2 * f0 + fx(-hi, 0)) / hi ** 2
            else:
                H[i, j] = H[j, i] = -(fx(hi, hj) - fx(hi, -hj) - fx(-hi, hj) + fx(-hi, -hj)) / (4 * hi * hj)
    return np.linalg.inv(H)


if __name__ == "__main__":
    data = tm.load()
    post = tm.Posterior(data, rho_prior=("uniform", 0.003, 5.0))
    E0 = 322
    x0 = np.array([tm.T0_INIT + E0 * tm.P_INIT, tm.P_INIT, 0.067, 0.85, 0.063, tm.U1_MU, tm.U2_MU] + [0.0] * 5 +
                  [np.log(data[s]["sigma0"]) for s in tm.SECTORS])
    best = ml_fit(post, x0)
    x = best.x
    print("ML:", dict(zip(post.names(), np.round(x, 6))), "logpost", -best.fun)
    cov = numeric_cov(post, x, [0, 1, 2, 3, 4], [2e-4, 2e-7, 5e-4, 5e-3, 2e-3])
    sd = np.sqrt(np.diag(cov))
    print("Laplace sd T0,P,k,b,rho:", sd)
    # epoch shift that removes the T0-P covariance: T0' = T0 + dE P, cov' = cov_TP + dE var_P = 0
    dE = -cov[0, 1] / cov[1, 1]
    E_ref = int(np.round(E0 + dE))
    print("covariance-minimising epoch offset %.2f -> E_ref = %d (from T0_INIT)" % (dE, E_ref))
    x[0] = x[0] + (E_ref - E0) * x[1]
    out = dict(E_ref=E_ref, x=x.tolist(), names=post.names(), laplace_sd=dict(zip(["T0", "P", "k", "b", "rho"], sd.tolist())),
               logpost=float(-best.fun))
    json.dump(out, open(os.path.join(HERE, "ml_start.json"), "w"), indent=1)
