"""emcee + batman MCMC for TIC 87055304.

Usage: python run_mcmc.py NAME [--rho uniform,lo,hi] [--dil D] [--bmode free|grazing|nongrazing]
                               [--drop S] [--only S] [--lc path] [--maxsteps N] [--nwalkers W] [--noldprior]
Writes mcmc_NAME.npz (flat post-burn-in chain + log-prob) and mcmc_NAME.json (summary, convergence diagnostics).
Convergence: run in blocks of 2000 steps until N_steps > 100 * max(tau) and tau changed < 5 % (free and dwarf runs: 2 %) (emcee
recommendation), or maxsteps. Burn-in = 5 * max(tau) (at least 20 % of the chain); thin = max(1, tau/2).
"""
import os, sys, json, time, argparse
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import emcee
import tmodel as tm

HERE = os.path.dirname(os.path.abspath(__file__))


def summarize(x):
    q = np.percentile(x, [2.275, 15.865, 50, 84.135, 97.725])
    return dict(med=float(q[2]), lo=float(q[2] - q[1]), hi=float(q[3] - q[2]), p2=float(q[0]), p97=float(q[4]),
                mean=float(np.mean(x)), sd=float(np.std(x)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("name"); ap.add_argument("--rho", default="uniform,0.003,5.0"); ap.add_argument("--dil", type=float, default=0.0)
    ap.add_argument("--bmode", default="free"); ap.add_argument("--drop", type=int, default=None); ap.add_argument("--only", type=int, default=None)
    ap.add_argument("--lc", default=None); ap.add_argument("--maxsteps", type=int, default=40000); ap.add_argument("--nwalkers", type=int, default=48)
    ap.add_argument("--noldprior", action="store_true"); ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--start", default="ml_start.json"); ap.add_argument("--minsteps", type=int, default=6000)
    ap.add_argument("--oddeven", action="store_true")
    a = ap.parse_args()
    kind, lo, hi = a.rho.split(","); rho_prior = (kind, float(lo), float(hi))
    st = json.load(open(os.path.join(HERE, a.start)))
    sectors = [s for s in tm.SECTORS if s != a.drop]
    if a.only is not None:
        sectors = [a.only]
    data = tm.load(a.lc, sectors=sectors)
    if a.oddeven:
        post = tm.PosteriorOE(data, st["E_ref"], rho_prior=rho_prior, dil=a.dil, b_mode=a.bmode, ld_prior=not a.noldprior, sectors=sectors)
    else:
        post = tm.Posterior(data, rho_prior=rho_prior, dil=a.dil, b_mode=a.bmode, ld_prior=not a.noldprior, sectors=sectors)
    x0full = dict(zip(st["names"], st["x"]))
    x0full["k_even"] = x0full["k"]
    x0 = np.array([x0full[n] for n in post.names()])
    rng = np.random.default_rng(a.seed)
    # make the start point valid for the chosen prior
    if kind in ("uniform", "loguniform") and not (rho_prior[1] < x0[4] < rho_prior[2]):
        x0[4] = np.sqrt(rho_prior[1] * rho_prior[2])
    if kind == "gauss":
        x0[4] = rho_prior[1]
    # adjust b so the start duration is roughly right: try a few b values, keep the best posterior
    best = None
    for bb in np.linspace(0.0, 1.05, 43):
        y = x0.copy(); y[3] = bb
        for kk in (0.055, 0.06, 0.065, 0.07, 0.08, 0.09):
            y[2] = kk
            if a.oddeven:
                y[-1] = kk
            lp = post(y)
            if best is None or lp > best[0]:
                best = (lp, y.copy())
    x0 = best[1]
    from scipy.optimize import minimize
    f = lambda x: -post(x) if np.isfinite(post(x)) else 1e25
    for it in range(3):
        x0 = minimize(f, x0, method="Nelder-Mead", options=dict(maxiter=8000, maxfev=8000, adaptive=True)).x
    print("start", dict(zip(post.names(), np.round(x0, 6))), post(x0), flush=True)
    nd = len(x0); nw = max(a.nwalkers, 2 * nd + 2)
    scale = np.array([3e-4, 3e-7, 1e-3, 0.01, 0.002, 0.02, 0.02] + [1e-4] * len(sectors) + [0.01] * len(sectors) + ([1e-3] if a.oddeven else []))
    p0 = []
    while len(p0) < nw:
        y = x0 + scale * rng.standard_normal(nd)
        if np.isfinite(post(y)):
            p0.append(y)
    sampler = emcee.EnsembleSampler(nw, nd, post, moves=[(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)])
    t0 = time.time(); old_tau = np.inf; state = np.array(p0); n = 0; taus = []
    while n < a.maxsteps:
        state = sampler.run_mcmc(state, 2000, progress=False)
        n += 2000
        tau = sampler.get_autocorr_time(tol=0)
        taus.append((n, float(np.max(tau))))
        conv = (n > 50 * np.max(tau)) and np.all(np.abs(old_tau - tau) / tau < 0.05) and n >= a.minsteps and n > 100 * np.max(tau)
        print(f"{n} steps, max tau {np.max(tau):.0f} ({post.names()[int(np.argmax(tau))]}), acc {np.mean(sampler.acceptance_fraction):.3f}, "
              f"{time.time() - t0:.0f}s", flush=True)
        if conv:
            break
        old_tau = tau
    tau = sampler.get_autocorr_time(tol=0)
    burn = int(max(5 * np.max(tau), 0.2 * n)); thin = max(1, int(np.max(tau) / 2))
    chain = sampler.get_chain(discard=burn, thin=thin, flat=True)
    lp = sampler.get_log_prob(discard=burn, thin=thin, flat=True)
    names = post.names()
    # derived
    der = tm.derived(chain)
    sub = chain[rng.choice(len(chain), size=min(4000, len(chain)), replace=False)]
    depth = np.array([tm.depth_of(th[:7], a.dil) for th in sub])
    res = dict(name=a.name, rho_prior=rho_prior, dil=a.dil, b_mode=a.bmode, sectors=sectors, ld_prior=not a.noldprior,
               E_ref=st["E_ref"], nwalkers=nw, nsteps=n, burn=burn, thin=thin, n_samples=len(chain),
               acceptance=float(np.mean(sampler.acceptance_fraction)), tau=dict(zip(names, map(float, tau))),
               n_over_tau=float(n / np.max(tau)), converged=bool(n > 50 * np.max(tau)), runtime_s=time.time() - t0,
               tau_history=taus, params={})
    for i, nm in enumerate(names):
        res["params"][nm] = summarize(chain[:, i])
    for nm, v in der.items():
        res["params"][nm] = summarize(v)
    res["params"]["depth_ppm"] = summarize(depth)
    res["params"]["k2_ppm"] = summarize(chain[:, 2] ** 2 * 1e6)
    if a.oddeven:
        dk = chain[:, 2] - chain[:, -1]
        res["params"]["k_odd_minus_even"] = summarize(dk)
        res["params"]["depth_ratio_odd_even"] = summarize((chain[:, 2] / chain[:, -1]) ** 2)
        res["oddeven_diff_sigma"] = float(np.median(dk) / np.std(dk))
        res["oddeven_frac_dk_gt0"] = float(np.mean(dk > 0))
    for i, s in enumerate(sectors):
        res["params"][f"sigma{s}_ppm"] = summarize(np.exp(chain[:, 7 + len(sectors) + i]) * 1e6)
    # max-likelihood sample
    ll = np.array([post.log_like(th) for th in chain[np.argsort(lp)[-200:]]])
    ib = np.argsort(lp)[-200:][int(np.argmax(ll))]
    res["best"] = dict(zip(names, map(float, chain[ib])))
    llb, chi2b = post.log_like(chain[ib], return_chi2=True)
    res["best_lnL"] = float(llb); res["best_chi2"] = float(chi2b); res["ndata"] = int(sum(len(d["t"]) for d in data.values()))
    np.savez_compressed(os.path.join(HERE, f"mcmc_{a.name}.npz"), chain=chain, lp=lp, names=names)
    json.dump(res, open(os.path.join(HERE, f"mcmc_{a.name}.json"), "w"), indent=1)
    p = res["params"]
    print("DONE", a.name, {k: (round(p[k]["med"], 5), round(p[k]["lo"], 5), round(p[k]["hi"], 5)) for k in ["k", "b", "rho", "a", "inc", "depth_ppm", "T14_h", "T23_h", "u1", "u2"]})


if __name__ == "__main__":
    main()
