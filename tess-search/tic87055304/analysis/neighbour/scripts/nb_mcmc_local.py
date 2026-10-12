"""Neighbour-host transit MCMC on the PSF pair light curve with LOCAL baselines (no global biweight detrending).

Data: data/lc_pair.npy f_raw (normalised, undetrended PSF pair flux F_t + F_n), cadences within +-0.45 d of each
transit; a transit is used if it has >= 3 in-transit (|tau| < 0.08 d) and >= 6 out-of-transit (|tau| > 0.15 d) points.
Model: f = [transit(T0, P, k, b, rho, u1, u2; dilution D)] x (c0 + c1 x + c2 x^2)_epoch, x = t - <t>_epoch.
The 3 baseline coefficients of every transit are profiled (solved by weighted linear least squares at each step,
batched 3x3 normal equations), so the posterior is the profile posterior over the baselines.
Free: T0, P, k, b, rho, u1, u2 (Claret prior as nbmodel), ln sigma_s (5 sectors), D ~ N(5.32, 0.06).
Exposure integration, priors and convergence rules as nb_mcmc.py (nbmodel.py).
Usage: python nb_mcmc_local.py NAME [--lc path] [--dilprior m,s] [--sectors 13,39,...]
Output: mcmc/mcmc_NAME.json, mcmc/mcmc_NAME.npz
"""
import os, sys, json, time, argparse
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import emcee
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nbmodel as tm

HERE = "/home/user/work/tic87/neighbour/mcmc"
ap = argparse.ArgumentParser()
ap.add_argument("name"); ap.add_argument("--lc", default="/home/user/work/tic87/neighbour/data/lc_pair.npy")
ap.add_argument("--dilprior", default="5.32,0.06"); ap.add_argument("--sectors", default="13,39,66,93,104")
ap.add_argument("--maxsteps", type=int, default=20000); ap.add_argument("--seed", type=int, default=4)
a = ap.parse_args()
SECT = [int(s) for s in a.sectors.split(",")]
Dm, Ds = (float(v) for v in a.dilprior.split(","))
P0, T00 = 4.1362006, 3271.5326
lc = np.load(a.lc, allow_pickle=True).item()
data = {}
for s in SECT:
    v = lc[s]; t = v["t"]; y = v["f_raw"]
    n = np.round((t - T00) / P0); tau = t - T00 - n * P0
    keep = np.zeros(len(t), bool); ep = np.full(len(t), -1)
    for i, e in enumerate(np.unique(n[np.abs(tau) < 0.45])):
        m = (n == e) & (np.abs(tau) < 0.45)
        if np.sum(m & (np.abs(tau) < 0.08)) < 3 or np.sum(m & (np.abs(tau) > 0.15)) < 6: continue
        keep |= m; ep[m] = e
    t, y, ep = t[keep], y[keep], ep[keep]
    o = np.lexsort((t, ep)); t, y, ep = t[o], y[o], ep[o]
    starts = np.r_[0, np.where(np.diff(ep) != 0)[0] + 1]
    xm = np.zeros_like(t)
    for j, s0 in enumerate(starts):
        s1 = starts[j + 1] if j + 1 < len(starts) else len(t)
        xm[s0:s1] = t[s0:s1] - t[s0:s1].mean()
        y[s0:s1] = y[s0:s1] / np.median(y[s0:s1])
    data[s] = dict(t=t, f=y, x=xm, starts=starts, exptime=tm.EXPT[s], sigma0=v["sigma"], n_ep=len(starts))
print({s: (len(d["t"]), d["n_ep"]) for s, d in data.items()}, flush=True)
NAMES = ["T0", "P", "k", "b", "rho", "u1", "u2"] + [f"lns{s}" for s in SECT] + ["D"]


def profile_chi2(d, mod):
    """min over per-epoch quadratic baselines of sum (y - mod*poly)^2; returns (rss, fitted model)."""
    A0, A1, A2 = mod, mod * d["x"], mod * d["x"] ** 2
    cols = [A0, A1, A2]
    st = d["starts"]
    G = np.empty((len(st), 3, 3)); r = np.empty((len(st), 3))
    for i in range(3):
        r[:, i] = np.add.reduceat(cols[i] * d["f"], st)
        for j in range(i, 3):
            G[:, i, j] = G[:, j, i] = np.add.reduceat(cols[i] * cols[j], st)
    c = np.linalg.solve(G, r[..., None])[..., 0]
    idx = np.repeat(np.arange(len(st)), np.diff(np.r_[st, len(d["t"])]))
    fit = A0 * c[idx, 0] + A1 * c[idx, 1] + A2 * c[idx, 2]
    return np.sum((d["f"] - fit) ** 2), fit


def logpost(th, return_fit=False):
    T0, P, k, b, rho, u1, u2 = th[:7]
    lns = th[7:7 + len(SECT)]; D = th[-1]
    if not (0.005 < k < 0.5 and 0 <= b < 1 + k and 4.0 < P < 4.3 and 0.003 < rho < 5 and tm.ld_ok(u1, u2)):
        return -np.inf
    if np.any(lns < -10) or np.any(lns > -3) or not (0 <= D < 100) or tm.a_from_rho(rho, P) < 1.5:
        return -np.inf
    lp = -0.5 * ((u1 - tm.U1_MU) / tm.U1_SD) ** 2 - 0.5 * ((u2 - tm.U2_MU) / tm.U2_SD) ** 2 - 0.5 * ((D - Dm) / Ds) ** 2
    mods = tm.phys_model(th[:7], data, D)
    ll = 0.0; fits = {}
    for i, s in enumerate(SECT):
        rss, fit = profile_chi2(data[s], mods[s])
        s2 = np.exp(2 * lns[i])
        ll += -0.5 * (rss / s2 + len(data[s]["t"]) * np.log(2 * np.pi * s2))
        fits[s] = fit
    if return_fit:
        return lp + ll, fits
    return lp + ll


x0 = np.array([T00, P0, 0.141, 0.56, 0.28, tm.U1_MU, tm.U2_MU] + [np.log(data[s]["sigma0"]) for s in SECT] + [Dm])
from scipy.optimize import minimize
f = lambda x: -logpost(x) if np.isfinite(logpost(x)) else 1e25
for it in range(3):
    x0 = minimize(f, x0, method="Nelder-Mead", options=dict(maxiter=8000, maxfev=8000, adaptive=True)).x
print("start", dict(zip(NAMES, np.round(x0, 6))), logpost(x0), flush=True)
rng = np.random.default_rng(a.seed)
nd = len(x0); nw = 48
scale = np.array([3e-4, 3e-7, 1e-3, 0.01, 0.002, 0.02, 0.02] + [0.01] * len(SECT) + [0.02])
p0 = []
while len(p0) < nw:
    y = x0 + scale * rng.standard_normal(nd)
    if np.isfinite(logpost(y)): p0.append(y)
sam = emcee.EnsembleSampler(nw, nd, logpost, moves=[(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)])
t0 = time.time(); state = np.array(p0); n = 0; old = np.inf; taus = []
while n < a.maxsteps:
    state = sam.run_mcmc(state, 2000, progress=False); n += 2000
    tau = sam.get_autocorr_time(tol=0); taus.append((n, float(np.max(tau))))
    print(f"{n} steps, max tau {np.max(tau):.0f} ({NAMES[int(np.argmax(tau))]}), acc {np.mean(sam.acceptance_fraction):.3f}, {time.time() - t0:.0f}s", flush=True)
    if n >= 6000 and n > 100 * np.max(tau) and np.all(np.abs(old - tau) / tau < 0.05): break
    old = tau
tau = sam.get_autocorr_time(tol=0)
burn = int(max(5 * np.max(tau), 0.2 * n)); thin = max(1, int(np.max(tau) / 2))
ch = sam.get_chain(discard=burn, thin=thin, flat=True); lp = sam.get_log_prob(discard=burn, thin=thin, flat=True)
der = tm.derived(ch)
sub = ch[rng.choice(len(ch), size=min(3000, len(ch)), replace=False)]
dep_u = np.array([tm.depth_of(th[:7], 0.0) for th in sub]); dep = np.array([tm.depth_of(th[:7], th[-1]) for th in sub])


def summ(x):
    q = np.percentile(x, [2.275, 15.865, 50, 84.135, 97.725])
    return dict(med=float(q[2]), lo=float(q[2] - q[1]), hi=float(q[3] - q[2]), p2=float(q[0]), p97=float(q[4]))


res = dict(name=a.name, lc=a.lc, sectors=SECT, dilprior=[Dm, Ds], baseline="per-transit quadratic, profiled",
           n_transits={s: int(data[s]["n_ep"]) for s in SECT}, ndata=int(sum(len(d["t"]) for d in data.values())),
           nwalkers=nw, nsteps=n, burn=burn, thin=thin, n_samples=len(ch), acceptance=float(np.mean(sam.acceptance_fraction)),
           tau=dict(zip(NAMES, map(float, tau))), n_over_tau=float(n / np.max(tau)), tau_history=taus,
           runtime_s=time.time() - t0, E_ref=390, params={})
for i, nm in enumerate(NAMES): res["params"][nm] = summ(ch[:, i])
for nm, v in der.items(): res["params"][nm] = summ(v)
res["params"]["depth_ppm"] = summ(dep); res["params"]["depth_undiluted_ppm"] = summ(dep_u)
ib = int(np.argmax(lp)); res["best"] = dict(zip(NAMES, map(float, ch[ib]))); res["best_logpost"] = float(lp[ib])
np.savez_compressed(os.path.join(HERE, f"mcmc_{a.name}.npz"), chain=ch, lp=lp, names=NAMES)
json.dump(res, open(os.path.join(HERE, f"mcmc_{a.name}.json"), "w"), indent=1)
p = res["params"]
print("DONE", a.name, {k: (round(p[k]["med"], 5), round(p[k]["lo"], 5), round(p[k]["hi"], 5)) for k in
                       ["k", "b", "rho", "a", "depth_undiluted_ppm", "T14_h", "T23_h", "D"]})
