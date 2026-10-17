"""r02: emcee transit fits (independent re-derivation of ../shape s03).
usage: r02_planet_mcmc.py DATASET(inv|ref) MODEL(free|single|twin) [beta: fake|white] [nsteps]
Parameters: T0, P, k (undiluted for twin), b in [0, 1+k), log10 rho_circ (free: U(-1, 2.3); single: N(9.41, 0.88);
twin: N(15.5, 1.45) with the light diluted by 2), u1, u2 with N(0.125, 0.10), N(0.533, 0.10) priors (Claret 2017,
as ../shape; twin uses the same) and physical LD bounds. Per-window offsets profiled. Max-likelihood start from a
multi-start Powell optimisation. Output: out/mcmc_{ds}_{model}_{beta}.npz/.json"""
import sys, json, numpy as np, time
sys.path.insert(0, 'tic214311908/analysis/referee2_shape/scripts')
from rlib import *
import emcee
from scipy.optimize import minimize
ds, model = sys.argv[1], sys.argv[2]
beta = sys.argv[3] if len(sys.argv) > 3 else 'fake'
nsteps = int(sys.argv[4]) if len(sys.argv) > 4 else 4000
D = load(ds, beta); lc = LC(D)
U1, U2, SU = 0.125, 0.533, 0.10
twin = model == 'twin'
def lnprior(th):
    T0, P, k, b, lr, u1, u2 = th
    if not (0.01 < k < 0.8 and 0 <= b < 1 + k and -1 < lr < 2.3 and abs(T0 - T0REF) < 0.02 and abs(P - PREF) < 3e-4): return -np.inf
    if not (u1 + u2 < 1 and u1 > 0 and u1 + 2*u2 > 0): return -np.inf
    lp = -0.5*((u1 - U1)/SU)**2 - 0.5*((u2 - U2)/SU)**2
    rho = 10**lr
    if model == 'single': lp += -0.5*((rho - 9.41)/0.88)**2 + np.log(rho)   # Gaussian in rho, sampled in log rho
    if model == 'twin': lp += -0.5*((rho - 15.5)/1.45)**2 + np.log(rho)
    return lp
def model_flux(th):
    T0, P, k, b, lr, u1, u2 = th
    m = lc.transit(T0, P, k, aR_from_rho(10**lr, P), b, u1, u2)
    return 1 + (m - 1)/2 if twin else m
def lnprob(th):
    lp = lnprior(th)
    if not np.isfinite(lp): return -np.inf, np.nan
    c = lc.chi2(model_flux(th))
    return lp - 0.5*c, c
t0 = time.time(); [lnprob([T0REF, PREF, 0.12, 0.5, 1.0, 0.125, 0.533]) for _ in range(20)]; dtm = (time.time() - t0)/20
print('time per call', dtm, 'n', lc.n, 'near', lc.near.sum(), D['info'], flush=True)
best = None
for b0 in (0.1, 0.5, 0.8):
    for lr0 in (0.7, 1.0, 1.2):
        x0 = [T0REF, PREF, 0.12*(1.41 if twin else 1), b0, lr0, U1, U2]
        r = minimize(lambda x: -lnprob(x)[0] if np.isfinite(lnprob(x)[0]) else 1e12, x0, method='Nelder-Mead', options=dict(maxiter=6000, xatol=1e-7, fatol=1e-4))
        if best is None or r.fun < best.fun: best = r
print('ML', best.x, best.fun, flush=True)
nw = 32; ndim = 7
p0 = best.x + np.array([3e-4, 2e-6, 2e-3, 0.03, 0.03, 0.02, 0.02])*np.random.default_rng(1).normal(size=(nw, ndim))
p0[:, 3] = np.abs(p0[:, 3])
sampler = emcee.EnsembleSampler(nw, ndim, lnprob)
sampler.run_mcmc(p0, nsteps, progress=False)
ch = sampler.get_chain(discard=nsteps//3, thin=5, flat=True); c2 = sampler.get_blobs(discard=nsteps//3, thin=5, flat=True)
try: tau = sampler.get_autocorr_time(discard=nsteps//3, quiet=True).tolist()
except Exception: tau = None
T0, P, k, b, lr, u1, u2 = ch.T
rho = 10**lr; aR = aR_from_rho(rho, P)
T14, T23 = t14(P, k, aR, b)
inc = np.degrees(np.arccos(b/aR))
# observed depth of the model at mid-transit (approximate via batman at t=T0 for a 200-s exposure is not needed: use analytic small-planet LD)
graz = b > 1 - k
np.savez(f'{RD}/out/mcmc_{ds}_{model}_{beta}.npz', chain=ch, chi2=c2, T14=T14, T23=T23, aR=aR, inc=inc)
def q(x):
    p = np.percentile(x, [0.135, 2.275, 15.865, 50, 84.135, 97.725, 99.865]); return [float(v) for v in p]
out = dict(dataset=ds, model=model, beta=beta, n=int(lc.n), info=D['info'], chi2_min_chain=float(np.nanmin(c2)), chi2_ML=float(lc.chi2(model_flux(best.x))), x_ML=best.x.tolist(),
           chi2_flat=float(lc.chi2(np.ones(lc.n))), accept=float(np.mean(sampler.acceptance_fraction)), tau=tau, frac_grazing=float(np.mean(graz)),
           pct_note='percentiles 0.135, 2.275, 15.865, 50, 84.135, 97.725, 99.865',
           post={nm: q(v) for nm, v in dict(T0=T0, P=P, k=k, b=b, rho=rho, u1=u1, u2=u2, aR=aR, inc=inc, T14_h=T14*24, T23_h=T23*24, flat=T23/np.where(T14 > 0, T14, np.nan), ingress_min=(T14 - T23)/2*1440).items()})
json.dump(out, open(f'{RD}/out/mcmc_{ds}_{model}_{beta}.json', 'w'), indent=1)
for kk, v in out['post'].items(): print(kk, [round(x, 5) for x in v[2:5]])
print('chi2 ML', out['chi2_ML'], 'flat', out['chi2_flat'], 'grazing frac', out['frac_grazing'], 'tau', tau)
