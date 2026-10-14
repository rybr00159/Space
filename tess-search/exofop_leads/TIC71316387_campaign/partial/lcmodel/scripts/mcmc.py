"""emcee sampling of the EB posterior.
Usage: mcmc.py <kind: geo|ms|msell> <nsteps> <tag>
geo  : dT0,dP,rsum,k,bn,J,L3,ecosw,esinw,q1a,q2a,q1b,q2b,s95,s105,s106   (uniform priors; isotropic cos i)
ms   : dT0,dP,M1,q,fR,dlnJ,bn,L3,ecosw,esinw,q1a,q2a,q1b,q2b,s95,s105,s106  (Pecaut&Mamajek MS relations)
msell: ms + ellipsoidal-amplitude likelihood term A_ell(obs) = -12 +- 42 ppm (sector-to-sector scatter, ooe_var)"""
import os
os.environ['NUMBA_NUM_THREADS'] = '1'; os.environ['OMP_NUM_THREADS'] = '1'; os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np, sys, json, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/scripts')
import fit3, fit2, lclib as L, ms_table as MS
import emcee
from multiprocessing import get_context
KIND = sys.argv[1]; NSTEP = int(sys.argv[2]); TAG = sys.argv[3] if len(sys.argv) > 3 else KIND
d = fit3.load_data(L.HERE + '/results/data_default.pkl')
SC = ['s95', 's105', 's106']
if KIND == 'geo':
    NAMES = ['dT0', 'dP', 'rsum', 'k', 'bn', 'J', 'L3', 'ecosw', 'esinw', 'q1a', 'q2a', 'q1b', 'q2b'] + SC
else:
    NAMES = ['dT0', 'dP', 'M1', 'q', 'fR', 'dlnJ', 'bn', 'L3', 'ecosw', 'esinw', 'q1a', 'q2a', 'q1b', 'q2b'] + SC
AELL_OBS, AELL_SIG = -12e-6, 42e-6
def alpha_ell(u=0.6, tau=0.4):
    return 0.15 * (15 + u) * (1 + tau) / (3 - u)
def full(x):
    p = dict(zip(NAMES, x))
    p['T0'] = L.T00 + p['dT0'] / 1440.0; p['P'] = L.P0 + p['dP'] * 1e-6
    if KIND != 'geo': p.update(fit3.ms_derived(p))
    p['cosi'] = p['bn'] * p['rsum']
    for s in fit3.SECS: p.setdefault(f's{s}', 1.0)
    return p
def inb(p):
    if KIND == 'geo':
        if not (0.03 < p['rsum'] < 0.7 and 0.05 < p['k'] < 2.0 and 0 < p['J'] < 1.6): return False
    else:
        if not (0.1 < p['M1'] < 1.4 and 0.09 / p['M1'] <= p['q'] <= 1.0 and 0.85 < p['fR'] < 1.3 and abs(p['dlnJ']) < 0.5): return False
    if not (0 <= p['bn'] < 1 and 0 <= p['L3'] < 0.999 and abs(p['ecosw']) < 0.05 and abs(p['esinw']) < 0.3 and abs(p['dT0']) < 30 and abs(p['dP']) < 50): return False
    for a in ('q1a', 'q2a', 'q1b', 'q2b'):
        if not (0 < p[a] < 1): return False
    for s in SC:
        if not (0.6 < p[s] < 1.4): return False
    return True
def aell(p):
    """predicted ellipsoidal semi-amplitude (fraction of total flux) [Morris & Naftilan 1993], both stars, diluted"""
    k2J = p['k'] ** 2 * p['J']; L2 = k2J / (1 + k2J); L1 = 1 - L2
    r1 = p['rsum'] / (1 + p['k']); r2 = p['rsum'] * p['k'] / (1 + p['k'])
    s2i = 1 - p['cosi'] ** 2
    a = alpha_ell()
    return (1 - p['L3']) * a * s2i * (L1 * p['q'] * r1 ** 3 + L2 / p['q'] * r2 ** 3)
def lnpost(x):
    p = full(x)
    if not inb(p): return -np.inf, np.nan
    pr = fit3.prior_resid(p, 'ms' if KIND != 'geo' else 'geo')
    lp = -0.5 * pr @ pr + np.log(p['rsum'])
    if KIND == 'ms' or KIND == 'msell':
        # uniform prior in M1 and in M2 (=q*M1): Jacobian dM2/dq = M1
        lp += np.log(p['M1'])
    ae = aell(p) if KIND != 'geo' else np.nan
    if KIND == 'msell':
        lp += -0.5 * ((ae - AELL_OBS) / AELL_SIG) ** 2
    c2 = fit3.chi2(d, p)
    return lp - 0.5 * c2, ae
if __name__ == '__main__':
    rng = np.random.default_rng(11)
    if KIND == 'geo':
        c0 = dict(dT0=0.0, dP=0.0, rsum=0.1685, k=0.419, bn=0.267, J=0.502, L3=0.966, ecosw=0.00088, esinw=0.02, q1a=0.39, q2a=0.17, q1b=0.39, q2b=0.9, s95=0.94, s105=0.93, s106=0.97)
        sd = dict(dT0=0.05, dP=0.2, rsum=0.005, k=0.03, bn=0.05, J=0.02, L3=0.005, ecosw=1e-4, esinw=0.01, q1a=0.05, q2a=0.05, q1b=0.05, q2b=0.05, s95=0.01, s105=0.01, s106=0.01)
    else:
        c0 = dict(dT0=0.0, dP=0.0, M1=0.55, q=0.4, fR=1.04, dlnJ=0.0, bn=0.3, L3=0.96, ecosw=0.00088, esinw=0.0, q1a=0.39, q2a=0.25, q1b=0.39, q2b=0.3, s95=0.94, s105=0.93, s106=0.97)
        sd = dict(dT0=0.05, dP=0.2, M1=0.1, q=0.1, fR=0.03, dlnJ=0.05, bn=0.1, L3=0.01, ecosw=1e-4, esinw=0.01, q1a=0.05, q2a=0.05, q1b=0.05, q2b=0.05, s95=0.01, s105=0.01, s106=0.01)
    if len(sys.argv) > 4:   # optional start from json of best fit
        c0.update({k: v for k, v in json.load(open(sys.argv[4])).items() if k in c0})
    nw = 64
    p0 = []
    while len(p0) < nw:
        x = np.array([c0[n] + sd[n] * rng.normal() for n in NAMES])
        if np.isfinite(lnpost(x)[0]): p0.append(x)
    p0 = np.array(p0)
    back = emcee.backends.HDFBackend(L.HERE + f'/results/chain_{TAG}.h5') if False else None
    with get_context("spawn").Pool(4) as pool:
        S = emcee.EnsembleSampler(nw, len(NAMES), lnpost, pool=pool, moves=[(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)],
                                  blobs_dtype=[('aell', float)])
        t0 = time.time()
        for i, st in enumerate(S.sample(p0, iterations=NSTEP)):
            if (i + 1) % 500 == 0:
                lp = S.get_log_prob()[-1]
                try: tau = S.get_autocorr_time(tol=0, discard=(i + 1) // 3)
                except Exception: tau = np.full(len(NAMES), np.nan)
                print(f'[{time.time()-t0:6.0f}s] step {i+1} lnp max {lp.max():.2f} med {np.median(lp):.2f} acc {S.acceptance_fraction.mean():.3f} tau max {np.nanmax(tau):.0f} ({NAMES[int(np.nanargmax(tau))]})', flush=True)
                np.savez_compressed(L.HERE + f'/results/chain_{TAG}.npz', chain=S.get_chain(), lnp=S.get_log_prob(), aell=S.get_blobs()['aell'], names=NAMES)
    np.savez_compressed(L.HERE + f'/results/chain_{TAG}.npz', chain=S.get_chain(), lnp=S.get_log_prob(), aell=S.get_blobs()['aell'], names=NAMES)
    print('done')
