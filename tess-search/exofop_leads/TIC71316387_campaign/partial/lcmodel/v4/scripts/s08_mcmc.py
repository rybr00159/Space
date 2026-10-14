"""emcee sampling of the EB posterior (default data set).
Usage: s08_mcmc.py <geo|ms> <nsteps> <tag> <nproc> [SIGJ FR_MU FR_SD]
geo: dT0,dP,rsum,k,bn,J,L3,ecosw,esinw,q1a,q2a,q1b,q2b,s95,s105,s106
     priors: uniform in rsum (0.03-0.7), k (0.05-2.5), J (0-2), L3 (0-0.998), e cos w, e sin w (|.|<0.05, 0.3), sector scales
     (0.6-1.4); isotropic orbit orientation (uniform cos i  ->  +ln rsum when sampling bn = cos i / rsum); u1+u2 ~ N(0.63,0.10).
ms : dT0,dP,M1,q,fR,dlnJ,bn,L3,ecosw,esinw,q1a,q2a,q1b,q2b,s95,s105,s106  (lib/msfit; uniform in M1 and M2 -> +ln M1)."""
import os
os.environ['NUMBA_NUM_THREADS'] = '1'; os.environ['OMP_NUM_THREADS'] = '1'; os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np, sys, json, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4/lib')
import fit3, lclib as L, msfit as M
import emcee
from multiprocessing import get_context
R = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4/results'
KIND = sys.argv[1]; NSTEP = int(sys.argv[2]); TAG = sys.argv[3]; NPROC = int(sys.argv[4])
if KIND == 'ms' and len(sys.argv) > 7:
    M.CFG.update(SIGJ=float(sys.argv[5]), FR_MU=float(sys.argv[6]), FR_SD=float(sys.argv[7]))
d = fit3.load_data(f'{R}/data_default.pkl')
SC = ['s95', 's105', 's106']
if KIND == 'geo':
    NAMES = ['dT0', 'dP', 'rsum', 'k', 'bn', 'J', 'L3', 'ecosw', 'esinw', 'q1a', 'q2a', 'q1b', 'q2b'] + SC
else:
    NAMES = ['dT0', 'dP', 'M1', 'q', 'fR', 'dlnJ', 'bn', 'L3', 'ecosw', 'esinw', 'q1a', 'q2a', 'q1b', 'q2b'] + SC
def full(x):
    p = dict(zip(NAMES, x))
    p['T0'] = L.T00 + p['dT0'] / 1440.0; p['P'] = L.P0 + p['dP'] * 1e-6
    if KIND == 'ms': p.update(M.ms_derived(p))
    p['cosi'] = p['bn'] * p['rsum']
    for s in fit3.SECS: p.setdefault(f's{s}', 1.0)
    return p
def inb(p):
    if KIND == 'geo':
        if not (0.03 < p['rsum'] < 0.7 and 0.05 < p['k'] < 2.5 and 0 < p['J'] < 2.0): return False
    else:
        if not (0.09 < p['M1'] < 1.25 and 0.08 <= p['M1'] * p['q'] and p['q'] <= 1.0 and 0.9 < p['fR'] < 1.3 and abs(p['dlnJ']) < 0.8): return False
    if not (0 <= p['bn'] < 1 and 0 <= p['L3'] < 0.998 and abs(p['ecosw']) < 0.05 and abs(p['esinw']) < 0.3 and abs(p['dT0']) < 30 and abs(p['dP']) < 50): return False
    for a in ('q1a', 'q2a', 'q1b', 'q2b'):
        if not (0 < p[a] < 1): return False
    for s in SC:
        if not (0.6 < p[s] < 1.4): return False
    return True
def lnpost(x):
    p = full(x)
    if not inb(p): return -np.inf
    if KIND == 'geo':
        pr = fit3.prior_resid(p, 'geo'); lp = -0.5 * pr @ pr + np.log(p['rsum'])
    else:
        pr = M.prior_resid(p); lp = -0.5 * pr @ pr + np.log(p['rsum']) + np.log(p['M1'])
    return lp - 0.5 * fit3.chi2(d, p)
if __name__ == '__main__':
    rng = np.random.default_rng(11 if KIND == 'geo' else 12)
    if KIND == 'geo':
        c0 = dict(dT0=0.9, dP=1.0, rsum=0.168, k=0.42, bn=0.26, J=0.50, L3=0.966, ecosw=0.00088, esinw=0.019, q1a=0.39, q2a=0.2, q1b=0.39, q2b=0.4, s95=0.94, s105=0.93, s106=0.97)
        sd = dict(dT0=0.1, dP=0.3, rsum=0.01, k=0.08, bn=0.15, J=0.02, L3=0.02, ecosw=1e-4, esinw=0.01, q1a=0.05, q2a=0.05, q1b=0.05, q2b=0.05, s95=0.01, s105=0.01, s106=0.01)
    else:
        c0 = dict(dT0=0.9, dP=1.0, M1=0.5, q=0.4, fR=1.04, dlnJ=0.0, bn=0.3, L3=0.965, ecosw=0.00088, esinw=0.01, q1a=0.39, q2a=0.25, q1b=0.39, q2b=0.3, s95=0.94, s105=0.93, s106=0.97)
        sd = dict(dT0=0.1, dP=0.3, M1=0.08, q=0.08, fR=0.03, dlnJ=0.05, bn=0.15, L3=0.01, ecosw=1e-4, esinw=0.01, q1a=0.05, q2a=0.05, q1b=0.05, q2b=0.05, s95=0.01, s105=0.01, s106=0.01)
    nw = int(os.environ.get("NW", "48")); p0 = []
    while len(p0) < nw:
        x = np.array([c0[n] + sd[n] * rng.normal() for n in NAMES])
        if np.isfinite(lnpost(x)): p0.append(x)
    p0 = np.array(p0)
    with get_context('spawn').Pool(NPROC) as pool:
        S = emcee.EnsembleSampler(nw, len(NAMES), lnpost, pool=pool, moves=[(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)])
        t0 = time.time()
        for i, st in enumerate(S.sample(p0, iterations=NSTEP)):
            if (i + 1) % 500 == 0:
                lp = S.get_log_prob()[-1]
                try: tau = S.get_autocorr_time(tol=0, discard=(i + 1) // 3)
                except Exception: tau = np.full(len(NAMES), np.nan)
                print(f'[{time.time()-t0:6.0f}s] step {i+1} lnp max {lp.max():.2f} med {np.median(lp):.2f} acc {S.acceptance_fraction.mean():.3f} '
                      f'tau max {np.nanmax(tau):.0f} ({NAMES[int(np.nanargmax(tau))]})', flush=True)
                np.savez_compressed(f'{R}/chain_{TAG}.npz', chain=S.get_chain(), lnp=S.get_log_prob(), names=NAMES, cfg=json.dumps(M.CFG))
    np.savez_compressed(f'{R}/chain_{TAG}.npz', chain=S.get_chain(), lnp=S.get_log_prob(), names=NAMES, cfg=json.dumps(M.CFG))
    print('done')
