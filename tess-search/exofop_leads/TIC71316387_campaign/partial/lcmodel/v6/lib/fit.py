"""Likelihood, priors and optimisation helpers for the v5 physical EB fit.
Per-event multiplicative quadratic baseline is profiled analytically (weighted linear least squares of y on
[m, m dt, m dt^2]), so no prior detrending touches the eclipses.  Per-product white noise = OOE robust rms x
max(1, beta).  Limb darkening: u1+u2 ~ N(0.64, 0.06) and 0.05 < u1 < 0.60 for each star, which spans the
Claret (2017) TESS quadratic coefficients of 3200-4800 K dwarfs (u1+u2 = 0.60-0.68, u1 = 0.13-0.51)."""
import numpy as np, sys
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6/lib')
import ebm
from numba import njit

BASE = ['T0', 'P', 'rsum', 'k', 'b0', 'J', 'L3', 'ecosw', 'esinw', 'usa', 'u1a', 'usb', 'u1b']
LD_MU, LD_SD = 0.64, 0.06
T0_REF = 3165.6697   # reference epoch for the T0 sanity bound (+-0.05 d); reset by Model() for injection data sets

@njit(cache=True)
def profile_chi2(y, m, dt, w, starts, ends, order):
    chi = 0.0
    nE = starts.size
    coefs = np.zeros((nE, 3))
    nb = order + 1
    for e in range(nE):
        a = starts[e]; b = ends[e]
        A = np.zeros((nb, nb)); v = np.zeros(nb)
        for i in range(a, b):
            x = 1.0
            bas = np.empty(nb)
            for j in range(nb):
                bas[j] = m[i] * x; x *= dt[i]
            for j in range(nb):
                v[j] += w[i] * bas[j] * y[i]
                for l in range(nb):
                    A[j, l] += w[i] * bas[j] * bas[l]
        c = np.linalg.solve(A, v)
        for j in range(nb): coefs[e, j] = c[j]
        for i in range(a, b):
            x = 1.0; mod = 0.0
            for j in range(nb):
                mod += c[j] * m[i] * x; x *= dt[i]
            chi += w[i] * (y[i] - mod) ** 2
    return chi, coefs

class Model:
    def __init__(self, d, scale_keys=(), order=2, P0=1.3036712, T00=3165.6697):
        self.d = d; self.order = order
        global T0_REF
        if 'truth' in d: T0_REF = float(d['truth']['T0'])
        self.gm = ebm.GridModel(P0)
        self.scale_keys = list(scale_keys)
        self.names = BASE + [f's_{k}' for k in self.scale_keys]
        self.keyidx = {k: np.where(d['key'] == k)[0] for k in self.scale_keys}
        self.w = 1.0 / d['sig'] ** 2
        self.T00 = T00

    def eb(self, p, t=None, texp=None):
        t = self.d['t'] if t is None else t; texp = self.d['texp'] if texp is None else texp
        P = p['P']
        tph = ((t - p['T0'] + 0.25 * P) % P) - 0.25 * P
        u1a = p['u1a']; u2a = p['usa'] - u1a; u1b = p['u1b']; u2b = p['usb'] - u1b
        return self.gm.eval(tph, texp, P, p['rsum'], p['k'], p['b0'], p['J'], p['ecosw'], p['esinw'], u1a, u2a, u1b, u2b)

    def model(self, p):
        f = self.eb(p)
        dep = (1.0 - p['L3']) * (1.0 - f)
        for k, idx in self.keyidx.items():
            dep[idx] *= p[f's_{k}']
        return 1.0 - dep

    def chi2(self, p, return_coefs=False):
        m = self.model(p)
        chi, c = profile_chi2(self.d['y'], m, self.d['dt'], self.w, self.d['starts'], self.d['ends'], self.order)
        return (chi, c, m) if return_coefs else chi

def in_bounds(p, kmax=3.0):
    if not (abs(p['P'] - 1.3036712) < 2e-4 and abs(p['T0'] - T0_REF) < 0.05):
        return False
    if not (0.03 < p['rsum'] < 0.8 and 0.08 < p['k'] < kmax and 0.0 <= p['b0'] < 1.0 + 0.0 and 0.0 < p['J'] < 3.0
            and 0.0 <= p['L3'] < 0.999 and abs(p['ecosw']) < 0.05 and abs(p['esinw']) < 0.4 and p['b0'] * p['rsum'] < 1.0):
        return False
    for a in 'ab':
        us, u1 = p['us' + a], p['u1' + a]
        if not (0.05 < u1 < 0.60 and 0.3 < us < 1.0): return False
    if p['ecosw'] ** 2 + p['esinw'] ** 2 > 0.16: return False
    for k, v in p.items():
        if k.startswith('s_') and not (0.6 < v < 1.5): return False
    return True

def log_prior(p, scale_sd=None, kmax=3.0):
    if not in_bounds(p, kmax): return -np.inf
    lp = -0.5 * ((p['usa'] - LD_MU) / LD_SD) ** 2 - 0.5 * ((p['usb'] - LD_MU) / LD_SD) ** 2
    if scale_sd is not None:
        for k, v in p.items():
            if k.startswith('s_'): lp += -0.5 * (np.log(v) / scale_sd) ** 2
    return lp

class Posterior:
    """log posterior over the free parameters (others fixed)."""
    def __init__(self, model, fixed, free, scale_sd=None, kmax=3.0, extra_prior=None):
        self.m = model; self.fixed = dict(fixed); self.free = list(free); self.scale_sd = scale_sd; self.kmax = kmax
        self.extra = extra_prior
    def full(self, x):
        p = dict(self.fixed); p.update(zip(self.free, x)); return p
    def __call__(self, x):
        p = self.full(x)
        lp = log_prior(p, self.scale_sd, self.kmax)
        if not np.isfinite(lp): return -np.inf
        if self.extra is not None:
            lp += self.extra(p)
            if not np.isfinite(lp): return -np.inf
        return lp - 0.5 * self.m.chi2(p)

BOUNDS = dict(rsum=(0.0301, 0.799), k=(0.0801, 2.999), b0=(0.0, 0.9999), J=(1e-4, 2.999), L3=(0.0, 0.9989), ecosw=(-0.0499, 0.0499),
              esinw=(-0.399, 0.399), usa=(0.3001, 0.9999), u1a=(0.0501, 0.5999), usb=(0.3001, 0.9999), u1b=(0.0501, 0.5999),
              M1=(0.1201, 0.8999), q=(0.1501, 1.0), dR1=(-0.299, 0.299), dR2=(-0.299, 0.299), dJ=(-0.399, 0.399), dmag=(-1.0, 1.0))
STEPS = dict(T0=2e-4, P=5e-7, rsum=0.01, k=0.05, b0=0.03, J=0.03, L3=0.01, ecosw=0.0003, esinw=0.01, usa=0.03, u1a=0.05, usb=0.03,
             u1b=0.05, M1=0.03, q=0.05, dR1=0.02, dR2=0.02, dJ=0.03, dmag=0.05)
def bounds_for(names):
    return [BOUNDS.get(n, (0.6001, 1.4999) if n.startswith('s_') else (None, None)) for n in names]
def steps_for(names):
    return [STEPS.get(n, 0.01) for n in names]
