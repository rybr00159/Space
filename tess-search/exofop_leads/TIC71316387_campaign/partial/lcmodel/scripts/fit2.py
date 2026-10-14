"""Physical EB fit library (v2).  Two-sphere model (ebmodel.eb_flux, 16-node Gauss-Legendre, error <2e-8),
per-event multiplicative quadratic baseline profiled analytically, per-sector white noise from OOE residuals.
Geometric model ('geo'): rsum, k, cosi, J, L3 free (isotropic prior on cos i, uniform on others).
Main-sequence model ('ms'): M1, q, fR (radius inflation), dlnJ (surface-brightness scatter), cosi, L3; rsum, k, J derived
from Pecaut & Mamajek relations + Kepler (distance independent)."""
import numpy as np, sys
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/scripts')
import lclib as L, fitlib as F
from ebmodel import eb_flux
import ms_table as MS
XG, WG = np.polynomial.legendre.leggauss(16)
SECS = [1, 28, 68, 95, 105, 106]
GM_SUN = 1.32712440018e20; RSUN = 6.957e8; DAY = 86400.0

def q2u(q1, q2):
    s = np.sqrt(q1); return 2 * s * q2, s * (1 - 2 * q2)

def a_over_R(Mtot, P):
    return (GM_SUN * Mtot * (P * DAY) ** 2 / (4 * np.pi ** 2)) ** (1 / 3) / RSUN

def model_all(d, p):
    u1a, u2a = q2u(p['q1a'], p['q2a']); u1b, u2b = q2u(p['q1b'], p['q2b'])
    m = np.empty(d['t'].size)
    for sec, idx, texp, nsub in d['groups']:
        f = eb_flux(d['t'][idx], p['T0'], p['P'], p['rsum'], p['k'], p['cosi'], p['J'], p['L3'], p['ecosw'], p['esinw'],
                    u1a, u2a, u1b, u2b, texp, nsub, XG, WG)
        s = p.get(f's{sec}', 1.0)
        m[idx] = 1 - s * (1 - f)
    return m

def chi2(d, p, order=2):
    m = model_all(d, p)
    c, coefs = F.profile_chi2(d['y'], m, d['dt'], 1.0 / d['sig'] ** 2, d['starts'], d['ends'], order)
    return c

LD_SUM_MU, LD_SUM_SD = 0.63, 0.10

def ms_derived(p):
    """map MS parameters -> geometric ones"""
    M1 = p['M1']; M2 = p['M1'] * p['q']
    R1 = MS.radius(M1) * p['fR']; R2 = MS.radius(M2) * p['fR']
    a = a_over_R(M1 + M2, p['P'])
    J = MS.SB_T(M2) / MS.SB_T(M1) * np.exp(p['dlnJ'])
    return dict(rsum=(R1 + R2) / a, k=R2 / R1, J=J, R1=R1, R2=R2, a=a, M2=M2)

class Post:
    """kind='geo' or 'ms'.  free: list of free names; fixed: dict."""
    def __init__(self, d, free, fixed, kind='geo', ldprior=True, sprior=None):
        self.d = d; self.free = list(free); self.fixed = dict(fixed); self.kind = kind; self.ldprior = ldprior
        self.sprior = sprior   # None (uniform 0.6-1.4) or sigma of Gaussian prior around 1
    def full(self, x):
        p = dict(self.fixed); p.update(dict(zip(self.free, x)))
        if self.kind == 'ms':
            p.update(ms_derived(p))
        return p
    def lnprior(self, p):
        if self.kind == 'ms':
            if not (0.09 < p['M1'] < 1.5 and 0.09 / p['M1'] <= p['q'] <= 1.0 and 0.9 < p['fR'] < 1.3 and abs(p['dlnJ']) < 0.6):
                return -np.inf
            lp = -0.5 * (p['dlnJ'] / 0.1) ** 2 - 0.5 * ((p['fR'] - 1.04) / 0.05) ** 2
        else:
            if not (0.03 < p['rsum'] < 0.7 and 0.1 < p['k'] < 1.6 and 0 < p['J'] < 1.6):
                return -np.inf
            lp = 0.0
        if not (0 <= p['cosi'] < p['rsum'] and 0 <= p['L3'] < 0.999 and abs(p['ecosw']) < 0.05 and abs(p['esinw']) < 0.3
                and 0 < p['q1a'] < 1 and 0 < p['q2a'] < 1 and 0 < p['q1b'] < 1 and 0 < p['q2b'] < 1):
            return -np.inf
        if p['ecosw'] ** 2 + p['esinw'] ** 2 > 0.09: return -np.inf
        for s in SECS:
            v = p.get(f's{s}', 1.0)
            if not (0.6 < v < 1.4): return -np.inf
            if self.sprior and f's{s}' in self.free: lp += -0.5 * ((v - 1) / self.sprior) ** 2
        if self.ldprior:
            for a in 'ab':
                u1, u2 = q2u(p['q1' + a], p['q2' + a])
                lp += -0.5 * ((u1 + u2 - LD_SUM_MU) / LD_SUM_SD) ** 2
        return lp
    def __call__(self, x):
        p = self.full(x)
        lp = self.lnprior(p)
        if not np.isfinite(lp): return -np.inf
        return lp - 0.5 * chi2(self.d, p)

def default_data(hw=0.20, verbose=True, lcs=None):
    return F.prepare(lcs if lcs is not None else L.default_set(), hw=hw, verbose=verbose)
