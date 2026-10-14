"""Fit library v3: same physics as fit2 (two-sphere model, quadratic LD, third light, per-event profiled quadratic
baseline), with (i) supersampling NSUB = 9 (1800 s), 3 (600 s), 2 (200 s) [model error < 12 ppm, test_nsub.py],
(ii) impact parameter in units of the sum of radii, bn = cos i / rsum (eclipse requires bn < 1),
(iii) residual-vector interface for scipy least_squares, (iv) isotropic-orientation prior (uniform in cos i)
implemented as +ln(rsum) when sampling in bn."""
import numpy as np, sys, pickle
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4/lib')
import lclib as L, fitlib as F, fit2
from ebmodel import eb_flux
from numba import njit
import ms_table as MS
XG, WG = fit2.XG, fit2.WG
NSUB = {1: 9, 28: 5, 68: 2, 95: 2, 105: 2, 106: 2}
SECS = [1, 28, 68, 95, 105, 106]
LD_SUM_MU, LD_SUM_SD = 0.63, 0.10
q2u = fit2.q2u

def load_data(path=None, hw=0.20, lcs=None):
    if path is not None:
        d = pickle.load(open(path, 'rb'))
    else:
        d = F.prepare(lcs if lcs is not None else L.default_set(), hw=hw, verbose=False)
    d['groups'] = [(sec, idx, texp, NSUB[sec]) for (sec, idx, texp, nsub) in d['groups']]
    return d

def model_all(d, p):
    u1a, u2a = q2u(p['q1a'], p['q2a']); u1b, u2b = q2u(p['q1b'], p['q2b'])
    m = np.empty(d['t'].size)
    cosi = p['bn'] * p['rsum'] if 'bn' in p else p['cosi']
    for sec, idx, texp, nsub in d['groups']:
        f = eb_flux(d['t'][idx], p['T0'], p['P'], p['rsum'], p['k'], cosi, p['J'], p['L3'], p['ecosw'], p['esinw'],
                    u1a, u2a, u1b, u2b, texp, nsub, XG, WG)
        s = p.get(f's{sec}', 1.0)
        m[idx] = 1 - s * (1 - f)
    return m

@njit(cache=True)
def profile_resid(y, m, dt, w, starts, ends):
    """weighted residuals after per-event LSQ fit of multiplicative quadratic baseline"""
    r = np.empty(y.size)
    for e in range(starts.size):
        a = starts[e]; b = ends[e]
        A = np.zeros((3, 3)); v = np.zeros(3)
        for i in range(a, b):
            b0 = m[i]; b1 = m[i] * dt[i]; b2 = m[i] * dt[i] * dt[i]
            bas = (b0, b1, b2)
            for j in range(3):
                v[j] += w[i] * w[i] * bas[j] * y[i]
                for l in range(3):
                    A[j, l] += w[i] * w[i] * bas[j] * bas[l]
        c = np.linalg.solve(A, v)
        for i in range(a, b):
            mod = m[i] * (c[0] + c[1] * dt[i] + c[2] * dt[i] * dt[i])
            r[i] = w[i] * (y[i] - mod)
    return r

def resid(d, p):
    m = model_all(d, p)
    return profile_resid(d['y'], m, d['dt'], 1.0 / d['sig'], d['starts'], d['ends'])

def chi2(d, p):
    r = resid(d, p); return float(r @ r)

def ms_derived(p):
    M1 = p['M1']; M2 = p['M1'] * p['q']
    R1 = MS.radius(M1) * p['fR']; R2 = MS.radius(M2) * p['fR']
    a = fit2.a_over_R(M1 + M2, p['P'])
    J = MS.SB_T(M2) / MS.SB_T(M1) * np.exp(p['dlnJ'])
    return dict(rsum=(R1 + R2) / a, k=R2 / R1, J=J, R1=R1, R2=R2, a=a, M2=M2)

def prior_resid(p, kind):
    """Gaussian priors expressed as residuals (for least_squares): LD sums; for 'ms' the fR and dlnJ priors"""
    out = []
    for a in 'ab':
        u1, u2 = q2u(p['q1' + a], p['q2' + a]); out.append((u1 + u2 - LD_SUM_MU) / LD_SUM_SD)
    if kind == 'ms':
        out.append(p['dlnJ'] / 0.1); out.append((p['fR'] - 1.04) / 0.05)
    return np.array(out)

class Model:
    def __init__(self, d, free, fixed, kind='geo', isotropic=True):
        self.d = d; self.free = list(free); self.fixed = dict(fixed); self.kind = kind; self.isotropic = isotropic
    def full(self, x):
        p = dict(self.fixed); p.update(dict(zip(self.free, x)))
        if self.kind == 'ms': p.update(ms_derived(p))
        p['cosi'] = p['bn'] * p['rsum']
        return p
    def inbounds(self, p):
        if self.kind == 'ms':
            if not (0.09 < p['M1'] < 1.5 and 0.09 / p['M1'] <= p['q'] <= 1.0 and 0.85 < p['fR'] < 1.3 and abs(p['dlnJ']) < 0.5):
                return False
        else:
            if not (0.03 < p['rsum'] < 0.7 and 0.1 < p['k'] < 1.6 and 0 < p['J'] < 1.6): return False
        if not (0 <= p['bn'] < 1 and 0 <= p['L3'] < 0.999 and abs(p['ecosw']) < 0.05 and abs(p['esinw']) < 0.3): return False
        if p['ecosw'] ** 2 + p['esinw'] ** 2 > 0.09: return False
        for a in ('q1a', 'q2a', 'q1b', 'q2b'):
            if not (0 < p[a] < 1): return False
        for s in SECS:
            if not (0.6 < p.get(f's{s}', 1.0) < 1.4): return False
        return True
    def lnpost(self, x):
        p = self.full(x)
        if not self.inbounds(p): return -np.inf
        pr = prior_resid(p, self.kind)
        lp = -0.5 * pr @ pr
        if self.isotropic: lp += np.log(p['rsum'])   # uniform in cos i  (cos i = bn * rsum)
        return lp - 0.5 * chi2(self.d, p)
    __call__ = lnpost
    def lsq_resid(self, x):
        p = self.full(x)
        return np.concatenate([resid(self.d, p), prior_resid(p, self.kind)])
