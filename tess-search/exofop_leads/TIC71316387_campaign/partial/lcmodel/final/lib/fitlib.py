"""Data preparation, likelihood and priors for the physical EB fit.
Per-event multiplicative quadratic baseline is profiled out analytically (weighted linear least squares
of y on [m, m*dt, m*dt^2]); per-sector white-noise sigma is estimated from out-of-eclipse residuals and
inflated by the time-averaging beta factor (red noise on 30-90 min scales)."""
import numpy as np, sys
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/final/lib')
import lclib as L
from ebmodel import eb_flux, _xg, _wg
from numba import njit

NAMES = ['T0', 'P', 'rsum', 'k', 'cosi', 'J', 'L3', 'ecosw', 'esinw', 'q1a', 'q2a', 'q1b', 'q2b', 's95', 's105', 's106']

def beta_factor(r, cad, scales=(30, 45, 60, 90)):
    """Pont et al. time-averaging beta: ratio of binned rms to white expectation, median over scales (min)."""
    b = []
    s1 = np.std(r)
    for sc in scales:
        n = int(round(sc / (cad * 1440)))
        if n < 2: continue
        m = r.size // n
        if m < 10: continue
        rb = r[:m * n].reshape(m, n).mean(1)
        b.append(np.std(rb) / (s1 / np.sqrt(n) * np.sqrt(m / (m - 1))))
    return max(1.0, float(np.median(b))) if b else 1.0

def prepare(lcs, hw=0.20, P=L.P0, T0=L.T00, inflate=True, verbose=True):
    """Concatenate eclipse windows. Returns dict of arrays: t, y, sig, ev (event id), dt, sector, texp, nsub, and event table."""
    T, Y, S, EV, DT, SEC, TEXP, NSUB, PH = [], [], [], [], [], [], [], [], []
    evtab = []
    eid = 0
    info = {}
    for lc in lcs:
        lc = L.clip_outliers(lc, P=P, T0=T0)
        ws = L.windows(lc, P=P, T0=T0, hw=hw)
        # noise from OOE parts of windows after quadratic baseline
        res = []; resc = []
        for w in ws:
            i = w['idx']; dt = lc['t'][i] - w['tc']; y = lc['f'][i]; oot = np.abs(dt) > 0.75 * L.T14
            c = np.polyfit(dt[oot], y[oot], 2); r = y[oot] / np.polyval(c, dt[oot]) - 1
            res.append(r)
        r_all = np.concatenate(res)
        sig = 1.4826 * np.median(np.abs(r_all - np.median(r_all)))
        # beta from the continuous OOE stretches of each window (concatenate per side)
        bet = []
        for r in res:
            if r.size > 20: bet.append(r)
        if not bet: bet = res
        beta = beta_factor(np.concatenate(bet), lc['texp']) if inflate else 1.0
        info[lc['sector']] = dict(prod=lc['prod'], sigma=sig, beta=beta, nev=len(ws), npts=int(sum(len(w['idx']) for w in ws)))
        if verbose:
            print(f"  S{lc['sector']} {lc['prod']}: {len(ws)} eclipse windows, sigma {sig*1e6:.0f} ppm/pt, beta {beta:.2f}")
        for w in ws:
            i = w['idx']
            T.append(lc['t'][i]); Y.append(lc['f'][i]); S.append(np.full(i.size, sig * beta)); EV.append(np.full(i.size, eid))
            DT.append(lc['t'][i] - w['tc']); SEC.append(np.full(i.size, lc['sector'])); TEXP.append(np.full(i.size, lc['texp']))
            NSUB.append(np.full(i.size, lc['nsub'])); PH.append(np.full(i.size, w['ph']))
            evtab.append(dict(eid=eid, n=w['n'], ph=w['ph'], tc=w['tc'], sector=lc['sector'], prod=lc['prod']))
            eid += 1
    d = dict(t=np.concatenate(T), y=np.concatenate(Y), sig=np.concatenate(S), ev=np.concatenate(EV), dt=np.concatenate(DT),
             sector=np.concatenate(SEC), texp=np.concatenate(TEXP), nsub=np.concatenate(NSUB), ph=np.concatenate(PH))
    d['evtab'] = evtab; d['info'] = info
    # event boundaries (contiguous by construction)
    ev = d['ev']; starts = np.r_[0, np.where(np.diff(ev) != 0)[0] + 1]; ends = np.r_[starts[1:], ev.size]
    d['starts'] = starts; d['ends'] = ends
    # groups of identical (texp,nsub) for model evaluation
    d['groups'] = []
    for sec in np.unique(d['sector']):
        m = np.where(d['sector'] == sec)[0]
        d['groups'].append((int(sec), m, float(d['texp'][m[0]]), int(d['nsub'][m[0]])))
    return d

@njit(cache=True)
def profile_chi2(y, m, dt, w, starts, ends, order):
    """Sum over events of min_c sum w (y - m*poly_c(dt))^2 ; also returns per-event baseline coefficients."""
    chi = 0.0
    nE = starts.size
    coefs = np.zeros((nE, 3))
    for e in range(nE):
        a = starts[e]; b = ends[e]
        nb = order + 1
        A = np.zeros((nb, nb)); v = np.zeros(nb)
        for i in range(a, b):
            basis = np.empty(nb)
            x = 1.0
            for j in range(nb):
                basis[j] = m[i] * x
                x *= dt[i]
            for j in range(nb):
                v[j] += w[i] * basis[j] * y[i]
                for l in range(nb):
                    A[j, l] += w[i] * basis[j] * basis[l]
        c = np.linalg.solve(A, v)
        for j in range(nb): coefs[e, j] = c[j]
        for i in range(a, b):
            x = 1.0; mod = 0.0
            for j in range(nb):
                mod += c[j] * m[i] * x
                x *= dt[i]
            chi += w[i] * (y[i] - mod) ** 2
    return chi, coefs

def q2u(q1, q2):
    s = np.sqrt(q1); return 2 * s * q2, s * (1 - 2 * q2)

def unpack(th):
    return dict(zip(NAMES, th))

def model_all(d, p, scales=None):
    u1a, u2a = q2u(p['q1a'], p['q2a']); u1b, u2b = q2u(p['q1b'], p['q2b'])
    m = np.empty(d['t'].size)
    for sec, idx, texp, nsub in d['groups']:
        f = eb_flux(d['t'][idx], p['T0'], p['P'], p['rsum'], p['k'], p['cosi'], p['J'], p['L3'], p['ecosw'], p['esinw'],
                    u1a, u2a, u1b, u2b, texp, nsub, _xg, _wg)
        s = p.get(f's{sec}', 1.0)
        m[idx] = 1 - s * (1 - f)
    return m

LD_SUM_MU, LD_SUM_SD = 0.63, 0.10   # Claret (2017) TESS quadratic, 3200-5800 K dwarfs: u1+u2 = 0.57-0.69

def log_prior(p, ldprior='claret'):
    if not (0.03 < p['rsum'] < 0.7 and 0.1 < p['k'] < 1.6 and 0 <= p['cosi'] < p['rsum'] and 0 < p['J'] < 1.6
            and 0 <= p['L3'] < 0.995 and abs(p['ecosw']) < 0.05 and abs(p['esinw']) < 0.3
            and 0 < p['q1a'] < 1 and 0 < p['q2a'] < 1 and 0 < p['q1b'] < 1 and 0 < p['q2b'] < 1
            and 0.6 < p['s95'] < 1.4 and 0.6 < p['s105'] < 1.4 and 0.6 < p['s106'] < 1.4):
        return -np.inf
    if p['ecosw'] ** 2 + p['esinw'] ** 2 > 0.09: return -np.inf
    lp = 0.0
    if ldprior == 'claret':
        for a in 'ab':
            u1, u2 = q2u(p['q1' + a], p['q2' + a])
            lp += -0.5 * ((u1 + u2 - LD_SUM_MU) / LD_SUM_SD) ** 2
    return lp

def log_like(d, p, order=2):
    m = model_all(d, p)
    w = 1.0 / d['sig'] ** 2
    chi, _ = profile_chi2(d['y'], m, d['dt'], w, d['starts'], d['ends'], order)
    return -0.5 * chi

class Posterior:
    def __init__(self, d, fixed=None, ldprior='claret', free=None):
        self.d = d; self.fixed = fixed or {}; self.ldprior = ldprior
        self.free = free or [n for n in NAMES if n not in self.fixed]
    def full(self, x):
        p = dict(self.fixed); p.update(dict(zip(self.free, x))); return p
    def __call__(self, x):
        p = self.full(x)
        lp = log_prior(p, self.ldprior)
        if not np.isfinite(lp): return -np.inf
        return lp + log_like(self.d, p)
