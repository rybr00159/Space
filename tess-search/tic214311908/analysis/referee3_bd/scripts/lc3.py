"""lc3: own light-curve likelihood machinery (referee 3).
Data: the CLEANED curves of ../clean (lc_M3V_REF_clean.csv: flux_norm, 21 strict windows; lc_M3V_INV_clean.csv:
flux_raw_norm, quality_ok & no flare, |dt| < 0.30 d of the final ephemeris), errors flux_err x beta_transit (clean's
fake-epoch calibration), per-window LINEAR baseline profiled (weighted projection).
Model: dark body, radius ratio k, impact parameter b, circular-equivalent density rho_c (a/R from rho_c), quadratic LD,
diluted by D:  m = 1 + D (lc - 1). The instantaneous curve is computed with batman on a symmetric fine dt grid
(1e-4 d) and box-averaged over each cadence's exposure with a cumulative trapezoid integral (own implementation).
Fitting: brute-force template grid in (k, b, rho_c) at nominal LD and dT0 = 0 (vectorised chi2), then Nelder-Mead
refinement of (dT0, k, b, ln rho_c, u1, u2[, ln D]) with Gaussian LD priors (sigma 0.1) and scenario bounds."""
import numpy as np, pandas as pd, batman
from scipy.optimize import minimize
from rlib3 import *
DTF = 1e-4
TFULL = np.arange(-0.34, 0.34 + DTF/2, DTF)
_pp = batman.TransitParams(); _pp.t0 = 0.; _pp.per = PC; _pp.rp = 0.1; _pp.a = 25.; _pp.inc = 89.; _pp.ecc = 0.; _pp.w = 90.
_pp.limb_dark = 'quadratic'; _pp.u = [U1, U2]
_BM = batman.TransitModel(_pp, TFULL)
CADS = (200., 600., 1800.)

def inst(k, b, aR, u1=U1, u2=U2):
    if k <= 0 or b >= aR: return np.zeros(len(TFULL))
    _pp.rp = float(k); _pp.a = float(aR); _pp.inc = float(np.degrees(np.arccos(b/aR))); _pp.u = [float(u1), float(u2)]
    return _BM.light_curve(_pp) - 1.0

def binned(k, b, aR, u1=U1, u2=U2):
    lc = inst(k, b, aR, u1, u2)
    ci = np.concatenate([[0.], np.cumsum(0.5*(lc[1:] + lc[:-1])*DTF)])
    out = {}
    for c in CADS:
        hw = c/DAY/2
        x = TFULL[(TFULL >= -0.32) & (TFULL <= 0.32)]
        out[c] = (x, (np.interp(x + hw, TFULL, ci) - np.interp(x - hw, TFULL, ci))/(2*hw))
    return out

class LC:
    def __init__(self, t, f, e, sec, cad, win, tc):
        self.t, self.f, self.e = np.asarray(t, float), np.asarray(f, float), np.asarray(e, float)
        self.sec, self.cad, self.tc = np.asarray(sec, int), np.asarray(cad, float), np.asarray(tc, float)
        _, self.wi = np.unique(np.asarray(win), return_inverse=True); self.nw = self.wi.max() + 1
        self.dt = self.t - self.tc
        self.w = 1/self.e**2
        self.y = self.f - 1
        self.X = np.vstack([np.ones_like(self.dt), self.dt]).T
        A = np.zeros((self.nw, 2, 2))
        for i in range(2):
            for j in range(2):
                A[:, i, j] = np.bincount(self.wi, self.w*self.X[:, i]*self.X[:, j], self.nw)
        self.Ainv = np.linalg.inv(A)
        self.ck = np.array([min(CADS, key=lambda c: abs(c - x)) for x in self.cad])
        self.cidx = {c: np.where(self.ck == c)[0] for c in CADS if np.any(self.ck == c)}
        self.yp = self.perp(self.y)
    def perp(self, a):
        """weighted residual of a after per-window linear fit (a: (N,) or (N, m))"""
        a2 = a if a.ndim == 2 else a[:, None]
        p0 = np.vstack([np.bincount(self.wi, self.w*self.X[:, 0]*a2[:, j], self.nw) for j in range(a2.shape[1])]).T   # nw x m
        p1 = np.vstack([np.bincount(self.wi, self.w*self.X[:, 1]*a2[:, j], self.nw) for j in range(a2.shape[1])]).T
        c0 = self.Ainv[:, 0, 0][:, None]*p0 + self.Ainv[:, 0, 1][:, None]*p1
        c1 = self.Ainv[:, 1, 0][:, None]*p0 + self.Ainv[:, 1, 1][:, None]*p1
        r = a2 - (c0[self.wi]*self.X[:, :1] + c1[self.wi]*self.X[:, 1:2])
        return r if a.ndim == 2 else r[:, 0]
    def template(self, curves, dT0=0.):
        s = np.zeros(len(self.t))
        for c, idx in self.cidx.items():
            x, v = curves[c]
            s[idx] = np.interp(self.dt[idx] - dT0, x, v, left=0., right=0.)
        return s
    def chi2_s(self, s):
        r = self.perp(self.y - s)
        return float(np.sum(self.w*r*r))
    def subset(self, m):
        n = LC.__new__(LC)
        return LC(self.t[m], self.f[m], self.e[m], self.sec[m], self.cad[m], self.wi[m], self.tc[m])

def load_ref(beta=True, fpath=None):
    r = pd.read_csv(fpath or f'{B}/clean/lc_M3V_REF_clean.csv')
    r = r[np.isfinite(r.flux_norm) & np.isfinite(r.flux_err)]
    tc = T0C + r.epoch.values*PC
    e = r.flux_err.values*(r.beta_transit.values if beta else 1.0)
    return LC(r.time_btjd_tdb.values, r.flux_norm.values, e, r.sector.values, r.cadence_s.values, r.sector.values*1000 + r.epoch.values, tc)

def load_inv(beta=True, half=0.30, T0=T0C, P=PC, avoid=None, df=None, min_in=2, min_side=10):
    d = df if df is not None else inv_frame()
    ep = np.round((d.time_btjd_tdb.values - T0)/P).astype(int)
    tc = T0 + ep*P; dt = d.time_btjd_tdb.values - tc
    m = np.abs(dt) < half
    if avoid is not None: m &= avoid
    win = d.sector.values*100000 + ep
    keep = np.zeros(len(d), bool)
    for wv in np.unique(win[m]):
        mm = m & (win == wv)
        if np.sum(mm & (np.abs(dt) < 0.03)) >= min_in and np.sum(mm & (dt < -0.07)) >= min_side and np.sum(mm & (dt > 0.07)) >= min_side:
            keep |= mm
    e = d.flux_err.values*(d.beta_transit.values if beta else 1.0)
    return LC(d.time_btjd_tdb.values[keep], d.flux_raw_norm.values[keep], e[keep], d.sector.values[keep], d.cadence_s.values[keep], win[keep], tc[keep])

_INV = {}
def inv_frame():
    if 'd' not in _INV:
        d = pd.read_csv(f'{B}/clean/lc_M3V_INV_clean.csv', usecols=['time_btjd_tdb', 'sector', 'cadence_s', 'quality_ok', 'flux_raw_norm', 'flux_err', 'beta_transit', 'flare_flag'])
        d = d[(d.quality_ok == 1) & (d.flare_flag == 0) & np.isfinite(d.flux_raw_norm) & np.isfinite(d.flux_err)].reset_index(drop=True)
        _INV['d'] = d
    return _INV['d']

# ------------------------------------------------------------------ fitting
SCEN = dict(planet=dict(k=(0.02, 0.35), D=None), grazBD=dict(k=(0.19, 0.35), D=None), grazBD20=dict(k=(0.20, 0.35), D=None),
            comp=dict(k=(0.25, 1.0), D=(0.01, 0.6)))
class TemplateBank:
    """precomputed binned templates on a (k, b, rho_c) grid at nominal LD"""
    def __init__(self, ks, nb=34, rhos=None):
        self.rhos = rhos if rhos is not None else np.exp(np.linspace(np.log(0.6), np.log(45.), 30))
        self.items = []
        for k in ks:
            for bb in np.linspace(0, 0.995, nb):
                b = bb*(1 + k)
                for rho in self.rhos:
                    self.items.append((k, b, rho, binned(k, b, float(aR_of_rho(rho)))))
    def scan(self, lc, D=None, kr=None, chunk=400):
        """return list of (chi2, k, b, rho, D) for all templates (D profiled linearly when D == 'free')"""
        res = []
        for i0 in range(0, len(self.items), chunk):
            its = self.items[i0:i0 + chunk]
            if kr is not None: its = [it for it in its if kr[0] <= it[0] <= kr[1]]
            if not its: continue
            S = np.column_stack([lc.template(it[3]) for it in its])
            Sp = lc.perp(S)
            if D == 'free':
                num = np.sum(lc.w[:, None]*Sp*lc.yp[:, None], axis=0); den = np.sum(lc.w[:, None]*Sp*Sp, axis=0)
                Dh = np.clip(num/np.maximum(den, 1e-30), 0.0, 1.0)
                c2 = np.sum(lc.w[:, None]*(lc.yp[:, None] - Dh[None, :]*Sp)**2, axis=0)
            else:
                Dh = np.ones(len(its))
                c2 = np.sum(lc.w[:, None]*(lc.yp[:, None] - Sp)**2, axis=0)
            res += [(float(c2[j]), its[j][0], its[j][1], its[j][2], float(Dh[j])) for j in range(len(its))]
        return res

def chi2_par(lc, x, scen):
    lim = SCEN[scen]
    dT0, lk, bb, lr, u1, u2 = x[:6]
    k = np.exp(lk); rho = np.exp(lr)
    D = np.exp(x[6]) if lim['D'] is not None else 1.0
    if not (lim['k'][0] <= k <= lim['k'][1]) or not (0 <= bb < 1) or not (0.3 <= rho <= 200) or abs(dT0) > 0.012: return 1e12
    if lim['D'] is not None and not (lim['D'][0] <= D <= lim['D'][1]): return 1e12
    if not (u1 > 0 and u1 + u2 < 1 and u1 + 2*u2 > 0): return 1e12
    s = D*lc.template(binned(k, bb*(1 + k), float(aR_of_rho(rho)), u1, u2), dT0)
    return lc.chi2_s(s) + ((u1 - U1)/SU)**2 + ((u2 - U2)/SU)**2

def refine(lc, scen, starts, maxfev=1500):
    best = None
    for s0 in starts:
        f = lambda x: chi2_par(lc, x, scen)
        r = minimize(f, np.array(s0, float), method='Nelder-Mead', options=dict(maxfev=maxfev, xatol=1e-5, fatol=1e-3, adaptive=True))
        r = minimize(f, r.x, method='Nelder-Mead', options=dict(maxfev=maxfev//2, xatol=1e-6, fatol=1e-4, adaptive=True))
        if best is None or r.fun < best.fun: best = r
    x = best.x
    k = float(np.exp(x[1])); b = float(x[2]*(1 + k)); rho = float(np.exp(x[3]))
    D = float(np.exp(x[6])) if len(x) > 6 else 1.0
    aR = float(aR_of_rho(rho)); inc = np.arccos(min(b/aR, 1))
    x14 = np.sqrt(max((1 + k)**2 - b**2, 0))/(aR*np.sin(inc)); x23 = np.sqrt(max((1 - k)**2 - b**2, 0))/(aR*np.sin(inc))
    T14 = PC/np.pi*np.arcsin(min(x14, 1)); T23 = PC/np.pi*np.arcsin(min(x23, 1)) if b < 1 - k else 0.
    return dict(chi2=float(best.fun), x=[float(v) for v in x], k=k, b=b, rho_c=rho, D=D, dT0_min=float(x[0])*1440, u1=float(x[4]), u2=float(x[5]),
                T14_h=T14*24, flat=T23/T14 if T14 > 0 else 0.)

def fit_scen(lc, bank, scen, nstart=3):
    lim = SCEN[scen]
    sc = bank.scan(lc, D='free' if lim['D'] is not None else None, kr=lim['k'])
    sc.sort(key=lambda r: r[0])
    starts, seen = [], []
    for c2, k, b, rho, D in sc:
        key = (round(np.log(rho), 1), round(b, 2))
        if any(abs(key[0] - s[0]) < 0.25 and abs(key[1] - s[1]) < 0.08 for s in seen): continue
        seen.append(key)
        bb = min(b/(1 + k), 0.99)
        st = [0., np.log(k), bb, np.log(rho), U1, U2]
        if lim['D'] is not None: st.append(np.log(max(D, lim['D'][0]*1.01)))
        starts.append(st)
        if len(starts) >= nstart: break
    r = refine(lc, scen, starts)
    r['grid_best_chi2'] = sc[0][0]
    return r
