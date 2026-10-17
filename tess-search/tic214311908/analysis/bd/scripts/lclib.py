"""lclib: fast transit-shape machinery for bd/.
Model: an opaque body of radius ratio k crossing a quadratic-limb-darkened host (batman, circular-equivalent orbit with
a/R from rho_circ), diluted by D (fraction of the deblended M3V-system flux supplied by the host):
    m(t) = 1 + D [lc_k,b,aR(t - tc - dT0) - 1].
The instantaneous curve is computed once per parameter set on a fine |dt| grid (2e-4 d), box-averaged per cadence
(200/600/1800 s, exposure = cadence) with cumulative sums, and interpolated at |dt - dT0| (the circular-equivalent
transit is symmetric). Per-window polynomial baselines (order 1 for REF, 2 for INV) are profiled analytically.
Noise: sigma = err x per-sector robust OOT scale x beta_sector."""
import numpy as np, batman
from blib import *
DTF = 2e-4
TGF = np.arange(0, 0.36 + DTF/2, DTF)          # |dt| grid (d)
_pp = batman.TransitParams(); _pp.t0 = 0.; _pp.per = P; _pp.rp = 0.1; _pp.a = 25.; _pp.inc = 89.; _pp.ecc = 0.; _pp.w = 90.
_pp.limb_dark = 'quadratic'; _pp.u = [U1, U2]
# symmetric full grid for batman (needed for exact box averaging near dt = 0)
_TFULL = np.arange(-0.36 - 0.011, 0.36 + 0.011 + DTF/2, DTF)
_BM = batman.TransitModel(_pp, _TFULL)
CADS = (200., 600., 1800.)

def inst_curve(k, b, aR, u1=U1, u2=U2):
    cosi = b/aR
    if cosi >= 1 or k <= 0: return np.ones(len(_TFULL))
    _pp.rp = k; _pp.a = aR; _pp.inc = float(np.degrees(np.arccos(cosi))); _pp.u = [u1, u2]
    return _BM.light_curve(_pp)

def cad_curves(k, b, aR, u1=U1, u2=U2):
    """dict cadence -> exposure-averaged (curve - 1) on TGF (|dt| >= 0): exact box average of the trapezoid-integrated
    instantaneous curve, CI(x + c/2) - CI(x - c/2) over c."""
    lc = inst_curve(k, b, aR, u1, u2) - 1
    ci = np.concatenate([[0.], np.cumsum(0.5*(lc[1:] + lc[:-1])*DTF)])
    out = {}
    for c in CADS:
        hw = c/DAY/2
        out[c] = (np.interp(TGF + hw, _TFULL, ci) - np.interp(TGF - hw, _TFULL, ci))/(2*hw)
    return out

class Data:
    def __init__(self, t, f, e, sec, cad, win, tc, oot, npoly=2, beta=None, scale=None):
        t = np.asarray(t, float); win = np.asarray(win)
        # drop windows with too few points (singular baselines)
        uw, cnt = np.unique(win, return_counts=True); keep = np.isin(win, uw[cnt >= 12])
        t, f, e, sec, cad, win, tc, oot = [np.asarray(v)[keep] for v in (t, f, e, sec, cad, win, tc, oot)]
        self.t = t; self.f = np.asarray(f, float); self.sec = np.asarray(sec).astype(int)
        self.cad = np.asarray(cad, float); self.tc = np.asarray(tc, float); self.oot = np.asarray(oot, bool)
        self.dt = self.t - self.tc
        _, self.wi = np.unique(np.asarray(win), return_inverse=True); self.nw = self.wi.max() + 1
        e = np.asarray(e, float)
        if scale is None:
            scale = {}
            for s in np.unique(self.sec):
                m = (self.sec == s) & self.oot
                scale[int(s)] = rstd(self.f[m] - 1)/np.median(e[m])
        self.scale = scale
        self.beta = beta if beta is not None else {int(s): 1.0 for s in np.unique(self.sec)}
        self.sig = e*np.array([scale[int(s)]*self.beta.get(int(s), 1.0) for s in self.sec])
        self.w = 1/self.sig**2
        self.npoly = npoly
        # per-window baseline design (centred dt)
        self.X = np.vstack([self.dt**p for p in range(npoly)]).T
        A = np.zeros((self.nw, npoly, npoly))
        for a_ in range(npoly):
            for b_ in range(npoly):
                A[:, a_, b_] = np.bincount(self.wi, self.w*self.X[:, a_]*self.X[:, b_], self.nw)
        self.Ainv = np.linalg.pinv(A)
        self.y = self.f - 1
        self.cadkey = np.array([min(CADS, key=lambda c: abs(c - x)) for x in self.cad])
        self.cidx = {c: np.where(self.cadkey == c)[0] for c in CADS if np.any(self.cadkey == c)}
    def proj(self, a):
        return np.vstack([np.bincount(self.wi, self.w*a*self.X[:, p], self.nw) for p in range(self.npoly)]).T
    def Q(self, a, b_):
        pa = self.proj(a); pb = pa if b_ is a else self.proj(b_)
        return float(np.sum(self.w*a*b_) - np.einsum('wi,wij,wj->', pa, self.Ainv, pb))
    def template(self, curves, dT0=0.):
        s = np.zeros(len(self.t))
        x = np.abs(self.dt - dT0)
        for c, idx in self.cidx.items():
            s[idx] = np.interp(x[idx], TGF, curves[c], right=0.)
        return s
    def chi2_D(self, s, D):
        r = self.y - D*s
        return self.Q(r, r)
    def best_D(self, s):
        qys = self.Q(self.y, s); qss = self.Q(s, s)
        return qys/qss if qss > 0 else 0., qys, qss

def load_npz(path, sel=None, **kw):
    z = np.load(path)
    m = np.ones(len(z['t']), bool) if sel is None else sel(z)
    return Data(z['t'][m], z['f'][m], z['e'][m], z['sec'][m], z['cad'][m], z['win'][m], z['tc'][m], z['oot'][m], **kw)
