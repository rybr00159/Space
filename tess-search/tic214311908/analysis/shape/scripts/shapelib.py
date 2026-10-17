"""shapelib: shared model code for the shape analysis of TIC 214311908.
* LC: a light-curve data set (INV or REF windows) with per-cadence batman models (exposure integration:
  1800 s -> 31 sub-samples, 600 s -> 11, 200 s -> 5, i.e. sub-steps of 40-60 s).
* transit(): batman quadratic-LD transit of a dark body (any radius ratio k; validated for k up to 5 in s02).
* eb_flux(): two luminous spheres: eclipse of A by B (k = R_B/R_A) and of B by A, light ratio L_B/L_A, third light.
* chi2 with per-window constant offsets profiled analytically (optional).
Constants: G, Rsun, Msun, P0/T0 reference ephemeris."""
import numpy as np, batman
G = 6.674e-8; RSUN = 6.957e10; MSUN = 1.989e33; AU = 1.496e13; DAY = 86400.0; REARTH = 6.371e8
T0REF, PREF = 3840.91414, 5.6915087
SS = {1800: 31, 600: 11, 200: 5}
SD = 'tic214311908/analysis/shape'

def rho_from_aR(aR, P_d):
    """stellar density (g/cm3) from a/R* and P for a circular orbit, M_p << M*"""
    return 3*np.pi*aR**3/(G*(P_d*DAY)**2)
def aR_from_rho(rho, P_d):
    return (rho*G*(P_d*DAY)**2/(3*np.pi))**(1/3.)

class LC:
    def __init__(self, path, use_offsets=True, sig_key='sig', mask=None):
        z = np.load(path) if isinstance(path, str) else path   # path to an npz, or a dict of arrays
        files = z.files if hasattr(z, 'files') else list(z.keys())
        sel = np.ones(len(z['t']), bool) if mask is None else mask
        self.t = z['t'][sel]; self.f = z['f'][sel]; self.sig = z[sig_key][sel]; self.cad = z['cad'][sel].astype(int)
        self.sector = z['sector'][sel]; self.epoch = z['epoch'][sel]
        wid = z['win'][sel] if 'win' in files else self.epoch
        _, self.wi = np.unique(self.sector*100000 + wid*0 + self.epoch, return_inverse=True)
        self.nw = self.wi.max()+1
        self.use_offsets = use_offsets
        self.w = 1/self.sig**2
        self.sw = np.bincount(self.wi, self.w, self.nw)
        self.groups = []
        p = self._params(T0REF, PREF, 0.1, 25., 89., [0.12, 0.53])
        for c in np.unique(self.cad):
            idx = np.where(self.cad == c)[0]
            m = batman.TransitModel(p, self.t[idx], supersample_factor=SS[int(c)], exp_time=c/DAY)
            self.groups.append((idx, m, int(c)))
        self.n = len(self.t)
    @staticmethod
    def _params(t0, per, rp, a, inc, u, ecc=0., w=90.):
        p = batman.TransitParams(); p.t0 = t0; p.per = per; p.rp = rp; p.a = a; p.inc = inc; p.ecc = ecc; p.w = w
        p.limb_dark = 'quadratic'; p.u = list(u); return p
    def transit(self, t0, per, k, aR, b, u1, u2, ecc=0., w=90.):
        """fraction of the eclipsed star's light visible (dark occulter of radius ratio k). b = impact parameter at
        conjunction in units of R_eclipsed: cos i = b/aR * (1+e sin w)/(1-e^2)"""
        fac = (1 + ecc*np.sin(np.radians(w)))/(1 - ecc**2)
        cosi = b/aR*fac
        if cosi >= 1: return np.ones(self.n)
        inc = np.degrees(np.arccos(cosi))
        p = self._params(t0, per, k, aR, inc, [u1, u2], ecc, w)
        out = np.empty(self.n)
        for idx, m, c in self.groups:
            out[idx] = m.light_curve(p)
        return out
    def chi2(self, model, return_offsets=False):
        r = self.f - model
        if self.use_offsets:
            o = np.bincount(self.wi, r*self.w, self.nw)/self.sw
            r = r - o[self.wi]
        c = np.sum(r*r*self.w)
        return (c, o) if return_offsets else c
    def loglike(self, model):
        return -0.5*self.chi2(model)

def eb_flux(lc, t0, per, k, aR, b, uA, uB, LB_LA, f_src, ecc=0., w=90.):
    """flux of a blended system normalised to the total out-of-eclipse flux of the aperture (=1):
    A = eclipsed star at the primary eclipse (time t0), radius R_A; B = companion, radius k R_A, light ratio
    LB_LA = L_B/L_A in the TESS band; f_src = fraction of the total (aperture) TESS flux that comes from A+B.
    Primary: A eclipsed by B (batman with rp=k, a/R_A). Secondary: B eclipsed by A (rp = 1/k, a/R_B = aR/k),
    at the opposite conjunction (w+180 deg for eccentric orbits)."""
    fA = lc.transit(t0, per, k, aR, b, uA[0], uA[1], ecc, w)
    # secondary: B behind A. impact parameter in units of R_B: b_sec = b * (1-e sin w)/(1+e sin w) / k (approx: same inclination)
    fac1 = (1 + ecc*np.sin(np.radians(w)))/(1 - ecc**2)
    cosi = b/aR*fac1
    if LB_LA > 0:
        # time of secondary eclipse for eccentric orbits handled by batman via w+180 and a phase shift
        tsec = t0 + sec_offset(per, ecc, w)
        aRB = aR/k
        fac2 = (1 + ecc*np.sin(np.radians(w+180)))/(1 - ecc**2)
        bB = cosi*aRB/fac2
        fB = lc.transit(tsec, per, 1.0/k, aRB, bB, uB[0], uB[1], ecc, (w+180) % 360)
    else:
        fB = np.ones(lc.n)
    sys = (fA + LB_LA*fB)/(1 + LB_LA)
    return 1 + f_src*(sys - 1)

def sec_offset(per, ecc, w):
    """time from primary to secondary conjunction (Sterne 1940 exact via anomalies)."""
    if ecc == 0: return per/2
    wr = np.radians(w)
    def tp(f):  # time of true anomaly f relative to periastron
        E = 2*np.arctan(np.sqrt((1-ecc)/(1+ecc))*np.tan(f/2))
        M = E - ecc*np.sin(E)
        return M/(2*np.pi)*per
    fprim = np.pi/2 - wr; fsec = 3*np.pi/2 - wr
    return (tp(fsec) - tp(fprim)) % per

def t14_t23(per, k, aR, b, ecc=0., w=90.):
    """Winn (2010) durations (days) with the eccentric factor sqrt(1-e^2)/(1+e sin w)."""
    fac = (1 + ecc*np.sin(np.radians(w)))/(1 - ecc**2)
    inc = np.arccos(np.clip(b/aR*fac, 0, 1)); si = np.sin(inc)
    ef = np.sqrt(1-ecc**2)/(1+ecc*np.sin(np.radians(w)))
    a14 = np.clip(((1+k)**2 - b**2), 0, None); a23 = np.clip(((1-k)**2 - b**2), 0, None)
    T14 = per/np.pi*np.arcsin(np.clip(np.sqrt(a14)/(aR*si), 0, 1))*ef
    T23 = per/np.pi*np.arcsin(np.clip(np.sqrt(a23)/(aR*si), 0, 1))*ef
    T23 = np.where(b < np.abs(1-k), T23, 0.)  # flat part: inside (k<1) or total occultation (k>1)
    return T14, T23
