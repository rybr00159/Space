"""Core transit model and likelihood for TIC 87055304 (shared by all fit scripts in this folder).

Model: batman quadratic-LD transit (circular orbit), computed once per call on a fine grid in time-from-mid-transit
(tau, 10 s spacing) and integrated over each sector's exposure time exactly (box average via the cumulative
integral of the grid model), then interpolated at the data. This is equivalent to batman supersampling with
10 s sub-exposures for every cadence (1800, 600, 200 s), but ~20x faster.

Parameters theta (17):
  T0 [BTJD, reference epoch E_REF], P [d], k = Rp/R*, b, rho* [rho_sun], u1, u2,
  c_s (per-sector multiplicative baseline offset, 5), ln sigma_s (per-sector white-noise level, 5)
Optional dilution D (contaminating / target flux, fixed): observed = (F + D) / (1 + D).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import json
import numpy as np
import batman

HERE = os.path.dirname(os.path.abspath(__file__))
G = 6.67430e-11; RSUN = 6.957e8; MSUN = 1.98892e30
RHO_SUN = MSUN / (4 / 3 * np.pi * RSUN ** 3)
SECTORS = [13, 39, 66, 93, 104]
EXPT = {13: 1800.0, 39: 600.0, 66: 200.0, 93: 200.0, 104: 200.0}
P_INIT, T0_INIT = 4.136195, 1658.41325
WIN = 0.45          # data half-window around each transit [d]
GRID_DT = 10.0 / 86400
LD = json.load(open(os.path.join(HERE, "ld_prior.json")))
U1_MU, U1_SD = LD["u1_mean"], LD["u1_prior_sigma"]
U2_MU, U2_SD = LD["u2_mean"], LD["u2_prior_sigma"]
NAMES = ["T0", "P", "k", "b", "rho", "u1", "u2"] + [f"c{s}" for s in SECTORS] + [f"lns{s}" for s in SECTORS]


def a_from_rho(rho, P):
    """a/R* from stellar density (solar units) and period (d), circular orbit, Mp << M*."""
    return (rho * RHO_SUN * G * (P * 86400.0) ** 2 / (3 * np.pi)) ** (1.0 / 3.0)


def rho_from_a(a, P):
    return 3 * np.pi * a ** 3 / (G * (P * 86400.0) ** 2) / RHO_SUN


def load(path=None, sectors=SECTORS, win=WIN, P=P_INIT, T0=T0_INIT, key="f"):
    """Data within +-win of each predicted transit. Returns dict sector -> (t, f, epoch)."""
    path = path or os.path.join(HERE, "lc_prepared.npy")
    lc = np.load(path, allow_pickle=True).item()
    out = {}
    for s in sectors:
        v = lc[s]
        n = np.round((v["t"] - T0) / P)
        tau = v["t"] - T0 - n * P
        m = np.abs(tau) < win
        out[s] = dict(t=v["t"][m], f=v[key][m], n=n[m].astype(int), sigma0=v["sigma"], exptime=EXPT[s],
                      f_raw=v["f_raw"][m], trend=v["trend"][m])
    return out


class GridModel:
    """Exposure-integrated transit profile as a function of tau (time from mid-transit)."""

    def __init__(self, half=WIN + 0.03, dt=GRID_DT):
        self.g = np.arange(-half, half + dt / 2, dt)
        self.p = batman.TransitParams()
        p = self.p
        p.t0, p.per, p.rp, p.a, p.inc, p.ecc, p.w, p.u, p.limb_dark = 0.0, 4.0, 0.06, 5.0, 85.0, 0.0, 90.0, [0.3, 0.2], "quadratic"
        self.m = batman.TransitModel(p, self.g)
        self.dt = dt

    def profile(self, P, k, a, b, u1, u2):
        p = self.p
        p.per, p.rp, p.a = P, k, a
        cosi = b / a
        if cosi >= 1:
            return np.ones_like(self.g)
        p.inc = np.degrees(np.arccos(cosi)); p.u = [u1, u2]
        return self.m.light_curve(p)

    def integrate(self, f, exptime, tau):
        """Box average of profile f over exposure exptime [s], evaluated at tau."""
        if exptime <= 0:
            return np.interp(tau, self.g, f)
        C = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * self.dt)])
        h = exptime / 86400 / 2
        return (np.interp(tau + h, self.g, C) - np.interp(tau - h, self.g, C)) / (2 * h)


GM = None


def gm():
    global GM
    if GM is None:
        GM = GridModel()
    return GM


def phys_model(theta7, data, dil=0.0):
    """Return dict sector -> transit model (no baseline) for theta7 = (T0, P, k, b, rho, u1, u2)."""
    T0, P, k, b, rho, u1, u2 = theta7
    a = a_from_rho(rho, P)
    g = gm()
    prof = g.profile(P, k, a, b, u1, u2)
    # exposure-integrated profiles on the grid; exposure half-widths are integer multiples of the 10 s grid step
    C = np.concatenate([[0.0], np.cumsum(0.5 * (prof[1:] + prof[:-1]) * g.dt)])  # C[i] = integral up to g[i]
    cache = {}
    out = {}
    for s, d in data.items():
        E = d["exptime"]
        if E not in cache:
            m = int(round(E / 86400 / 2 / g.dt))
            Fe = np.ones_like(prof)
            Fe[m:-m] = (C[2 * m:] - C[:-2 * m]) / (2 * m * g.dt)
            cache[E] = Fe
        tau = d["t"] - T0 - np.round((d["t"] - T0) / P) * P
        mod = np.interp(tau, g.g, cache[E])
        if dil:
            mod = (mod + dil) / (1 + dil)
        out[s] = mod
    return out


def ld_ok(u1, u2):
    return (u1 > 0) and (u1 + u2 < 1) and (u1 + 2 * u2 > 0)


class Posterior:
    """log posterior; rho_prior = ('uniform', lo, hi). b_range = (lo, hi) in units of b; allow grazing up to 1+k."""

    def __init__(self, data, rho_prior=("uniform", 0.003, 5.0), dil=0.0, b_mode="free", ld_prior=True,
                 sectors=None, t0_ref=None, p_ref=None):
        self.data = data
        self.sectors = sectors or list(data.keys())
        self.rho_prior = rho_prior
        self.dil = dil
        self.b_mode = b_mode  # 'free' (0 <= b < 1+k), 'grazing' (1-k < b < 1+k), 'nongrazing' (b < 1-k)
        self.ld_prior = ld_prior
        self.npar = 7 + 2 * len(self.sectors)

    def names(self):
        return ["T0", "P", "k", "b", "rho", "u1", "u2"] + [f"c{s}" for s in self.sectors] + [f"lns{s}" for s in self.sectors]

    def log_prior(self, th):
        T0, P, k, b, rho, u1, u2 = th[:7]
        if not (0.005 < k < 0.5 and 0 <= b < 1 + k and 3.9 < P < 4.4):
            return -np.inf
        if self.b_mode == "grazing" and not (b > 1 - k):
            return -np.inf
        if self.b_mode == "nongrazing" and not (b < 1 - k):
            return -np.inf
        kind, lo, hi = self.rho_prior
        if kind == "uniform":
            if not (lo < rho < hi):
                return -np.inf
            lp = 0.0
        elif kind == "loguniform":
            if not (lo < rho < hi):
                return -np.inf
            lp = -np.log(rho)
        elif kind == "gauss":
            if rho <= 0:
                return -np.inf
            lp = -0.5 * ((rho - lo) / hi) ** 2
        if not ld_ok(u1, u2):
            return -np.inf
        if self.ld_prior:
            lp += -0.5 * ((u1 - U1_MU) / U1_SD) ** 2 - 0.5 * ((u2 - U2_MU) / U2_SD) ** 2
        ns = len(self.sectors)
        c = th[7:7 + ns]; lns = th[7 + ns:7 + 2 * ns]
        if np.any(np.abs(c) > 0.05) or np.any(lns < -10) or np.any(lns > -3):
            return -np.inf
        # a/R* must exceed 1+k (star not inside orbit)
        if a_from_rho(rho, P) < 1.5:
            return -np.inf
        return lp

    def log_like(self, th, return_chi2=False):
        ns = len(self.sectors)
        mods = phys_model(th[:7], {s: self.data[s] for s in self.sectors}, self.dil)
        ll = 0.0; chi2 = 0.0
        for i, s in enumerate(self.sectors):
            d = self.data[s]
            mod = mods[s] * (1 + th[7 + i])
            sig2 = np.exp(2 * th[7 + ns + i])
            r2 = (d["f"] - mod) ** 2
            ll += -0.5 * np.sum(r2 / sig2 + np.log(2 * np.pi * sig2))
            chi2 += np.sum(r2) / sig2
        if return_chi2:
            return ll, chi2
        return ll

    def __call__(self, th):
        lp = self.log_prior(th)
        if not np.isfinite(lp):
            return -np.inf
        ll = self.log_like(th)
        return lp + ll if np.isfinite(ll) else -np.inf


def derived(chain, P_col=1):
    """Derived quantities from a flat chain with columns T0,P,k,b,rho,u1,u2,..."""
    T0, P, k, b, rho, u1, u2 = [chain[:, i] for i in range(7)]
    a = a_from_rho(rho, P)
    inc = np.degrees(np.arccos(np.clip(b / a, -1, 1)))
    sini = np.sin(np.radians(inc))
    arg14 = np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) / (a * sini)
    arg23 = np.sqrt(np.clip((1 - k) ** 2 - b ** 2, 0, None)) / (a * sini)
    T14 = P / np.pi * np.arcsin(np.clip(arg14, 0, 1)) * 24
    T23 = P / np.pi * np.arcsin(np.clip(arg23, 0, 1)) * 24
    return dict(a=a, inc=inc, T14_h=T14, T23_h=T23)


def depth_of(th7, dil=0.0):
    """Model depth at mid-transit (unbinned, no exposure smearing), in ppm."""
    T0, P, k, b, rho, u1, u2 = th7
    g = gm()
    prof = g.profile(P, k, a_from_rho(rho, P), b, u1, u2)
    d = 1 - prof[np.argmin(np.abs(g.g))]
    return d / (1 + dil) * 1e6


class PosteriorOE(Posterior):
    """Odd/even variant: theta = base 17 parameters (k = k_odd) + k_even appended at the end.
    Parity is that of the absolute epoch number counted from T0_INIT (epoch 0 = BTJD 1658.413)."""

    def __init__(self, data, e_ref, **kw):
        super().__init__(data, **kw)
        self.e_ref = e_ref
        self.npar += 1
        for s, d in self.data.items():
            d["parity"] = None

    def names(self):
        return super().names() + ["k_even"]

    def log_prior(self, th):
        lp = super().log_prior(th[:-1])
        if not np.isfinite(lp) or not (0.005 < th[-1] < 0.5):
            return -np.inf
        if self.b_mode == "free" and not (th[3] < 1 + th[-1]):
            return -np.inf
        return lp

    def log_like(self, th, return_chi2=False):
        ns = len(self.sectors)
        base = th[:-1]
        T0, P = base[0], base[1]
        th_e = np.array(base[:7]); th_e[2] = th[-1]
        mo = phys_model(base[:7], {s: self.data[s] for s in self.sectors}, self.dil)
        me = phys_model(th_e, {s: self.data[s] for s in self.sectors}, self.dil)
        ll = 0.0; chi2 = 0.0
        for i, s in enumerate(self.sectors):
            d = self.data[s]
            ep = np.round((d["t"] - T0) / P).astype(int) + self.e_ref
            odd = (ep % 2) == 1
            mod = np.where(odd, mo[s], me[s]) * (1 + base[7 + i])
            sig2 = np.exp(2 * base[7 + ns + i])
            r2 = (d["f"] - mod) ** 2
            ll += -0.5 * np.sum(r2 / sig2 + np.log(2 * np.pi * sig2))
            chi2 += np.sum(r2) / sig2
        if return_chi2:
            return ll, chi2
        return ll
