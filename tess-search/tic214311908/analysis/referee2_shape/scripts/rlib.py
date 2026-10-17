"""rlib: referee-2 shape library (independent of ../shape/scripts).
Data:
  INV = ../phot/lc_M3V_all.csv, quality_ok==1 & flare_flag==0 & finite; fractional flux rescaled to 979 e-/s:
        f = 1 + (flux_norm-1)*935/979; windows |dt| < 0.30 d around the phase-1 ephemeris.
  REF = ../referee_phot/lc_M3V_referee_windows.csv (all finite rows; already 1 + dF/979).
Noise: sigma_i = err_i * s_sector * beta_sector, s = robust std(OOT)/median(err). beta options:
  'white' (1), 'fake' (per-sector depth-error calibration from fake epochs; INV: own calibration on the full
  detrended series, file out/beta_inv.json; REF: the referee's fake-window ratios ../referee_phot/out/per_sector.csv).
Model: batman quadratic-LD transit with exposure integration (sub-steps <= 60 s), evaluated only where
|t - t_mid| < 0.15 d (exactly 1 elsewhere for the planet / box models used here). Per-window constant offsets
are profiled analytically (optionally offset + slope)."""
import numpy as np, pandas as pd, json, batman
B = 'tic214311908/analysis'
RD = f'{B}/referee2_shape'
T0REF, PREF = 3840.91414, 5.6915087
G = 6.674e-8; RSUN = 6.957e10; MSUN = 1.989e33; DAY = 86400.0; AU = 1.496e13; C = 2.998e10
SUB = {1800: 31, 600: 11, 200: 5}

def rstd(x):
    x = x[np.isfinite(x)]; return 1.4826*np.median(np.abs(x - np.median(x)))

def aR_from_rho(rho, P):
    return (rho*G*(P*DAY)**2/(3*np.pi))**(1/3.)
def rho_from_aR(aR, P):
    return 3*np.pi*aR**3/(G*(P*DAY)**2)

def load(ds, beta='fake', win=0.30):
    if ds == 'inv':
        d = pd.read_csv(f'{B}/phot/lc_M3V_all.csv')
        d = d[(d.quality_ok == 1) & (d.flare_flag == 0) & np.isfinite(d.flux_norm) & np.isfinite(d.flux_err)].copy()
        sc = 935/979.
        t = d.time_btjd_tdb.values; f = 1 + (d.flux_norm.values - 1)*sc; e = d.flux_err.values*sc
        sec = d.sector.values; cad = d.cadence_s.values.astype(int)
        ep = np.round((t - T0REF)/PREF).astype(int); dt = t - (T0REF + ep*PREF)
        m = np.abs(dt) < win
        oot = (np.abs(dt) > 0.10)
        wid = sec*1000 + ep
    else:
        d = pd.read_csv(f'{B}/referee_phot/lc_M3V_referee_windows.csv')
        d = d[np.isfinite(d.flux_norm) & np.isfinite(d.flux_err)].copy()
        t = d.time_btjd_tdb.values; f = d.flux_norm.values; e = d.flux_err.values
        sec = d.sector.values; cad = d.cadence_s.values.astype(int); ep = d.epoch.values.astype(int)
        dt = t - (T0REF + ep*PREF); m = np.ones(len(t), bool); oot = d.oot.values == 1
        wid = sec*1000 + ep
    sig = np.empty(len(t)); info = {}
    if beta == 'fake':
        if ds == 'inv':
            bj = json.load(open(f'{RD}/out/beta_inv.json'))
            bt = {int(k): v['beta'] for k, v in bj.items() if not k.startswith('_')}
        else:
            ps = pd.read_csv(f'{B}/referee_phot/out/per_sector.csv').set_index('sector')
            bt = {int(s): max(1., 0.5*(ps.loc[s, 'err_fake'] + ps.loc[s, 'err_fake_robust'])/ps.loc[s, 'err_white']) for s in ps.index}
    for s in np.unique(sec):
        k = sec == s
        o = k & oot & m
        sc_ = rstd(f[o] - 1)/np.median(e[o])
        b_ = 1.0 if beta == 'white' else (bt[int(s)] if isinstance(beta, str) else beta[int(s)])
        sig[k] = e[k]*sc_*b_
        info[int(s)] = dict(scale=float(sc_), beta=float(b_))
    out = dict(t=t[m], f=f[m], sig=sig[m], sector=sec[m], cad=cad[m], epoch=ep[m], dt=dt[m], wid=wid[m], oot=oot[m])
    out['info'] = info
    return out

class LC:
    def __init__(self, d, slope=False, near=0.15):
        self.t = d['t']; self.f = d['f']; self.sig = d['sig']; self.cad = d['cad']; self.sector = d['sector']; self.epoch = d['epoch']
        _, self.wi = np.unique(d['wid'], return_inverse=True); self.nw = self.wi.max() + 1
        self.w = 1/self.sig**2; self.n = len(self.t); self.slope = slope
        self.near = np.abs(d['dt']) < near
        # exposure integration with the MIDPOINT rule (offsets (j+0.5)/n - 0.5 of the exposure). NB batman's own
        # supersample_factor uses linspace(-0.5, 0.5, n) offsets, i.e. an effective exposure n/(n-1) x too long
        # (1.25 x for 200-s data with n = 5); it is therefore not used here.
        idx = np.where(self.near)[0]
        ts, own = [], []
        for c in np.unique(self.cad[idx]):
            ii = idx[self.cad[idx] == c]; ns = SUB[int(c)]
            off = ((np.arange(ns) + 0.5)/ns - 0.5)*c/DAY
            ts.append((self.t[ii][:, None] + off[None, :]).ravel()); own.append(np.repeat(ii, ns))
        self.ts = np.concatenate(ts); self.own = np.concatenate(own)
        self.cnt = np.bincount(self.own, minlength=self.n if hasattr(self, 'n') else len(self.t)).astype(float)
        p = self.params(T0REF, PREF, 0.1, 25., 89., [0.2, 0.4])
        self.bm = batman.TransitModel(p, self.ts)
        # window design for offsets (+ slopes)
        self.sw = np.bincount(self.wi, self.w, self.nw)
        if slope:
            tm = np.bincount(self.wi, self.t*self.w, self.nw)/self.sw
            self.x = self.t - tm[self.wi]
            self.sxx = np.bincount(self.wi, self.w*self.x**2, self.nw)
    @staticmethod
    def params(t0, per, k, aR, inc, u, e=0., w=90.):
        p = batman.TransitParams(); p.t0 = t0; p.per = per; p.rp = k; p.a = aR; p.inc = inc; p.ecc = e; p.w = w
        p.limb_dark = 'quadratic'; p.u = list(u); return p
    def transit(self, t0, per, k, aR, b, u1, u2, e=0., w=90.):
        fac = (1 + e*np.sin(np.radians(w)))/(1 - e**2)
        cosi = b/aR*fac
        out = np.ones(self.n)
        if cosi >= 1: return out
        p = self.params(t0, per, k, aR, np.degrees(np.arccos(cosi)), [u1, u2], e, w)
        v = self.bm.light_curve(p)
        s = np.bincount(self.own, v, minlength=self.n); m = self.cnt > 0
        out[m] = s[m]/self.cnt[m]
        return out
    def chi2(self, model, ret=False):
        r = self.f - model
        o = np.bincount(self.wi, r*self.w, self.nw)/self.sw
        r = r - o[self.wi]
        if self.slope:
            s = np.bincount(self.wi, r*self.w*self.x, self.nw)/self.sxx
            r = r - s[self.wi]*self.x
        c = float(np.sum(r*r*self.w))
        return (c, r) if ret else c

def r_nauenberg(M):
    mu = M/1.454
    return 0.0112*np.sqrt(mu**(-2/3.) - mu**(2/3.))

def kepler_a_rsun(Mtot, P):
    return (G*Mtot*MSUN*(P*DAY)**2/(4*np.pi**2))**(1/3.)/RSUN

def t14(P, k, aR, b):
    inc = np.arccos(np.clip(b/aR, 0, 1))
    x = np.clip(np.sqrt(np.clip((1 + k)**2 - b**2, 0, None))/(aR*np.sin(inc)), 0, 1)
    y = np.clip(np.sqrt(np.clip((1 - k)**2 - b**2, 0, None))/(aR*np.sin(inc)), 0, 1)
    T14 = P/np.pi*np.arcsin(x); T23 = np.where(b < 1 - k, P/np.pi*np.arcsin(y), 0.)
    return T14, T23
