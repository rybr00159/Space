"""Sensitivity of the transit parameters to the detrending / data-cleaning choices (maximum-likelihood fits).

Variants (all with the same likelihood as the MCMC: per-sector offset and white noise, Claret LD prior, free density):
  base        biweight 1.0 d, mask 0.30 d, +-1 h around momentum dumps removed (= the MCMC data set)
  win0.5/0.75/1.5/2.0   biweight window
  mask0.25/0.40         transit-mask full width
  dump0 / dump2h        no dump padding / +-2 h
  context     /home/user/work/tic87/data/lc_sap_clean.npy (the earlier data set behind the CONTEXT.md numbers)
  localpoly   raw SAP flux (no detrending), per-transit quadratic baseline profiled analytically at every
              likelihood call (segment-restricted windows of +-0.45 d)
  gp          raw SAP flux, celerite2 SHO (Q = 1/sqrt 2) GP + transit per sector (ML hyper-parameters),
              windows of +-0.45 d
For each: k, b, rho*, depth (mid-transit), T14, chi2 at the ML. Writes detrend_sys.json.
"""
import os, sys, json
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.optimize import minimize
import tmodel as tm
import prep

HERE = os.path.dirname(os.path.abspath(__file__))
NAMES7 = ["T0", "P", "k", "b", "rho", "u1", "u2"]


def ml(post, x0, restarts=2, seed=0):
    rng = np.random.default_rng(seed)
    f = lambda x: -post(x) if np.isfinite(post(x)) else 1e25
    best = None
    for r in range(restarts):
        x = np.array(x0, float)
        if r:
            x[3] = np.clip(x[3] + 0.08 * rng.standard_normal(), 0, 0.98); x[4] *= np.exp(0.3 * rng.standard_normal())
        if not np.isfinite(post(x)):
            x = np.array(x0, float)
        for it in range(4):
            res = minimize(f, x, method="Nelder-Mead", options=dict(maxiter=15000, maxfev=15000, xatol=1e-9, fatol=1e-5, adaptive=True))
            x = res.x
        if best is None or res.fun < best.fun:
            best = res
    return best


def summarize(x, extra=None):
    th7 = x[:7]
    der = tm.derived(np.array([th7]))
    d = dict(zip(NAMES7, map(float, th7)))
    d.update(depth_ppm=float(tm.depth_of(th7)), T14_h=float(der["T14_h"][0]), T23_h=float(der["T23_h"][0]), a=float(der["a"][0]))
    if extra:
        d.update(extra)
    return d


class LocalPolyPosterior(tm.Posterior):
    """Raw flux; per-transit quadratic baseline (multiplicative) solved by linear least squares at each call."""

    def __init__(self, data, **kw):
        super().__init__(data, **kw)
        lc = np.load(os.path.join(HERE, "lc_prepared.npy"), allow_pickle=True).item()
        self.win = {}
        for s in self.sectors:
            v = lc[s]
            n = np.round((v["t"] - tm.T0_INIT) / tm.P_INIT).astype(int)
            tau = v["t"] - tm.T0_INIT - n * tm.P_INIT
            m = np.abs(tau) < tm.WIN
            t, f, seg, nn = v["t"][m], v["f_raw"][m], v["seg"][m], n[m]
            groups = []
            for e in np.unique(nn):
                k = nn == e
                # restrict to the segment with most points near mid-transit
                kk = k & (np.abs(tau[m]) < 0.1)
                if kk.sum() == 0:
                    vals, cnt = np.unique(seg[k], return_counts=True)
                else:
                    vals, cnt = np.unique(seg[kk], return_counts=True)
                k = k & (seg == vals[np.argmax(cnt)])
                if k.sum() > 10:
                    groups.append(np.where(k)[0])
            self.data[s] = dict(t=t, f=f, exptime=tm.EXPT[s], groups=groups)

    def log_like(self, th, return_chi2=False):
        ns = len(self.sectors)
        mods = tm.phys_model(th[:7], {s: self.data[s] for s in self.sectors}, self.dil)
        ll = 0.0; chi2 = 0.0
        for i, s in enumerate(self.sectors):
            d = self.data[s]; sig2 = np.exp(2 * th[7 + ns + i]); n = 0; r2 = 0.0
            for g in d["groups"]:
                t = d["t"][g]; tc = t.mean(); m = mods[s][g]
                X = np.vstack([m, m * (t - tc), m * (t - tc) ** 2]).T
                c, *_ = np.linalg.lstsq(X, d["f"][g], rcond=None)
                r2 += np.sum((d["f"][g] - X @ c) ** 2); n += len(g)
            ll += -0.5 * (r2 / sig2 + n * np.log(2 * np.pi * sig2))
            chi2 += r2 / sig2
        return (ll, chi2) if return_chi2 else ll


def gp_fit(x0):
    """celerite2 SHO GP on raw flux in transit windows (one GP per sector, shared hyper-parameters per sector)."""
    import celerite2
    from celerite2 import terms
    lc = np.load(os.path.join(HERE, "lc_prepared.npy"), allow_pickle=True).item()
    data = {}
    for s in tm.SECTORS:
        v = lc[s]
        n = np.round((v["t"] - tm.T0_INIT) / tm.P_INIT).astype(int)
        m = np.abs(v["t"] - tm.T0_INIT - n * tm.P_INIT) < tm.WIN
        # normalise each segment piece by its median to remove dump jumps
        f = v["f_raw"][m].copy(); seg = v["seg"][m]
        for k in np.unique(seg):
            q = seg == k; f[q] /= np.median(f[q])
        data[s] = dict(t=v["t"][m], f=f, exptime=tm.EXPT[s], seg=seg)
    ns = len(tm.SECTORS)

    def unpack(x):
        th7 = x[:7]; hp = x[7:].reshape(ns, 3)  # per sector: ln sigma_gp, ln rho_gp, ln jitter
        return th7, hp

    gps = {}
    def nll(x):
        th7, hp = unpack(x)
        T0, P, k, b, rho, u1, u2 = th7
        if not (0.01 < k < 0.3 and 0 <= b < 1 + k and 0.003 < rho < 5 and tm.ld_ok(u1, u2) and tm.a_from_rho(rho, P) > 1.5):
            return 1e25
        if np.any(hp[:, 1] < np.log(0.1)) or np.any(hp[:, 1] > np.log(10)) or np.any(np.abs(hp[:, 0] + 7) > 4) or np.any(np.abs(hp[:, 2] + 6) > 4):
            return 1e25
        mods = tm.phys_model(th7, data)
        out = ((u1 - tm.U1_MU) / tm.U1_SD) ** 2 / 2 + ((u2 - tm.U2_MU) / tm.U2_SD) ** 2 / 2
        for i, s in enumerate(tm.SECTORS):
            d = data[s]
            kern = terms.SHOTerm(sigma=np.exp(hp[i, 0]), rho=np.exp(hp[i, 1]), Q=1 / np.sqrt(2))
            gp = celerite2.GaussianProcess(kern, mean=0.0)
            gp.compute(d["t"], diag=np.exp(2 * hp[i, 2]) * np.ones(len(d["t"])), quiet=True)
            r = d["f"] - mods[s]
            # per-segment constant offsets absorbed approximately by the GP (rho up to 10 d); remove mean first
            out -= gp.log_likelihood(r - np.mean(r))
        return out if np.isfinite(out) else 1e25

    hp0 = []
    for s in tm.SECTORS:
        hp0 += [np.log(1e-3), np.log(1.0), np.log(data[s]["f"].std())]
    x = np.r_[x0[:7], hp0]
    best = None
    for it in range(5):
        res = minimize(nll, x, method="Nelder-Mead", options=dict(maxiter=30000, maxfev=30000, xatol=1e-9, fatol=1e-5, adaptive=True))
        x = res.x
        if best is None or res.fun < best.fun:
            best = res
    th7, hp = unpack(best.x)
    return best.x, dict(gp_hyper={str(s): dict(sigma_ppm=float(np.exp(h[0]) * 1e6), rho_d=float(np.exp(h[1])), jitter_ppm=float(np.exp(h[2]) * 1e6)) for s, h in zip(tm.SECTORS, hp)},
                        nll=float(best.fun))


if __name__ == "__main__":
    st = json.load(open(os.path.join(HERE, "ml_start.json")))
    x0 = np.array(st["x"])
    variants = dict(base=dict(), **{f"win{w}": dict(win=w) for w in (0.5, 0.75, 1.5, 2.0)},
                    **{f"mask{m}": dict(maskw=m) for m in (0.25, 0.40)}, dump0=dict(dump_pad=0.0), dump2h=dict(dump_pad=2 / 24))
    which = sys.argv[1:] or (list(variants) + ["context", "localpoly", "gp"])
    out_path = os.path.join(HERE, "detrend_sys.json")
    out = json.load(open(out_path)) if os.path.exists(out_path) else {}
    for v in which:
        if v in variants:
            lc = prep.build(**variants[v])
            path = os.path.join(HERE, f"lc_variant_{v}.npy"); np.save(path, lc, allow_pickle=True)
            data = tm.load(path)
            post = tm.Posterior(data)
            r = ml(post, x0)
            out[v] = summarize(r.x, dict(lnpost=float(-r.fun), sigma_ppm={str(s): float(np.exp(r.x[12 + i]) * 1e6) for i, s in enumerate(tm.SECTORS)}))
            os.remove(path)
        elif v == "context":
            src = np.load("/home/user/work/tic87/data/lc_sap_clean.npy", allow_pickle=True).item()
            for s in src:
                src[s]["sigma"] = src[s]["sigma"]
            path = os.path.join(HERE, "lc_variant_context.npy"); np.save(path, src, allow_pickle=True)
            data = tm.load(path)
            post = tm.Posterior(data)
            r = ml(post, x0)
            out[v] = summarize(r.x, dict(lnpost=float(-r.fun)))
            os.remove(path)
        elif v == "localpoly":
            data = tm.load()
            post = LocalPolyPosterior(data)
            r = ml(post, x0)
            out[v] = summarize(r.x, dict(lnpost=float(-r.fun)))
        elif v == "gp":
            x, extra = gp_fit(x0)
            out[v] = summarize(x, extra)
        print(v, {k: (round(val, 5) if isinstance(val, float) else val) for k, val in out[v].items() if k in ("k", "b", "rho", "depth_ppm", "T14_h", "T23_h", "P")}, flush=True)
        json.dump(out, open(out_path, "w"), indent=1)
