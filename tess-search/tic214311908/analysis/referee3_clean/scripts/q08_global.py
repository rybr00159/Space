"""q08: independent global transit fits of the analyst's delivered cleaned light curves (own windows, own noise
calibration, own likelihood and sampler set-up).

Data: INV = lc_M3V_INV_clean.csv flux_raw_norm (quality_ok, flare_flag == 0), +-0.30 d windows at the analyst's
final ephemeris, acceptance as analyst; REF = lc_M3V_REF_clean.csv windows as delivered (optionally without E58).
Noise: per-window OOT robust std x beta_s; beta_s = std(fake-window depths) / median(formal error) with the fixed
shape, from fake windows (INV: phases 0.12-0.88 of the delivered LC; REF: clean/data/win_ref_fake.npz).
Model: batman quadratic LD, circular, exposure-integrated (5/9/21 sub-steps); per-window ADDITIVE linear baseline
profiled (analyst: multiplicative; difference negligible at 1.5%); parameters T0 (E=0), P, k, b, log10 rho, u1, u2.
Priors: rho ~ N(9.3, 0.9) g/cm3 ('single'); u ~ N((0.125, 0.533), 0.1) (as analyst) and a WIDE variant
(q1, q2 uniform in [0,1], Kipping 2013) to show the LD-prior dependence of the depth.
Sampler: emcee 32 walkers x 3000 steps, burn 1000.
Output: out/q08_global.json
"""
import json, sys, numpy as np, pandas as pd, batman, emcee
from scipy.optimize import minimize
sys.path.insert(0, "tic214311908/analysis/referee3_clean/scripts")
from q02_bigap import shape, accept
R3 = "tic214311908/analysis/referee3_clean"
CL = "tic214311908/analysis/clean"
T0A, PA = 3840.91488, 5.6915020
G = 6.674e-8


def rstd(x): return 1.4826 * np.median(np.abs(x - np.median(x)))


def linfit_depth(t, f, tc, ex):
    s = shape(t, tc, ex); tau = t - tc
    oot = np.abs(tau) > 0.07
    sg = rstd(f[oot] - np.polyval(np.polyfit(tau[oot], f[oot], 1), tau[oot]))
    M = np.column_stack([-s, np.ones_like(t), tau])
    c, *_ = np.linalg.lstsq(M, f, rcond=None)
    return c[0], sg * np.sqrt(np.linalg.inv(M.T @ M)[0, 0]), sg


def windows(pipe, drop58=False):
    W = []; betas = {}
    if pipe == "inv":
        d = pd.read_csv(f"{CL}/lc_M3V_INV_clean.csv", usecols=["time_btjd_tdb", "sector", "cadence_s", "quality_ok", "flux_raw_norm", "flare_flag"])
        d = d[(d.quality_ok == 1) & np.isfinite(d.flux_raw_norm) & (d.flare_flag == 0)]
        for sec, g in d.groupby("sector"):
            t = g.time_btjd_tdb.values; f = g.flux_raw_norm.values; ex = g.cadence_s.values[0] / 86400
            e0 = int(np.floor((t.min() - T0A) / PA)); e1 = int(np.ceil((t.max() - T0A) / PA))
            fk = []
            for e in range(e0, e1 + 1):
                for fph in np.r_[0.0, np.arange(0.12, 0.881, 0.04)]:
                    tc = T0A + (e + fph) * PA; m = np.abs(t - tc) < 0.30
                    if m.sum() < 8 or not accept(t[m], tc, ex): continue
                    dd, de, sg = linfit_depth(t[m], f[m], tc, ex)
                    if fph == 0: W.append(dict(sector=sec, epoch=e, t=t[m], f=f[m], ex=ex, sg=sg))
                    else: fk.append((dd, de))
            fk = np.array(fk); betas[int(sec)] = float(max(np.std(fk[:, 0]) / np.median(fk[:, 1]), 1.0))
    else:
        d = pd.read_csv(f"{CL}/lc_M3V_REF_clean.csv")
        for (sec, ep), g in d.groupby(["sector", "epoch"]):
            if drop58 and sec == 103 and ep == 58: continue
            t = g.time_btjd_tdb.values; f = g.flux_norm.values; ex = g.cadence_s.values[0] / 86400
            oot = np.abs(t - (T0A + ep * PA)) > 0.07
            sg = rstd(f[oot] - np.polyval(np.polyfit(t[oot] - t.mean(), f[oot], 1), t[oot] - t.mean()))
            W.append(dict(sector=sec, epoch=int(ep), t=t, f=f, ex=ex, sg=sg))
        z = np.load(f"{CL}/data/win_ref_fake.npz")
        for sec in np.unique(z["sector"]):
            fk = []
            for wi in np.where(z["sector"] == sec)[0]:
                m = z["w"] == wi; t = z["t"][m]; f = z["f"][m]; ex = z["ex"][m][0]; tc = z["tc0"][wi]
                if np.sum(np.abs(t - tc) < 0.028) < 2 or np.sum(np.abs(t - tc) > 0.07) < 8: continue
                fk.append(linfit_depth(t, f, tc, ex)[:2])
            fk = np.array(fk); betas[int(sec)] = float(max(np.std(fk[:, 0]) / np.median(fk[:, 1]), 1.0))
    return W, betas


class Like:
    def __init__(self, W, betas, ld="gauss"):
        self.t = np.concatenate([w["t"] for w in W]); self.f = np.concatenate([w["f"] for w in W])
        self.s = np.concatenate([np.full(len(w["t"]), w["sg"] * betas[int(w["sector"])]) for w in W])
        self.idx = np.concatenate([np.full(len(w["t"]), i) for i, w in enumerate(W)]); self.nw = len(W)
        self.tc0 = np.array([T0A + w["epoch"] * PA for w in W]); self.tau = self.t - self.tc0[self.idx]
        ex = np.concatenate([np.full(len(w["t"]), w["ex"]) for w in W])
        self.groups = []
        for e in np.unique(ex):
            m = ex == e; ns = 5 if e * 86400 < 300 else (9 if e * 86400 < 900 else 21)
            off = ((np.arange(ns) + 0.5) / ns - 0.5) * e
            tt = (self.t[m][:, None] + off[None]).ravel()
            p = batman.TransitParams(); p.t0, p.per, p.rp, p.a, p.inc, p.ecc, p.w, p.u, p.limb_dark = T0A, PA, 0.12, 25, 89, 0, 90, [0.1, 0.5], "quadratic"
            self.groups.append((m, ns, p, batman.TransitModel(p, tt)))
        self.w = 1 / self.s ** 2; self.ld = ld

    def model(self, T0, P, k, b, aR, u1, u2):
        out = np.ones(len(self.t))
        for m, ns, p, bm in self.groups:
            p.t0, p.per, p.rp, p.a, p.inc, p.u = T0, P, k, aR, np.degrees(np.arccos(b / aR)), [u1, u2]
            out[m] = bm.light_curve(p).reshape(-1, ns).mean(1)
        return out

    def chi2(self, mod):
        r = self.f - mod; w = self.w; i = self.idx; tau = self.tau; n = self.nw
        S0 = np.bincount(i, w, n); S1 = np.bincount(i, w * tau, n); S2 = np.bincount(i, w * tau ** 2, n)
        b0 = np.bincount(i, w * r, n); b1 = np.bincount(i, w * tau * r, n)
        det = S0 * S2 - S1 ** 2
        q = (S2 * b0 ** 2 - 2 * S1 * b0 * b1 + S0 * b1 ** 2) / det
        return float(np.sum(w * r * r) - q.sum())

    def lnp(self, p):
        T0, P, k, b, lr, a1, a2 = p
        if not (0.03 < k < 0.3 and 0 <= b < 1 + k and abs(P - PA) < 1e-3 and abs(T0 - T0A) < 0.02 and -1 < lr < 2): return -np.inf
        if self.ld == "gauss":
            u1, u2 = a1, a2
            if not (u1 > 0 and u1 + u2 < 1 and u1 + 2 * u2 > 0): return -np.inf
            lp = -0.5 * ((u1 - 0.125) / 0.1) ** 2 - 0.5 * ((u2 - 0.533) / 0.1) ** 2
        else:
            if not (0 < a1 < 1 and 0 < a2 < 1): return -np.inf
            sq = np.sqrt(a1); u1 = 2 * sq * a2; u2 = sq * (1 - 2 * a2); lp = 0.0
        rho = 10 ** lr
        lp += -0.5 * ((rho - 9.3) / 0.9) ** 2 + np.log(rho)
        aR = (G * rho * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
        return lp - 0.5 * self.chi2(self.model(T0, P, k, b, aR, u1, u2))


def middepth(k, b, aR, u1, u2):
    p = batman.TransitParams(); p.t0, p.per, p.rp, p.a, p.inc, p.ecc, p.w, p.u, p.limb_dark = 0, PA, k, aR, np.degrees(np.arccos(b / aR)), 0, 90, [u1, u2], "quadratic"
    return 1 - batman.TransitModel(p, np.array([0.0])).light_curve(p)[0]


def run(pipe, ld, drop58=False, seed=3):
    W, betas = windows(pipe, drop58)
    L = Like(W, betas, ld)
    x0 = np.array([T0A, PA, 0.118, 0.5, np.log10(9.3)] + ([0.125, 0.533] if ld == "gauss" else [0.43, 0.095]))
    best = minimize(lambda p: -L.lnp(p) if np.isfinite(L.lnp(p)) else 1e25, x0, method="Nelder-Mead", options=dict(maxiter=4000, xatol=1e-7, fatol=1e-3))
    rng = np.random.default_rng(seed); nw = 32
    sc = np.array([5e-4, 5e-6, 3e-3, 0.05, 0.02, 0.02, 0.02])
    pos = []
    while len(pos) < nw:
        x = best.x + sc * rng.normal(size=7)
        if np.isfinite(L.lnp(x)): pos.append(x)
    S = emcee.EnsembleSampler(nw, 7, L.lnp); S.run_mcmc(np.array(pos), 3000, progress=False)
    ch = S.get_chain(discard=1000, thin=5, flat=True)
    T0, P, k, b, lr, a1, a2 = ch.T
    if ld == "gauss": u1, u2 = a1, a2
    else: sq = np.sqrt(a1); u1 = 2 * sq * a2; u2 = sq * (1 - 2 * a2)
    rho = 10 ** lr; aR = (G * rho * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
    sub = np.arange(0, len(k), max(1, len(k) // 2000))
    dep = np.array([middepth(k[i], b[i], aR[i], u1[i], u2[i]) for i in sub]) * 1e3
    q = lambda x: [float(np.percentile(x, 50)), float(np.percentile(x, 50) - np.percentile(x, 16)), float(np.percentile(x, 84) - np.percentile(x, 50))]
    out = dict(pipe=pipe, ld=ld, drop58=drop58, n_windows=len(W), betas=betas, k=q(k), b=q(b), rho=q(rho), depth_ppt=q(dep), T0=q(T0), P=q(P),
               u1=q(u1), u2=q(u2), Rp_RE=q(k * np.random.default_rng(1).normal(0.383, 0.012, len(k)) / 0.0091577),
               acc=float(np.mean(S.acceptance_fraction)), tau=[float(x) for x in S.get_autocorr_time(discard=1000, quiet=True)])
    print(pipe, ld, drop58, "k", np.round(out["k"], 4), "b", np.round(out["b"], 2), "depth", np.round(out["depth_ppt"], 2), "betas", {k_: round(v, 2) for k_, v in betas.items()}, "tau", np.round(out["tau"], 0), flush=True)
    return out


if __name__ == "__main__":
    res = []
    for pipe, ld, d58 in [("inv", "gauss", False), ("ref", "gauss", False), ("ref", "gauss", True), ("inv", "wide", False), ("ref", "wide", False)]:
        res.append(run(pipe, ld, d58))
        json.dump(res, open(f"{R3}/out/q08_global.json", "w"), indent=1)
