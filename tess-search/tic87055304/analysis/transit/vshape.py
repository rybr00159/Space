"""U-shape (physical, flat-bottomed) versus V-shape (grazing / no flat bottom) models for TIC 87055304.

Models (all exposure-integrated, per-sector multiplicative offsets profiled analytically, per-sector white noise
fixed at the free-density MCMC best-fit values):
  U  : batman, free k, b (0..1+k), rho*, u1, u2 (Claret Gaussian prior added to chi^2)    -> best physical model
  G  : batman forced grazing, b > 1 - k (same LD prior)
  V  : symmetric trapezoid with no flat bottom (triangle): depth, T14
  T  : general trapezoid: depth, T14, ingress fraction f = T12 / (T14/2) in (0, 1]  (f = 1 is V)
  Uc : batman forced non-grazing, b < 1 - k  (to show the flat bottom is required/allowed)
Red noise: beta per sector from the time-averaging method (Pont et al. 2006; Winn et al. 2008) applied to the U
residuals (all cadences in the transit windows), bin sizes 10-120 min, beta = mean ratio of observed to white-noise
expected rms (floored at 1). Delta chi^2_red = sum_s Delta chi^2_s / beta_s^2.
Also a profile-likelihood scan of the trapezoid ingress fraction f.
Writes vshape.json and vshape.png.
"""
import os, json
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.optimize import minimize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import tmodel as tm

HERE = os.path.dirname(os.path.abspath(__file__))
data = tm.load()
SECT = tm.SECTORS
mc = json.load(open(os.path.join(HERE, "mcmc_free.json")))
SIG = {s: np.exp(mc["best"][f"lns{s}"]) for s in SECT}
g = tm.gm()


def integ(prof):
    C = np.concatenate([[0.0], np.cumsum(0.5 * (prof[1:] + prof[:-1]) * g.dt)])
    F = {}
    for E in (200.0, 600.0, 1800.0):
        m = int(round(E / 86400 / 2 / g.dt)); Fe = np.ones_like(prof)
        Fe[m:-m] = (C[2 * m:] - C[:-2 * m]) / (2 * m * g.dt); F[E] = Fe
    return F


def eval_models(F, T0, P):
    out = {}
    for s in SECT:
        d = data[s]
        tau = d["t"] - T0 - np.round((d["t"] - T0) / P) * P
        out[s] = np.interp(tau, g.g, F[d["exptime"]])
    return out


def chi2_from_models(mods, per_sector=False):
    cs = {}
    for s in SECT:
        m, f = mods[s], data[s]["f"]
        c = np.sum(f * m) / np.sum(m * m)          # profiled multiplicative offset
        cs[s] = np.sum((f - c * m) ** 2) / SIG[s] ** 2
    return cs if per_sector else sum(cs.values())


def batman_prof(x):
    T0, P, k, b, rho, u1, u2 = x
    return g.profile(P, k, tm.a_from_rho(rho, P), b, u1, u2)


def chi2_batman(x, mode="free", per_sector=False):
    T0, P, k, b, rho, u1, u2 = x
    if not (0.01 < k < 0.6 and 0 <= b < 1 + k and 0.002 < rho < 10 and tm.ld_ok(u1, u2) and 3.9 < P < 4.4):
        return 1e30
    if mode == "grazing" and b <= 1 - k:
        return 1e30
    if mode == "nongrazing" and b >= 1 - k:
        return 1e30
    if tm.a_from_rho(rho, P) < 1.5:
        return 1e30
    mods = eval_models(integ(batman_prof(x)), T0, P)
    c = chi2_from_models(mods, per_sector)
    pen = ((u1 - tm.U1_MU) / tm.U1_SD) ** 2 + ((u2 - tm.U2_MU) / tm.U2_SD) ** 2
    return c if per_sector else c + pen


def trap_prof(depth, T14, f):
    a = np.abs(g.g); half = T14 / 2; T12 = f * half
    sh = np.clip((half - a) / max(T12, 1e-6), 0, 1)
    return 1 - depth * sh


def chi2_trap(x, f_fixed=None, per_sector=False):
    if f_fixed is None:
        T0, P, depth, T14, f = x
    else:
        T0, P, depth, T14 = x; f = f_fixed
    if not (0 < depth < 0.05 and 0.05 < T14 < 0.4 and 0 < f <= 1 and 3.9 < P < 4.4):
        return 1e30
    mods = eval_models(integ(trap_prof(depth, T14, f)), T0, P)
    return chi2_from_models(mods, per_sector)


def nm(fun, x0, scales, n=4, seed=3):
    rng = np.random.default_rng(seed)
    best = None
    for r in range(n):
        x = np.array(x0, float) + (0 if r == 0 else 1) * np.array(scales) * rng.standard_normal(len(x0))
        if fun(x) >= 1e29:
            x = np.array(x0, float)
        for it in range(3):
            res = minimize(fun, x, method="Nelder-Mead", options=dict(maxiter=6000, maxfev=6000, xatol=1e-9, fatol=1e-6, adaptive=True))
            x = res.x
        if best is None or res.fun < best.fun:
            best = res
    return best


def beta_time_avg(resid, t, exptime, bins_min=(10, 15, 20, 30, 45, 60, 90, 120)):
    s1 = np.std(resid)
    betas = []
    # contiguous chunks (split at gaps > 3 cadences)
    gaps = np.where(np.diff(t) > 3 * exptime / 86400)[0]
    chunks = np.split(np.arange(len(t)), gaps + 1)
    for bm in bins_min:
        N = int(round(bm * 60 / exptime))
        if N < 2:
            continue
        means = []
        for ch in chunks:
            M = len(ch) // N
            if M < 1:
                continue
            r = resid[ch][:M * N].reshape(M, N).mean(axis=1)
            means.append(r)
        means = np.concatenate(means)
        M = len(means)
        if M < 8:
            continue
        sN = np.std(means)
        exp = s1 / np.sqrt(N) * np.sqrt(M / (M - 1))
        betas.append(sN / exp)
    return max(1.0, float(np.mean(betas))), [float(b) for b in betas]


if __name__ == "__main__":
    xb = np.array([mc["best"][n] for n in ["T0", "P", "k", "b", "rho", "u1", "u2"]])
    res = {}
    U = nm(lambda x: chi2_batman(x, "free"), xb, [1e-4, 1e-7, 1e-3, 0.02, 0.005, 0.02, 0.02])
    G0 = xb.copy(); G0[3] = 1.0; G0[2] = 0.1
    Gs = [nm(lambda x: chi2_batman(x, "grazing"), y, [1e-4, 1e-7, 3e-3, 0.02, 0.01, 0.02, 0.02]) for y in
          (G0, np.r_[xb[:2], 0.2, 0.95, xb[4], xb[5:]], np.r_[xb[:2], 0.07, 0.94, xb[4], xb[5:]], np.r_[xb[:2], 0.4, 1.2, 0.05, xb[5:]])]
    G = min(Gs, key=lambda r: r.fun)
    Uc0 = xb.copy(); Uc0[3] = min(xb[3], 0.9 * (1 - xb[2]))
    Uc = nm(lambda x: chi2_batman(x, "nongrazing"), Uc0, [1e-4, 1e-7, 1e-3, 0.02, 0.005, 0.02, 0.02])
    dep0 = 1 - batman_prof(U.x).min(); T14_0 = tm.derived(np.array([U.x]))["T14_h"][0] / 24
    V = nm(lambda x: chi2_trap(x, f_fixed=1.0), [xb[0], xb[1], dep0 * 1.3, T14_0], [1e-4, 1e-7, 2e-4, 0.005])
    T = nm(lambda x: chi2_trap(x), [xb[0], xb[1], dep0, T14_0, 0.4], [1e-4, 1e-7, 2e-4, 0.005, 0.1])
    fits = dict(U=(U, lambda x: chi2_batman(x, "free", True), 5 + 2), G=(G, lambda x: chi2_batman(x, "grazing", True), 7),
                Uc=(Uc, lambda x: chi2_batman(x, "nongrazing", True), 7), V=(V, lambda x: chi2_trap(x, 1.0, True), 4),
                T=(T, lambda x: chi2_trap(x, None, True), 5))
    # per-sector residuals of U for beta
    mods = eval_models(integ(batman_prof(U.x)), U.x[0], U.x[1])
    beta = {}
    for s in SECT:
        m, f = mods[s], data[s]["f"]; c = np.sum(f * m) / np.sum(m * m)
        r = f - c * m
        beta[s], bl = beta_time_avg(r, data[s]["t"], data[s]["exptime"])
        res[f"beta_detail_{s}"] = bl
    persec = {k: v[1](v[0].x) for k, v in fits.items()}
    chi_tot = {k: float(sum(persec[k].values())) for k in fits}
    N = sum(len(data[s]["t"]) for s in SECT)
    out = dict(sigma_fixed_ppm={s: SIG[s] * 1e6 for s in SECT}, beta=beta, N=N, beta_detail=res,
               models={k: dict(x=v[0].x.tolist(), chi2=chi_tot[k], chi2_with_ldpen=float(v[0].fun), npar_shape=v[2],
                               chi2_per_sector={str(s): float(persec[k][s]) for s in SECT}) for k, v in fits.items()})
    for k in ("G", "V", "T", "Uc"):
        dc = {s: persec[k][s] - persec["U"][s] for s in SECT}
        out[f"dchi2_{k}_minus_U_white"] = float(sum(dc.values()))
        out[f"dchi2_{k}_minus_U_red"] = float(sum(dc[s] / beta[s] ** 2 for s in SECT))
        out[f"dchi2_{k}_minus_U_per_sector"] = {str(s): float(dc[s]) for s in SECT}
        out[f"dBIC_{k}_minus_U"] = float(out[f"dchi2_{k}_minus_U_white"] + (fits[k][2] - fits["U"][2]) * np.log(N))
    # derived for U and G
    for k in ("U", "G", "Uc"):
        x = fits[k][0].x
        der = tm.derived(np.array([x]))
        out["models"][k].update(dict(depth_ppm=float((1 - batman_prof(x).min()) * 1e6), T14_h=float(der["T14_h"][0]), T23_h=float(der["T23_h"][0]),
                                     a=float(der["a"][0]), inc=float(der["inc"][0])))
    xT = T.x
    out["models"]["T"].update(dict(depth_ppm=xT[2] * 1e6, T14_h=xT[3] * 24, f_ingress=xT[4], T23_over_T14=1 - xT[4]))
    out["models"]["V"].update(dict(depth_ppm=V.x[2] * 1e6, T14_h=V.x[3] * 24))
    # profile likelihood in f (trapezoid)
    fgrid = np.linspace(0.1, 1.0, 19)
    prof = []
    x = T.x[:4].copy()
    for fv in fgrid:
        r = nm(lambda y: chi2_trap(y, f_fixed=fv), x, [5e-5, 5e-8, 1e-4, 0.002], n=2)
        prof.append(r.fun)
    prof = np.array(prof)
    out["trap_profile"] = dict(f=fgrid.tolist(), chi2=prof.tolist())
    # 1-sigma / 2-sigma upper and lower bounds on f from the profile (white and red-scaled with the mean beta^2 weighting)
    beff = sum(out["models"]["T"]["chi2_per_sector"][str(s)] for s in SECT) / sum(out["models"]["T"]["chi2_per_sector"][str(s)] / beta[s] ** 2 for s in SECT)
    out["beta_eff_sq"] = float(beff)
    json.dump(out, open(os.path.join(HERE, "vshape.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k.startswith("dchi2") or k.startswith("dBIC") or k == "beta"}, indent=1))
    print({k: (round(v["chi2"], 2), {kk: vv for kk, vv in v.items() if kk in ("depth_ppm", "T14_h", "T23_h", "f_ingress", "x")}) for k, v in out["models"].items()})

    # figure: folded binned data (all sectors, 15-min bins) with U, G, V models
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    T0, P = U.x[0], U.x[1]
    tt = np.concatenate([data[s]["t"] for s in SECT]); ff = np.concatenate([data[s]["f"] / (np.sum(data[s]["f"] * mods[s]) / np.sum(mods[s] ** 2)) for s in SECT])
    ww = np.concatenate([np.full(len(data[s]["t"]), 1 / SIG[s] ** 2) for s in SECT])
    tau = (tt - T0 - np.round((tt - T0) / P) * P) * 24
    be = np.arange(-0.45 * 24, 0.45 * 24 + 0.01, 0.25)
    idx = np.digitize(tau, be)
    bc, bv, bs = [], [], []
    for i in range(1, len(be)):
        m = idx == i
        if m.sum() < 5: continue
        bc.append(0.5 * (be[i - 1] + be[i])); w = ww[m]
        bv.append(np.sum(ff[m] * w) / np.sum(w)); bs.append(1 / np.sqrt(np.sum(w)))
    bc, bv, bs = map(np.array, (bc, bv, bs))
    for a_ in ax:
        a_.errorbar(bc, (bv - 1) * 1e6, bs * 1e6, fmt="o", color="k", ms=3, lw=0.8, label="all sectors, 15-min bins")
    gg = g.g * 24
    ax[0].plot(gg, (integ(batman_prof(U.x))[200.0] - 1) * 1e6, "C0", label=f"U (batman, b={U.x[3]:.2f})  chi2={chi_tot['U']:.1f}")
    ax[0].plot(gg, (integ(batman_prof(G.x))[200.0] - 1) * 1e6, "C3", label=f"G (grazing, b={G.x[3]:.2f}, k={G.x[2]:.2f})  dchi2={chi_tot['G'] - chi_tot['U']:+.1f}")
    ax[0].plot(gg, (integ(trap_prof(V.x[2], V.x[3], 1.0))[200.0] - 1) * 1e6, "C2", label=f"V (triangle)  dchi2={chi_tot['V'] - chi_tot['U']:+.1f}")
    ax[0].set_xlim(-5, 5); ax[0].legend(fontsize=7, loc="lower left"); ax[0].set_xlabel("hours from mid-transit"); ax[0].set_ylabel("ppm")
    ax[1].plot(gg, (integ(batman_prof(U.x))[200.0] - integ(batman_prof(U.x))[200.0]) * 1e6, "C0")
    ax[1].errorbar(bc, (bv - np.interp(bc / 24, g.g, integ(batman_prof(U.x))[200.0])) * 1e6, bs * 1e6, fmt="o", color="k", ms=3, lw=0.8)
    ax[1].plot(gg, (integ(batman_prof(G.x))[200.0] - integ(batman_prof(U.x))[200.0]) * 1e6, "C3")
    ax[1].plot(gg, (integ(trap_prof(V.x[2], V.x[3], 1.0))[200.0] - integ(batman_prof(U.x))[200.0]) * 1e6, "C2")
    ax[1].set_xlim(-5, 5); ax[1].set_title("residuals from U and model differences (200-s integration shown)", fontsize=9)
    ax[1].set_xlabel("hours from mid-transit")
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "vshape.png"), dpi=130)
