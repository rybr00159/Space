"""Equal-mass (twin) eclipsing binary at 2P = 8.27 d as the source of the dips on TIC 87055315?

Data and noise model identical to the neighbour-host transit MCMC (nbmodel.load of data/lc_pair.npy, +-0.45 d windows,
per-sector baseline offsets c_s and white-noise levels sigma_s fixed at the transit fit's best values, so chi2 values
are directly comparable).
Twin model: star A (radius R, quadratic LD u1, u2 with the same Gaussian prior as the transit fit) eclipsed by an
identical star B (k = R_B/R_A in [0.9, 1.1], surface brightness equal -> flux ratio F_B/F_A = k^2), orbital period 2P,
circular orbit, so every dip at P is identical (primary = secondary). The EB contributes a fraction f_eb (0-1] of the
neighbour's light (f_eb < 1 = third light from an unresolved companion), and the neighbour is diluted in the pair by
D = 1/r (fixed at 5.32). Free: T0, b, a/R_A, k, f_eb, u1, u2, c_s. Maximum likelihood from 40 Nelder-Mead starts.
Also: the same comparison for the best transit model (k free) with b forced grazing (b > 1 - k) - the V-shape test.
Output: results/eb2p.json, fig/eb2p_fold.png
"""
import os, sys, json
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.optimize import minimize
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nbmodel as tm
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = tm.NBROOT
mc = json.load(open(os.path.join(ROOT, "mcmc/mcmc_pair_Dfree.json")))
best = mc["best"]
st = json.load(open(os.path.join(ROOT, "results/stellar.json")))["adopted"]
D = st["D_pair"]
data = tm.load(os.path.join(ROOT, "data/lc_pair.npy"))
S = tm.SECTORS
sig = {s: np.exp(best[f"lns{s}"]) for s in S}
P = best["P"]


def chi2_of(models):
    return float(sum(np.sum((data[s]["f"] - models[s]) ** 2) / sig[s] ** 2 for s in S))


# reference: best transit model, evaluated here with the same sigma_s
th7 = [best[n] for n in ["T0", "P", "k", "b", "rho", "u1", "u2"]]
mods = tm.phys_model(th7, data, best["D"])
ref = {s: mods[s] * (1 + best[f"c{s}"]) for s in S}
chi2_tr = chi2_of(ref)
ld = tm.LD


def ld_pen(u1, u2):
    return ((u1 - ld["u1_mean"]) / ld["u1_prior_sigma"]) ** 2 + ((u2 - ld["u2_mean"]) / ld["u2_prior_sigma"]) ** 2


g = tm.gm()


def twin_models(x):
    T0, b, aR, k, feb, u1, u2 = x[:7]
    c = x[7:]
    prof = g.profile(2 * P, k, aR, b, u1, u2)  # flux of star A alone (normalised to A)
    lam_sys = 1 - (1 - prof) * 1 / (1 + k ** 2)  # A + B, equal surface brightness
    lam_nb = 1 - feb * (1 - lam_sys)
    out = {}
    C = np.concatenate([[0.0], np.cumsum(0.5 * (lam_nb[1:] + lam_nb[:-1]) * g.dt)])
    for i, s in enumerate(S):
        d = data[s]; E = d["exptime"]
        m = int(round(E / 86400 / 2 / g.dt))
        Fe = np.ones_like(lam_nb); Fe[m:-m] = (C[2 * m:] - C[:-2 * m]) / (2 * m * g.dt)
        tau = d["t"] - T0 - np.round((d["t"] - T0) / P) * P
        mod = np.interp(tau, g.g, Fe)
        out[s] = ((mod + D) / (1 + D)) * (1 + c[i])
    return out


def obj_twin(x):
    T0, b, aR, k, feb, u1, u2 = x[:7]
    if not (0.9 <= k <= 1.1 and 0 <= b < 1 + k and 1.5 < aR < 60 and 0 < feb <= 1 and tm.ld_ok(u1, u2)) or np.any(np.abs(x[7:]) > 0.01):
        return 1e12
    return chi2_of(twin_models(x)) + ld_pen(u1, u2)


rng = np.random.default_rng(11)
fits = []
for i in range(40):
    x0 = np.array([best["T0"] + rng.normal(0, 0.002), rng.uniform(0.3, 1.9), rng.uniform(4, 30), rng.uniform(0.95, 1.05),
                   rng.uniform(0.1, 1.0), ld["u1_mean"], ld["u2_mean"]] + [best[f"c{s}"] for s in S])
    if obj_twin(x0) >= 1e12: continue
    r = minimize(obj_twin, x0, method="Nelder-Mead", options=dict(maxiter=6000, maxfev=6000, xatol=1e-6, fatol=1e-3, adaptive=True))
    r = minimize(obj_twin, r.x, method="Nelder-Mead", options=dict(maxiter=6000, maxfev=6000, xatol=1e-6, fatol=1e-3, adaptive=True))
    fits.append((r.fun, r.x))
fits.sort(key=lambda z: z[0])
bt = fits[0]
x = bt[1]
ndata = sum(len(data[s]["t"]) for s in S)
# implied stellar density of A from a/R at period 2P
rhoA = tm.rho_from_a(x[2], 2 * P)
res = dict(chi2_transit_best=chi2_tr + ld_pen(best["u1"], best["u2"]), chi2_twin_best=float(bt[0]),
           delta_chi2_twin_minus_transit=float(bt[0] - chi2_tr - ld_pen(best["u1"], best["u2"])), ndata=int(ndata),
           twin_best=dict(T0=x[0], b=x[1], aR=x[2], k=x[3], f_eb=x[4], u1=x[5], u2=x[6], rho_A_sun=float(rhoA)),
           n_starts=len(fits), top5_chi2=[float(f[0]) for f in fits[:5]],
           note="Same data and sigma_s as the transit fit; transit model has 7 shape parameters + D, twin model 7.")


# grazing transit (V-shape) test with the transit model at P
def obj_graz(y):
    T0, k, b, rho, u1, u2 = y[:6]; c = y[6:]
    if not (0.005 < k < 0.5 and 1 - k < b < 1 + k and 0.003 < rho < 5 and tm.ld_ok(u1, u2)) or tm.a_from_rho(rho, P) < 1.5:
        return 1e12
    mm = tm.phys_model([T0, P, k, b, rho, u1, u2], data, D)
    return chi2_of({s: mm[s] * (1 + c[i]) for i, s in enumerate(S)}) + ld_pen(u1, u2)


gf = []
for i in range(30):
    k0 = rng.uniform(0.1, 0.45)
    y0 = np.array([best["T0"], k0, 1 - k0 + rng.uniform(0.01, 0.5) * k0, rng.uniform(0.05, 1.0), ld["u1_mean"], ld["u2_mean"]] + [best[f"c{s}"] for s in S])
    if obj_graz(y0) >= 1e12: continue
    r = minimize(obj_graz, y0, method="Nelder-Mead", options=dict(maxiter=6000, maxfev=6000, adaptive=True))
    gf.append((r.fun, r.x))
gf.sort(key=lambda z: z[0])
res["grazing_transit_best"] = dict(chi2=float(gf[0][0]), delta_chi2_vs_free=float(gf[0][0] - res["chi2_transit_best"]),
                                   k=float(gf[0][1][1]), b=float(gf[0][1][2]), rho=float(gf[0][1][3]))
print(json.dumps(res, indent=1, default=float))
json.dump(res, open(os.path.join(ROOT, "results/eb2p.json"), "w"), indent=1, default=float)

# figure: binned residual comparison (all sectors, neighbour scale)
tw = twin_models(x)
fig, ax = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
for lab, mm, col in (("transit (best)", ref, "C1"), ("twin EB at 2P (best)", tw, "C2")):
    tt = np.concatenate([data[s]["t"] for s in S]); ff = np.concatenate([data[s]["f"] for s in S])
    mo = np.concatenate([mm[s] for s in S])
    tau = (tt - best["T0"] - np.round((tt - best["T0"]) / P) * P) * 24
    o = np.argsort(tau)
    be = np.arange(-10.8, 10.81, 0.3); idx = np.digitize(tau, be)
    bx = 0.5 * (be[1:] + be[:-1])
    by = [np.mean(ff[idx == i]) for i in range(1, len(be))]; bm = [np.mean(mo[idx == i]) for i in range(1, len(be))]
    sc = (1 + D) * 1e3
    if lab.startswith("transit"):
        ax[0].plot(bx, (np.array(by) - 1) * sc, "k.", label="pair data x (1+D), 18-min bins")
    ax[0].plot(bx, (np.array(bm) - 1) * sc, "-", color=col, label=lab)
    ax[1].plot(bx, (np.array(by) - np.array(bm)) * sc, ".-", color=col, label=lab)
ax[0].set_ylabel("neighbour flux - 1 (ppt)"); ax[0].legend(fontsize=7)
ax[1].axhline(0, color="k", lw=0.5); ax[1].set_ylabel("residual (ppt)"); ax[1].set_xlabel("hours from mid-eclipse")
ax[1].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "fig/eb2p_fold.png"), dpi=110)
print("DONE")
