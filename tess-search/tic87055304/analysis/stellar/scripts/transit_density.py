"""Re-derive the transit-implied circular-orbit stellar density (rho_circ) of TIC 87055304 with an MCMC fit of the
QLP SAP light curves (data/lc_sap_clean.npy, all 5 sectors, quality==0), with correct exposure-time integration per
sector (S13 1800 s, S39 600 s, S66/93/104 200 s), quadratic limb darkening (Kipping q1,q2) and a fixed dilution D
(fraction of in-aperture flux from other stars).

Usage: python transit_density.py <D> <ld_mode: free|prior> <nsteps> <tag>
Output: out/transit/<tag>_chain.npz, <tag>_summary.json, fig/transit_<tag>.png
Model: f_obs = 1 - (1 - D) * (1 - m_batman(t; P, t0, k, a/R*, b, u1, u2))
Parameters: t0, P, k (=Rp/R* undiluted), b, log10(rho_circ/rho_sun), q1, q2 ; a/R* from rho_circ and P.
"""
import os, sys, json, time
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import batman, emcee
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = float(sys.argv[1]); LD = sys.argv[2]; NSTEP = int(sys.argv[3]); TAG = sys.argv[4]
SECTORS = [int(x) for x in sys.argv[5].split(",")] if len(sys.argv) > 5 and sys.argv[5] != "all" else None
PSD = len(sys.argv) > 6 and sys.argv[6] == "psd"   # per-sector depth scaling (free extra dilution per sector)
OUT = "/home/user/work/tic87/stellar/out/transit"; os.makedirs(OUT, exist_ok=True)
FIG = "/home/user/work/tic87/stellar/fig"
G = 6.674e-11; RSUN = 6.957e8; MSUN = 1.989e30
RHOSUN = MSUN / (4 / 3 * np.pi * RSUN ** 3)
P0, T00, DUR = 4.136195, 1658.41325, 0.1875

lc = np.load("/home/user/work/tic87/data/lc_sap_clean.npy", allow_pickle=True).item()
secs = []
for s in sorted(lc):
    if SECTORS is not None and s not in SECTORS: continue
    v = lc[s]; t, f = v["t"], v["f"]
    m = np.isfinite(t) & np.isfinite(f)
    # mask the S104 scattered-light event noted in CONTEXT.md (BTJD ~4178.97, start of S104)
    if s == 104: m &= ~((t > 4178.5) & (t < 4179.5))
    ph = ((t - T00) / P0 + 0.5) % 1 - 0.5
    m &= np.abs(ph * P0) < 3 * DUR
    exptime = float(v["exptime"])
    ss = {1800.0: 15, 600.0: 5, 200.0: 3}.get(exptime, max(1, int(exptime / 120)))
    pr = batman.TransitParams()
    pr.t0, pr.per, pr.rp, pr.a, pr.inc, pr.ecc, pr.w, pr.u, pr.limb_dark = T00, P0, 0.07, 5., 85., 0., 90., [0.3, 0.2], "quadratic"
    mod = batman.TransitModel(pr, t[m], supersample_factor=ss, exp_time=exptime / 86400.)
    secs.append(dict(sector=s, t=t[m], f=f[m], sig=np.full(m.sum(), float(v["sigma"])), exptime=exptime, mod=mod, pr=pr))
print("points per sector:", [(s["sector"], len(s["t"])) for s in secs], flush=True)


def a_from_rho(rho, P):
    return ((rho * RHOSUN) * G * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)


NS = len(secs)


def lc_model(theta, sec):
    t0, P, k, b, lrho, q1, q2 = theta[:7]
    a = a_from_rho(10 ** lrho, P)
    inc = np.degrees(np.arccos(b / a))
    u1 = 2 * np.sqrt(q1) * q2; u2 = np.sqrt(q1) * (1 - 2 * q2)
    pr = sec["pr"]; pr.t0, pr.per, pr.rp, pr.a, pr.inc, pr.u = t0, P, k, a, inc, [u1, u2]
    m = sec["mod"].light_curve(pr)
    sc = 1.0
    if PSD:
        i = secs.index(sec)
        sc = 1.0 if i == 0 else theta[7 + i - 1]
    return 1 - (1 - D) * sc * (1 - m)

# Claret (2017) TESS LD for Teff~5900 K, logg~4.3, [Fe/H]~-0.3 (approx.): u1~0.32, u2~0.22 ; prior sigma 0.1
U1P, U2P = 0.32, 0.22


def lnprob(theta):
    t0, P, k, b, lrho, q1, q2 = theta[:7]
    if PSD and np.any((np.asarray(theta[7:]) < 0.3) | (np.asarray(theta[7:]) > 2.0)): return -np.inf
    if not (0 < q1 < 1 and 0 < q2 < 1 and 0.01 < k < 0.5 and 0 <= b < 1 + k and -3 < lrho < 1.5): return -np.inf
    if abs(P - P0) > 0.001 or abs(t0 - T00) > 0.05: return -np.inf
    a = a_from_rho(10 ** lrho, P)
    if a <= 1.5 or b >= a: return -np.inf
    lp = 0.0
    if LD == "prior":
        u1 = 2 * np.sqrt(q1) * q2; u2 = np.sqrt(q1) * (1 - 2 * q2)
        lp += -0.5 * ((u1 - U1P) / 0.1) ** 2 - 0.5 * ((u2 - U2P) / 0.1) ** 2
    chi2 = 0.0
    for s in secs:
        chi2 += np.sum(((s["f"] - lc_model(theta, s)) / s["sig"]) ** 2)
    return lp - 0.5 * chi2


# start: earlier least-squares solution, k rescaled for dilution
k0 = 0.0668 / np.sqrt(1 - D)
x0 = np.array([T00, P0, k0, 0.8, np.log10(0.063), 0.36, 0.30] + ([1.0] * (NS - 1) if PSD else []))
from scipy.optimize import minimize
best = None
for b0, lr0 in [(0.1, -0.6), (0.5, -0.9), (0.85, -1.2), (0.95, -1.4)]:
    x = x0.copy(); x[3] = b0; x[4] = lr0
    r = minimize(lambda th: -lnprob(th), x, method="Nelder-Mead", options=dict(maxiter=6000, xatol=1e-7, fatol=1e-3))
    print("start b0", b0, "->", r.x, -r.fun, flush=True)
    if best is None or r.fun < best.fun: best = r
ndim = len(x0); nw = max(36, 2 * ndim + 2)
scale = np.array([5e-4, 2e-6, 2e-3, 0.02, 0.03, 0.05, 0.05] + ([0.03] * (NS - 1) if PSD else []))
p0 = []
while len(p0) < nw:
    p = best.x + scale * np.random.randn(ndim)
    if np.isfinite(lnprob(p)): p0.append(p)
sampler = emcee.EnsembleSampler(nw, ndim, lnprob)
t1 = time.time()
sampler.run_mcmc(p0, NSTEP, progress=False)
print("mcmc time", time.time() - t1, flush=True)
chain = sampler.get_chain(discard=NSTEP // 3, thin=5, flat=True)
lnp = sampler.get_log_prob(discard=NSTEP // 3, thin=5, flat=True)
try:
    tau = sampler.get_autocorr_time(discard=NSTEP // 3, quiet=True).tolist()
except Exception:
    tau = None
t0, P, k, b, lrho, q1, q2 = chain.T[:7]
rho = 10 ** lrho
a = a_from_rho(rho, P)
inc = np.arccos(b / a)
T14 = P / np.pi * np.arcsin(np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) / (a * np.sin(inc)))
T23 = P / np.pi * np.arcsin(np.sqrt(np.clip((1 - k) ** 2 - b ** 2, 0, None)) / (a * np.sin(inc)))
u1 = 2 * np.sqrt(q1) * q2; u2 = np.sqrt(q1) * (1 - 2 * q2)
# max density allowed by T14 for a central transit (b=0), circular orbit, for each sample
amax = (1 + k) / np.sin(np.pi * T14 / P)
rho_max_b0 = 3 * np.pi * (amax * 1.0) ** 3 / (G * (P * 86400) ** 2) / RHOSUN


def q(x):
    p = np.percentile(x, [2.275, 15.87, 50, 84.13, 97.725])
    return dict(med=float(p[2]), lo1=float(p[2] - p[1]), hi1=float(p[3] - p[2]), p2=float(p[0]), p97=float(p[4]))


ibest = np.argmax(lnp)
summ = dict(D=D, ld=LD, nstep=NSTEP, nwalkers=nw, autocorr=tau, acc=float(np.mean(sampler.acceptance_fraction)),
            npts=int(sum(len(s["t"]) for s in secs)), chi2_best=float(-2 * (lnp[ibest] - (0 if LD == "free" else 0))),
            best=dict(zip(["t0", "P", "k", "b", "lrho", "q1", "q2"], chain[ibest].tolist())),
            t0=q(t0), P=q(P), k=q(k), b=q(b), rho_circ=q(rho), log10rho=q(lrho), aRs=q(a), T14_h=q(T14 * 24),
            T23_h=q(T23 * 24), u1=q(u1), u2=q(u2), rho_max_b0=q(rho_max_b0), depth_k2=q(k ** 2),
            sectors=[s_["sector"] for s_ in secs], psd=PSD,
            sector_scale={str(secs[i + 1]["sector"]): q(chain[:, 7 + i]) for i in range(NS - 1)} if PSD else None)
json.dump(summ, open(os.path.join(OUT, f"{TAG}_summary.json"), "w"), indent=1)
np.savez_compressed(os.path.join(OUT, f"{TAG}_chain.npz"), chain=chain, lnp=lnp, names=["t0", "P", "k", "b", "lrho", "q1", "q2"])
print(json.dumps({k_: summ[k_] for k_ in ["rho_circ", "aRs", "b", "k", "T14_h", "T23_h", "rho_max_b0", "acc", "autocorr"]}, indent=1))

# figure: phase-folded data (binned) with best model per cadence class
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
th = chain[ibest]
for s, c in zip(secs, ["C0", "C1", "C2", "C3", "C4"][:len(secs)]):
    ph = ((s["t"] - th[0]) / th[1] + 0.5) % 1 - 0.5
    o = np.argsort(ph)
    bins = np.arange(-0.6, 0.6, 0.02)
    hb = ph * th[1] * 24
    idx = np.digitize(hb, bins)
    xb = [hb[idx == i].mean() for i in range(1, len(bins)) if (idx == i).sum() > 3]
    yb = [s["f"][idx == i].mean() for i in range(1, len(bins)) if (idx == i).sum() > 3]
    ax[0].plot(xb, np.array(yb) - 0.004 * secs.index(s), ".", color=c, label=f"S{s['sector']}")
    ax[0].plot(hb[o], lc_model(th, s)[o] - 0.004 * secs.index(s), "-", color="k", lw=0.8)
ax[0].set_xlim(-12, 12); ax[0].set_xlabel("hours from mid-transit"); ax[0].set_ylabel("rel. flux + offset"); ax[0].legend(fontsize=7)
ax[0].set_title(f"D={D}, LD={LD}")
ax[1].hist2d(b, np.log10(rho), bins=60, cmap="Greys")
ax[1].set_xlabel("b"); ax[1].set_ylabel("log10 rho_circ [rho_sun]")
fig.tight_layout(); fig.savefig(os.path.join(FIG, f"transit_{TAG}.png"), dpi=110)
