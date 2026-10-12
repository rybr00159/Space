"""Single-star SED fit of TIC 87055304 with Castelli-Kurucz (2004) atmospheres + Fitzpatrick (1999) R_V=3.1 dust.
Free: Teff, logg, [M/H], A_V, angular diameter theta (mas), extra jitter sigma_j (mag, added in quadrature).
Variants (argv[1]):
  all      : all bands except W3 (W3 zero point inconsistent, see validate_synphot.py)
  nouv     : without SkyMapper u, v (the gravity/metallicity-sensitive bands)
  core     : Gaia G/BP/RP + 2MASS + WISE W1/W2 only
  g35      : all bands, logg fixed at 3.5 (GSP-Phot MARCS solution)
  g43      : all bands, logg fixed at 4.3 (dwarf solutions)
  mhprior  : all bands, [M/H] ~ N(-0.25, 0.15) (catalogue consensus)
  noapass  : all bands except APASS and Tycho (APASS/Tycho are excluded from every variant except 'all')
  avprior  : A_V <= total-column A_V (SF11) prior ; avmh: avprior + [M/H] prior ; g35av / g43av: logg fixed + A_V prior
Output: out/sed/<variant>.json, chain npz, fig/sed_<variant>.png
"""
import os, sys, json
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np, pandas as pd, emcee
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
import sedlib

VAR = sys.argv[1] if len(sys.argv) > 1 else "all"
NSTEP = int(sys.argv[2]) if len(sys.argv) > 2 else 5000
OUT = "/home/user/work/tic87/stellar/out/sed"; os.makedirs(OUT, exist_ok=True)
FIG = "/home/user/work/tic87/stellar/fig"
ph = pd.read_csv("/home/user/work/tic87/stellar/out/sed_photometry.csv")
ph = ph[ph.band != "W3"]
if VAR == "nouv": ph = ph[~ph.band.isin(["SM_u", "SM_v"])]
if VAR == "core": ph = ph[ph.band.isin(["G", "BP", "RP", "J", "H", "Ks", "W1", "W2"])]
# APASS (large apertures; the G=14.0 neighbour at 8.9" is 19% of the target flux) is 0.06-0.2 mag too bright and
# Tycho-2 has S/N ~ 4-5 -> both are excluded from every variant except 'all'
if VAR != "all": ph = ph[~ph.source.isin(["APASS DR9", "Tycho-2"])]
bands = list(ph.band); mobs = ph.mag.values; sig = ph.sigma.values
allb = list(pd.read_csv("/home/user/work/tic87/stellar/out/sed_photometry.csv").band) + ["TESS", "U", "R", "I"]
M = sedlib.SEDModel(allb)
idx = np.array([allb.index(b) for b in bands])
LOGG_FIX = {"g35": 3.5, "g43": 4.3, "g35av": 3.5, "g43av": 4.3}.get(VAR)


def unpack(p):
    if LOGG_FIX is None:
        T, lg, mh, av, lth, lsj = p
    else:
        T, mh, av, lth, lsj = p; lg = LOGG_FIX
    return T, lg, mh, av, lth, lsj


def lnprob(p):
    T, lg, mh, av, lth, lsj = unpack(p)
    if not (4000 < T < 8000 and 3.0 <= lg <= 5.0 and -1.5 < mh < 0.5 and 0 <= av < 1.5 and -7 < lsj < -1): return -np.inf
    mm = M.mags(T, lg, mh, av, 10 ** lth)[idx]
    if not np.all(np.isfinite(mm)): return -np.inf
    s2 = sig ** 2 + np.exp(2 * lsj)
    lp = -0.5 * np.sum((mobs - mm) ** 2 / s2 + np.log(s2))
    if VAR in ("mhprior", "avmh"): lp += -0.5 * ((mh + 0.25) / 0.15) ** 2
    if VAR in ("avprior", "avmh", "g35av", "g43av"):
        # A_V cannot exceed the total Galactic column (Schlafly & Finkbeiner 2011 E(B-V)=0.112 -> A_V=0.35);
        # half-Gaussian soft edge of 0.05 mag above 0.35
        if av > 0.35: lp += -0.5 * ((av - 0.35) / 0.05) ** 2
    return lp


p0c = [5950, 4.2, -0.3, 0.3, np.log10(0.030), -4] if LOGG_FIX is None else [5950, -0.3, 0.3, np.log10(0.030), -4]
ndim = len(p0c); nw = 40
sc = np.array([100, 0.2, 0.15, 0.08, 0.01, 0.3] if LOGG_FIX is None else [100, 0.15, 0.08, 0.01, 0.3])
p0 = []
rng = np.random.default_rng(1)
while len(p0) < nw:
    p = np.array(p0c) + sc * rng.standard_normal(ndim)
    if np.isfinite(lnprob(p)): p0.append(p)
S = emcee.EnsembleSampler(nw, ndim, lnprob)
S.run_mcmc(p0, NSTEP, progress=False)
ch = S.get_chain(discard=NSTEP // 3, thin=10, flat=True); lnp = S.get_log_prob(discard=NSTEP // 3, thin=10, flat=True)
pars = np.array([unpack(p) for p in ch])
T, lg, mh, av, lth, lsj = pars.T
th = 10 ** lth
fbol = np.array([M.fbol(a, b, c, d) for a, b, c, d in zip(T, lg, mh, th)])
ib = np.argmax(lnp); best = unpack(ch[ib])
mbest = M.mags(best[0], best[1], best[2], best[3], 10 ** best[4])
res = pd.DataFrame(dict(band=bands, obs=mobs, sig=sig, model=mbest[idx]))
res["resid"] = res.obs - res.model
chi2 = float(np.sum(res.resid ** 2 / res.sig ** 2))
chi2_j = float(np.sum(res.resid ** 2 / (res.sig ** 2 + np.exp(2 * best[5]))))


def q(x):
    p = np.percentile(x, [2.275, 15.87, 50, 84.13, 97.725])
    return dict(med=float(p[2]), lo1=float(p[2] - p[1]), hi1=float(p[3] - p[2]), p2=float(p[0]), p97=float(p[4]))


# radius per mas of angular diameter per pc: R[Rsun] = theta[rad]/2 * d / Rsun = theta[mas]*d[pc]*(1.495978707e11/6.957e8)/2/1000
K_R = 1.495978707e11 / 6.957e8 / 2 / 1000.
out = dict(variant=VAR, bands=bands, nstep=NSTEP, acc=float(np.mean(S.acceptance_fraction)),
           Teff=q(T), logg=q(lg), mh=q(mh), AV=q(av), theta_mas=q(th), sigma_jit=q(np.exp(lsj)), fbol_cgs=q(fbol),
           R_per_100pc=q(th * 100 * K_R), chi2_best_nojit=chi2, chi2_best_withjit=chi2_j, nbands=len(bands),
           lnp_max=float(lnp[ib]), best=dict(Teff=best[0], logg=best[1], mh=best[2], AV=best[3], theta=10 ** best[4], sj=np.exp(best[5])),
           resid=res.to_dict(orient="list"),
           synthetic_best={b: float(m) for b, m in zip(allb, mbest)})
json.dump(out, open(os.path.join(OUT, f"{VAR}.json"), "w"), indent=1)
np.savez_compressed(os.path.join(OUT, f"{VAR}_chain.npz"), T=T, logg=lg, mh=mh, av=av, theta=th, sj=np.exp(lsj), fbol=fbol, lnp=lnp)
print(VAR, "Teff %.0f +%.0f -%.0f" % (out["Teff"]["med"], out["Teff"]["hi1"], out["Teff"]["lo1"]),
      "logg %.2f +%.2f -%.2f" % (out["logg"]["med"], out["logg"]["hi1"], out["logg"]["lo1"]),
      "[M/H] %.2f +%.2f -%.2f" % (out["mh"]["med"], out["mh"]["hi1"], out["mh"]["lo1"]),
      "AV %.3f +%.3f -%.3f" % (out["AV"]["med"], out["AV"]["hi1"], out["AV"]["lo1"]),
      "theta %.5f +%.5f -%.5f mas" % (out["theta_mas"]["med"], out["theta_mas"]["hi1"], out["theta_mas"]["lo1"]),
      "Fbol %.4e" % out["fbol_cgs"]["med"], "jit %.3f" % out["sigma_jit"]["med"], "chi2(nojit) %.1f / %d" % (chi2, len(bands)),
      "lnpmax %.2f" % out["lnp_max"])
print(res.round(3).to_string())
# figure
fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw=dict(height_ratios=[3, 1]))
lp = np.array([sedlib.FILT[b]["lpiv"] for b in bands]) / 1e4
fobs = np.array([sedlib.FILT[b]["zp_flam"] for b in bands]) * 10 ** (-0.4 * (mobs - M.off[idx])) * lp * 1e4
fmod = np.array([sedlib.FILT[b]["zp_flam"] for b in bands]) * 10 ** (-0.4 * (mbest[idx] - M.off[idx])) * lp * 1e4
ax[0].errorbar(lp, fobs, yerr=fobs * sig / 1.0857, fmt="o", ms=4, label="observed")
ax[0].plot(lp, fmod, "s", mfc="none", label="model (band avg)")
fn = sedlib.ck04_files()
T0 = int(round(best[0] / 250) * 250); mh0 = min(sedlib.MHS, key=lambda x: abs(x - best[2]))
wl, fl = sedlib.load_spec(fn[(T0, mh0)], min(sedlib.LOGGS, key=lambda x: abs(x - best[1])))
import extinction
m_ = (wl > 2000) & (wl < 60000)
flx = fl[m_] * 10 ** (-0.4 * extinction.fitzpatrick99(wl[m_], best[3], 3.1)) * (10 ** best[4] / 206264806.247 / 2) ** 2
ax[0].plot(wl[m_] / 1e4, flx * wl[m_], "-", color="0.6", lw=0.6, zorder=0, label=f"CK04 {T0} K (nearest node)")
ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_ylabel(r"$\lambda F_\lambda$ (erg s$^{-1}$ cm$^{-2}$)")
ax[0].legend(fontsize=8); ax[0].set_title(f"SED fit '{VAR}': Teff={out['Teff']['med']:.0f} K, A_V={out['AV']['med']:.2f}, theta={out['theta_mas']['med']:.4f} mas")
ax[1].errorbar(lp, res.resid, yerr=sig, fmt="o", ms=4); ax[1].axhline(0, color="k", lw=0.5)
for x, y, b in zip(lp, res.resid, bands): ax[1].annotate(b, (x, y), fontsize=6)
ax[1].set_xlabel("wavelength (micron)"); ax[1].set_ylabel("obs - model (mag)")
fig.tight_layout(); fig.savefig(os.path.join(FIG, f"sed_{VAR}.png"), dpi=110)
