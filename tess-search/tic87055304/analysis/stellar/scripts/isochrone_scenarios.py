"""Isochrone (MIST v1.2) inference of M*, R*, rho*, logg, age for TIC 87055304 under each distance scenario,
and comparison with the transit-implied circular-orbit density (photo-eccentric analysis) and companion radius.

Inputs
  * SED posterior (out/sed/<SEDVAR>_chain.npz): angular diameter theta, Teff, (A_V, [M/H]).
  * Adopted spectroscopic-like constraints: Teff = 5900 +- 150 K (spread of SED + 10 catalogue estimates),
    [Fe/H] = -0.30 +- 0.15 (spread of XP/SkyMapper catalogue values and SED).
  * Distance scenarios (parallax- or catalogue-based), each a distribution of d.
  * Transit posteriors (out/transit/<tag>_chain.npz): rho_circ and k for dilution D = 0, 0.336, 0.51.
Method
  * MIST isochrones interpolated to [Fe/H] steps of 0.05 at fixed (age, EEP); weights = Salpeter IMF (dN/dM ~ M^-2.35)
    x mass spacing along the isochrone x age prior uniform in linear age (0.3-13.8 Gyr) x Gaussian likelihoods in
    Teff, [Fe/H] and R (R from theta * d, with d drawn from the scenario; R likelihood width from the d spread and theta).
  * For each scenario: posterior samples of (M, R, rho) by weighted resampling of isochrone points.
  * Photo-eccentric: g = (rho_circ/rho_true)^(1/3) = (1 + e sin w)/sqrt(1 - e^2). e_min = |g^2 - 1|/(g^2 + 1).
    Posterior on e with uniform priors on e in [0, 0.95) and w in [0, 2pi), likelihood of g from the sampled
    distribution (KDE in log g), following Dawson & Johnson (2012).
  * Companion radius Rp = k R* for D = 0, 0.336 (TIC contratio 0.507 = F_cont/F_target), 0.51.
Output: out/scenarios.json, out/scenarios_table.csv, fig/scenarios_density.png, fig/scenarios_hrd.png
"""
import os, sys, json
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

rng = np.random.default_rng(42)
ROOT = "/home/user/work/tic87/stellar"
SEDVAR = sys.argv[1] if len(sys.argv) > 1 else "noapass"
TEFF, ETEFF = 5900., 150.
FEH, EFEH = -0.30, 0.15
RJ_PER_RSUN = 9.731  # R_sun / R_Jup (equatorial, IAU nominal: 6.957e8 / 7.1492e7)
K_R = 1.495978707e11 / 6.957e8 / 2 / 1000.  # R[Rsun] = theta[mas] * d[pc] * K_R

# ---------------------------------------------------------------- MIST, interpolated in [Fe/H]
z = np.load(os.path.join(ROOT, "cache/mist_compact.npz"))
tab = pd.DataFrame(z["tab"], columns=list(z["cols"]))
fehs = np.sort(tab.feh.unique())
parts = []
for f1, f2 in zip(fehs[:-1], fehs[1:]):
    a = tab[tab.feh == f1].set_index(["logage", "eep"]); b = tab[tab.feh == f2].set_index(["logage", "eep"])
    j = a.join(b, lsuffix="_1", rsuffix="_2", how="inner")
    for x in np.arange(0, 1.0, 0.2):
        fe = f1 + x * (f2 - f1)
        p = pd.DataFrame({"feh": fe, "logage": j.index.get_level_values(0), "eep": j.index.get_level_values(1)})
        for c in ["mini", "mass", "logL", "logT", "logR", "logg", "phase"]:
            p[c] = (1 - x) * j[c + "_1"].values + x * j[c + "_2"].values
        parts.append(p)
iso = pd.concat(parts + [tab[tab.feh == fehs[-1]]], ignore_index=True)
iso = iso.sort_values(["feh", "logage", "eep"]).reset_index(drop=True)
# mass spacing along each isochrone (IMF weight x dM)
dm = np.zeros(len(iso))
for (fe, la), g in iso.groupby(["feh", "logage"]):
    mi = g.mini.values
    d_ = np.gradient(mi) if len(mi) > 1 else np.ones(1)
    dm[g.index] = np.abs(d_)
iso["w_prior"] = dm * iso.mini ** -2.35 * 10 ** iso.logage  # age uniform in linear age on a log-age grid
iso = iso[(iso.logage <= np.log10(13.8e9)) & (iso.logage >= np.log10(3e8))]
iso["T"] = 10 ** iso.logT; iso["R"] = 10 ** iso.logR; iso["rho"] = iso.mass / iso.R ** 3
print("isochrone points", len(iso))

# ---------------------------------------------------------------- SED angular diameter
sed = np.load(os.path.join(ROOT, f"out/sed/{SEDVAR}_chain.npz"))
theta_s = sed["theta"]
TH, ETH = float(np.median(theta_s)), float(np.std(theta_s))
print("theta = %.5f +- %.5f mas" % (TH, ETH))

# ---------------------------------------------------------------- distance scenarios
N = 200000
def plx_dist(p, e, zp=0.0):
    x = rng.normal(p - zp, e, N); x = x[x > 0.2]; return 1000 / x
def lognorm_from_ci(med, lo, hi):
    s = (np.log(hi) - np.log(lo)) / 2; return np.exp(rng.normal(np.log(med), s, N))
gzp = -0.025   # approx. Lindegren+21 zero point for G=12.2, nu_eff=1.548, beta=-20.6 deg (5-p); see README
scen = {
    "S1_DR3_plx": ("Gaia DR3 parallax 2.807+-0.357 mas (ZP -0.025)", plx_dist(2.807, 0.357, gzp)),
    "S1b_DR3_plx_inflated": ("DR3 parallax, error x2 for RUWE 18.8 (2.807+-0.714)", plx_dist(2.807, 0.714, gzp)),
    "S2_DR2_plx": ("Gaia DR2 parallax 3.512+-0.735 mas (ZP -0.03)", plx_dist(3.512, 0.735, -0.03)),
    "S3_BJ21_rpgeo": ("Bailer-Jones+21 photogeometric 359 (318-394) pc", lognorm_from_ci(358.7, 318.3, 394.3)),
    "S4_StarHorse21": ("StarHorse 2021: 395 (288-524) pc", lognorm_from_ci(394.8, 288.1, 524.1)),
    "S5_Zhang23_plx": ("Zhang+23 XP+Gaia revised parallax 2.888+-0.275 mas", plx_dist(2.888, 0.275)),
    "S6_MSC_458": ("Gaia MSC (binary model) 458 (344-575) pc", lognorm_from_ci(457.6, 344.3, 574.9)),
    "S7_TIC_464": ("TIC v8.2 distance 464 pc (+-10%)", lognorm_from_ci(464, 420, 510)),
    "S8_GSPPhot_A_843": ("GSP-Phot A-library 843 pc (+-5%)", lognorm_from_ci(843, 801, 885)),
    "S9_GSPPhot_MARCS_998": ("GSP-Phot MARCS (best) 998 pc (+-5%)", lognorm_from_ci(998, 948, 1048)),
}
extra = {}
if os.path.exists(os.path.join(ROOT, "out/xp/knn_result.json")):
    kr = json.load(open(os.path.join(ROOT, "out/xp/knn_result.json")))
    scen["S10_XP_kNN"] = (f"XP k-NN spectro-photometric distance {kr['d_med']:.0f} pc", lognorm_from_ci(kr["d_med"], kr["d_p16"], kr["d_p84"]))

# ---------------------------------------------------------------- transit posteriors
tr = {}
for D, tag in [(0.0, "D000_prior_psd"), (0.336, "D0336_prior"), (0.51, "D051_prior")]:
    f = os.path.join(ROOT, f"out/transit/{tag}_chain.npz")
    if not os.path.exists(f):
        f = os.path.join(ROOT, "out/transit/D000_prior_chain.npz") if D == 0 else None
    if f is None or not os.path.exists(f): continue
    c = np.load(f)["chain"]
    tr[D] = dict(rho=10 ** c[:, 4], k=c[:, 2], b=c[:, 3], tag=os.path.basename(f))
rho_circ = tr[0.0]["rho"]
print("transit chains:", {k: v["tag"] for k, v in tr.items()})


def q(x, p=(2.275, 15.87, 50, 84.13, 97.725)):
    v = np.percentile(x, p)
    return dict(med=float(v[2]), lo1=float(v[2] - v[1]), hi1=float(v[3] - v[2]), p2=float(v[0]), p97=float(v[4]))


def iso_posterior(R_obs, eR, T=TEFF, eT=ETEFF, fe=FEH, efe=EFEH, n=40000):
    lnl = -0.5 * ((iso["T"] - T) / eT) ** 2 - 0.5 * ((iso.feh - fe) / efe) ** 2 - 0.5 * ((iso.R - R_obs) / eR) ** 2
    w = iso.w_prior * np.exp(lnl - lnl.max())
    w = w / w.sum()
    idx = rng.choice(len(iso), n, p=w.values)
    s = iso.iloc[idx]
    return s, float(np.log(np.sum(iso.w_prior * np.exp(lnl)) / np.sum(iso.w_prior)))


# photo-eccentric posterior on e (Dawson & Johnson 2012): uniform e in [0,0.95), w uniform
E_GRID = np.linspace(0, 0.95, 96); W_GRID = np.linspace(0, 2 * np.pi, 181)[:-1]
EE, WW = np.meshgrid(E_GRID, W_GRID, indexing="ij")
G_EW = (1 + EE * np.sin(WW)) / np.sqrt(1 - EE ** 2)


def ecc_post(g_samples):
    lg = np.log(g_samples[np.isfinite(g_samples) & (g_samples > 0)])
    kde = gaussian_kde(lg[rng.choice(len(lg), min(len(lg), 20000), replace=False)])
    L = kde(np.log(G_EW).ravel()).reshape(G_EW.shape)
    pe = L.sum(axis=1); pe /= pe.sum()
    cdf = np.cumsum(pe)
    return dict(e_med=float(np.interp(0.5, cdf, E_GRID)), e_p2=float(np.interp(0.0228, cdf, E_GRID)),
                e_p16=float(np.interp(0.1587, cdf, E_GRID)), P_e_lt_0p1=float(pe[E_GRID < 0.1].sum()),
                P_e_lt_0p3=float(pe[E_GRID < 0.3].sum()), P_e_lt_0p5=float(pe[E_GRID < 0.5].sum()))


rows = []; out = {}
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for key, (lab, dsamp) in scen.items():
    dsamp = dsamp[np.isfinite(dsamp)]
    th = rng.normal(TH, ETH, len(dsamp))
    Rs = th * dsamp * K_R
    R_obs = float(np.median(Rs)); eR = float(0.5 * (np.percentile(Rs, 84.13) - np.percentile(Rs, 15.87)))
    s, lnZ = iso_posterior(R_obs, eR)
    M, R, rho = s.mass.values, s.R.values, s.rho.values
    # combine with transit (independent draws)
    nn = min(len(rho), len(rho_circ))
    rc = rng.choice(rho_circ, nn); rt = rng.choice(rho, nn)
    g = (rc / rt) ** (1 / 3)
    emin = np.abs(g ** 2 - 1) / (g ** 2 + 1)
    ep = ecc_post(g)
    comp = {}
    for D, t in tr.items():
        k = rng.choice(t["k"], nn); Rp = k * rng.choice(R, nn) * RJ_PER_RSUN
        comp[f"D={D}"] = q(Rp)
    o = dict(label=lab, d_pc=q(dsamp), R_from_theta=dict(med=R_obs, sig=eR), lnZ_iso=lnZ,
             M=q(M), R=q(R), rho=q(rho), logg=q(s.logg.values), age_Gyr=q(10 ** s.logage.values / 1e9), Teff_iso=q(s["T"].values),
             feh_iso=q(s.feh.values), frac_MS=float(np.mean(s.phase.values < 1)),
             g=q(g), e_min=q(emin), ecc_posterior=ep, Rp_RJ=comp,
             rho_ratio_circ_over_true=q(rc / rt))
    out[key] = o
    rows.append(dict(scenario=key, d_pc=f"{o['d_pc']['med']:.0f} (+{o['d_pc']['hi1']:.0f}/-{o['d_pc']['lo1']:.0f})",
                     R_theta=f"{R_obs:.2f}+-{eR:.2f}",
                     M=f"{o['M']['med']:.2f} (+{o['M']['hi1']:.2f}/-{o['M']['lo1']:.2f})",
                     R=f"{o['R']['med']:.2f} (+{o['R']['hi1']:.2f}/-{o['R']['lo1']:.2f})",
                     rho=f"{o['rho']['med']:.3f} (+{o['rho']['hi1']:.3f}/-{o['rho']['lo1']:.3f})",
                     logg=f"{o['logg']['med']:.2f}", age=f"{o['age_Gyr']['med']:.1f}",
                     g=f"{o['g']['med']:.2f} (+{o['g']['hi1']:.2f}/-{o['g']['lo1']:.2f})",
                     emin=f"{o['e_min']['med']:.2f} [{o['e_min']['p2']:.2f}]",
                     P_e_lt_0p1=f"{ep['P_e_lt_0p1']:.3f}", e_post=f"{ep['e_med']:.2f} (2.3%: {ep['e_p2']:.2f})",
                     **{f"Rp_RJ_{kk}": f"{vv['med']:.2f} (+{vv['hi1']:.2f}/-{vv['lo1']:.2f})" for kk, vv in comp.items()}))
    ax[0].hist(np.log10(rho), bins=80, histtype="step", density=True, label=f"{key}")
ax[0].hist(np.log10(rho_circ), bins=80, histtype="stepfilled", color="k", alpha=0.3, density=True, label="transit rho_circ (D=0)")
ax[0].set_xlabel("log10 rho* [rho_sun]"); ax[0].legend(fontsize=6); ax[0].set_title("Stellar density per scenario vs transit (circular)")
# R* vs distance curve
dd = np.linspace(200, 1200, 200)
ax[1].plot(dd, TH * dd * K_R, "k-", label="R* = theta d / 2 (SED)")
ax[1].fill_between(dd, (TH - ETH) * dd * K_R, (TH + ETH) * dd * K_R, color="k", alpha=0.2)
for key, o in out.items():
    ax[1].errorbar(o["d_pc"]["med"], o["R"]["med"], xerr=[[o["d_pc"]["lo1"]], [o["d_pc"]["hi1"]]], yerr=[[o["R"]["lo1"]], [o["R"]["hi1"]]], fmt="o", ms=4, label=key)
# radius at which a MIST star (Teff, [Fe/H]) has rho = rho_circ
sel = (np.abs(iso["T"] - TEFF) < ETEFF) & (np.abs(iso.feh - FEH) < EFEH)
rc_med = np.median(rho_circ)
ok = sel & (np.abs(np.log10(iso.rho / rc_med)) < 0.05)
if ok.sum():
    Rm = np.median(iso.R[ok]); ax[1].axhline(Rm, color="r", ls="--", label=f"R* with rho*=rho_circ ({Rm:.2f})")
    extra["R_at_rho_circ"] = dict(R=float(Rm), d_pc=float(Rm / (TH * K_R)), M=float(np.median(iso.mass[ok])), n=int(ok.sum()),
                                  R_p16=float(np.percentile(iso.R[ok], 16)), R_p84=float(np.percentile(iso.R[ok], 84)))
ax[1].set_xlabel("distance (pc)"); ax[1].set_ylabel("R* (R_sun)"); ax[1].legend(fontsize=6)
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "fig/scenarios_density.png"), dpi=110)

# maximum density allowed by the transit duration for b = 0 (circular) -> minimum R*
f0 = os.path.join(ROOT, "out/transit/D000_prior_psd_summary.json")
if not os.path.exists(f0): f0 = os.path.join(ROOT, "out/transit/D000_prior_summary.json")
ts = json.load(open(f0))
extra["transit_summary_used"] = os.path.basename(f0)
extra["rho_circ"] = ts["rho_circ"]; extra["rho_max_b0"] = ts["rho_max_b0"]
extra["theta_mas"] = dict(med=TH, sig=ETH); extra["sed_variant"] = SEDVAR
extra["adopted"] = dict(Teff=[TEFF, ETEFF], feh=[FEH, EFEH])
json.dump(dict(scenarios=out, extra=extra), open(os.path.join(ROOT, "out/scenarios.json"), "w"), indent=1)
pd.DataFrame(rows).to_csv(os.path.join(ROOT, "out/scenarios_table.csv"), index=False)
print(pd.DataFrame(rows).to_string())
print(json.dumps(extra, indent=1))
