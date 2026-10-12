"""Stellar parameters of the neighbour TIC 87055315 = Gaia DR3 6720986208895208960.

1. Catalogue compilation (cached CSVs from s01_catalogs.py): TIC v8.2, Gaia DR3 gaia_source / GSP-Phot / FLAME / MSC,
   StarHorse 2019 (I/349) and 2021 (I/354), SHBoost (V/160), Zhang+23, Huson+25, Ye+25, Huang+22 (SMSS), BJ21.
2. Distance: Gaia DR3 parallax corrected with the Lindegren+21 zero point (`zero_point` package); error inflated by 1.3
   (El-Badry+21 for G ~ 14) and 0.010 mas (systematic floor) added in quadrature; cross-checked with Bailer-Jones+21.
3. SED fit (local copy of stellar/scripts/sedlib.py: Castelli & Kurucz 2004 ATLAS9 + Fitzpatrick 1999 R_V = 3.1):
   bands Gaia G/BP/RP, 2MASS J/H/Ks, SkyMapper DR4 g/r/i/z (PSF), AllWISE W1/W2. Free: Teff, logg, [M/H], A_V,
   distance d, radius R (theta = 2R/d), jitter. Priors: parallax likelihood; [M/H] ~ N(-0.05, 0.20) (catalogue
   consensus, see below); A_V flat in [0, 0.36] with a 0.05-mag half-Gaussian edge above the total Galactic column
   (SFD E(B-V) = 0.117-0.131 x 0.86 x 3.1); logg flat 3.0-5.0. emcee, 32 walkers x 6000 steps.
   SkyMapper u/v are left out (strongly metallicity/gravity sensitive, model-dependent); WISE W1/W2 have a 0.05-mag
   floor (6" beam, 1.8-mag brighter star at 9.3").
4. Mass, age, rho*: grid-based MIST v1.2 inference (stellar/cache/mist_compact.npz, read-only): every isochrone point
   is weighted by exp(-chi2/2) for (log Teff, log R) from the SED posterior (2-D Gaussian with the posterior
   covariance) and [Fe/H] ~ N(feh_sed, sd), times the prior: Kroupa IMF x d(M_init) along the isochrone, uniform age
   (10^logage x d logage), uniform [Fe/H] over the grid nodes. Gives M*, age, logg, rho* = M/R^3 and the phase.
5. Flux ratio r = F_nb/F_tg in the TESS band (from s03 results/fluxratio.json) and the SED-synthetic T magnitude.
Output: results/stellar.json, fig/sed_nb.png, fig/hrd_nb.png
"""
import os, sys, json
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np, pandas as pd
import emcee
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sedlib

ROOT = "/home/user/work/tic87/neighbour"
CAT = os.path.join(ROOT, "data/cat")
NB = 6720986208895208960
rng = np.random.default_rng(7)
res = {}

# ------------------------------------------------------------------ 1. catalogue compilation
g = pd.read_csv("/home/user/work/tic87/data/gaia_dr3_90arcsec.csv").set_index("source_id").loc[NB]
ap = pd.read_csv(os.path.join(CAT, "gaiadr3_astrophysical_parameters_nb.csv")).iloc[0]
tic = pd.read_csv(os.path.join(CAT, "tic82_87055315.csv")).iloc[0]
rd = lambda n: pd.read_csv(os.path.join(CAT, f"viz_{n}.csv")).iloc[0]
sh19, sh21, shb, z23, h25, y25, h22, bj = (rd(n) for n in ["StarHorse19", "StarHorse21", "SHBoost24", "Zhang23", "Huson25",
                                                           "Ye25", "Huang22", "BJ21"])
comp = {
    "TIC_v8.2": dict(Teff=tic.Teff, e_Teff=tic.e_Teff, logg=tic.logg, R=tic.rad, M=tic.mass, rho=tic.rho, d_pc=tic.d,
                     ebv=tic.ebv, Tmag=tic.Tmag, plx=tic.plx, note="DR2-based (plx 0.761), BJ18 distance"),
    "Gaia_GSPPhot": dict(Teff=ap.teff_gspphot, Teff_lo=ap.teff_gspphot_lower, Teff_hi=ap.teff_gspphot_upper,
                         logg=ap.logg_gspphot, mh=ap.mh_gspphot, R=ap.radius_gspphot, R_lo=ap.radius_gspphot_lower,
                         R_hi=ap.radius_gspphot_upper, d_pc=ap.distance_gspphot, A0=ap.azero_gspphot, lib=ap.libname_gspphot),
    "Gaia_FLAME": dict(R=ap.radius_flame, R_lo=ap.radius_flame_lower, R_hi=ap.radius_flame_upper, L=ap.lum_flame,
                       M=ap.mass_flame, M_lo=ap.mass_flame_lower, M_hi=ap.mass_flame_upper, age=ap.age_flame,
                       age_lo=ap.age_flame_lower, age_hi=ap.age_flame_upper, evolstage=int(ap.evolstage_flame),
                       flags=str(ap.flags_flame)),
    "Gaia_MSC": dict(Teff1=ap.teff_msc1, Teff2=ap.teff_msc2, logg1=ap.logg_msc1, logg2=ap.logg_msc2, mh=ap.mh_msc,
                     d_pc=ap.distance_msc, A0=ap.azero_msc, logpost=ap.logposterior_msc,
                     note="binary-star model; distance 738 pc conflicts with the parallax (1300 pc)"),
    "Gaia_ESP": dict(spectraltype_esphs=ap.spectraltype_esphs, flags_esphs=int(ap.flags_esphs),
                     ew_halpha=ap.ew_espels_halpha, classprob_binarystar=ap.classprob_dsc_combmod_binarystar),
    "StarHorse19_I349": dict(Teff=sh19.teff50, logg=sh19.logg50, mh=sh19.met50, M=sh19.mass50, M16=sh19.mass16,
                             M84=sh19.mass84, d_kpc=sh19.dist50, AV=sh19.AV50, note="DR2-based"),
    "StarHorse21_I354": dict(Teff=sh21.teff50, logg=sh21.logg50, mh=sh21.met50, M=sh21.mass50, M16=sh21.mass16,
                             M84=sh21.mass84, d_kpc=sh21.dist50, AV=sh21.AV50),
    "SHBoost24_V160": dict(Teff=10 ** shb.logTeff, logg=shb.logg, mh=shb.Met, M=shb.Mass, sM=shb.s_Massmean,
                           AV=shb.AV, d_kpc=shb.Dist),
    "Zhang23": dict(Teff=z23.Teff * 1000, e_Teff=z23.e_Teff * 1000, logg=z23.logg, e_logg=z23.e_logg, mh=z23["[Fe/H]"],
                    E=z23.Ext),
    "Huson25": dict(Teff=10 ** h25.logTeff, logg=h25["log(g)"], mh=h25["[Fe/H]"]),
    "Ye25": dict(Teff=y25.Teff, logg=y25.logg, mh=y25["[Fe/H]"]),
    "Huang22_SMSS": dict(Teff=h22.Teff, e_Teff=h22.e_Teff, mh=h22["[Fe/H]"], e_mh=h22["e_[Fe/H]"], age=h22.Age,
                         e_age=h22.e_Age, type=f"{h22.Type}/{h22.subType}", ebv_sfd=h22["E(B-V)sfd"]),
    "BJ21": dict(rgeo=bj.rgeo, rgeo_lo=bj.b_rgeo, rgeo_hi=bj.B_rgeo, rpgeo=bj.rpgeo, rpgeo_lo=bj.b_rpgeo,
                 rpgeo_hi=bj.B_rpgeo),
}
comp = {k: {kk: (float(vv) if isinstance(vv, (int, float, np.floating, np.integer)) else vv) for kk, vv in v.items()}
        for k, v in comp.items()}
res["catalogues"] = comp
teffs = [comp[k]["Teff"] for k in ["TIC_v8.2", "Gaia_GSPPhot", "StarHorse21_I354", "SHBoost24_V160", "Zhang23", "Huson25",
                                   "Ye25", "Huang22_SMSS"]]
mhs = [comp[k]["mh"] for k in ["Gaia_GSPPhot", "Gaia_MSC", "StarHorse21_I354", "SHBoost24_V160", "Zhang23", "Huson25",
                               "Ye25", "Huang22_SMSS"]]
res["catalogue_consensus"] = dict(Teff_median=float(np.median(teffs)), Teff_sd=float(np.std(teffs, ddof=1)),
                                  mh_median=float(np.median(mhs)), mh_sd=float(np.std(mhs, ddof=1)),
                                  note="StarHorse19 (Teff 7131, A_V 1.14) excluded: DR2-based, A_V above the total column")
gaia_astrom = {k: (float(g[k]) if np.isscalar(g[k]) and not isinstance(g[k], str) else str(g[k])) for k in
               ["parallax", "parallax_error", "pmra", "pmdec", "ruwe", "astrometric_excess_noise",
                "astrometric_excess_noise_sig", "ipd_frac_multi_peak", "ipd_gof_harmonic_amplitude", "ipd_frac_odd_win",
                "visibility_periods_used", "non_single_star", "duplicated_source", "radial_velocity",
                "radial_velocity_error", "rv_method_used", "rv_nb_transits", "rv_nb_deblended_transits",
                "rv_expected_sig_to_noise", "rv_template_teff", "rv_template_logg", "grvs_mag", "phot_g_mean_mag",
                "phot_bp_mean_mag", "phot_rp_mean_mag", "bp_rp", "phot_bp_rp_excess_factor", "phot_variable_flag",
                "has_rvs", "has_epoch_photometry", "l", "b"]}
res["gaia_dr3_source"] = gaia_astrom

# ------------------------------------------------------------------ 2. distance
from zero_point import zpt
zpt.load_tables()
zp = float(zpt.get_zpt(np.array([g.phot_g_mean_mag]), np.array([g.nu_eff_used_in_astrometry]), np.array([np.nan]),
                     np.array([g.ecl_lat]), np.array([int(g.astrometric_params_solved)]))[0])
plx = float(g.parallax - zp)
plx_err = float(np.hypot(1.3 * g.parallax_error, 0.010))
res["distance"] = dict(parallax=float(g.parallax), zero_point_L21=zp, parallax_corr=plx, parallax_err_adopted=plx_err,
                       d_inv_plx_pc=1000 / plx, d_inv_plx_err=1000 / plx ** 2 * plx_err,
                       BJ21_rgeo=comp["BJ21"]["rgeo"], BJ21_rpgeo=comp["BJ21"]["rpgeo"])
print("distance", res["distance"], flush=True)

# ------------------------------------------------------------------ 3. SED fit
two = rd("2MASS"); sm = rd("SMSS_DR4"); wi = rd("AllWISE")
phot = [  # band, mag, sigma (catalogue error + floor)
    ("G", g.phot_g_mean_mag, np.hypot(1.0857 / g.phot_g_mean_flux_over_error, 0.010)),
    ("BP", g.phot_bp_mean_mag, np.hypot(1.0857 / g.phot_bp_mean_flux_over_error, 0.015)),
    ("RP", g.phot_rp_mean_mag, np.hypot(1.0857 / g.phot_rp_mean_flux_over_error, 0.015)),
    ("J", two.Jmag, np.hypot(two.e_Jmag, 0.02)), ("H", two.Hmag, np.hypot(two.e_Hmag, 0.02)),
    ("Ks", two.Kmag, np.hypot(two.e_Kmag, 0.02)),
    ("SM_g", sm.gPSF, np.hypot(sm.e_gPSF, 0.03)), ("SM_r", sm.rPSF, np.hypot(sm.e_rPSF, 0.03)),
    ("SM_i", sm.iPSF, np.hypot(sm.e_iPSF, 0.03)), ("SM_z", sm.zPSF, np.hypot(sm.e_zPSF, 0.03)),
    ("W1", wi.W1mag, np.hypot(wi.e_W1mag, 0.05)), ("W2", wi.W2mag, np.hypot(wi.e_W2mag, 0.05)),
]
bands = [p[0] for p in phot]; mobs = np.array([p[1] for p in phot], float); sig = np.array([p[2] for p in phot], float)
res["sed_photometry"] = {b: dict(mag=float(m), sigma=float(s)) for b, m, s in zip(bands, mobs, sig)}
allb = [str(b) for b in np.load(sedlib.GRID, allow_pickle=True)["bands"]]  # band order of the cached grid
M = sedlib.SEDModel(allb)
idx = np.array([allb.index(b) for b in bands])
iT = allb.index("TESS")
MH_MU, MH_SD = -0.05, 0.20
AV_MAX = 0.36
RSUN_PC = 2.25461e-8  # R_sun in pc


def theta_mas(R, d):
    return 2 * R * RSUN_PC / d * 206264806.247


def lnprob(p):
    T, lg, mh, av, d, R, lsj = p
    if not (4500 < T < 7500 and 3.0 <= lg <= 5.0 and -1.0 < mh < 0.5 and 0 <= av < 1.0 and 300 < d < 5000
            and 0.3 < R < 6 and -7 < lsj < -1):
        return -np.inf
    mm = M.mags(T, lg, mh, av, theta_mas(R, d))[idx]
    if not np.all(np.isfinite(mm)): return -np.inf
    s2 = sig ** 2 + np.exp(2 * lsj)
    lp = -0.5 * np.sum((mobs - mm) ** 2 / s2 + np.log(s2))
    lp += -0.5 * ((1000 / d - plx) / plx_err) ** 2
    lp += -0.5 * ((mh - MH_MU) / MH_SD) ** 2
    if av > AV_MAX: lp += -0.5 * ((av - AV_MAX) / 0.05) ** 2
    lp += 2 * np.log(d)  # uniform space density prior (r^2); weak for plx/err ~ 25
    return lp


x0 = np.array([5800, 4.0, -0.05, 0.3, 1000 / plx, 1.95, -4.0])
from scipy.optimize import minimize
x0 = minimize(lambda x: -lnprob(x) if np.isfinite(lnprob(x)) else 1e20, x0, method="Nelder-Mead",
              options=dict(maxiter=6000, maxfev=6000)).x
nw, nd = 32, 7
p0 = x0 + np.array([30, 0.05, 0.03, 0.02, 10, 0.02, 0.1]) * rng.standard_normal((nw, nd))
p0[:, 3] = np.abs(p0[:, 3])
sam = emcee.EnsembleSampler(nw, nd, lnprob)
sam.run_mcmc(p0, 6000, progress=False)
tau = sam.get_autocorr_time(tol=0)
ch = sam.get_chain(discard=2000, thin=10, flat=True)
T, lg, mh, av, d, R, lsj = ch.T
th = theta_mas(R, d)
fbol = np.array([M.fbol(a, b, c, e) for a, b, c, e in zip(T[::5], lg[::5], mh[::5], th[::5])])  # erg/s/cm2
L = 4 * np.pi * (d[::5] * 3.0857e18) ** 2 * fbol / 3.828e33
Tsyn = np.array([M.mags(a, b, c, e_, f)[iT] for a, b, c, e_, f in zip(T[::20], lg[::20], mh[::20], av[::20], th[::20])])


def summ(x):
    q = np.percentile(x, [15.865, 50, 84.135])
    return dict(med=float(q[1]), lo=float(q[1] - q[0]), hi=float(q[2] - q[1]))


best = ch[np.argmax(sam.get_log_prob(discard=2000, thin=10, flat=True))]
mbest = M.mags(*best[:4], theta_mas(best[5], best[4]))[idx]
chi2 = float(np.sum((mobs - mbest) ** 2 / sig ** 2))
res["sed_fit"] = dict(Teff=summ(T), logg=summ(lg), mh=summ(mh), A_V=summ(av), d_pc=summ(d), R=summ(R),
                      jitter_mag=summ(np.exp(lsj)), L=summ(L), Fbol_cgs=summ(fbol), T_synthetic=summ(Tsyn),
                      chi2_best_no_jitter=chi2, n_bands=len(bands), tau_max=float(np.max(tau)), n_samples=len(ch),
                      residuals_best={b: float(o - m) for b, o, m in zip(bands, mobs, mbest)},
                      corr_logT_logR=float(np.corrcoef(np.log10(T), np.log10(R))[0, 1]),
                      priors=dict(mh=[MH_MU, MH_SD], AV_max=AV_MAX, plx=[plx, plx_err]))
print("SED", json.dumps({k: v for k, v in res["sed_fit"].items() if k != "residuals_best"}, default=float), flush=True)

# ------------------------------------------------------------------ 4. MIST grid inference
z = np.load("/home/user/work/tic87/stellar/cache/mist_compact.npz", allow_pickle=True)
tab = z["tab"]; cols = list(z["cols"])
C = {c: tab[:, i] for i, c in enumerate(cols)}
mu = np.array([np.mean(np.log10(T)), np.mean(np.log10(R))])
cov = np.cov(np.vstack([np.log10(T), np.log10(R)]))
# model systematics: 0.005 dex in log Teff (~70 K) and 0.01 dex in log R added in quadrature
cov_s = cov + np.diag([0.005 ** 2, 0.010 ** 2])
icov = np.linalg.inv(cov_s)
dv = np.vstack([C["logT"] - mu[0], C["logR"] - mu[1]])
chi2 = np.einsum("in,ij,jn->n", dv, icov, dv)
feh_mu, feh_sd = float(np.mean(mh)), float(np.hypot(np.std(mh), 0.10))
chi2 += ((C["feh"] - feh_mu) / feh_sd) ** 2
# prior weights
w = np.zeros(len(tab))
key = np.round(C["feh"], 2) * 1000 + np.round(C["logage"], 2)
for kk in np.unique(key):
    m = np.where(key == kk)[0]
    m = m[np.argsort(C["eep"][m])]
    mi = C["mini"][m]
    dm = np.gradient(mi) if len(mi) > 1 else np.ones(1)
    imf = np.where(mi < 0.5, mi ** -1.3, 0.5 ** -1.3 * (mi / 0.5) ** -2.3)
    w[m] = np.abs(dm) * imf
w *= 10 ** C["logage"]  # uniform in age for uniform logage spacing
w *= np.exp(-0.5 * (chi2 - chi2.min()))
w /= w.sum()
dens = C["mass"] / (10 ** C["logR"]) ** 3


def wq(x, q=(0.15865, 0.5, 0.84135)):
    o = np.argsort(x); cw = np.cumsum(w[o])
    v = np.interp(q, cw, x[o])
    return dict(med=float(v[1]), lo=float(v[1] - v[0]), hi=float(v[2] - v[1]))


neff = float(1 / np.sum(w ** 2))
res["mist"] = dict(M=wq(C["mass"]), R=wq(10 ** C["logR"]), Teff=wq(10 ** C["logT"]), logg=wq(C["logg"]),
                   age_Gyr=wq(10 ** (C["logage"] - 9)), feh=wq(C["feh"]), rho_sun=wq(dens), L=wq(10 ** C["logL"]),
                   frac_MS=float(np.sum(w[C["phase"] == 0])), frac_post_MS=float(np.sum(w[C["phase"] == 2])),
                   n_eff=neff, inputs=dict(logT=[float(mu[0]), float(np.sqrt(cov_s[0, 0]))],
                                           logR=[float(mu[1]), float(np.sqrt(cov_s[1, 1]))], feh=[feh_mu, feh_sd]),
                   prior="Kroupa IMF x uniform age x uniform [Fe/H] grid nodes; MIST v1.2 vvcrit0.4, 0.5-3 Msun")
# rho* from the SED radius and the MIST mass (posterior samples drawn from the weights) -- equivalent cross-check
sel = rng.choice(len(w), size=20000, p=w)
res["mist"]["rho_from_samples"] = summ(dens[sel])
res["mist"]["M_samples_pct"] = [float(v) for v in np.percentile(C["mass"][sel], [2.275, 15.865, 50, 84.135, 97.725])]
# posterior samples for downstream scripts (companion radius / mass): SED samples and MIST-weighted draws
np.savez(os.path.join(ROOT, "data/stellar_samples.npz"), sed_T=T, sed_R=R, sed_logg=lg, sed_mh=mh, sed_AV=av, sed_d=d,
         mist_M=C["mass"][sel], mist_R=10 ** C["logR"][sel], mist_rho=dens[sel], mist_age=10 ** (C["logage"][sel] - 9),
         mist_logg=C["logg"][sel], mist_T=10 ** C["logT"][sel])
print("MIST", json.dumps(res["mist"], default=float), flush=True)

# ------------------------------------------------------------------ adopted values
fr = json.load(open(os.path.join(ROOT, "results/fluxratio.json")))
rv = [fr["A1_stassun_gaiadr3"]["r"], fr["A2_tic82"]["r"], fr["A3_xp_sampled"]["TESS"]["r"]]
r_ad = float(np.mean(rv)); r_err = float(max(np.std(rv, ddof=1), fr["A1_stassun_gaiadr3"]["r_err"]))
Rm, Mm, rhom = res["mist"]["R"], res["mist"]["M"], res["mist"]["rho_from_samples"]
res["adopted"] = dict(
    Teff=res["sed_fit"]["Teff"], R=res["sed_fit"]["R"], M=Mm, rho_sun=rhom, logg=res["mist"]["logg"],
    age_Gyr=res["mist"]["age_Gyr"], d_pc=res["sed_fit"]["d_pc"], A_V=res["sed_fit"]["A_V"], L=res["sed_fit"]["L"],
    mh=res["sed_fit"]["mh"], r_TESS=r_ad, r_TESS_err=r_err, D_pair=1 / r_ad, D_pair_err=r_err / r_ad ** 2,
    note="R, Teff, L from the SED+parallax fit; M, rho, age, logg from the MIST grid fed with (Teff, R, [M/H])")
json.dump(res, open(os.path.join(ROOT, "results/stellar.json"), "w"), indent=1, default=float)

# ------------------------------------------------------------------ figures
lam = np.array([sedlib.FILT[b]["lpiv"] for b in bands]) / 1e4
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].errorbar(lam, mobs, yerr=sig, fmt="o", label="observed")
ax[0].plot(lam, mbest, "s", mfc="none", label="best CK04 model")
ax[0].invert_yaxis(); ax[0].set_xscale("log"); ax[0].set_xlabel("pivot wavelength (um)"); ax[0].set_ylabel("mag")
ax[0].set_title(f"TIC 87055315 SED: Teff {res['sed_fit']['Teff']['med']:.0f} K, R {res['sed_fit']['R']['med']:.2f} Rsun")
ax[0].legend(fontsize=8)
sub = sel[:3000]
ax[1].scatter(10 ** C["logT"][sub], C["logg"][sub], s=2, c=C["mass"][sub], cmap="viridis")
for k in np.unique(np.round(C["logage"], 2)):
    if k not in (9.5, 9.7, 9.8, 9.9, 10.0): continue
    m = (np.round(C["logage"], 2) == k) & (np.abs(C["feh"] - 0.0) < 0.01)
    o = np.argsort(C["eep"][m])
    ax[1].plot(10 ** C["logT"][m][o], C["logg"][m][o], lw=0.7, label=f"MIST [Fe/H]=0, {10 ** (k - 9):.1f} Gyr")
ax[1].set_xlim(7000, 4800); ax[1].set_ylim(4.7, 3.3); ax[1].set_xlabel("Teff (K)"); ax[1].set_ylabel("log g")
ax[1].legend(fontsize=7); ax[1].set_title("MIST posterior samples (colour = mass)")
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "fig/sed_hrd_nb.png"), dpi=110)
print("DONE")
