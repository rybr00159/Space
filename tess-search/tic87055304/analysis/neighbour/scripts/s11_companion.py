"""Nature of the companion eclipsing TIC 87055315.

Inputs: neighbour-host transit MCMC on the PSF pair light curve (mcmc/mcmc_pair_local.npz: k, b, rho_transit, a/R*,
T14, D free with prior N(5.32, 0.06)); stellar posterior samples (data/stellar_samples.npz: SED radius, MIST mass and
density); pair-based secondary / odd-even / phase-curve limits (results/pair_tests.json); psf_photometry neighbour-only
secondary (0.30 +- 1.87 ppt) and odd/even; Gaia RV test (results/gaia_rv.json); twin-EB test (results/eb2p.json).
1. rho_transit vs rho_star (log-space z and posterior overlap).
2. Companion radius R2 = k R* (independent draws), in R_sun and R_Jup (R_Jup = 71492 km). Systematic variant: only
   84 % of the decrement on the neighbour (psf_photometry f_target = 0.16): k scales by ~sqrt(0.84).
3. Secondary-eclipse brightness-temperature limit: delta_occ = k^2 F_TESS(T2) / F_TESS,star with a blackbody companion
   and the CK04 surface flux of the star (sedlib grid, A_V = 0); 2-sigma limit -> T_b,max and L2,max = R2^2 (T/5772)^4.
   Reflected light is ignored (it would only lower the thermal limit).
4. Expected signals vs companion mass (M2 = 1 M_Jup ... 0.5 M_sun): occultation depth (M dwarfs: Teff from a simple
   main-sequence relation, blackbody), ellipsoidal A_ell = alpha (M2/M1)(R1/a)^3 sin^2 i with alpha = 0.15 (15+u)(1+tau)
   /(3-u), u = 0.55, tau = 0.40 (alpha ~ 1.33), beaming A_beam = alpha_b 4 K/c with alpha_b = 1.0, and the RV
   semi-amplitude K.  M-dwarf mass from R2 with the Boyajian+12 relation R = 0.32 M^2 + 0.6063 M + 0.0906 (+5 % for
   inflation seen in EBLM systems as an alternative).
Output: results/companion.json
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
ROOT = "/home/user/work/tic87/neighbour"
rng = np.random.default_rng(5)
R_J = 71492e3 / 6.957e8  # R_Jup in R_sun
M_J = 1.89813e27 / 1.98892e30
G = 6.67430e-11; MSUN = 1.98892e30; RSUN = 6.957e8; c = 2.99792458e8
z = np.load(os.path.join(ROOT, "mcmc/mcmc_pair_local.npz"))
ch = z["chain"]; names = list(z["names"])
col = lambda n: ch[:, names.index(n)]
k, b, rho_t, P = col("k"), col("b"), col("rho"), np.median(col("P"))
a_R = (rho_t * 1.40894e3 * G * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)  # rho_sun = 1.40894 g/cm3 -> kg/m3 * 1e3
inc = np.arccos(b / a_R)
ss = np.load(os.path.join(ROOT, "data/stellar_samples.npz"))
N = 20000
iR = rng.integers(0, len(ss["sed_R"]), N); iM = rng.integers(0, len(ss["mist_M"]), N); iC = rng.integers(0, len(k), N)
Rs, Ms, rhos = ss["sed_R"][iR], ss["mist_M"][iM], ss["mist_rho"][iM]
Ts = ss["sed_T"][iR]
kk, bb, aa, ii, rt = k[iC], b[iC], a_R[iC], inc[iC], rho_t[iC]
res = {}


def summ(x):
    q = np.percentile(x, [2.275, 15.865, 50, 84.135, 97.725])
    return dict(med=float(q[2]), lo=float(q[2] - q[1]), hi=float(q[3] - q[2]), p2=float(q[0]), p97=float(q[4]))


# 1. density
lt, ls = np.log(rho_t), np.log(ss["mist_rho"])
zrho = (np.median(lt) - np.median(ls)) / np.hypot(np.std(lt), np.std(ls))
# fraction of the transit-density posterior inside the stellar 95 % interval
lo, hi = np.percentile(ss["mist_rho"], [2.5, 97.5])
res["density"] = dict(rho_transit=summ(rho_t), rho_star_MIST=summ(ss["mist_rho"]), z_log=float(zrho),
                      ratio_transit_over_star=summ(rt / rhos), frac_transit_in_star95=float(np.mean((rho_t > lo) & (rho_t < hi))))
# 2. radius
R2 = kk * Rs
res["radius"] = dict(k=summ(k), b=summ(b), aR=summ(a_R), inc_deg=summ(np.degrees(inc)), R_star=summ(Rs),
                     R2_Rsun=summ(R2), R2_RJup=summ(R2 / R_J),
                     R2_RJup_if_84pct_on_nb=summ(R2 * np.sqrt(0.84) / R_J),
                     R2_Rsun_if_84pct_on_nb=summ(R2 * np.sqrt(0.84)))
# first-round cross-check (transit/mcmc_nbhost.json: QLP SAP, D = 4.35 fixed, S13-S93)
fr = json.load(open("/home/user/work/tic87/transit/mcmc_nbhost.json"))["params"]
k1 = rng.normal(fr["k"]["med"], 0.5 * (fr["k"]["lo"] + fr["k"]["hi"]), N)
res["radius"]["first_round_nbhost"] = dict(k=fr["k"], b=fr["b"], rho=fr["rho"], T14_h=fr["T14_h"],
                                           R2_RJup_with_this_Rstar=summ(k1 * Rs / R_J))
print(json.dumps(res, indent=1, default=float)[:3000], flush=True)

# 3. secondary limit
import sedlib
grid = np.load(sedlib.GRID, allow_pickle=True)
bands = [str(x) for x in grid["bands"]]
SM = sedlib.SEDModel(bands)
iT = bands.index("TESS")
f = sedlib.FILT["TESS"]
wl = f["w"]; tr = f["t"]
h, kB = 6.62607015e-34, 1.380649e-23


def F_bb(T):  # photon-weighted band mean of pi B_lambda, erg/s/cm2/A
    lam = wl * 1e-10
    B = 2 * h * c ** 2 / lam ** 5 / np.expm1(h * c / (lam * kB * T))  # W/m2/m/sr
    Fl = np.pi * B * 1e-7  # W/m2/m -> erg/s/cm2/A (1e7 erg/J / 1e4 cm2/m2 / 1e10 A/m)
    return np.trapezoid(Fl * tr * wl, wl) / np.trapezoid(tr * wl, wl)


Tst, lgst, mhst = np.median(ss["sed_T"]), np.median(ss["mist_logg"]), np.median(ss["sed_mh"])
F_star = 10 ** SM.interp([[Tst, lgst, mhst, 0.0]])[0][iT]
F_star_bb = F_bb(Tst)
pt = json.load(open(os.path.join(ROOT, "results/pair_tests.json")))
psf = json.load(open("/home/user/work/tic87/psf_photometry/lc/nominal/secondary_oddeven.json"))
k2 = np.median(k) ** 2
Tgrid = np.arange(500, 7001, 10.0)
dgrid = np.array([k2 * F_bb(T) / F_star for T in Tgrid])


def Tb_of(delta):
    return float(np.interp(delta, dgrid, Tgrid)) if delta > dgrid[0] else float(Tgrid[0])


lim_pair = pt["secondary_nb"]["upper2sig_ppt"] * 1e-3
lim_psf = (max(0.30, 0) + 2 * 1.87) * 1e-3
res["secondary"] = dict(pair_based_ppt=[pt["secondary_nb"]["depth_ppt"], pt["secondary_nb"]["err_ppt"]],
                        psf_nb_only_ppt=[0.30, 1.87], upper2sig_pair_ppt=lim_pair * 1e3, upper2sig_psf_ppt=lim_psf * 1e3,
                        Tb_max_pair_K=Tb_of(lim_pair), Tb_max_psf_K=Tb_of(lim_psf),
                        star_surface_flux_ratio_CK04_over_BB=float(F_star / F_star_bb),
                        R2_Rsun_med=float(np.median(R2)))
for key in ("pair", "psf"):
    T_ = res["secondary"][f"Tb_max_{key}_K"]
    res["secondary"][f"L2_max_{key}_Lsun"] = float(np.median(R2) ** 2 * (T_ / 5772) ** 4)
res["secondary"]["scan_any_phase_upper2sig_ppt"] = pt["secondary_scan"]["upper2sig_any_ppt"]
res["secondary"]["Tb_max_any_phase_K"] = Tb_of(pt["secondary_scan"]["upper2sig_any_ppt"] * 1e-3)
# equilibrium temperature of a non-luminous companion (zero albedo, full redistribution)
Teq = Ts * np.sqrt(1 / (2 * aa))
res["secondary"]["Teq_companion_K"] = summ(Teq)
print("secondary", res["secondary"], flush=True)


# 4. expectations vs mass
def ms_teff(M):  # crude M-dwarf / K-dwarf main-sequence Teff(M) (Pecaut & Mamajek 2013 table, interpolated)
    Mt = np.array([0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
    Tt = np.array([2650, 2900, 3100, 3200, 3300, 3400, 3550, 3750, 4000, 4300, 4700])
    return np.interp(M, Mt, Tt)


def R_boy(M):
    return 0.32 * M ** 2 + 0.6063 * M + 0.0906


M1 = np.median(Ms); R1 = np.median(Rs); aR_med = np.median(a_R); si2 = np.median(np.sin(inc) ** 2)
u_ld, tau_gd = 0.55, 0.40
alpha = 0.15 * (15 + u_ld) * (1 + tau_gd) / (3 - u_ld)
rows = []
for M2 in [1 * M_J, 5 * M_J, 13 * M_J, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]:
    K = (2 * np.pi * G / (P * 86400)) ** (1 / 3) * M2 * MSUN * np.sqrt(si2) / ((M1 + M2) * MSUN) ** (2 / 3)
    # a/R from Kepler with the actual masses (vs the fitted a/R; the companion mass changes a only slightly)
    A_ell = alpha * (M2 / M1) * (1 / aR_med) ** 3 * si2
    A_beam = 1.0 * 4 * K / c
    T2 = ms_teff(M2) if M2 >= 0.08 else float(np.median(Teq))
    occ = k2 * F_bb(T2) / F_star
    rows.append(dict(M2_Msun=float(M2), M2_MJup=float(M2 / M_J), K_kms=float(K / 1e3), A_ell_ppt=float(A_ell * 1e3),
                     A_beam_ppt=float(A_beam * 1e3), T2_K=float(T2), occultation_ppt=float(occ * 1e3),
                     R_ms_Boyajian12_Rsun=float(R_boy(M2)) if M2 >= 0.08 else None))
res["expectations_vs_mass"] = dict(alpha_ell=alpha, rows=rows, M1=float(M1), R1=float(R1), aR=float(aR_med))
pc = pt["phase_curve"]["fit"]
Alim = pc["A_ell"]["upper2sig_ppt_nb"] * 1e-3
q_lim = Alim / (alpha * (1 / aR_med) ** 3 * si2)
Blim = pc["A_beam"]["upper2sig_ppt_nb"] * 1e-3
K_lim = Blim * c / 4 / 1.0
res["mass_limits"] = dict(A_ell_meas_ppt=[pc["A_ell"]["ppt_nb"], pc["A_ell"]["err_ppt_nb"]],
                          A_beam_meas_ppt=[pc["A_beam"]["ppt_nb"], pc["A_beam"]["err_ppt_nb"]],
                          M2_max_from_ellipsoidal_Msun=float(q_lim * M1), K_max_from_beaming_kms=float(K_lim / 1e3))
# M-dwarf mass implied by the radius
Mg = np.linspace(0.08, 0.7, 2000)
Mfrom = lambda R: np.interp(R, R_boy(Mg), Mg)
res["mass_from_radius"] = dict(M2_Boyajian12=summ(Mfrom(R2)), M2_inflated5pct=summ(Mfrom(R2 / 1.05)),
                               note="valid only if the companion is a main-sequence M dwarf")
Mmd = np.median(Mfrom(R2))
Kmd = (2 * np.pi * G / (P * 86400)) ** (1 / 3) * Mmd * MSUN * np.sqrt(si2) / ((M1 + Mmd) * MSUN) ** (2 / 3)
res["mass_from_radius"]["K_expected_kms"] = float(Kmd / 1e3)
res["mass_from_radius"]["occultation_expected_ppt"] = float(k2 * F_bb(ms_teff(Mmd)) / F_star * 1e3)
res["mass_from_radius"]["A_ell_expected_ppt"] = float(alpha * (Mmd / M1) * (1 / aR_med) ** 3 * si2 * 1e3)
res["mass_from_radius"]["A_beam_expected_ppt"] = float(4 * Kmd / c * 1e3)
try:
    res["gaia_rv"] = json.load(open(os.path.join(ROOT, "results/gaia_rv.json")))
except Exception as e:
    res["gaia_rv"] = f"n/a: {e}"
try:
    res["eb2p"] = json.load(open(os.path.join(ROOT, "results/eb2p.json")))
except Exception as e:
    res["eb2p"] = f"n/a: {e}"
res["oddeven"] = dict(pair=pt["oddeven_nb"], psf_nb_only=dict(odd=psf["combined"]["odd_n"], even=psf["combined"]["even_n"],
                                                       diff=psf["combined"]["odd_minus_even_n"]))
json.dump(res, open(os.path.join(ROOT, "results/companion.json"), "w"), indent=1, default=float)
print(json.dumps({k_: v for k_, v in res.items() if k_ in ("mass_limits", "mass_from_radius")}, indent=1, default=float))
print("DONE")
