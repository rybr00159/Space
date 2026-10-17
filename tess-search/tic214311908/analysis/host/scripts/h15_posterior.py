"""h15: posterior odds 'metal-rich single' (S) vs 'unresolved near-equal companion' (T) for the M3V, and the
consequences for the planet radius.
Model. [Fe/H] = f on a grid; S: CMD offsets (dMKs, dMG) ~ N(a + b f, Sigma_int + Sigma_meas) with (a, b, scatter) of
the Gaussian core of the calibrator fits (h13, RUWE < 1.1) and the dMKs-dMG residual correlation of the calibrators;
T(q): the same plus the binary offsets (dK(q), dG(q)) computed at fixed combined BP-RP with the Mamajek main sequence
(q = 1 -> -0.753 in both bands). Luminosity-independent metallicity indicators with their calibrated forward relations
(value = A0 f + A1, scatter): XP-shape regression (h12b), colour residuals r(G-W1) and r(G-Ks) at fixed BP-RP (h13;
the measured target error is added), Gaia Net (h16; found uninformative for M dwarfs); the wide
companions' blend colour (h14) as a likelihood in f. Prior on f: GSP-Spec dwarfs with the M3V's kinematics (h10b).
Prior on T: multiplicity rate 0.268 (Winters+19); separation log-normal (Winters+19: mu = log10 25 au, sigma 1.16;
alternative Susemiehl & Meyer 22: 48 au, 0.97); q-distribution p(q) ~ q^0.5 on 0.3-1 (El-Badry+19, 0.1-0.4 Msun)
plus a twin excess F_twin = 0.25 (0.15-0.35) uniform in 0.95-1; survival of T under the data that do not involve
[Fe/H]: Gaia RUWE + DR3-DR2 PM (h09, per q and a), Gaia RV error (h07, per a), Gaia IPD (phase 1: L = 1 at
a < 5 au, 0.05 at 5-160 au, 0 beyond (resolved)).
Output: out/h15_posterior.json, figs/h15_posterior.png"""
import json, numpy as np, warnings
warnings.filterwarnings('ignore')
from scipy.stats import norm
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
B = 'tic214311908/analysis/host'
S_ = 'tic214311908/analysis/star'
J = lambda f: json.load(open(f'{B}/out/{f}'))
h13 = J('h13_colour_feh.json'); h12 = J('h12b_xp.json'); h16 = J('h16_lit_xp_feh.json'); h14 = J('h14_wide_pair_feh.json')
h10 = J('h10b_feh_prior.json'); h09 = J('h09_twin_astrometry.json'); h07 = J('h07_gaia_rv.json')
fg = np.round(np.arange(-1.2, 0.81, 0.005), 4)
# ---- Mamajek MS for binary offsets
tab = []
for line in open(f'{S_}/data/ext/EEM_dwarf.txt'):
    p = line.split()
    if len(p) > 30 and p[0][0] == 'M' and p[0].endswith('V'):
        try: tab.append((float(p[30]), float(p[11]), float(p[13]), float(p[22])))
        except Exception: pass
tab = np.array(sorted(tab))           # mass, BP-RP, M_G, M_Ks
mass = np.linspace(0.09, 0.6, 2000)
c_ = np.interp(mass, tab[:, 0], tab[:, 1]); g_ = np.interp(mass, tab[:, 0], tab[:, 2]); k_ = np.interp(mass, tab[:, 0], tab[:, 3])
# BP and RP from G and BP-RP using G-RP(BP-RP) from the table
def ms_at(m):
    return np.interp(m, mass, c_), np.interp(m, mass, g_), np.interp(m, mass, k_)
ridgeK = lambda col: np.interp(col, c_[::-1], k_[::-1]); ridgeG = lambda col: np.interp(col, c_[::-1], g_[::-1])
def binary_offsets(q):
    best = None
    for m1 in np.linspace(0.2, 0.45, 501):
        c1, G1, K1 = ms_at(m1); c2, G2, K2 = ms_at(q * m1)
        # BP, RP via G and colour assuming BP-G ~ 0.55*(BP-RP)+ ... use flux-weighted colour approximation in G
        Gc = -2.5 * np.log10(10 ** (-0.4 * G1) + 10 ** (-0.4 * G2)); Kc = -2.5 * np.log10(10 ** (-0.4 * K1) + 10 ** (-0.4 * K2))
        # combined BP-RP: weight colours by G flux in BP and RP separately (BP-G and G-RP split ~ 0.56/0.44 of BP-RP at these colours)
        BP1, RP1 = G1 + 0.562 * c1, G1 - 0.438 * c1; BP2, RP2 = G2 + 0.562 * c2, G2 - 0.438 * c2
        BPc = -2.5 * np.log10(10 ** (-0.4 * BP1) + 10 ** (-0.4 * BP2)); RPc = -2.5 * np.log10(10 ** (-0.4 * RP1) + 10 ** (-0.4 * RP2))
        cc = BPc - RPc
        if best is None or abs(cc - 2.756) < best[0]:
            best = (abs(cc - 2.756), m1, Kc - ridgeK(cc), Gc - ridgeG(cc), cc, K2 - K1)
    return dict(M1=best[1], dK=best[2], dG=best[3], comb_col=best[4], dK_components=best[5])
qg = np.array(h09['q'])
boff = {float(q): binary_offsets(q) for q in qg}
print('binary offsets:', {q: (round(v['dK'], 3), round(v['dG'], 3)) for q, v in boff.items()})
# ---- CMD likelihood
fK = h13['fits']['dMK_ruwe1.1']; fG = h13['fits']['dMG_ruwe1.1']
rho = h13.get('core_resid_corr_KG', 0.8)
tK, tG = h13['target']['dMK'], h13['target']['dMG']
sK = np.hypot(fK['scatter'], 0.035); sG = np.hypot(fG['scatter'], 0.02)     # + measurement and ridge-zero-point
CMD2D = False; SK_OVERRIDE = None
def L_cmd(f, dK=0.0, dG=0.0):
    mK = fK['a'] + fK['b'] * f + dK; mG = fG['a'] + fG['b'] * f + dG
    if not CMD2D:   # 1D in dMKs (the dMG - dMKs difference is ~ the G-Ks colour residual, used separately as an indicator)
        s1 = np.hypot(SK_OVERRIDE, 0.035) if SK_OVERRIDE else np.hypot(fK['scatter'], 0.035)
        return norm.pdf(tK, mK, s1)
    x = (tK - mK) / sK; y = (tG - mG) / sG
    return np.exp(-0.5 * (x ** 2 - 2 * rho * x * y + y ** 2) / (1 - rho ** 2)) / (2 * np.pi * sK * sG * np.sqrt(1 - rho ** 2))
# ---- indicator likelihoods (forward-calibrated: value = A0 f + A1, scatter s)
ind = {}
CAL = 'main'
def set_cal(which):
    global fK, fG, CAL
    CAL = which
    src = h13['fits'] if which == 'main' else h13['fits_highprec']
    fK = src['dMK_ruwe1.1'] if which == 'main' else src['dMK']; fG = src['dMG_ruwe1.1'] if which == 'main' else src['dMG']
    ind.clear()
    xp = h12['xp_feh']; ind['XP_shape'] = (xp['target_pred'], xp['cv_pred_vs_true_slope'], xp['cv_pred_vs_true_intercept'], np.hypot(xp['cv_robust_scatter_local_colour'], xp['target_pred_mc_sd']))
    g = src['rGK_ruwe1.1'] if which == 'main' else src['rGK']; ind['rGK'] = (h13['target']['rGK'], g['b'], g['a'], np.hypot(g['scatter'], 0.027 * np.sqrt(2)))
    g = h13['fits']['rGW1_ruwe1.1']; ind['rGW1'] = (h13['target']['rGW1'], g['b'], g['a'], np.hypot(g['scatter'], 0.03 * np.sqrt(2)))
    gn = h16['gaianet']; ind['GaiaNet'] = (gn['target_value'], gn['slope_pred_vs_true'], gn['intercept'], gn['local_robust_sd'] or gn['robust_sd'])
set_cal('main')
def L_ind(name, f):
    v, a0, a1, s = ind[name]; return norm.pdf(v, a0 * f + a1, s)
wp = np.array([[r['feh'], r['L']] for r in h14['grid']])
L_wide = lambda f: np.interp(f, wp[:, 0], wp[:, 1] / wp[:, 1].max())
# ---- [Fe/H] prior
def gmm_prior(key):
    g = h10[key]['gmm']; p = sum(w * norm.pdf(fg, m, s) for w, m, s in zip(g['w'], g['mu'], g['sd'])); return p / np.trapezoid(p, fg)
# ---- twin survival per q (astrometry x RV x IPD), separation prior
ag = np.array(h09['a_au']); LJ = np.array(h09['L_joint'])
rv = h07['simulation']['sw15']; a_rv = np.array([r['a_au'] for r in rv]); L_rv = np.array([r['L_ratio_floor0p5'] for r in rv])
def survival(q, mu, sig):
    la = np.linspace(-2.0, 4.0, 1201); a = 10 ** la
    pa = norm.pdf(la, mu, sig)
    iq = np.argmin(np.abs(qg - q))
    Lastro = np.interp(la, np.log10(ag), LJ[iq], left=LJ[iq][0], right=LJ[iq][-1])
    Lastro = np.where(a > ag[-1], 1.0, Lastro)
    Lr = np.interp(la, np.log10(a_rv), L_rv, left=L_rv[0], right=1.0)
    Lipd = np.where(a < 5.0, 1.0, np.where(a < 160, 0.05, 0.0))
    return float(np.trapezoid(pa * Lastro * Lr * Lipd, la)), float(np.trapezoid(pa * (a < 5.0), la))
def q_prior(gamma=0.5, Ftw=0.25, qmin=0.3):
    qq = np.linspace(qmin, 1, 1401); p = (1 - Ftw) * qq ** gamma / np.trapezoid(qq ** gamma, qq) + Ftw * ((qq >= 0.95) / 0.05)
    return qq, p
def odds(prior_key='analogues_absW_30_65', inds=('XP_shape', 'rGW1'), wide=True, mu=np.log10(25), sig=1.16, Ftw=0.25, MR=0.268, gamma=0.5, verbose=False):
    pf = gmm_prior(prior_key)
    Lind = np.ones_like(fg)
    for n in inds: Lind *= L_ind(n, fg)
    if wide: Lind *= L_wide(fg)
    # single
    ES = np.trapezoid(pf * Lind * L_cmd(fg), fg)
    qq, pq = q_prior(gamma, Ftw)
    # bin q into the h09 grid cells (0.8-1.0); below 0.8 the binary offset is too small to matter: counted as 'single-like'
    edges = np.r_[0.775, (qg[1:] + qg[:-1]) / 2, 1.0001]
    ET = 0.0; detail = []
    for i, q in enumerate(qg):
        w = np.trapezoid(pq * ((qq >= edges[i]) & (qq < edges[i + 1])), qq)
        surv, fa5 = survival(q, mu, sig)
        Lf = np.trapezoid(pf * Lind * L_cmd(fg, boff[float(q)]['dK'], boff[float(q)]['dG']), fg)
        ET += MR * w * surv * Lf
        detail.append(dict(q=float(q), w_q=float(w), survival=surv, frac_a_lt5=fa5, ev=float(MR * w * surv * Lf)))
    PS = 1 - MR * np.trapezoid(pq * (qq >= 0.775), qq)
    O = ET / (PS * ES)
    # posterior on f under each hypothesis
    postS = pf * Lind * L_cmd(fg); postS /= np.trapezoid(postS, fg)
    postT = pf * Lind * L_cmd(fg, boff[1.0]['dK'], boff[1.0]['dG']); postT /= np.trapezoid(postT, fg)
    return dict(odds_T_over_S=float(O), P_T=float(O / (1 + O)), feh_post_S_mean=float(np.trapezoid(fg * postS, fg)), feh_post_T_mean=float(np.trapezoid(fg * postT, fg)),
                detail=detail if verbose else None)
out = dict(binary_offsets=boff, indicators={k: dict(value=v[0], slope=v[1], intercept=v[2], sigma=v[3]) for k, v in ind.items()}, rv_used='sw15 L_ratio_floor0p5',
           cmd=dict(target_dMK=tK, target_dMG=tG, fit_K=fK, fit_G=fG, rho=rho))
runs = {}
runs['baseline (CMD dMKs + XP + r(G-W1) + wide pair; analogue prior |W| 30-65; Winters sep; Ftwin .25)'] = odds(verbose=True)
SK_OVERRIDE = 0.15
runs['baseline with intrinsic dMKs scatter 0.15 mag (label noise removed)'] = odds()
runs['CMD only, intrinsic scatter 0.15'] = odds(inds=(), wide=False)
SK_OVERRIDE = None
CMD2D = True
runs['2D CMD (dMKs+dMG, implicitly G-Ks) + XP + wide'] = odds(inds=('XP_shape',))
CMD2D = False
runs['CMD only (no metallicity indicators)'] = odds(inds=(), wide=False)
runs['CMD + XP shape'] = odds(inds=('XP_shape',), wide=False)
runs['CMD + r(G-W1)'] = odds(inds=('rGW1',), wide=False)
runs['CMD + r(G-Ks)'] = odds(inds=('rGK',), wide=False)
runs['CMD + wide pair'] = odds(inds=(), wide=True)
runs['CMD + Gaia Net (uninformative check)'] = odds(inds=('GaiaNet',), wide=False)
runs['CMD + XP + r(G-Ks) + wide (2MASS colour instead of WISE)'] = odds(inds=('XP_shape', 'rGK'))
runs['CMD + XP + wide (no photometric colour)'] = odds(inds=('XP_shape',))
runs['baseline, prior all local dwarfs'] = odds(prior_key='all_dwarfs')
runs['baseline, prior tight analogues (n=82)'] = odds(prior_key='analogues_W_vtot')
runs['baseline, Susemiehl&Meyer separations'] = odds(mu=1.68, sig=0.97)
runs['baseline, Ftwin 0.15'] = odds(Ftw=0.15)
runs['baseline, Ftwin 0.35'] = odds(Ftw=0.35)
runs['baseline, MR x1.5 (triple-system enhancement)'] = odds(MR=0.268 * 1.5)
runs['baseline, prior all dwarfs, MR x1.5, Ftwin .35, S&M sep (twin-friendly)'] = odds(prior_key='all_dwarfs', MR=0.268 * 1.5, Ftw=0.35, mu=1.68, sig=0.97)
set_cal('highprec')
runs['high-precision (non-Birky) calibration of dMK/dMG/rGK: CMD + XP + r(G-Ks) + wide'] = odds(inds=('XP_shape', 'rGK'))
runs['high-precision calibration: CMD only'] = odds(inds=(), wide=False)
set_cal('main')
for k, v in runs.items():
    print(f'{k:70s} odds T:S = {v["odds_T_over_S"]:.3g}  P(T) = {v["P_T"]:.3f}  <[Fe/H]>|S = {v["feh_post_S_mean"]:+.2f}  <[Fe/H]>|T = {v["feh_post_T_mean"]:+.2f}')
out['runs'] = runs
# ---- planet radius under each hypothesis (depth 15.7 +- 0.8 +- 1.0 ppt on the combined light; b ~ 0.5)
depth = 15.7e-3; ed = np.hypot(0.8, 1.0) * 1e-3
ld = 1.04   # limb-darkening boost of depth for b ~ 0.5 (quadratic LD 0.2/0.38): Rp/R* = sqrt(depth/LDfac)... approx factor on k^2
def rp(Rs, dilution):
    k2 = depth * dilution / ld * 1.0
    return np.sqrt(k2) * Rs * 109.08, np.sqrt(k2) * Rs * 109.08 * 0.5 * ed / depth
res_rp = {}
res_rp['single_solar_R0.381'] = rp(0.381, 1.0)
res_rp['single_metalrich_R0.385'] = rp(0.385, 1.0)
for q in (1.0, 0.95, 0.9):
    b_ = boff[q]; m1 = b_['M1']
    # TESS-band flux ratio ~ G-band ratio from MS; radii from Mann+15 M_Ks relation
    c1, G1, K1 = ms_at(m1); c2, G2, K2 = ms_at(q * m1)
    fr = 10 ** (-0.4 * (0.5 * (G2 - G1) + 0.5 * (K2 - K1)))   # T band between G and Ks
    def R_mann(MK): return 1.9515 - 0.3520 * MK + 0.01680 * MK ** 2
    MKobs = 6.4406
    MK1 = MKobs + 2.5 * np.log10(1 + 10 ** (-0.4 * (K2 - K1))); MK2 = MK1 + (K2 - K1)
    R1, R2 = R_mann(MK1), R_mann(MK2)
    res_rp[f'twin_q{q}_host_primary'] = rp(R1, 1 + fr)
    res_rp[f'twin_q{q}_host_secondary'] = rp(R2, 1 + 1 / fr)
out['planet_radius_Rearth'] = {k: dict(Rp=float(v[0]), err_depth_only=float(v[1])) for k, v in res_rp.items()}
print('Rp:', {k: round(v[0], 2) for k, v in res_rp.items()})
json.dump(out, open(f'{B}/out/h15_posterior.json', 'w'), indent=1, default=float)
