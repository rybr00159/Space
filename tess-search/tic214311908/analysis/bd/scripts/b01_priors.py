"""b01: literature priors for the brown-dwarf (BD) channels.
(a) Table of transiting BDs (13-80 M_J) around M dwarfs (M* <= 0.65 Msun) with mass, radius, period, eccentricity and
    host mass. Values are taken from the discovery-paper abstracts fetched from arXiv on 2026-10-17 (bibcodes from the
    journal/DOI metadata via Crossref) or, where the abstract gives none, from the Carmichael catalogue page
    (theroncarmichael.com/browndwarfs, fetched 2026-10-17; 2023MNRAS.519.5177C) and the uniform RV re-fit table of
    Ferreira et al. (2026, arXiv:2607.02682, Table A1). Every row lists its source.
(b) Occurrence of BDs at P = 3-10 d around early/mid M dwarfs: f_BD = f_giant-size(TESS, M dwarfs) x (BD share among
    confirmed transiting giant-size companions of M dwarfs at P < 10 d, weighted by 1/p_transit ~ P^(2/3)) x (share
    of the BD rate at 3-10 d, same weighting). Giant planets: NASA Exoplanet Archive pscomppars (queried 2026-10-17,
    data/nea/nea_mdwarf_giants.csv; Teff < 4000 K, 0.08 < M* < 0.65, R > 0.55 R_J, M < 13 M_J, P < 10 d).
(c) Eccentricity distributions: Ferreira+26 short-period (P < 16 d) and long-period BD Beta laws (with their quoted
    hyper-uncertainties), Kipping (2013) and Van Eylen+19 planet laws, and a hierarchical Beta fit to the M-dwarf BDs
    with 2.5 < P < 16 d of table (a). Gives P(e > x) and the probability of the configurations a grazing BD needs.
(d) Tidal circularisation (equilibrium tide, constant Q', Jackson et al. 2008, 2008ApJ...678.1396J, their eq. 2):
    1/tau_e = [63/4 sqrt(G M*^3) R_bd^5/(Q'_bd M_bd) + 171/16 sqrt(G/M*) R*^5 M_bd/Q'_*] a^(-13/2)
    for M* 0.370, R* 0.383, P 5.6915 d, M_bd 13-80 M_J, R_bd 0.8-1.1 R_J, Q'_bd 1e4-1e9, Q'_* 1e5-1e8.
Outputs: out/b01_mdwarf_bds.csv, out/b01_priors.json"""
import sys; sys.path.insert(0, 'tic214311908/analysis/bd/scripts')
from blib import *
import pandas as pd
from scipy import stats
from scipy.special import betaln
rng = np.random.default_rng(1)
T = [
 # name, M_bd(MJ), sM, R_bd(RJ), sR, P(d), e, e_note, Mstar, Teff, age_note, grazing, source
 ('NLTT 41135 b', 33.7, 2.8, 1.13, 0.27, 2.889, 0.0, '<0.02 (circular fit)', 0.19, 3230, 'field', True, '2010ApJ...718.1353I (abstract: 31-34+-3 MJ, circular, grazing); R, M from Carmichael cat.; M* Ferreira+26 Tab. A1'),
 ('LP 261-75 C (TOI-1779)', 68.1, 2.1, 0.90, 0.02, 1.882, 0.0, '<0.01', 0.30, 3100, 'young/field (debated)', False, '2018AJ....156..140I; values Carmichael cat.; obliquity arXiv:2411.10402'),
 ('AD 3116 b', 54.6, 6.9, 0.95, 0.07, 1.983, 0.14, '0.14+-0.05 (cat.); 0.133+-0.025 (Ferreira refit)', 0.28, 3165, 'Praesepe (sub-Gyr)', False, '2017ApJ...849...11G; values Carmichael cat.'),
 ('NGTS-7A b', 75.5, 13.7, 1.38, 0.14, 0.676, 0.0, '0 (fixed)', 0.48, 3359, '~55 Myr', False, '2019MNRAS.489.5146J (abstract: 16.2 h, ~55 Myr); values Carmichael cat.'),
 ('TOI-263 b', 61.6, 4.0, 0.91, 0.07, 0.557, 0.0, '<0.02 (cat.); 0.150+-0.016 (Ferreira refit)', 0.44, 3471, 'field', False, '2021A&A...650A..55P; values Carmichael cat.; M* Ferreira+26'),
 ('TOI-2119 b', 64.4, 2.3, 1.08, 0.03, 7.2009, 0.337, '0.337+-0.002 (Carmichael+22), 0.3362+-0.0005 (Canas+22); Ferreira refit 0.885 (rejected: contradicts secondary-eclipse timing)', 0.53, 3621, 'field', False, '2022MNRAS.514.4944C, 2022AJ....163...89C (abstracts; TESS secondary: T_b 2100+-80 K / Teff 2030+-84 K)'),
 ('TOI-5375 B', 77.0, 8.0, 0.99, 0.16, 1.7216, 0.02, '0.021+-0.012 (Ferreira refit)', 0.64, 3885, 'active', False, '2023A&A...674A.132M (abstract), 2023AJ....165..218L'),
 ('NGTS-28A b', 69.0, 5.0, 0.95, 0.05, 1.254, 0.10, '0.095+-0.098 (Ferreira refit)', 0.56, 3626, '>0.5 Gyr', False, '2024MNRAS.530..318H (arXiv:2402.09943 abstract); M* Ferreira+26'),
 ('TOI-5610 b', 40.4, 1.0, 0.887, 0.031, 7.95346, 0.354, '0.354+-0.012', 0.53, 3618, '~1.5 Gyr (Sonora)', False, '2025AJ....169..246L (arXiv:2501.16554 abstract; M* from paper text)'),
 ('TOI-5389A b', 68.0, 2.2, 0.824, 0.032, 10.40046, 0.096, '0.096 (+0.003/-0.005)', 0.43, 3569, '>~2 Gyr', False, '2025AJ....169..246L (arXiv:2501.16554 abstract; M* from paper text)'),
 ('LHS 6343 C', 62.1, 1.2, 0.783, 0.011, 12.71, 0.056, 'e cos w = 0.0228 (Montet+15); 0.056 (Johnson+11); 0.088 (Ferreira refit)', 0.37, 3130, '>1 Gyr', False, '2011ApJ...730...79J, 2015ApJ...800..134M (abstracts); M* Ferreira+26 (Johnson+11 quote 0.45/0.36 for the AB pair)'),
 ('TOI-1278 B', 18.5, 0.5, 1.09, 0.24, 14.476, 0.013, '0.013 (slight but non-zero)', 0.55, 3799, 'field', True, '2021AJ....162..144A (abstract: grazing transit); R, e Carmichael cat./NEA'),
 ('TOI-6508 b', 72.5, 6.4, 1.03, 0.03, 19.0, 0.28, '0.28 (+0.09/-0.08)', 0.174, 2930, 'field', False, 'arXiv:2502.19940 (abstract)'),
 ('TOI-5575 b', 72.4, 4.1, 0.84, 0.07, 32.0, 0.187, '0.187+-0.002', 0.21, None, 'field', False, '2025ApJ...988L..78G (arXiv:2507.09461 abstract)'),
 ('ZTF J2020+5033 b', 80.1, 1.6, 1.05, 0.02, 0.0792, 0.0, '0', None, 2856, 'old (thick disk)', False, 'arXiv:2307.15729 (abstract; OJAp 2023); at the H-burning limit'),
]
cols = ['name', 'M_MJ', 'sM_MJ', 'R_RJ', 'sR_RJ', 'P_d', 'e', 'e_note', 'Mstar', 'Teff', 'age_note', 'grazing', 'source']
tab = pd.DataFrame(T, columns=cols)
tab.to_csv(f'{BD}/out/b01_mdwarf_bds.csv', index=False)
out = {'table_file': 'out/b01_mdwarf_bds.csv'}
# ---------------- (b) occurrence
nea = pd.read_csv(f'{BD}/data/nea/nea_mdwarf_giants.csv')
gp = nea[(nea.st_teff < 4000) & (nea.pl_bmassj < 13) & (nea.pl_orbper < 10) & (nea.st_mass > 0.08)]
bdP10 = tab[(tab.P_d < 10) & (tab.M_MJ < 79) & (tab.Mstar.fillna(0.15) <= 0.65)]
w = lambda Pd: np.asarray(Pd, float)**(2/3.)
WG, WB = w(gp.pl_orbper).sum(), w(bdP10.P_d).sum()
WB310 = w(bdP10.P_d[(bdP10.P_d >= 3) & (bdP10.P_d < 10)]).sum()
share = WB/(WB + WG)
share310 = WB310/WB
f_giant = dict(Bryant23_all=(0.194e-2, 0.072e-2), Bryant23_026_042=(0.108e-2, 0.083e-2), Gan23_early=(0.27e-2, 0.09e-2), Glusman25_mid=(0.153e-2, 0.069e-2))
fg = 0.194e-2
# Monte Carlo of the estimate: Poisson counts for BDs (n_BD at P<10, n at 3-10), giant rate uncertainty (Bryant all)
nB, nB310, nG = len(bdP10), int(((bdP10.P_d >= 3) & (bdP10.P_d < 10)).sum()), len(gp)
draws = []
for i in range(20000):
    lamB = rng.gamma(nB + 0.5); lamG = rng.gamma(nG + 0.5); lam310 = rng.gamma(nB310 + 0.5)
    sh = lamB*WB/nB/(lamB*WB/nB + lamG*WG/nG)
    s310 = min(1, lam310*(WB310/max(nB310, 1))/(lamB*WB/nB))
    fgi = max(1e-5, rng.normal(fg, 0.072e-2))
    draws.append(fgi*sh*s310)
draws = np.array(draws)
LNP = np.log(10/3.)
out['occurrence'] = dict(
    n_giant_planets_P_lt_10=int(nG), n_BD_P_lt_10=int(nB), n_BD_P_3_10=int(nB310), BD_names_P_lt_10=list(bdP10.name),
    BD_share_weighted=float(share), share_of_BD_rate_at_3_10d=float(share310),
    f_giant_literature=f_giant,
    f_BD_3_10d_per_star=dict(median=float(np.median(draws)), p16=float(np.percentile(draws, 16)), p84=float(np.percentile(draws, 84)),
                             p05=float(np.percentile(draws, 5)), p95=float(np.percentile(draws, 95))),
    f_BD_per_lnP_at_5p7d=float(np.median(draws)/LNP),
    FGK_comparison=dict(Csizmadia15='~0.2% at P < 10 d (2015A&A...584A..13C), ~6x below hot Jupiters',
                        Santerne16='0.29 +- 0.17% within 400 d (2016A&A...587A..64S)', Sahlmann11='< 0.6% (2011A&A...525A..95S)'),
    referee2_fpp_value='1e-3 per star at 3-10 d (own order-of-magnitude choice)',
    note='Transit-probability weighting only; detection/confirmation completeness assumed equal for BDs and giant planets. '
         'The Bryant+23 giant-size rate is a candidate-based rate for 0.6-2 R_J objects (BDs fall in it), so the product '
         'is approximately the BD rate. Uncertainty dominated by the 2 BDs known at 3-10 d.')
print('occurrence', json.dumps(out['occurrence']['f_BD_3_10d_per_star']), 'share', share, share310, nG, nB, nB310, flush=True)
# ---------------- (c) eccentricity distributions
def beta_draws(a, sa, b, sb, n):
    A = np.clip(rng.normal(a, sa, n), 0.05, None); Bb = np.clip(rng.normal(b, sb, n), 0.1, None)
    return rng.beta(A, Bb)
N = 400000
E = {}
E['ferreira_SP'] = beta_draws(0.426, 0.092, 2.354, 0.73, N)
E['ferreira_LP'] = beta_draws(1.879, 0.92, 2.470, 1.25, N)
E['kipping13'] = rng.beta(0.867, 3.03, N)
E['vaneylen19_halfnormal'] = np.abs(rng.normal(0, 0.32, N))
# hierarchical Beta fit to the M-dwarf BDs with 2.5 < P < 16 d (upper limits -> 0.01)
sel = tab[(tab.P_d > 2.5) & (tab.P_d < 16) & (tab.M_MJ < 79)]
ev = np.clip(sel.e.values, 0.005, 0.99)
ga, gb = np.meshgrid(np.logspace(-1.2, 1.3, 120), np.logspace(-1.0, 1.6, 120))
ll = np.sum([(ga - 1)*np.log(x) + (gb - 1)*np.log(1 - x) - betaln(ga, gb) for x in ev], axis=0)
lp = ll - 0.0          # log-uniform hyperpriors on alpha, beta (grid is log-spaced)
pw = np.exp(lp - lp.max()); pw /= pw.sum()
idx = rng.choice(pw.size, size=4000, p=pw.ravel())
Ed = np.concatenate([rng.beta(ga.ravel()[i], gb.ravel()[i], 100) for i in idx])
E['mdwarf_BD_hier'] = Ed
out['ecc_mdwarf_sample'] = dict(names=list(sel.name), e=list(map(float, sel.e)),
                                post_mean_alpha=float(np.sum(pw*ga)), post_mean_beta=float(np.sum(pw*gb)))
for k in list(E): E[k] = np.clip(E[k], 0, 0.95)
out['ecc'] = {k: dict(mean=float(v.mean()), P_gt_0p3=float(np.mean(v > 0.3)), P_gt_0p4=float(np.mean(v > 0.4)), P_gt_0p5=float(np.mean(v > 0.5)),
                      P_gt_0p6=float(np.mean(v > 0.6)), P_gt_0p7=float(np.mean(v > 0.7))) for k, v in E.items()}
np.savez_compressed(f'{BD}/data/ecc_draws.npz', **{k: v[:200000].astype(np.float32) for k, v in E.items()})
for k, v in out['ecc'].items(): print('ecc', k, v, flush=True)
# probability of rho_circ/rho_true <= x (needed by a grazing BD) incl. transit-probability weighting
wv = rng.uniform(0, 2*np.pi, N)
out['rho_ratio_tail'] = {}
for k, ee in E.items():
    ee = ee[:N]; g, geo = g_ecc(ee, wv)
    r = g**3
    wt = geo/np.mean(geo)
    out['rho_ratio_tail'][k] = {f'P(rho_c/rho <= {x})': float(np.mean(wt*(r <= x))) for x in (0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7)}
    print('rho tail', k, out['rho_ratio_tail'][k], flush=True)
# ---------------- (d) tides
Ms = 0.370*MSUN; Rs = 0.383*RSUN; a = (G*(0.370 + 0.05)*MSUN*(P*DAY)**2/(4*np.pi**2))**(1/3.)
tide = []
for Mbd in (13, 30, 50, 65, 80):
    for Rbd in (0.8, 1.0, 1.1):
        for lQb in (4, 5, 6, 7, 8, 9):
            for lQs in (5, 6, 7, 8):
                Mb = Mbd*MJUP; Rb = Rbd*RJUP
                r1 = 63/4*np.sqrt(G*Ms**3)*Rb**5/(10**lQb*Mb)*a**-6.5
                r2 = 171/16*np.sqrt(G/Ms)*Rs**5*Mb/(10**lQs)*a**-6.5
                tau = 1/(r1 + r2)/YR/1e9
                tide.append(dict(M_bd=Mbd, R_bd=Rbd, logQbd=lQb, logQs=lQs, tau_e_Gyr=float(tau), tau_bd_only_Gyr=float(1/r1/YR/1e9),
                                 tau_star_only_Gyr=float(1/r2/YR/1e9), e_surv_frac_1Gyr=float(np.exp(-1/tau)), e_surv_frac_5Gyr=float(np.exp(-5/tau)),
                                 e_surv_frac_10Gyr=float(np.exp(-10/tau))))
tdf = pd.DataFrame(tide)
out['tides'] = dict(a_au=float(a/AU), formula='Jackson+2008 eq. (2), small-e limit; e(t) = e0 exp(-t/tau)',
                    ferreira26_Q=dict(Q_bd_only='10^(8.1+-1.0)', Q_bd_with_star='10^(7.1+-0.3)', Q_star='10^(6.0+-0.1)'),
                    table=tide,
                    summary={f'M{m}_R1.0_Qbd7.1_Qs6': float(tdf[(tdf.M_bd == m) & (tdf.R_bd == 1.0) & (tdf.logQbd == 7) & (tdf.logQs == 6)].tau_e_Gyr.iloc[0]) for m in (13, 30, 50, 65, 80)})
print('tides summary (Gyr)', out['tides']['summary'], flush=True)
print(tdf[(tdf.R_bd == 1.0) & (tdf.M_bd == 65)].to_string(), flush=True)
jdump(dict(method=__doc__, **out), f'{BD}/out/b01_priors.json')
