"""f04: explicit Bayesian scenario comparison for the 5.69-d signal on TIC 214311908, with sensitivity analysis.

Evidence of each scenario S:  E_S = (prior probability per target star that S exists, per unit ln P at P = 5.69 d)
                                     x (likelihood of all data under S, integrated over S's internal parameters).
All E_S use the same period normalisation (per ln P) and the same light-curve reference chi2 (best planet fit on the grid),
so FPP = 1 - E(planet on the M3V or its twin) / sum E, NFPP = E(resolved neighbours) / sum E.
Three FPP definitions are reported:
  FPP        (task definition): planet on the M3V (single or diluted) or on a near-equal companion (Delta T <= dT_twin = 0.75)
  FPP_strict (TRICERATOPS-like): only a planet on the brighter component counts; a planet on the twin is an FP
  FPP_system : any planet-sized body transiting any star of the unresolved bound system counts (PTP radii are 5.7-6.1 R_e)

PLANET scenarios (TP single; TP diluted by a companion; TP-twin = planet on a near-equal companion; PTP = planet on a
fainter bound companion = false positive), from the f02 grid:
   Z = < (R/a)_true (1+e sin w)/(1-e^2) * sum_b db sqrt(2 pi) sigma_k g(k_hat R_h) / k_hat exp(-(chi2 - chi2_ref)/2) >_{e,w}
   evaluated at rho_circ = rho_host * [(1+e sin w)/sqrt(1-e^2)]^-3 (photo-eccentric effect), with g = planet occurrence per
   ln R per ln P for the host (classes: 2-4, 4-7, 7-10, 10-22 R_earth; f_class integrated over P = 3-10 d, uniform in ln R
   and ln P inside each class; host-mass scaling (M/0.37)^beta_class). e prior: Kipping (2013) beta(0.867, 3.03).
   Binary configurations from f03 (Winters+19 population; IPD, astrometry, DECaPS2, CMD likelihoods), planet suppression
   S(a) (Moe & Kratter 2021 / Kraus+16: S_close for a <= 10 au, log-linear to 1 at 200 au) and Holman & Wiegert (1999)
   stability for the S-type orbit.
EB scenarios use the maximum-likelihood Delta chi2 of the physical models fitted by the shape agent (shape/out/summary.json;
   relative to the free planet fit, including secondary eclipse, photometry, parallax and RUWE terms where stated),
   times an Occam factor f_occ (fraction of the class's parameter space that reproduces depth and duration; 0.1, generous
   compared with the planet's own ~0.05), times literature rates and eclipse probabilities.
BACKGROUND / NEARBY scenarios add the injection-calibrated localisation likelihoods (f01) and archival-imaging limits.

Outputs: out/fpp_baseline.json, out/fpp_variants.csv, out/fpp_random.csv, out/fpp_summary.json"""
import numpy as np, json, sys, itertools, pandas as pd
from scipy.interpolate import RegularGridInterpolator
sys.path.insert(0, 'tic214311908/analysis/fpp/scripts')
from fpplib import *
rng = np.random.default_rng(7)
LNP = np.log(10/3.)
CLASSES = {'sub': (2.0, 4.0), 'N': (4.0, 7.0), 'S': (7.0, 10.0), 'J': (10.0, 22.0)}
def aR_rho(rho): return (rho*G*(P*DAY)**2/(3*np.pi))**(1/3.)

# ------------------------------------------------------------------ light-curve grid -> J tables
def load_grid(ds):
    z = np.load(f'{F}/out/lc_grid_{ds}.npz')
    C = z['chi2']; c = C - np.nanmin(C)
    return dict(D=z['D'], rho=z['rho'], b=z['b'], c=c, k=z['k'], sk=np.nan_to_num(z['sigk']), chi2min=float(np.nanmin(C)))
RH = np.logspace(np.log10(0.09), np.log10(0.75), 34)        # host radius grid (R_sun); 0.75 covers pair B's R_eff
def J_tables(gr):
    """J_class[iD, irho_c, iR] = sum_b db sqrt(2pi) sk * [k R in class] / k * exp(-c/2) / (ln(Rhi/Rlo) LNP)  (unit f_class)"""
    db = gr['b'][1] - gr['b'][0]
    L = np.exp(-0.5*np.nan_to_num(gr['c'], nan=1e9))*np.sqrt(2*np.pi)*gr['sk']/gr['k']*db
    out = {}
    for cl, (lo, hi) in CLASSES.items():
        Rp = gr['k'][..., None]*RH[None, None, None, :]*RSUN/REARTH
        ind = (Rp >= lo) & (Rp < hi)
        out[cl] = np.nansum(L[..., None]*ind, axis=2)/(np.log(hi/lo)*LNP)    # sum over b
    return out
def ecc_samples(kind, n=3000):
    if kind == 'kipping':
        e = rng.beta(0.867, 3.03, n)
    elif kind == 'vaneylen':
        e = np.abs(rng.normal(0, 0.32, n))
    else:
        e = np.zeros(n)
    e = np.clip(e, 0, 0.9); w = rng.uniform(0, 2*np.pi, n)
    ge = (1 + e*np.sin(w))/np.sqrt(1 - e**2)
    geo = (1 + e*np.sin(w))/(1 - e**2)
    return ge, geo
RHOH = np.logspace(np.log10(3.0), np.log10(200.), 44)
def Zbar_tables(gr, J, kind):
    """Zbar_class[iD, irho_host, iR] = < (R/a) geo J_class(D, rho_host ge^-3, R) >_e"""
    ge, geo = ecc_samples(kind)
    lr = np.log(gr['rho'])
    out = {}
    for cl in CLASSES:
        T = np.zeros((len(gr['D']), len(RHOH), len(RH)))
        lJ = np.log(np.clip(J[cl], 1e-300, None))
        for ih, rh in enumerate(RHOH):
            rc = np.log(rh*ge**-3)
            ok = (rc >= lr[0]) & (rc <= lr[-1])
            ra = 1/aR_rho(rh)
            for iD in range(len(gr['D'])):
                for iR in range(len(RH)):
                    v = np.interp(rc[ok], lr, lJ[iD, :, iR])
                    T[iD, ih, iR] = ra*np.sum(geo[ok]*np.exp(v))/len(ge)
        out[cl] = T
    return out
def interp_Z(gr, Zt, cl, D, rho_h, R_h):
    """trilinear interpolation of log Zbar in (log D, log rho_h, log R_h); D < min grid D -> 0"""
    f = RegularGridInterpolator((np.log(gr['D'][::-1]), np.log(RHOH), np.log(RH)),
                                np.log(np.clip(Zt[cl][::-1], 1e-300, None)), bounds_error=False, fill_value=None)
    D = np.atleast_1d(D).astype(float); rho_h = np.atleast_1d(rho_h).astype(float); R_h = np.atleast_1d(R_h).astype(float)
    pts = np.column_stack([np.log(np.clip(D, gr['D'].min(), 1.0)), np.log(np.clip(rho_h, RHOH[0], RHOH[-1])),
                           np.log(np.clip(R_h, RH[0], RH[-1]))])
    v = np.exp(f(pts))
    v[D < gr['D'].min()*0.999] = 0.0
    return v

PRE = {}
for ds in ('inv', 'ref'):
    gr = load_grid(ds); J = J_tables(gr)
    PRE[ds] = dict(gr=gr, J=J)
    for kind in ('kipping', 'vaneylen', 'circ'):
        PRE[ds][kind] = Zbar_tables(gr, J, kind)
    print('tables done', ds, flush=True)

# ------------------------------------------------------------------ companions (f03)
cz = np.load(f'{F}/out/companions_mc.npz')
CP = {k: cz[k] for k in cz.files}
NC = len(CP['a'])
def hw_crit(mu, e):   # Holman & Wiegert (1999) S-type critical semi-major axis / a_bin
    return 0.464 - 0.380*mu - 0.631*e + 0.586*mu*e + 0.150*e**2 - 0.198*mu*e**2
ap1 = (G*CP['m1']*MSUN*(P*DAY)**2/(4*np.pi**2))**(1/3.)/AU
ap2 = (G*CP['m2']*MSUN*(P*DAY)**2/(4*np.pi**2))**(1/3.)/AU
stab1 = (ap1 < hw_crit(CP['m2']/(CP['m1'] + CP['m2']), CP['e'])*CP['a']).astype(float)
stab2 = (ap2 < hw_crit(CP['m1']/(CP['m1'] + CP['m2']), CP['e'])*CP['a']).astype(float)
ZC = {}
for ds in ('inv', 'ref'):
    gr = PRE[ds]['gr']
    for kind in ('kipping', 'vaneylen', 'circ'):
        Zt = PRE[ds][kind]
        for cl in CLASSES:
            ZC[(ds, kind, cl, 1)] = interp_Z(gr, Zt, cl, 1 - CP['Dc'], CP['rho1'], CP['R1'])
            ZC[(ds, kind, cl, 2)] = interp_Z(gr, Zt, cl, CP['Dc'], CP['rho2'], CP['R2'])
print('companion Z interpolated', flush=True)
# single-star TP, twin reference values and neighbours
RHO_TP = rng.normal(9.42, 0.88, 400)
def Z_single(ds, kind, cl, rho=RHO_TP, R=0.381, D=1.0):
    return float(np.mean(interp_Z(PRE[ds]['gr'], PRE[ds][kind], cl, np.full(len(np.atleast_1d(rho)), D), rho,
                                  np.full(len(np.atleast_1d(rho)), R))))
def mann_R(MKs): return 1.9515 - 0.3520*MKs + 0.01680*MKs**2
# pair A (TIC 214311910; G 14.56, M_G 10.61 -> ~0.30 Msun) and pair B (Gaia ...253632, M_G 11.70 -> ~0.20 Msun)
PAIR = {'A': dict(M=0.30, R=0.31, T=13.375), 'B': dict(M=0.20, R=0.22, T=15.273)}
for k, v in PAIR.items():
    v['rho'] = v['M']*MSUN/(4/3*np.pi*(v['R']*RSUN)**3)
    v['F'] = 979*10**(-0.4*(v['T'] - M3V['Tmag']))
    v['R_eff'] = v['R']*np.sqrt(979/v['F'])     # radius at which the D = 1 grid reproduces the planet radius on the neighbour
ZN = {}
for ds in ('inv', 'ref'):
    for kind in ('kipping', 'vaneylen', 'circ'):
        for cl in CLASSES:
            ZN[(ds, kind, cl, 'TP')] = Z_single(ds, kind, cl)
            for k, v in PAIR.items():
                ZN[(ds, kind, cl, 'pair' + k)] = Z_single(ds, kind, cl, rho=np.full(50, v['rho']), R=v['R_eff'])

# sanity checks against the shape agent: light-curve evidence of the single-star vs the equal-twin host (N class only)
CHK = {}
for ds in ('inv', 'ref'):
    gr = PRE[ds]['gr']
    zs = Z_single(ds, 'kipping', 'N'); zt = Z_single(ds, 'kipping', 'N', rho=rng.normal(15.45, 1.45, 400), R=0.289, D=0.5)
    zsc = Z_single(ds, 'circ', 'N'); ztc = Z_single(ds, 'circ', 'N', rho=rng.normal(15.45, 1.45, 400), R=0.289, D=0.5)
    i1 = 0; j1 = int(np.argmin(np.abs(np.log(gr['rho']/9.42)))); i5 = int(np.argmin(np.abs(gr['D'] - 0.5))); j5 = int(np.argmin(np.abs(np.log(gr['rho']/15.45))))
    CHK[ds] = dict(Z_single_kipping=zs, Z_twin_kipping=zt, lnB_single_twin_kipping=float(np.log(zs/zt)),
                   lnB_single_twin_circ=float(np.log(zsc/ztc)),
                   dchi2_min_single_rho9p4=float(np.nanmin(gr['c'][i1, j1])), dchi2_min_twin_rho15p5=float(np.nanmin(gr['c'][i5, j5])),
                   best_b_single=float(gr['b'][np.nanargmin(gr['c'][i1, j1])]), k_single=float(gr['k'][i1, j1][np.nanargmin(gr['c'][i1, j1])]))
    print('CHECK', ds, CHK[ds], flush=True)
LOC = jload(f'{F}/out/loc_likelihoods.json')
SH = {r['scenario'] + '_' + r['dataset']: r for r in jload('shape/out/summary.json')['scenarios']}
BD = jload('star/out/blend_density.json')['gaia_density_per_arcsec2']

# ------------------------------------------------------------------ parameters
BASE = dict(
    ds='inv', ecc='kipping',
    f_sub=0.15, f_N=0.003, f_S=0.0005, f_J=0.0015,          # per star, P = 3-10 d (ctx/README synthesis; Gan+23, Bryant+23)
    beta_sub=0.0, beta_N=0.0, beta_S=1.0, beta_J=1.5,          # host-mass scaling (M/0.37)^beta of each class
    f_mult=0.27, mu_loga=1.30, sig_loga=1.16, gamma_q=0.0, F_twin=0.10,   # Winters+19; twin excess (q > 0.95)
    S_close=0.25, sig_s=0.30, ipd='LR', ast='base', dec=True, dT_twin=0.75,
    r_close=0.024, p_ecl_ms=0.05, dchi2_bdms=45.0,          # Shan+15 0.11 close companions (P 1-90 d) per M dwarf
    r_wd=2e-5, f_Twd=0.3, p_ecl_wd=0.029, dchi2_wd=13.2, dchi2_wd_twin=11.1, f_occ=0.1,
    dchi2_heb=27.5, p_ecl_heb=0.06,
    r_eb_field=0.0017, f_compact=0.01, dchi2_bg3_ms=65.5, dchi2_bg3_c=19.6, dchi2_bg5_ms=76.3, dchi2_bg5_c=27.4,
    lr_bg3=10**LOC['catalogued']['G175']['log10LR_t'], lr_bg5=10**LOC['catalogued']['G177']['log10LR_t'],
    lr_pairA=10**LOC['catalogued']['pairA']['log10LR_t'], lr_pairB=10**LOC['catalogued']['pairB']['log10LR_t'],
    lr_mira=10**LOC['catalogued']['Mira']['log10LR_t'],
    kappa=2.2, eps_dss=0.02, f_Mdw_bg=0.2, dchi2_unc_ms=30.0,
    S_pair=None)
# dataset-dependent physical-model Delta chi2 (shape/out/summary.json): WD occulted by one twin component 17.9 (INV) /
# 11.1 (REF); bound hierarchical EB 27.5 / 33.7; BD or M dwarf eclipsing the M3V 45 (INV k-profile) / 84 (REF grid).
DCHI2_DS = {'inv': dict(dchi2_wd_twin=SH['twin_inv']['dchi2_total'], dchi2_heb=27.5, dchi2_bdms=45.0),
            'ref': dict(dchi2_wd_twin=SH['twin_ref']['dchi2_total'], dchi2_heb=33.7, dchi2_bdms=84.0)}
for k in ('dchi2_wd_twin', 'dchi2_heb', 'dchi2_bdms'): BASE[k] = None
# hypothetical follow-up (what-if) switches
BASE.update(hri=None, chroma=None, sb2=False)

def a_weight(par):
    lg = CP['loga']
    p_new = np.exp(-0.5*((lg - par['mu_loga'])/par['sig_loga'])**2)/par['sig_loga']
    p_old = np.exp(-0.5*((lg - 1.30)/1.16)**2)/1.16
    return p_new/p_old
def q_weight(par):
    # drawn uniform in [qmin, 1]; target p(q) ~ (1 - F_twin) q^gamma / norm + F_twin U(0.95, 1)
    q = CP['q']; qmin = 0.075/CP['m1']; g = par['gamma_q']
    norm = (1 - qmin**(g + 1))/(g + 1) if abs(g + 1) > 1e-9 else -np.log(qmin)
    p_new = (1 - par['F_twin'])*q**g/norm + par['F_twin']*(q >= 0.95)/0.05
    p_old = 1/(1 - qmin)
    return p_new/p_old
def S_of_a(a, Sc):
    s = np.where(a <= 10, Sc, np.where(a >= 200, 1.0, Sc + (1 - Sc)*np.log(a/10)/np.log(20)))
    return s
def occ(par, cl, M):
    f = par['f_' + cl]; beta = par['beta_' + cl]
    return f*(np.asarray(M)/0.37)**beta

# companion light fractions in V and Ic (Mamajek) for the chromatic what-if: depth(V)/depth(I) for a planet on the
# companion is Dc_V/Dc_I, on the primary (1-Dc_V)/(1-Dc_I); an on-target planet gives 1.
sys.path.insert(0, f'{B}/shape/scripts'); import eem
MV1 = eem.ms(CP['m1'], 'Mv'); MV2 = eem.ms(CP['m2'], 'Mv')
MI1 = MV1 - eem.ms(CP['m1'], 'V-Ic'); MI2 = MV2 - eem.ms(CP['m2'], 'V-Ic')
DcV = 1/(1 + 10**(0.4*(MV2 - MV1))); DcI = 1/(1 + 10**(0.4*(MI2 - MI1)))
RATIO = {2: DcV/DcI, 1: (1 - DcV)/(1 - DcI)}
def chroma_L(par, host):
    if par['chroma'] is None: return 1.0
    return np.exp(-0.5*((RATIO[host] - 1)/par['chroma'])**2)

def evaluate(par, detail=False):
    ds, kind = par['ds'], par['ecc']
    par = dict(par)
    for k, v in DCHI2_DS[ds].items():
        if par.get(k) is None: par[k] = v
    # configuration weights
    Lc = np.ones(NC)
    Lc *= {'LR': CP['L_ipd'], 'LR84': CP['L_ipd84'], 'off': 1.0}[par['ipd']]
    Lc *= {'base': CP['L_ast'], 'soft': CP['L_ast']**0.4, 'off': 1.0}[par['ast']]
    if par['dec']: Lc = Lc*CP['L_dec']
    if par['hri'] is not None:      # hypothetical speckle/AO null detection: companions brighter than the curve -> 0.01
        hs, hv = np.array(par['hri'][0]), np.array(par['hri'][1])
        limv = np.where(CP['sep'] >= hs[0], np.interp(CP['sep'], hs, hv), -99)
        Lc = Lc*np.where(CP['dI'] < limv, 0.01, 1.0)
    if par['sb2']:                  # hypothetical high-resolution spectra: no SB2 / RV variation: companions with
        Lc = Lc*np.where((CP['dT'] < 2.0) & (CP['a'] < 3.0), 0.02, 1.0)   # dT < 2 inside 3 au excluded
    s = par['sig_s']
    pe = lambda x: np.exp(-0.5*(x/s)**2)/(np.sqrt(2*np.pi)*s)
    Lcmd = pe(CP['eps_needed']); Lcmd_single = float(pe(-0.753))
    w = q_weight(par)*a_weight(par)
    wl = w*Lc*Lcmd
    fm = par['f_mult']
    Zcfg_single = (1 - fm)*Lcmd_single
    Zcfg_bin = fm*np.mean(wl)
    Zcfg = Zcfg_single + Zcfg_bin
    Sa = S_of_a(CP['a'], par['S_close'])
    E = {}
    # planets
    Z1 = sum(occ(par, cl, CP['m1'])*ZC[(ds, kind, cl, 1)] for cl in CLASSES)
    Z2 = sum(occ(par, cl, CP['m2'])*ZC[(ds, kind, cl, 2)] for cl in CLASSES)
    Ztp = sum(occ(par, cl, 0.37)*ZN[(ds, kind, cl, 'TP')] for cl in CLASSES)
    E['TP_single'] = Zcfg_single*Ztp
    E['TP_diluted'] = fm*np.mean(wl*Sa*stab1*Z1)
    twin = CP['dT'] <= par['dT_twin']
    E['TP_on_twin'] = fm*np.mean(wl*Sa*stab2*Z2*twin)
    E['PTP_companion'] = fm*np.mean(wl*Sa*stab2*Z2*(~twin)*chroma_L(par, 2))
    if par['chroma'] is not None:   # chromatic penalty also applies to the diluted primary-hosted planet (small)
        E['TP_diluted'] = fm*np.mean(wl*Sa*stab1*Z1*chroma_L(par, 1))
    # EBs on the M3V system
    ex = lambda d: np.exp(-0.5*d)
    nt = CP['dT'] > 0.3
    E['EB_WD_on_M3V'] = (Zcfg_single + fm*np.mean(wl*nt))*par['r_wd']*par['f_Twd']*par['p_ecl_wd']*par['f_occ']*ex(par['dchi2_wd'])
    E['EB_WD_on_twin'] = fm*np.mean(wl*(~nt))*2*par['r_wd']*par['f_Twd']*0.033*par['f_occ']*ex(par['dchi2_wd_twin'])
    E['EB_MS_BD_on_M3V'] = Zcfg*par['r_close']*par['p_ecl_ms']*par['f_occ']*ex(par['dchi2_bdms'])
    E['HEB_bound'] = fm*np.mean(wl*(CP['dT'] < 4.46))*par['r_close']*par['p_ecl_heb']*par['f_occ']*ex(par['dchi2_heb'])
    # catalogued background stars (3.1" and 5.5"): MS+MS and compact-object EB classes; a planet cannot give 54-76 %
    rc = par['r_eb_field']
    E['BEB_bg3_3.1as'] = Zcfg*rc*((1 - par['f_compact'])*ex(par['dchi2_bg3_ms']) + par['f_compact']*ex(par['dchi2_bg3_c']))*par['f_occ']*par['lr_bg3']
    E['BEB_bg5_5.5as'] = Zcfg*rc*((1 - par['f_compact'])*ex(par['dchi2_bg5_ms']) + par['f_compact']*ex(par['dchi2_bg5_c']))*par['f_occ']*par['lr_bg5']
    # uncatalogued background stars: density x effective localisation area x probability of escaping the archival
    # imaging (DSS look-back for unbound stars within 3", DECaPS2 5-sigma limits beyond 1", Gaia/IPD for bright ones)
    Aeff = LOC['uncatalogued'][f'kappa_{par["kappa"]}']['A_eff_by_annulus']
    eps = par['eps_dss']
    miss_faint = {'0-1': eps, '1-1.25': eps*0.5, '1.25-3': eps*0.05, '3-5': 0.02*0.1, '5-8': 0.02*0.1, '8-15': 0.0, '15-40': 0.0}
    miss_bright = {k: v*0.1 for k, v in miss_faint.items()}      # T <= 15.8 also flagged by Gaia IPD (Delta I < 3)
    nb_faint = sum((BD['T<=17.47'] - BD['T<=16'])*Aeff[k]*miss_faint[k] for k in Aeff)
    nb_bright = sum(BD['T<=16']*Aeff[k]*miss_bright[k] for k in Aeff)
    # BTP: Jupiter on a background M dwarf (T <= 16, depth <= 25 %); BEB compact (M dwarf + hot WD, LC dchi2 = 0);
    # BEB MS+MS (penalised by the shape, dchi2 >= 30 even for M+M)
    rate_hj = occ(par, 'J', 0.4)/LNP*0.05*0.3     # per ln P, transit prob ~0.05, depth/size match 0.3
    E['BTP_uncatalogued'] = Zcfg*nb_bright*par['f_Mdw_bg']*rate_hj
    E['BEB_uncat_compact'] = Zcfg*(nb_faint + nb_bright)*par['f_Mdw_bg']*par['r_wd']*0.1*par['p_ecl_wd']*par['f_occ']
    E['BEB_uncat_MSMS'] = Zcfg*(nb_faint + nb_bright)*rc*par['f_occ']*ex(par['dchi2_unc_ms'])
    # resolved neighbours: planets / EBs on the pair (bound wide companions) and the Mira
    Spair = S_of_a(np.array([53.0]), par['S_close'])[0] if par['S_pair'] is None else par['S_pair']
    for k in ('A', 'B'):
        v = PAIR[k]
        Zn = sum(occ(par, cl, v['M'])*ZN[(ds, kind, cl, 'pair' + k)] for cl in CLASSES)
        E[f'NTP_pair{k}'] = Zcfg*Spair*Zn*par[f'lr_pair{k}']
        E[f'NEB_pair{k}'] = Zcfg*par['r_wd']*par['f_Twd']*par['p_ecl_wd']*par['f_occ']*ex(6.7)*par[f'lr_pair{k}']
    E['NEB_NTP_Mira'] = 0.0     # a 5.69-d orbit (a ~ 15 R_sun) lies inside the Mira (hundreds of R_sun); pixels: log10 LR -26
    tot = sum(E.values())
    planet = E['TP_single'] + E['TP_diluted'] + E['TP_on_twin']
    nfp = sum(v for k, v in E.items() if k.startswith('N'))
    res = dict(FPP=1 - planet/tot, NFPP=nfp/tot, post={k: v/tot for k, v in E.items()},
               FPP_strict=1 - (E['TP_single'] + E['TP_diluted'])/tot,
               FPP_system=1 - (planet + E['PTP_companion'])/tot)
    fps = {k: v for k, v in res['post'].items() if not k.startswith('TP')}
    res['dominant_FP'] = max(fps, key=fps.get)
    if detail:
        res['E'] = E
        res['config'] = dict(P_single=Zcfg_single/Zcfg, P_binary=Zcfg_bin/Zcfg,
                             P_twin_dT_le_0p3=fm*np.mean(wl*(CP['dT'] <= 0.3))/Zcfg,
                             P_companion_dT_0p3_2=fm*np.mean(wl*((CP['dT'] > 0.3) & (CP['dT'] <= 2)))/Zcfg,
                             P_companion_dT_2_4p5=fm*np.mean(wl*((CP['dT'] > 2) & (CP['dT'] <= 4.5)))/Zcfg)
        # PTP posterior breakdown
        ptp = wl*Sa*stab2*Z2*(~twin)
        bins_dT = [(0.75, 1.5), (1.5, 2.5), (2.5, 3.5), (3.5, 5)]
        bins_a = [(0, 0.5), (0.5, 4), (4, 6), (6, 65), (65, 1e5)]
        res['PTP_by_dT'] = {f'{lo}-{hi}': float(fm*np.mean(ptp*(CP['dT'] > lo)*(CP['dT'] <= hi))/tot) for lo, hi in bins_dT}
        res['PTP_by_a'] = {f'{lo}-{hi}': float(fm*np.mean(ptp*(CP['a'] > lo)*(CP['a'] <= hi))/tot) for lo, hi in bins_a}
        rp2 = 0.1189*CP['R2']/np.sqrt(np.clip(CP['Dc'], 1e-3, None))*RSUN/REARTH   # approx. planet radius on the companion
        res['PTP_Rp_pct_16_50_84'] = [float(x) for x in np.interp([0.16, 0.5, 0.84], np.cumsum(ptp[np.argsort(rp2)])/ptp.sum(), np.sort(rp2))] if ptp.sum() > 0 else None
        res['PTP_m2_pct_16_50_84'] = [float(x) for x in np.interp([0.16, 0.5, 0.84], np.cumsum(ptp[np.argsort(CP['m2'])])/ptp.sum(), np.sort(CP['m2']))] if ptp.sum() > 0 else None
        res['PTP_sep_arcsec_pct_16_50_84'] = [float(x) for x in np.interp([0.16, 0.5, 0.84], np.cumsum(ptp[np.argsort(CP['sep'])])/ptp.sum(), np.sort(CP['sep']))] if ptp.sum() > 0 else None
        res['PTP_depth_ratio_V_over_I_pct_16_50_84'] = [float(x) for x in np.interp([0.16, 0.5, 0.84], np.cumsum(ptp[np.argsort(RATIO[2])])/ptp.sum(), np.sort(RATIO[2]))] if ptp.sum() > 0 else None
        np.savez(f'{F}/out/ptp_weights.npz', ptp=(ptp/ptp.sum()).astype(np.float32),
                 twin=((wl*Sa*stab2*Z2*twin)/max((wl*Sa*stab2*Z2*twin).sum(), 1e-300)).astype(np.float32))
        res['PTP_mean_companion_mass'] = float(np.sum(ptp*CP['m2'])/np.sum(ptp)) if ptp.sum() > 0 else None
        res['PTP_mean_dT'] = float(np.sum(ptp*CP['dT'])/np.sum(ptp)) if ptp.sum() > 0 else None
        res['Z_TP_single'] = Ztp
        res['nb_uncatalogued'] = dict(faint=nb_faint, bright=nb_bright)
    return res

base = evaluate(BASE, detail=True)
print(json.dumps({k: v for k, v in base.items() if k not in ('E',)}, indent=1, default=float))
jdump(dict(parameters=BASE, result=base, grid_checks=CHK, method=__doc__), 'out/fpp_baseline.json')

# ------------------------------------------------------------------ named one-at-a-time variants
V = []
def var(name, **kw):
    p = dict(BASE); p.update(kw); r = evaluate(p)
    V.append(dict(variant=name, FPP=r['FPP'], NFPP=r['NFPP'], FPP_strict=r['FPP_strict'], FPP_system=r['FPP_system'], dominant=r['dominant_FP'],
                  P_PTP=r['post']['PTP_companion'], P_TP_single=r['post']['TP_single'], P_twin=r['post']['TP_on_twin'],
                  P_TP_dil=r['post']['TP_diluted'], **{k: v for k, v in kw.items() if not isinstance(v, (list, dict, tuple))}))
var('baseline')
var('data REF', ds='ref'); var('e prior Van Eylen', ecc='vaneylen'); var('e circular', ecc='circ')
var('REF + circular', ds='ref', ecc='circ')
for v in (0.001, 0.01, 0.03): var(f'f_N={v}', f_N=v)
for v in (0.0001, 0.002, 0.005): var(f'f_S={v}', f_S=v)
for v in (0.0005, 0.004, 0.01): var(f'f_J={v}', f_J=v)
for v in (-1.0, 1.0, 2.0): var(f'beta_N={v}', beta_N=v)
for v in (0.0, 3.0): var(f'beta_J={v}', beta_J=v)
for v in (0.0, 2.0): var(f'beta_S={v}', beta_S=v)
for v in (0.20, 0.35, 0.45): var(f'f_mult={v}', f_mult=v)
for v in (0.9, 1.6): var(f'mu_loga={v}', mu_loga=v)
for v in (0.8, 1.5): var(f'sig_loga={v}', sig_loga=v)
for v in (-0.5, 1.0): var(f'gamma_q={v}', gamma_q=v)
for v in (0.0, 0.3): var(f'F_twin={v}', F_twin=v)
for v in (0.05, 0.5, 1.0): var(f'S_close={v}', S_close=v)
for v in (0.25, 0.35, 0.40, 0.50): var(f'sig_s={v}', sig_s=v)
var('IPD LR84', ipd='LR84'); var('IPD off', ipd='off')
var('astrometry soft', ast='soft'); var('astrometry off', ast='off'); var('DECaPS2 off', dec=False)
var('IPD+astrometry+DECaPS2 off', ipd='off', ast='off', dec=False)
for v in (0.5, 1.0, 1.5): var(f'dT_twin={v}', dT_twin=v)
for v in (2e-6, 2e-4, 2e-3, 1e-2): var(f'r_wd={v}', r_wd=v)
for v in (6.7, 8.0, 15.2): var(f'dchi2_wd={v}', dchi2_wd=v)
var('dchi2_wd_twin=11.1 (REF value on INV)', dchi2_wd_twin=11.1)
for v in (0.03, 0.3, 1.0): var(f'f_occ={v}', f_occ=v)
for v in (20.0, 33.7): var(f'dchi2_heb={v}', dchi2_heb=v)
var('dchi2_bdms=30', dchi2_bdms=30.0)
for v in (10**-1.4, 10**-1.63, 10**-4.2): var(f'lr_bg3=10^{np.log10(v):.2f}', lr_bg3=v)
var('bg3 compact class dchi2=11.1 (unphysical object)', dchi2_bg3_c=11.1)
var('bg3 compact fraction 0.1', f_compact=0.1)
for v in (10**-2.0, 10**-6.74): var(f'lr_pairA=10^{np.log10(v):.2f}', lr_pairA=v)
var('lr_pairB=10^-3', lr_pairB=1e-3)
for v in (1.5, 4.8): var(f'kappa={v}', kappa=v)
for v in (0.005, 0.1, 0.3): var(f'eps_dss={v}', eps_dss=v)
var('f_Mdw_bg=0.5', f_Mdw_bg=0.5)
var('dchi2_unc_ms=15', dchi2_unc_ms=15.0)
var('S_pair=1', S_pair=1.0)
# combined pessimistic variants (every FP-favouring choice at once)
pess = dict(ds='ref', ecc='vaneylen', f_N=0.001, f_S=0.002, f_J=0.004, beta_N=-1.0, beta_S=0.0, beta_J=0.0, f_mult=0.45,
            F_twin=0.0, gamma_q=-0.5, S_close=1.0, sig_s=0.40, ipd='LR84', ast='soft', dT_twin=0.5, r_wd=2e-3, dchi2_wd=8.0,
            f_occ=1.0, dchi2_heb=20.0, lr_bg3=10**-1.4, lr_pairA=1e-2, kappa=4.8, eps_dss=0.1, f_Mdw_bg=0.5)
var('EXTREME stress (all FP-favouring ends at once, incl. r_wd=2e-3, f_occ=1)', **pess)
p2 = dict(pess); p2.update(ipd='off', ast='off', dec=False); var('EXTREME stress + no imaging/astrometry constraints', **p2)
opt = dict(f_N=0.01, f_S=0.0001, f_J=0.0005, beta_N=1.0, beta_J=3.0, f_mult=0.2, S_close=0.05, sig_s=0.40, r_wd=2e-6, f_occ=0.03)
var('OPTIMISTIC combined', **opt)
pr = dict(f_N=0.001, beta_N=-1.0, f_mult=0.35, F_twin=0.0, S_close=0.5, sig_s=0.25, ipd='LR84', ast='soft', r_wd=2e-4, f_occ=0.3,
          lr_pairA=10**-2.5, ecc='vaneylen')
var('PESSIMISTIC realistic (each input at its plausible FP-favouring end)', **pr)
# ---- what-if follow-up observations (null results), applied to the baseline and to the realistic-pessimistic case
HRI = ([0.03, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0], [2.0, 3.0, 4.0, 4.5, 5.0, 5.5, 6.0])   # typical 8-m speckle (Zorro 832 nm) 5-sigma
for tag, extra in (('', {}), (' | realistic-pessimistic', pr)):
    var('WHAT-IF speckle/AO imaging null' + tag, hri=HRI, **extra)
    var('WHAT-IF multicolour transit depths equal to 5%' + tag, chroma=0.05, **extra)
    var('WHAT-IF high-res spectra: no SB2/RV (dT<2, a<3 au)' + tag, sb2=True, **extra)
    var('WHAT-IF imaging + multicolour' + tag, hri=HRI, chroma=0.05, **extra)
    var('WHAT-IF imaging + multicolour + spectra' + tag, hri=HRI, chroma=0.05, sb2=True, **extra)
    # a seeing-limited ground transit that detects the dip on the M3V itself also resolves pair A (11.7") and the Mira
    e2 = dict(extra); e2.update(chroma=0.05, lr_pairA=1e-8, lr_pairB=1e-12)
    var('WHAT-IF ground multicolour transit on target (chromatic test + pair resolved)' + tag, **e2)
    e3 = dict(e2); e3.update(hri=HRI)
    var('WHAT-IF ground multicolour transit on target + speckle imaging' + tag, **e3)
dfv = pd.DataFrame(V)
dfv.to_csv(f'{F}/out/fpp_variants.csv', index=False)
pd.set_option('display.width', 250)
print(dfv[['variant', 'FPP', 'FPP_strict', 'FPP_system', 'NFPP', 'dominant', 'P_PTP', 'P_twin']].to_string())

# ------------------------------------------------------------------ random variants: all uncertain inputs drawn jointly
def lu(lo, hi): return float(np.exp(rng.uniform(np.log(lo), np.log(hi))))
R = []
for i in range(3000):
    p = dict(BASE)
    p.update(ds=rng.choice(['inv', 'ref']), ecc=rng.choice(['kipping', 'vaneylen', 'kipping', 'circ']),
             f_N=lu(0.001, 0.01), f_S=lu(0.0001, 0.003), f_J=lu(0.0005, 0.005), f_sub=lu(0.05, 0.3),
             beta_N=float(rng.uniform(-1, 2)), beta_S=float(rng.uniform(0, 2)), beta_J=float(rng.uniform(0, 3)),
             f_mult=float(rng.uniform(0.22, 0.40)), mu_loga=float(rng.uniform(1.0, 1.5)), sig_loga=float(rng.uniform(0.9, 1.4)),
             gamma_q=float(rng.uniform(-0.5, 1.0)), F_twin=float(rng.uniform(0, 0.3)), S_close=lu(0.05, 1.0),
             sig_s=float(rng.uniform(0.25, 0.40)), ipd=str(rng.choice(['LR', 'LR84'])), ast=str(rng.choice(['base', 'soft'])),
             dT_twin=float(rng.uniform(0.5, 1.0)), r_wd=lu(2e-6, 2e-4), f_Twd=float(rng.uniform(0.1, 0.5)),
             dchi2_wd=float(rng.uniform(8, 15.2)), f_occ=lu(0.03, 0.3), dchi2_heb=float(rng.uniform(20, 33.7)),
             r_eb_field=lu(0.001, 0.003), f_compact=lu(0.001, 0.1), lr_bg3=10**float(rng.uniform(-4.2, -1.4)),
             lr_pairA=10**float(rng.uniform(-6.7, -3.0)), kappa=float(rng.choice([1.5, 2.2, 3.0, 4.8])),
             eps_dss=lu(0.005, 0.1), f_Mdw_bg=float(rng.uniform(0.1, 0.5)))
    r = evaluate(p)
    R.append(dict(FPP=r['FPP'], NFPP=r['NFPP'], FPP_strict=r['FPP_strict'], FPP_system=r['FPP_system'], dominant=r['dominant_FP'], P_PTP=r['post']['PTP_companion'],
                  **{k: (v if not isinstance(v, (np.floating,)) else float(v)) for k, v in p.items() if k in
                     ('ds', 'ecc', 'f_N', 'f_S', 'f_J', 'beta_N', 'beta_J', 'f_mult', 'S_close', 'sig_s', 'ipd', 'ast', 'F_twin',
                      'gamma_q', 'dT_twin', 'f_occ', 'lr_bg3', 'lr_pairA', 'r_wd')}))
dfr = pd.DataFrame(R); dfr.to_csv(f'{F}/out/fpp_random.csv', index=False)
q = lambda x, p_: float(np.percentile(x, p_))
summ = dict(baseline=dict(FPP=base['FPP'], FPP_strict=base['FPP_strict'], FPP_system=base['FPP_system'], NFPP=base['NFPP'], dominant_FP=base['dominant_FP'], post=base['post'], config=base['config']),
            named_variants=dict(n=len(dfv), FPP_min=float(dfv.FPP.min()), FPP_max=float(dfv.FPP.max()),
                                NFPP_max=float(dfv.NFPP.max()), n_FPP_ge_0p015=int((dfv.FPP >= 0.015).sum()),
                                n_NFPP_ge_0p001=int((dfv.NFPP >= 0.001).sum())),
            random=dict(n=len(dfr), FPP_pct_5_16_50_84_95=[q(dfr.FPP, p_) for p_ in (5, 16, 50, 84, 95)],
                        NFPP_pct_5_50_95=[q(dfr.NFPP, p_) for p_ in (5, 50, 95)],
                        FPP_strict_pct_5_50_95=[q(dfr.FPP_strict, p_) for p_ in (5, 50, 95)],
                        FPP_system_pct_5_50_95=[q(dfr.FPP_system, p_) for p_ in (5, 50, 95)],
                        frac_FPP_system_lt_0p015=float(np.mean(dfr.FPP_system < 0.015)),
                        FPP_pct_by_dataset={d_: [q(dfr.FPP[dfr.ds == d_], p_) for p_ in (16, 50, 84)] for d_ in ('inv', 'ref')},
                        frac_FPP_lt_0p015=float(np.mean(dfr.FPP < 0.015)), frac_NFPP_lt_0p001=float(np.mean(dfr.NFPP < 0.001)),
                        dominant_counts=dfr.dominant.value_counts().to_dict()))
# rank-correlation of FPP with each random input
from scipy.stats import spearmanr
corr = {}
for c in ('f_N', 'f_S', 'f_J', 'beta_N', 'beta_J', 'f_mult', 'S_close', 'sig_s', 'F_twin', 'gamma_q', 'dT_twin', 'f_occ', 'lr_bg3', 'lr_pairA', 'r_wd'):
    corr[c] = float(spearmanr(dfr[c], dfr.FPP).correlation)
summ['random']['spearman_FPP'] = corr
jdump(summ, 'out/fpp_summary.json')
print(json.dumps(summ, indent=1, default=float))
