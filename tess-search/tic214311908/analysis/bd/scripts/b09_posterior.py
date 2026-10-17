"""b09: scenario posteriors with the brown-dwarf (BD) and companion-hosted channels quantified (bd/).
FPP definition (as the task / ../referee2_fpp): FPP = P(the occulter is NOT a planet transiting the M3V or its twin)
  = 1 - (TP_single + TP_diluted + TP_twin)/total.  P_system = P(the occulter is a planet on SOME star of the bound system)
  = (TP_single + TP_diluted + TP_twin + PTP + NTP_pairA)/total (pair A is the wide bound companion).
Framework: evidence E_S = rate per ln P at 5.69 d x likelihood of all data, integrated over the scenario's parameters.
* Configuration prior (single vs unresolved bound companion; companion population with IPD / astrometry / DECaPS2 /
  Gaia-wide / CMD+colour likelihoods; q, a re-weighting; S-type stability; suppression of planets in close binaries)
  and the EB / NTP / background terms are taken READ-ONLY from ../referee2_fpp/scripts/r05_fpp.py (refereed; same
  parameter values, BASE dict) - imported, not modified.
* Light-curve likelihood: own grid (b09a: chi2 for any dilution D from Q coefficients; REF or INV-clean; beta or white
  noise), integrated over b (uniform in [0, 1+k]) with the transit probability (1+k)/aR x (1 + e sin w)/(1 - e^2), and
  over e, w with rho_circ = rho_host * g^3, g = (1 + e sin w)/sqrt(1 - e^2) (CORRECT sign; ../fpp and ../referee2_fpp
  used g^-3, see b10).
* Planet classes and rates as r05 (sub 2-4, N 4-7, S 7-10, J 10-22 R_E; flat in ln R), planet e-prior Kipping (2013).
* BD class: radius from ATMO 2020 + inflation (b03 draws), rate f_BD per star at 3-10 d (b01: median 1.3e-4), BD
  e-prior (b01: hierarchical Beta fit to M-dwarf BDs at 2.5-16 d = baseline; variants Ferreira+26 short-period,
  Ferreira+26 long-period (tidally unprocessed, pessimistic), Kipping).
* BD on the M3V (single star): full Monte Carlo over (M_bd, age, R_bd, Teff, e, w, rho*, b) including (i) the occultation
  constraint at the eccentric phase/duration (b07 scan, INV-clean), (ii) the Gaia RVS error-smearing likelihood ratio
  (b08), (iii) the Gaia RUWE likelihood ratio (b08).
* Shape calibration (b06 injections at fake epochs): factor C = LR_calibrated / LR_formal for the statistic Lambda
  (grazing BD vs planet; companion-hosted k = 0.4 body vs planet), applied to the corresponding scenario evidence only
  when C > 1 (i.e. only to loosen the formal penalty).
usage: b09_posterior.py -> out/b09_posterior.json, out/b09_variants.csv"""
import sys, os
sys.path.insert(0, 'tic214311908/analysis/bd/scripts')
from blib import *
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
sys.path.insert(0, f'{B}/referee2_fpp/scripts')
import r05_fpp as R5                         # read-only import (configuration prior, companion population, EB terms)
rng = np.random.default_rng(21)
LNP = np.log(10/3.)
GR = np.load(os.environ.get('B09_GRID', f'{BD}/data/b09a_grid.npz'))
KS, BBS, RHOS = GR['ks'], GR['bbs'], GR['rhos']
LK, LR_ = np.log(KS), np.log(RHOS)
wbb = np.gradient(BBS)                         # trapezoid-like weights in beta_b (sum ~ 1)
wbb = wbb/np.sum(wbb)*(BBS[-1] - BBS[0] + 0.002)
DG = np.r_[np.logspace(np.log10(0.004), np.log10(0.999), 44), 1.0]
ECC = np.load(f'{BD}/data/ecc_draws.npz')
BDD = np.load(f'{BD}/data/bd_draws.npz')
SEC = np.load(f'{BD}/out/b07_secondary_scan.npz')
RV = np.load(f'{BD}/out/b08_rv_lr.npz')
GJ = jload(f'{BD}/out/b08_gaia.json')
CLASSES = {'sub': (2.0, 4.0), 'N': (4.0, 7.0), 'S': (7.0, 10.0), 'J': (10.0, 22.0)}

def chi2_grid(combo, D):
    return GR[f'QYY_{combo}'] - 2*D*GR[f'QYS_{combo}'] + D*D*GR[f'QSS_{combo}'] - GR[f'chi2ML_{combo}']

_JC = {}
def Jtab(combo):
    """J[D, k, rho_c] = (1 + k) * int_0^1 dbeta exp(-dchi2/2)  (b uniform in [0, 1+k])"""
    if combo in _JC: return _JC[combo]
    J = np.zeros((len(DG), len(KS), len(RHOS)))
    for i, D in enumerate(DG):
        c = chi2_grid(combo, D)
        J[i] = (1 + KS)[:, None]*np.einsum('kbr,b->kr', np.exp(-0.5*np.clip(c, -50, 1400)), wbb)
    _JC[combo] = J
    return J

def ecc_draws(kind, n):
    if kind in ECC.files: e = rng.choice(ECC[kind], n)
    elif kind == 'circ': e = np.zeros(n)
    else: raise ValueError(kind)
    e = np.clip(e, 0, 0.9); w = rng.uniform(0, 2*np.pi, n)
    g, geo = g_ecc(e, w)
    return e, w, g, geo

def Z_table(combo, ekind, wfun, RH, RHOH, ne=1500):
    """Z[D, rho_host, R_host] = (1/aR) < geo J(D, k, rho_h g^3) >_{e,w} integrated over ln k with weight wfun(R = k R_h, in R_E)
    (wfun = occurrence density per ln R, unit normalised)."""
    J = Jtab(combo)
    e, w, g, geo = ecc_draws(ekind, ne)
    dlk = LK[1] - LK[0]
    T = np.zeros((len(DG), len(RHOH), len(RH)))
    for ih, rh in enumerate(RHOH):
        lrc = np.log(rh*g**3)
        ok = (lrc >= LR_[0]) & (lrc <= LR_[-1])
        x = np.interp(lrc[ok], LR_, np.arange(len(LR_)))
        i0 = np.clip(np.floor(x).astype(int), 0, len(LR_) - 2); f = x - i0
        Jw = (J[:, :, i0]*(1 - f) + J[:, :, i0 + 1]*f)                       # [D, k, ne_ok]
        Jm = np.einsum('dkn,n->dk', Jw, geo[ok])/len(e)                      # <geo J>
        ra = 1/aR_of_rho(rh)
        for ir, Rh in enumerate(RH):
            Rp = KS*Rh*RSUN/REARTH
            T[:, ih, ir] = ra*np.sum(Jm*wfun(Rp)[None, :], axis=1)*dlk
    return T

def w_flat(lo, hi):
    return lambda R: ((R >= lo) & (R < hi))/np.log(hi/lo)
_rbd = BDD['R_rsun']*RSUN/REARTH
_hist, _edges = np.histogram(np.log(_rbd), bins=60, density=True)
def w_bd(R):
    lr = np.log(np.clip(R, 1e-3, None)); i = np.searchsorted(_edges, lr) - 1
    return np.where((i >= 0) & (i < len(_hist)), _hist[np.clip(i, 0, len(_hist) - 1)], 0.0)

RH = np.logspace(np.log10(0.09), np.log10(0.75), 26)
RHOH = np.logspace(np.log10(3.0), np.log10(150.), 30)

def interp_Z(T, D, rho, Rh):
    f = RegularGridInterpolator((np.log(DG), np.log(RHOH), np.log(RH)), np.log(np.clip(T, 1e-300, None)), bounds_error=False, fill_value=None)
    D = np.clip(D, DG[0], 1.0); rho = np.clip(rho, RHOH[0], RHOH[-1]); Rh = np.clip(Rh, RH[0], RH[-1])
    return np.exp(f(np.column_stack([np.log(D), np.log(rho), np.log(Rh)])))

# ---------------- BD on the M3V: detailed Monte Carlo
SECP, SECD = SEC['phases'], SEC['durs_d']
def lr_secondary(phi, dur, depth):
    i = np.clip(np.round((phi - SECP[0])/(SECP[1] - SECP[0])).astype(int), 0, len(SECP) - 1)
    j = np.argmin(np.abs(np.log(np.clip(dur, 1e-3, None))[:, None] - np.log(SECD)[None, :]), axis=1)
    dm = SEC['depth'][i, j]; se = SEC['err_emp'][i, j]
    ok = np.isfinite(dm) & np.isfinite(se)
    out = np.ones(len(phi))
    out[ok] = np.exp(-0.5*((dm[ok] - depth[ok])**2 - dm[ok]**2)/se[ok]**2)
    return out
def mean_anom(f, e):
    E = 2*np.arctan2(np.sqrt(1 - e)*np.sin(f/2), np.sqrt(1 + e)*np.cos(f/2))
    return E - e*np.sin(E)
RVK, RVE = RV['Ks'], RV['es']
_rvkeys = [k_ for k_ in RV.files if k_.startswith('LR_w')]
RVLR = np.mean([RV[k_] for k_ in _rvkeys], axis=0)   # average over CCF width (12/16/20 km/s) and peak (0.3/0.6)
_frv = RegularGridInterpolator((RVE, RVK), RVLR, bounds_error=False, fill_value=None)
_ru = GJ['ruwe']['table']['sigma_ccd_0.22']
_ruM = np.array([r['M_bd'] for r in _ru]); _ruL = np.array([r['LR_mean'] for r in _ru])

def bd_on_m3v(combo, ekind, N=1_500_000, use_sec=True, use_rv=True, use_ruwe=True, chunk=250_000):
    fQys = RegularGridInterpolator((LK, BBS, LR_), GR[f'QYS_{combo}'], bounds_error=False, fill_value=0.)
    fQss = RegularGridInterpolator((LK, BBS, LR_), GR[f'QSS_{combo}'], bounds_error=False, fill_value=0.)
    qyy, c0 = float(GR[f'QYY_{combo}']), float(GR[f'chi2ML_{combo}'])
    tot = 0.; parts = dict(lc=0., sec=0., rv=0., ruwe=0.); wsum_e = []; acc = []
    n_done = 0
    while n_done < N:
        n = min(chunk, N - n_done); n_done += n
        idx = rng.integers(0, len(BDD['k']), n)
        k = BDD['k'][idx].astype(float); Mbd = BDD['M_MJ'][idx].astype(float); fr = BDD['flux_ratio'][idx].astype(float)
        fr = fr*np.exp(rng.normal(0, 0.35, n))                   # surface-brightness systematics (BB vs empirical, b03 checks)
        rho = rng.normal(RHO_STAR, RHO_ERR, n)
        e, w, g, geo = ecc_draws(ekind, n)
        rc = rho*g**3
        bb = rng.uniform(0, 1, n)
        aR = aR_of_rho(rho)
        pts = np.column_stack([np.log(k), bb, np.log(np.clip(rc, RHOS[0], RHOS[-1]))])
        inside = (rc >= RHOS[0]) & (rc <= RHOS[-1]) & (k >= KS[0]) & (k <= KS[-1])
        c = qyy - 2*fQys(pts) + fQss(pts) - c0
        L = np.where(inside, np.exp(-0.5*np.clip(c, -50, 1400)), 0.)
        wt = (1 + k)/aR*geo*L
        # occultation
        b = bb*(1 + k)
        f_tr = np.pi/2 - w; f_oc = -np.pi/2 - w
        dphi = ((mean_anom(f_oc, e) - mean_anom(f_tr, e))/(2*np.pi)) % 1
        b_oc = b*(1 + e*np.sin(w))/(1 - e*np.sin(w))
        inc = np.arccos(np.clip(b/aR*(1 + e*np.sin(w))/(1 - e**2), 0, 1))
        x = np.sqrt(np.clip((1 + k)**2 - b_oc**2, 0, None))/(aR*np.sin(inc))
        Toc = P/np.pi*np.arcsin(np.clip(x, 0, 1))*np.sqrt(1 - e**2)/(1 - e*np.sin(w))
        full = np.clip(((1 + k) - b_oc)/(2*k), 0, 1)                # crude partial-occultation fraction
        dep = np.where(b_oc < 1 + k, fr*full, 0.)
        Ls = lr_secondary(dphi, np.clip(Toc, 0.4/24, 2.4/24), dep) if use_sec else np.ones(n)
        Mtot = (0.370 + Mbd/1047.57)*MSUN
        K = (2*np.pi*G/(P*DAY))**(1/3.)*(Mbd*MJUP)/Mtot**(2/3.)/np.sqrt(1 - e**2)/1e5   # km/s
        Lrv = _frv(np.column_stack([np.clip(e, 0, 0.7), np.clip(K, 0, 40)])) if use_rv else np.ones(n)
        Lru = np.interp(Mbd, _ruM, _ruL) if use_ruwe else np.ones(n)
        tot += np.sum(wt*Ls*Lrv*Lru)
        parts['lc'] += np.sum(wt); parts['sec'] += np.sum(wt*Ls); parts['rv'] += np.sum(wt*Lrv); parts['ruwe'] += np.sum(wt*Lru)
        acc.append(np.column_stack([e, w, k, rc, Mbd, wt*Ls*Lrv*Lru, K, b]))
    A = np.concatenate(acc)
    W = A[:, 5]; sw = W.sum()
    stats = {}
    if sw > 0:
        q = lambda col, p: float(np.interp(p, np.cumsum(W[np.argsort(A[:, col])])/sw, np.sort(A[:, col])))
        stats = dict(e_med=q(0, 0.5), e_p05=q(0, 0.05), e_p95=q(0, 0.95), sinw_med=float(np.sum(W*np.sin(A[:, 1]))/sw),
                     k_med=q(2, 0.5), rhoc_med=q(3, 0.5), rhoc_p05=q(3, 0.05), rhoc_p95=q(3, 0.95), M_med=q(4, 0.5), M_p05=q(4, 0.05), M_p95=q(4, 0.95),
                     K_med=q(6, 0.5), K_p05=q(6, 0.05), K_p95=q(6, 0.95), b_med=q(7, 0.5), n_eff=float(sw**2/np.sum(W**2)))
    ratios = {kk: (v/parts['lc'] if parts['lc'] > 0 else np.nan) for kk, v in parts.items()}
    sub = A[rng.choice(len(A), 20000, p=W/sw)] if sw > 0 else A[:0]
    return tot/N, ratios, stats, sub

# ---------------- configuration machinery from r05 (read-only)
CP = R5.CP
def config(par):
    Lc = np.ones(R5.NC)
    if par['ipd'] != 'off': Lc *= CP[par['ipd']]
    Lc *= {'base': CP['L_ast'], 'soft': CP['L_ast']**0.4, 'off': 1.0}[par['ast']]
    if par['dec'] == 'graded': Lc *= CP['L_dec']
    if par['wide']: Lc *= CP['L_wide']
    q = CP['q']; qmin = 0.075/CP['m1']; gq = par['gamma_q']
    normq = (1 - qmin**(gq + 1))/(gq + 1)
    wq = ((1 - par['F_twin'])*q**gq/normq + par['F_twin']*(q >= 0.95)/0.05)*(1 - qmin)
    lg = CP['loga']
    wa = (np.exp(-0.5*((lg - par['mu_loga'])/par['sig_loga'])**2)/par['sig_loga'])/(np.exp(-0.5*((lg - 1.30)/1.16)**2)/1.16)
    Lcmd = R5.cmd_like(par, R5.DX_C, R5.DR_C); Lcmd_single = float(np.atleast_1d(R5.cmd_like(par, np.zeros(1), np.zeros(1)))[0])
    fm = par['f_mult']
    wl = wq*wa*Lc*Lcmd
    Zs = (1 - fm)*Lcmd_single
    Sa = np.where(CP['a'] <= 10, par['S_close'], np.where(CP['a'] >= 200, 1.0, par['S_close'] + (1 - par['S_close'])*np.log(CP['a']/10)/np.log(20)))
    return wl, Zs, fm, Sa

_ZT = {}
def tables(combo, ekind_pl, ekind_bd):
    key = (combo, ekind_pl, ekind_bd)
    if key in _ZT: return _ZT[key]
    T = {cl: Z_table(combo, ekind_pl, w_flat(*CLASSES[cl]), RH, RHOH) for cl in CLASSES}
    T['BD'] = Z_table(combo, ekind_bd, w_bd, RH, RHOH)
    _ZT[key] = T
    return T

def evaluate(par):
    combo = par['combo']
    T = tables(combo, par['ecc_pl'], par['ecc_bd'])
    wl, Zs, fm, Sa = config(par)
    occ = lambda cl, M: par['f_' + cl]/LNP*(np.asarray(M)/0.37)**par['beta_' + cl]
    PL = ('sub', 'N', 'S', 'J')
    rho_tp = rng.normal(9.42, 0.88, 400)
    Ztp = {cl: float(np.mean(interp_Z(T[cl], np.ones(400), rho_tp, np.full(400, 0.381)))) for cl in PL}
    E = {}
    E['TP_single'] = Zs*sum(occ(cl, 0.37)*Ztp[cl] for cl in PL)
    key = ('bdm3v', combo, par['ecc_bd'], par['bd_sec'], par['bd_rv'], par['bd_ruwe'])
    if key not in _ZT:
        _ZT[key] = bd_on_m3v(combo, par['ecc_bd'], use_sec=par['bd_sec'], use_rv=par['bd_rv'], use_ruwe=par['bd_ruwe'])
    zbd, bd_ratios, bd_stats, bd_sub = _ZT[key]
    E['BD_on_M3V_single'] = Zs*occ('BD', 0.37)*zbd*par['C_graz']
    Z1 = {cl: interp_Z(T[cl], 1 - CP['Dc'], CP['rho1'], CP['R1']) for cl in list(PL) + ['BD']}
    Z2 = {cl: interp_Z(T[cl], CP['Dc'], CP['rho2'], CP['R2']) for cl in list(PL) + ['BD']}
    p1 = sum(occ(cl, CP['m1'])*Z1[cl] for cl in PL)*Sa*CP['stab1']
    p2 = sum(occ(cl, CP['m2'])*Z2[cl] for cl in PL)*Sa*CP['stab2']
    SBD = Sa if par['S_BD'] else 1.0
    b1 = occ('BD', CP['m1'])*Z1['BD']*SBD*CP['stab1']
    b2 = occ('BD', CP['m2'])*Z2['BD']*SBD*CP['stab2']*par['C_comp']
    twin = CP['dT'] <= par['dT_twin']
    E['TP_diluted'] = fm*np.mean(wl*p1)
    E['TP_on_twin'] = fm*np.mean(wl*p2*twin)
    E['PTP_companion'] = fm*np.mean(wl*p2*~twin)
    E['BD_on_primary_of_binary'] = fm*np.mean(wl*b1)*par['C_graz']
    E['BD_on_companion'] = fm*np.mean(wl*b2)
    ex = lambda d: np.exp(-0.5*d)
    nt = CP['dT'] > 0.3
    E['EB_WD_on_M3V'] = (Zs + fm*np.mean(wl*nt))*par['r_wd']*par['f_Twd']*par['p_ecl_wd']*par['f_occ']*ex(par['dchi2_wd'])
    E['HEB_MSMS'] = fm*np.mean(wl*(CP['dT'] < 4.5))*par['r_close']*par['p_ecl_heb']*par['f_occ']*ex(par['dchi2_heb'])
    Spair = par['S_close'] + (1 - par['S_close'])*np.log(53/10)/np.log(20)
    PA = R5.PAIR_A
    ZA = sum(occ(cl, PA['M'])*float(interp_Z(T[cl], np.ones(1), np.array([PA['rho']]), np.array([PA['Reff']]))[0]) for cl in PL)
    Zcfg = Zs + fm*np.mean(wl)
    E['NTP_pairA'] = Zcfg*Spair*ZA*par['lr_pairA']
    tot0 = sum(E.values())
    E['BG_all'] = par['bg_post']*tot0
    tot = sum(E.values())
    post = {k: float(v/tot) for k, v in E.items()}
    planet = post['TP_single'] + post['TP_diluted'] + post['TP_on_twin']
    res = dict(FPP=1 - planet, P_planet_M3V_or_twin=planet, P_planet_any_star_of_system=planet + post['PTP_companion'] + post['NTP_pairA'],
               FPP_strict=1 - post['TP_single'] - post['TP_diluted'], BD_total=post['BD_on_M3V_single'] + post['BD_on_primary_of_binary'] + post['BD_on_companion'],
               post=post, bd_m3v_ratios=bd_ratios, bd_m3v_stats=bd_stats)
    if par.get('detail'):
        ptp = wl*p2*~twin; bd2 = wl*b2
        wsum = lambda w_, x_: float(np.sum(w_*x_)/np.sum(w_))
        wq_ = lambda w_, x_, p_: float(np.interp(p_, np.cumsum(w_[np.argsort(x_)])/np.sum(w_), np.sort(x_)))
        res['PTP_props'] = dict(m2_mean=wsum(ptp, CP['m2']), dT_med=wq_(ptp, CP['dT'], 0.5), dT_p05=wq_(ptp, CP['dT'], 0.05), dT_p95=wq_(ptp, CP['dT'], 0.95),
                                sep_med=wq_(ptp, CP['sep'], 0.5), sep_p05=wq_(ptp, CP['sep'], 0.05), sep_p95=wq_(ptp, CP['sep'], 0.95),
                                a_med=wq_(ptp, CP['a'], 0.5), dI_med=wq_(ptp, CP['dI'], 0.5), dI_p95=wq_(ptp, CP['dI'], 0.95),
                                frac_sep_gt_0p05=wsum(ptp, CP['sep'] > 0.05), frac_sep_gt_0p1=wsum(ptp, CP['sep'] > 0.1), frac_sep_gt_0p2=wsum(ptp, CP['sep'] > 0.2),
                                VI_depth_ratio_med=wq_(ptp, CP['DcV']/CP['DcI'], 0.5))
        res['BDcomp_props'] = dict(m2_mean=wsum(bd2, CP['m2']), dT_med=wq_(bd2, CP['dT'], 0.5), sep_med=wq_(bd2, CP['sep'], 0.5),
                                   frac_sep_gt_0p1=wsum(bd2, CP['sep'] > 0.1), k_med=float('nan'))
        res['_arrays'] = dict(ptp=ptp, bd2=bd2, fm=fm, tot=tot)
    return res

BASE = dict(R5.BASE)
BASE.update(combo='ref_beta', ecc_pl='kipping13', ecc_bd='mdwarf_BD_hier', f_BD=1.3e-4, beta_BD=0.0, S_BD=False,
            bd_sec=True, bd_rv=True, bd_ruwe=True, C_graz=1.0, C_comp=1.0)
if __name__ == '__main__':
    calib = jload(f'{BD}/out/b06_calibration.json') if os.path.exists(f'{BD}/out/b06_calibration.json') else None
    out = dict(method=__doc__, base_params={k: v for k, v in BASE.items() if not isinstance(v, (np.ndarray,))})
    rows = []
    def run(name, **kw):
        p = dict(BASE); p.update(kw)
        r = evaluate(p)
        po = r['post']
        row = dict(variant=name, FPP=r['FPP'], P_system=r['P_planet_any_star_of_system'], BD_total=r['BD_total'], BD_M3V=po['BD_on_M3V_single'],
                   BD_primary=po['BD_on_primary_of_binary'], BD_comp=po['BD_on_companion'], PTP=po['PTP_companion'], TP_single=po['TP_single'],
                   TP_diluted=po['TP_diluted'], TP_twin=po['TP_on_twin'], NTP_pairA=po['NTP_pairA'], WD=po['EB_WD_on_M3V'], HEB=po['HEB_MSMS'],
                   **{k: v for k, v in kw.items() if not isinstance(v, (list, dict, np.ndarray))})
        rows.append(row)
        print(f"{name:45s} FPP {r['FPP']:.4f} Psys {r['P_planet_any_star_of_system']:.4f} BD {r['BD_total']:.2e} [M3V {po['BD_on_M3V_single']:.1e} prim {po['BD_on_primary_of_binary']:.1e} comp {po['BD_on_companion']:.1e}] PTP {po['PTP_companion']:.2e} NTP {po['NTP_pairA']:.1e}", flush=True)
        return r
    # calibration factors from b06 (if available)
    Cg = {}; Cc = {}
    if calib:
        for ds in ('ref', 'inv'):
            Cg[ds] = max(1.0, calib[ds]['grazBD']['C_factor']) if 'grazBD' in calib[ds] else 1.0
            Cc[ds] = max(1.0, calib[ds]['compBD40']['C_factor']) if 'compBD40' in calib[ds] else 1.0
    base = run('baseline (REF beta, mdwarf-hier BD e, f_BD 1.3e-4)', detail=True)
    out['baseline'] = {k: v for k, v in base.items() if k != '_arrays'}
    for kk_, vv_ in _ZT.items():
        if isinstance(kk_, tuple) and kk_[0] == 'bdm3v' and kk_[1] == 'ref_beta':
            np.savez_compressed(f'{BD}/data/b09_bdm3v_post_draws.npz', cols=np.array(['e', 'w', 'k', 'rho_c', 'M_MJ', 'weight', 'K_kms', 'b']), A=vv_[3])
    for combo in ('inv_beta', 'ref_white', 'inv_white'):
        r = run(f'combo={combo}', combo=combo, detail=(combo == 'inv_beta'))
        if combo == 'inv_beta': out['inv_beta'] = {k: v for k, v in r.items() if k != '_arrays'}
    if calib:
        for ds in ('ref', 'inv'):
            r = run(f'{ds}_beta + injection calibration (C_graz {Cg[ds]:.2f}, C_comp {Cc[ds]:.2f})', combo=f'{ds}_beta', C_graz=Cg[ds], C_comp=Cc[ds], detail=True)
            out[f'{ds}_calibrated'] = {k: v for k, v in r.items() if k != '_arrays'}
    for ek in ('ferreira_SP', 'ferreira_LP', 'kipping13', 'vaneylen19_halfnormal', 'circ'):
        for combo in ('ref_beta', 'inv_beta'):
            run(f'BD e-prior {ek} ({combo})', ecc_bd=ek, combo=combo)
    for fb in (5e-5, 2.5e-4, 3.5e-4, 1e-3):
        for combo in ('ref_beta', 'inv_beta'):
            run(f'f_BD={fb} ({combo})', f_BD=fb, combo=combo)
    run('BD rate scales with host mass (beta_BD=1)', beta_BD=1.0)
    run('BD suppressed in close binaries like planets', S_BD=True)
    run('no occultation term', bd_sec=False)
    run('no Gaia RV term', bd_rv=False)
    run('no Gaia RUWE term', bd_ruwe=False)
    run('no secondary/RV/RUWE terms', bd_sec=False, bd_rv=False, bd_ruwe=False)
    run('planet e-prior Van Eylen', ecc_pl='vaneylen19_halfnormal')
    run('planet e circular', ecc_pl='circ')
    run('pessimistic: INV, LP e, f_BD 1e-3, no Gaia', combo='inv_beta', ecc_bd='ferreira_LP', f_BD=1e-3, bd_rv=False, bd_ruwe=False)
    # follow-up what-ifs (REF baseline)
    p_ = dict(BASE); p_['detail'] = True
    rb = evaluate(p_); A = rb['_arrays']; po = rb['post']
    def whatif(name, ptp_keep, bd2_keep, bdm3v_keep, bdprim_keep=None):
        E = {k: v for k, v in po.items()}
        E['PTP_companion'] = po['PTP_companion']*ptp_keep
        E['BD_on_companion'] = po['BD_on_companion']*bd2_keep
        E['BD_on_M3V_single'] = po['BD_on_M3V_single']*bdm3v_keep
        E['BD_on_primary_of_binary'] = po['BD_on_primary_of_binary']*(bdm3v_keep if bdprim_keep is None else bdprim_keep)
        tot = sum(E.values()); pl = (E['TP_single'] + E['TP_diluted'] + E['TP_on_twin'])/tot
        print(f'what-if {name:60s} FPP {1 - pl:.5f}', flush=True)
        return dict(FPP=1 - pl, post={k: v/tot for k, v in E.items()})
    sep = CP['sep']; dI = CP['dI']
    # 8-m speckle (assumed 5-sigma contrast, 832 nm ~ I): dI < 3.0 at 0.05", 4.5 at 0.1", 5.5 at 0.2", 6.5 at >= 0.5"; detected -> excluded
    sp_s = np.array([0.05, 0.1, 0.2, 0.5, 1.2]); sp_v = np.array([3.0, 4.5, 5.5, 6.5, 7.0])
    lim = np.where(sep >= 0.05, np.interp(sep, sp_s, sp_v), -99)
    keep = (dI >= lim)
    fk = lambda arr: float(np.sum(arr*keep)/np.sum(arr)) if np.sum(arr) > 0 else 1.0
    WI = {}
    WI['speckle_8m'] = whatif('8-m speckle null (dI 3.0@0.05", 4.5@0.1", 5.5@0.2")', fk(A['ptp']), fk(A['bd2']), 1.0)
    WI['rv_2epochs_1kms'] = whatif('two RV epochs at quadrature, sigma 0.5-1 km/s (BD on M3V/primary excluded)', 1.0, 1.0, 0.0, 0.0)
    WI['twoband_5pct'] = whatif('two-band transit depth (V/I or g/z) at 5% -> companion-hosted excluded', 0.02, 0.02, 1.0)
    WI['speckle+rv'] = whatif('speckle + RV', fk(A['ptp']), fk(A['bd2']), 0.0, 0.0)
    WI['speckle+rv+twoband'] = whatif('speckle + RV + two-band', fk(A['ptp'])*0.02, fk(A['bd2'])*0.02, 0.0, 0.0)
    out['whatif'] = WI; out['speckle_keep'] = dict(PTP=fk(A['ptp']), BD_comp=fk(A['bd2']))
    pd.DataFrame(rows).to_csv(f'{BD}/out/b09_variants.csv', index=False)
    jdump(out, f'{BD}/out/b09_posterior.json')
