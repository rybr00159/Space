"""r05: own Bayesian scenario comparison for the 5.69-d signal (referee-2 re-derivation of fpp/f04).
Evidence E_S = prior per target per ln P at 5.69 d x likelihood of all data (integrated over internal parameters).
Planet-like bodies (planets of four radius classes AND brown dwarfs, a class the analyst omitted) on any host use the
r03 light-curve grid: Z = < (R/a) geo sum_b db I(D, rho_c, b) w(k_best R_host) >_{e, w}, rho_c = rho_host g(e,w)^-3,
with I = int exp(-(chi2 - chi2_ref)/2) dlnk and w = occurrence density per ln R (flat in ln R within a class, or a
power law R^-alpha across 4-10 R_e, normalised to the same 4-7 R_e rate).
Configurations: single star (1 - f_mult) S(x_obs) or a bound companion from r04 (weights: q and a re-weighting, IPD,
astrometry, DECaPS2, Gaia-wide, CMD S(x_obs - delta)), with S the single-star offset density from r02 (mixture fit)
or a Gaussian (analyst's form). Planet suppression in binaries S(a) and S-type stability.
EBs: WD occulted by the M3V (referee2_shape revised dchi2 ~9-11), MS hierarchical EB (shape dchi2 27.5), NTP/NEB pair A
with the r01 pixel LRs. Background-star scenarios are < 1e-4 in the analyst's own bookkeeping and are carried as a fixed
3e-5 posterior-equivalent term (checked, not re-derived).
Usage: r05_fpp.py   -> out/r05_baseline.json, out/r05_variants.csv, out/r05_random.csv"""
import numpy as np, json, sys, pandas as pd
from scipy.interpolate import RegularGridInterpolator
sys.path.insert(0, 'tic214311908/analysis/referee2_fpp/scripts')
from rlib import *
rng = np.random.default_rng(11)
LNP = np.log(10/3.)
def aR_of_rho(rho): return (np.asarray(rho)*G*(P*DAY)**2/(3*np.pi))**(1/3.)
RH = np.logspace(np.log10(0.09), np.log10(0.75), 30)
RHOH = np.logspace(np.log10(3.0), np.log10(150.), 36)
CLASSES = {'sub': (2.0, 4.0), 'N': (4.0, 7.0), 'S': (7.0, 10.0), 'J': (10.0, 22.0), 'BD': (8.2, 12.3)}
def w_class(cl, R, alpha):
    lo, hi = CLASSES[cl]
    ins = (R >= lo) & (R < hi)
    if cl in ('N', 'S') and alpha > 0:
        # power law across 4-10 R_e, normalised so the 4-7 class integrates to 1 (f_N) and 7-10 keeps its own rate
        norm = (lo**-alpha - hi**-alpha)/alpha
        return ins*R**-alpha/norm
    return ins/np.log(hi/lo)
def ecc(kind, n=2500):
    if kind == 'kipping': e = rng.beta(0.867, 3.03, n)
    elif kind == 'vaneylen': e = np.abs(rng.normal(0, 0.32, n))
    else: e = np.zeros(n)
    e = np.clip(e, 0, 0.9); w = rng.uniform(0, 2*np.pi, n)
    return (1 + e*np.sin(w))/np.sqrt(1 - e**2), (1 + e*np.sin(w))/(1 - e**2)
GR = {}; ZT = {}; ZCACHE = {}
def build(ds, kind, alpha):
    key = (ds, kind, alpha)
    if key in ZT: return ZT[key]
    if ds not in GR: GR[ds] = dict(np.load(f'{R2}/out/r03_grid_{ds}.npz'))
    g = GR[ds]
    db = g['b'][1] - g['b'][0]
    I = np.nan_to_num(g['I']); K = np.nan_to_num(g['k'], nan=1.0)
    Rp = K[..., None]*RH[None, None, None, :]*RSUN/REARTH
    ge, geo = ecc(kind)
    lr = np.log(g['rho'])
    out = {}
    for cl in CLASSES:
        J = np.sum(I[..., None]*w_class(cl, Rp, alpha), axis=2)*db          # [D, rho_c, R]
        lJ = np.log(np.clip(J, 1e-300, None))
        T = np.zeros((len(g['D']), len(RHOH), len(RH)))
        for ih, rh in enumerate(RHOH):
            rc = np.log(rh*ge**-3); ok = (rc >= lr[0]) & (rc <= lr[-1])
            ra = 1/aR_of_rho(rh)
            for iD in range(len(g['D'])):
                v = np.array([np.interp(rc[ok], lr, lJ[iD, :, iR]) for iR in range(len(RH))])   # [R, ne]
                T[iD, ih, :] = ra*np.sum(geo[ok][None, :]*np.exp(v), axis=1)/len(ge)
        out[cl] = RegularGridInterpolator((np.log(g['D'][::-1]), np.log(RHOH), np.log(RH)), np.log(np.clip(T[::-1], 1e-300, None)),
                                          bounds_error=False, fill_value=None)
    out['Dmin'] = float(g['D'].min())
    ZT[key] = out
    return out
def Zof(tab, cl, D, rho, R):
    D = np.atleast_1d(D).astype(float); rho = np.atleast_1d(rho).astype(float); R = np.atleast_1d(R).astype(float)
    n = max(len(D), len(rho), len(R))
    D, rho, R = [np.broadcast_to(v, (n,)) for v in (D, rho, R)]
    v = np.exp(tab[cl](np.column_stack([np.log(np.clip(D, tab['Dmin'], 1.0)), np.log(np.clip(rho, RHOH[0], RHOH[-1])),
                                        np.log(np.clip(R, RH[0], RH[-1]))])))
    return np.where(D < tab['Dmin']*0.999, 0.0, v)

CP = dict(np.load(f'{R2}/out/r04_companions.npz'))
NC = len(CP['a'])
CMDJ = json.load(open(f'{R2}/out/r02_cmd_colour.json'))
SLOPE = CMDJ['ridge']['slope']; CGS = CMDJ['GK_fit'][0]
XOBS, ROBS = CMDJ['target']['x'], CMDJ['target']['r']
fl_ = lambda mm: 10**(-0.4*np.asarray(mm))
DX_C = CP['dMbin'] - SLOPE*CP['dcol']                                     # companion shift of the M_Ks offset
DR_C = (-2.5*np.log10(1 + fl_(CP['dG'])) - CP['dMbin']) - CGS*CP['dcol']   # companion shift of the G-Ks residual
def _px(x, p):
    a = 2/(np.sqrt(2*np.pi)*(p['sb'] + p['sf']))
    c = np.where(x < p['mu'], a*np.exp(-0.5*((x - p['mu'])/p['sb'])**2), a*np.exp(-0.5*((x - p['mu'])/p['sf'])**2))
    return (1 - p['t'])*c + p['t']*np.exp(-0.5*((x - p['mt'])/p['st'])**2)/(np.sqrt(2*np.pi)*p['st'])
def _pr(r, x, p):
    return np.exp(-0.5*((r - p['a'] - p['c']*x - p['c2']*x*x)/p['sr'])**2)/(np.sqrt(2*np.pi)*p['sr'])
def cmd_like(par, dx, dr):
    """predictive density of the target's (x, r) given the companion shifts (dx, dr); dx = dr = 0 for a single star"""
    xo = XOBS if par['xobs'] is None else par['xobs']
    if par['cmd'] == 'gauss':      # analyst's form: Gaussian single-star offset, no colour term
        s = par['sig_s']; return np.exp(-0.5*((xo - dx)/s)**2)/(np.sqrt(2*np.pi)*s)
    if par['cmd'] == 'flat':
        return np.ones_like(np.asarray(dx, float))
    prof = CMDJ['profiles'][f"F_twin={par['cmd_Ft']}"]
    rows = prof['rows']
    if par['cmd_sb'] == 'marg': wts = prof['profile_weights']
    else: wts = [1.0 if abs(r_['sig_b'] - par['cmd_sb']) < 1e-6 else 0.0 for r_ in rows]
    out = 0.0
    for w_, r_ in zip(wts, rows):
        if w_ < 1e-6: continue
        p = r_['params']
        L = _px(xo - dx, p)
        if par['cmd'] == 'colour': L = L*_pr(ROBS - dr, xo - dx, p)
        out = out + w_*L
    return out
PAIR_A = dict(M=0.30, R=0.31, T=13.375)
PAIR_A['rho'] = PAIR_A['M']*MSUN/(4/3*np.pi*(PAIR_A['R']*RSUN)**3)
PAIR_A['Reff'] = PAIR_A['R']*np.sqrt(979/(979*10**(-0.4*(PAIR_A['T'] - 13.0135))))
LR = jload(f'{R2}/out/r01_loc_lr.json')['summary']['pairA']
BASE = dict(ds='ref_beta', ecc='kipping', alpha=0.0,
            f_sub=0.15, f_N=0.003, f_S=0.0005, f_J=0.0015, f_BD=0.001,
            beta_sub=0.0, beta_N=0.0, beta_S=1.0, beta_J=1.5, beta_BD=0.0,
            f_mult=0.27, mu_loga=1.30, sig_loga=1.16, gamma_q=0.0, F_twin=0.10, S_close=0.25, S_BD=False,
            cmd='colour', cmd_Ft=0.1, cmd_sb='marg', sig_s=0.30, xobs=None,
            ipd='L_ipd', ast='base', dec='graded', wide=True, dT_twin=0.75,
            r_wd=2e-5, f_Twd=0.3, p_ecl_wd=0.029, dchi2_wd=10.0, f_occ=0.1, dchi2_heb=27.5, r_close=0.024, p_ecl_heb=0.06,
            lr_pairA=10**LR['log10LR_t']['base'], bg_post=3e-5, hri=None, chroma=None)
def evaluate(par, detail=False):
    par = dict(par)
    tab = build(par['ds'], par['ecc'], par['alpha'])
    # configuration weights
    Lc = np.ones(NC)
    if par['ipd'] != 'off': Lc *= CP[par['ipd']]
    Lc *= {'base': CP['L_ast'], 'soft': CP['L_ast']**0.4, 'off': 1.0}[par['ast']]
    if par['dec'] == 'graded': Lc *= CP['L_dec']
    elif par['dec'] == 'flat0.02': Lc *= np.where(CP['L_dec'] < 1, 0.02, 1.0)
    if par['wide']: Lc *= CP['L_wide']
    if par['hri'] is not None:
        hs, hv = map(np.array, par['hri']); lim = np.where(CP['sep'] >= hs[0], np.interp(CP['sep'], hs, hv), -99)
        Lc *= np.where(CP['dI'] < lim, 0.01, 1.0)
    # q and a re-weighting
    q = CP['q']; qmin = 0.075/CP['m1']; gq = par['gamma_q']
    normq = (1 - qmin**(gq + 1))/(gq + 1)
    wq = ((1 - par['F_twin'])*q**gq/normq + par['F_twin']*(q >= 0.95)/0.05)*(1 - qmin)
    lg = CP['loga']
    wa = (np.exp(-0.5*((lg - par['mu_loga'])/par['sig_loga'])**2)/par['sig_loga'])/(np.exp(-0.5*((lg - 1.30)/1.16)**2)/1.16)
    Lcmd = cmd_like(par, DX_C, DR_C); Lcmd_single = float(np.atleast_1d(cmd_like(par, np.zeros(1), np.zeros(1)))[0])
    fm = par['f_mult']
    wl = wq*wa*Lc*Lcmd
    Zs = (1 - fm)*Lcmd_single; Zb = fm*np.mean(wl); Zcfg = Zs + Zb
    Sa = np.where(CP['a'] <= 10, par['S_close'], np.where(CP['a'] >= 200, 1.0, par['S_close'] + (1 - par['S_close'])*np.log(CP['a']/10)/np.log(20)))
    occ = lambda cl, M: par['f_' + cl]/LNP*(np.asarray(M)/0.37)**par['beta_' + cl]
    PL = ('sub', 'N', 'S', 'J')
    rho_tp = rng.normal(9.42, 0.88, 300)
    Ztp = {cl: float(np.mean(Zof(tab, cl, 1.0, rho_tp, 0.381))) for cl in CLASSES}
    E = {}
    E['TP_single'] = Zs*sum(occ(cl, 0.37)*Ztp[cl] for cl in PL)
    E['BD_on_M3V_single'] = Zs*occ('BD', 0.37)*Ztp['BD']
    key = (par['ds'], par['ecc'], par['alpha'])
    if key not in ZCACHE:
        ZCACHE[key] = ({cl: Zof(tab, cl, 1 - CP['Dc'], CP['rho1'], CP['R1']).astype(np.float32) for cl in CLASSES},
                       {cl: Zof(tab, cl, CP['Dc'], CP['rho2'], CP['R2']).astype(np.float32) for cl in CLASSES})
    Z1, Z2 = ZCACHE[key]
    chrom = lambda ratio: 1.0 if par['chroma'] is None else np.exp(-0.5*((ratio - 1)/par['chroma'])**2)
    r1 = (1 - CP['DcV'])/(1 - CP['DcI']); r2 = CP['DcV']/CP['DcI']
    p1 = sum(occ(cl, CP['m1'])*Z1[cl] for cl in PL)*Sa*CP['stab1']*chrom(r1)
    p2 = sum(occ(cl, CP['m2'])*Z2[cl] for cl in PL)*Sa*CP['stab2']*chrom(r2)
    SBD = Sa if par['S_BD'] else 1.0
    b1 = occ('BD', CP['m1'])*Z1['BD']*SBD*CP['stab1']*chrom(r1); b2 = occ('BD', CP['m2'])*Z2['BD']*SBD*CP['stab2']*chrom(r2)
    twin = CP['dT'] <= par['dT_twin']
    E['TP_diluted'] = fm*np.mean(wl*p1)
    E['TP_on_twin'] = fm*np.mean(wl*p2*twin)
    E['PTP_companion'] = fm*np.mean(wl*p2*~twin)
    E['BD_on_primary_of_binary'] = fm*np.mean(wl*b1)
    E['BDHEB_on_companion'] = fm*np.mean(wl*b2)
    ex = lambda d: np.exp(-0.5*d)
    nt = CP['dT'] > 0.3
    # EB terms in the same units: rate per ln P x p_ecl x Occam x exp(-dchi2/2); the planet evidences carry
    # exp(-(chi2 - chi2_ref)/2) with chi2_ref = best planet fit, so dchi2 (relative to the free planet fit) is consistent
    E['EB_WD_on_M3V'] = (Zs + fm*np.mean(wl*nt))*par['r_wd']*par['f_Twd']*par['p_ecl_wd']*par['f_occ']*ex(par['dchi2_wd'])
    E['HEB_MSMS'] = fm*np.mean(wl*(CP['dT'] < 4.5))*par['r_close']*par['p_ecl_heb']*par['f_occ']*ex(par['dchi2_heb'])
    Spair = par['S_close'] + (1 - par['S_close'])*np.log(53/10)/np.log(20)
    ZA = sum(occ(cl, PAIR_A['M'])*float(np.mean(Zof(tab, cl, 1.0, PAIR_A['rho'], PAIR_A['Reff']))) for cl in PL)
    E['NTP_pairA'] = Zcfg*Spair*ZA*par['lr_pairA']
    tot0 = sum(E.values())
    E['BG_all'] = par['bg_post']*tot0       # background-star scenarios (analyst: 3e-5 posterior; negligible)
    tot = sum(E.values())
    planet = E['TP_single'] + E['TP_diluted'] + E['TP_on_twin']
    res = dict(FPP=1 - planet/tot, NFPP=E['NTP_pairA']/tot, FPP_strict=1 - (E['TP_single'] + E['TP_diluted'])/tot,
               FPP_system=1 - (planet + E['PTP_companion'])/tot, post={k: v/tot for k, v in E.items()})
    fps = {k: v for k, v in res['post'].items() if not k.startswith('TP')}
    res['dominant_FP'] = max(fps, key=fps.get)
    if detail:
        res['config'] = dict(P_single=Zs/Zcfg, P_binary=Zb/Zcfg, P_dT_le_0p3=fm*np.mean(wl*(CP['dT'] <= 0.3))/Zcfg,
                             P_dT_0p3_0p75=fm*np.mean(wl*(CP['dT'] > 0.3)*(CP['dT'] <= 0.75))/Zcfg,
                             P_dT_0p75_2=fm*np.mean(wl*(CP['dT'] > 0.75)*(CP['dT'] <= 2))/Zcfg,
                             P_dT_gt2=fm*np.mean(wl*(CP['dT'] > 2))/Zcfg,
                             CMD_single_over_twin=float(Lcmd_single/np.atleast_1d(cmd_like(par, np.array([-0.753]), np.array([0.0])))[0]))
        ptp = wl*p2*~twin
        res['PTP_by_a'] = {f'{lo}-{hi}': float(fm*np.mean(ptp*(CP['a'] > lo)*(CP['a'] <= hi))/tot) for lo, hi in [(0, 4), (4, 6), (6, 65), (65, 1e7)]}
        res['PTP_by_dT'] = {f'{lo}-{hi}': float(fm*np.mean(ptp*(CP['dT'] > lo)*(CP['dT'] <= hi))/tot) for lo, hi in [(0.75, 1.5), (1.5, 2.5), (2.5, 9)]}
        bd = wl*b2
        res['BDHEB_by_dT'] = {f'{lo}-{hi}': float(fm*np.mean(bd*(CP['dT'] > lo)*(CP['dT'] <= hi))/tot) for lo, hi in [(0, 1), (1, 2), (2, 3), (3, 9)]}
        wsum = lambda w_, x_: float(np.sum(w_*x_)/np.sum(w_))
        res['PTP_mean'] = dict(m2=wsum(ptp, CP['m2']), dT=wsum(ptp, CP['dT']), sep=wsum(ptp, CP['sep']), a=wsum(ptp, CP['a']))
        res['BDHEB_mean'] = dict(m2=wsum(bd, CP['m2']), dT=wsum(bd, CP['dT']), sep=wsum(bd, CP['sep']))
        res['Ztp_by_class'] = Ztp
    return res

if __name__ == '__main__':
    base = evaluate(BASE, detail=True)
    print(json.dumps(base, indent=1, default=float))
    jdump(dict(parameters={k: v for k, v in BASE.items()}, result=base, method=__doc__), 'out/r05_baseline.json')
    V = []
    def var(name, **kw):
        p = dict(BASE); p.update(kw); r = evaluate(p)
        V.append(dict(variant=name, FPP=r['FPP'], NFPP=r['NFPP'], FPP_strict=r['FPP_strict'], FPP_system=r['FPP_system'],
                      dominant=r['dominant_FP'], P_PTP=r['post']['PTP_companion'], P_BDHEB=r['post']['BDHEB_on_companion'],
                      P_TP_single=r['post']['TP_single'], P_twin=r['post']['TP_on_twin']))
        print(f"{name:55s} FPP {r['FPP']:.4f} strict {r['FPP_strict']:.3f} NFPP {r['NFPP']:.1e} PTP {r['post']['PTP_companion']:.4f} BDHEB {r['post']['BDHEB_on_companion']:.1e} dom {r['dominant_FP']}", flush=True)
    var('baseline')
    for d_ in ('inv_beta', 'ref_white', 'inv_white'): var(f'data {d_}', ds=d_)
    var('e Van Eylen', ecc='vaneylen'); var('e circular', ecc='circ')
    for a_ in (2.0, 3.0, 4.0): var(f'radius prior R^-{a_:g} (4-10 R_e)', alpha=a_)
    for v in (-1.0, 1.0, 2.0): var(f'beta_N={v}', beta_N=v)
    for v in (0.2, 0.35, 0.45): var(f'f_mult={v}', f_mult=v)
    for v in (1.0, 1.5): var(f'mu_loga={v}', mu_loga=v)
    for v in (0.0, 0.25): var(f'F_twin={v}', F_twin=v)
    for v in (-0.5, 1.0): var(f'gamma_q={v}', gamma_q=v)
    for v in (0.05, 0.15, 0.5, 1.0): var(f'S_close={v}', S_close=v)
    var('CMD x only (no G-Ks colour term), profile-marginalised', cmd='xonly')
    for ft in (0.0, 0.25): var(f'CMD colour model fitted with F_twin={ft}', cmd_Ft=ft)
    for sb in (0.20, 0.30, 0.40): var(f'CMD colour model at fixed sig_b={sb}', cmd_sb=sb)
    for sb in (0.20, 0.30): var(f'CMD x only at fixed sig_b={sb}', cmd='xonly', cmd_sb=sb)
    for s_ in (0.25, 0.30, 0.40): var(f'CMD Gaussian sig_s={s_} (analyst form)', cmd='gauss', sig_s=s_)
    var('CMD analyst form sig 0.30 at x_obs=-0.753 (analyst baseline CMD)', cmd='gauss', sig_s=0.30, xobs=-0.753)
    var('CMD ignored (flat)', cmd='flat')
    var('x_obs = -0.753 (referee ridge)', xobs=-0.753)
    var('IPD LR84', ipd='L_ipd84'); var('IPD red-star cells', ipd='L_ipd_red'); var('IPD only beyond 0.25"', ipd='L_ipd_far'); var('IPD off', ipd='off')
    var('astrometry soft', ast='soft'); var('astrometry off', ast='off')
    var('DECaPS2 flat 0.02 (analyst)', dec='flat0.02'); var('DECaPS2 off', dec='off'); var('Gaia-wide off', wide=False)
    var('all archival constraints off', ipd='off', ast='off', dec='off', wide=False)
    for v in (0.5, 1.0): var(f'dT_twin={v}', dT_twin=v)
    for v in (0.0, 3e-4, 3e-3): var(f'f_BD={v}', f_BD=v)
    var('BD suppressed in close binaries like planets', S_BD=True)
    for v in (6.0, 13.2): var(f'dchi2_wd={v}', dchi2_wd=v)
    for v in (2e-6, 2e-4): var(f'r_wd={v}', r_wd=v)
    var('dchi2_heb=20', dchi2_heb=20.0)
    for nm_, v in (('Gauss base', LR['log10LR_gauss']['base']), ('t median of 19 variants', LR['log10LR_t']['median']),
                   ('sys-inflated, worst variant', LR['sysinfl']['log10LR_sysinfl_at_min_variant'])):
        var(f'pair A LR {nm_} (10^{v:.1f})', lr_pairA=10**v)
    var('bg scenarios x10', bg_post=3e-4)
    pess = dict(beta_N=-1.0, f_mult=0.40, F_twin=0.0, S_close=0.5, ipd='L_ipd84', ast='soft', ecc='vaneylen', f_BD=3e-3,
                lr_pairA=10**-2.8, dchi2_wd=6.0, r_wd=2e-4)
    var('PESSIMISTIC realistic', **pess)
    var('PESSIMISTIC + IPD only beyond 0.25"', **dict(pess, ipd='L_ipd_far'))
    opt = dict(beta_N=1.0, alpha=3.0, S_close=0.1, f_BD=3e-4, lr_pairA=10**-6.7)
    var('OPTIMISTIC realistic', **opt)
    HRI = ([0.03, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0], [2.0, 3.0, 4.0, 4.5, 5.0, 5.5, 6.0])
    var('WHAT-IF speckle null (8-m 832 nm curve)', hri=HRI)
    var('WHAT-IF two-band transit depth ratio = 1 +- 5% (+ pair resolved)', chroma=0.05, lr_pairA=1e-8)
    var('WHAT-IF both', hri=HRI, chroma=0.05, lr_pairA=1e-8)
    dfv = pd.DataFrame(V); dfv.to_csv(f'{R2}/out/r05_variants.csv', index=False)
    # joint random draws
    def lu(lo, hi): return float(np.exp(rng.uniform(np.log(lo), np.log(hi))))
    R = []
    for i in range(1500):
        p = dict(BASE)
        p.update(ds=str(rng.choice(['ref_beta', 'inv_beta', 'ref_white', 'inv_white'], p=[0.4, 0.4, 0.1, 0.1])),
                 ecc=str(rng.choice(['kipping', 'vaneylen', 'circ'])), alpha=float(rng.choice([0.0, 2.0, 3.0])),
                 f_N=lu(0.001, 0.01), f_S=lu(0.0001, 0.003), f_J=lu(0.0005, 0.005), f_BD=lu(3e-4, 3e-3),
                 beta_N=float(rng.uniform(-1, 2)), beta_S=float(rng.uniform(0, 2)), beta_J=float(rng.uniform(0, 3)),
                 f_mult=float(rng.uniform(0.22, 0.42)), mu_loga=float(rng.uniform(1.0, 1.5)), sig_loga=float(rng.uniform(0.9, 1.4)),
                 gamma_q=float(rng.uniform(-0.5, 1.0)), F_twin=float(rng.uniform(0, 0.3)), S_close=lu(0.05, 1.0),
                 cmd=str(rng.choice(['colour', 'colour', 'xonly'])), cmd_Ft=float(rng.choice([0.0, 0.1, 0.25])),
                 cmd_sb=['marg', 'marg', 0.30, 0.35, 0.40, 0.45][int(rng.integers(6))], xobs=float(rng.uniform(-0.85, -0.78)),
                 ipd=str(rng.choice(['L_ipd', 'L_ipd84', 'L_ipd_red', 'L_ipd_far'])), ast=str(rng.choice(['base', 'soft'])),
                 dT_twin=float(rng.uniform(0.5, 1.0)), r_wd=lu(2e-6, 2e-4), dchi2_wd=float(rng.uniform(6, 13)),
                 f_occ=lu(0.03, 0.3), dchi2_heb=float(rng.uniform(20, 33.7)),
                 lr_pairA=10**float(rng.uniform(-9.0, -2.8)), S_BD=bool(rng.random() < 0.5))
        r = evaluate(p)
        R.append(dict(FPP=r['FPP'], NFPP=r['NFPP'], FPP_strict=r['FPP_strict'], FPP_system=r['FPP_system'], dominant=r['dominant_FP'],
                      P_PTP=r['post']['PTP_companion'], P_BDHEB=r['post']['BDHEB_on_companion'],
                      **{k: p[k] for k in ('ds', 'ecc', 'alpha', 'beta_N', 'f_mult', 'S_close', 'cmd', 'cmd_Ft', 'ipd', 'ast', 'dT_twin', 'f_BD', 'F_twin', 'xobs')},
                      log10_lr_pairA=float(np.log10(p['lr_pairA']))))
    dfr = pd.DataFrame(R); dfr.to_csv(f'{R2}/out/r05_random.csv', index=False)
    qq = lambda x, p_: float(np.percentile(x, p_))
    summ = dict(n=len(dfr), FPP_pct_5_16_50_84_95=[qq(dfr.FPP, p_) for p_ in (5, 16, 50, 84, 95)],
                frac_FPP_lt_0p015=float(np.mean(dfr.FPP < 0.015)), NFPP_pct_5_50_95=[qq(dfr.NFPP, p_) for p_ in (5, 50, 95)],
                frac_NFPP_lt_1em3=float(np.mean(dfr.NFPP < 1e-3)), FPP_strict_pct_5_50_95=[qq(dfr.FPP_strict, p_) for p_ in (5, 50, 95)],
                FPP_system_pct_5_50_95=[qq(dfr.FPP_system, p_) for p_ in (5, 50, 95)], dominant=dfr.dominant.value_counts().to_dict())
    from scipy.stats import spearmanr
    summ['spearman_FPP'] = {c: float(spearmanr(dfr[c], dfr.FPP).correlation) for c in ('alpha', 'beta_N', 'f_mult', 'S_close', 'dT_twin', 'f_BD', 'F_twin', 'xobs', 'log10_lr_pairA')}
    summ['FPP_median_by_ipd'] = dfr.groupby('ipd').FPP.median().to_dict()
    summ['FPP_median_by_ds'] = dfr.groupby('ds').FPP.median().to_dict()
    summ['FPP_median_by_cmd'] = dfr.groupby('cmd').FPP.median().to_dict()
    print(json.dumps(summ, indent=1, default=float))
    jdump(dict(named_variants=dict(n=len(dfv), n_FPP_ge_0p015=int((dfv.FPP >= 0.015).sum()), min=float(dfv.FPP.min()), max=float(dfv.FPP.max())),
               random=summ), 'out/r05_summary.json')
