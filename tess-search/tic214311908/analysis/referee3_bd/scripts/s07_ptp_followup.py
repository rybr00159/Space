"""s07: the PTP residual (planet on a fainter bound companion) in the sign-corrected referee2_fpp framework, and what
follow-up removes it.
(1) Imports r05_fpp.py READ-ONLY with my sign patch (rho_c = rho_host g^+3; 10,000 e draws), for planet e priors
    Kipping13 and Sagear & Ballard 23 Rayleigh 0.19; f_BD 1.3e-4 with the Ferreira+26 short-period e law for BDs.
(2) Posterior-weighted distribution of the PTP companions (r04 Monte Carlo draws x r05 weights): separation, a, mass,
    dT, dI (TESS/I contrast), dKs, and the radius of the planet needed.
(3) Fraction of the PTP posterior that each imaging set-up would detect (companion brighter than the 5-sigma
    contrast curve at its separation): Gemini/Zorro-like 8-m speckle at 832 nm (curve used by referee2_fpp r05:
    0.03/0.05/0.1/0.2/0.5/1/2" -> 2/3/4/4.5/5/5.5/6 mag), a 4-m SOAR/HRCam-like I-band speckle (0.05/0.1/0.15/0.5/1" ->
    0.5/1.5/2.5/4/4.5 mag; my generic assumption), an 8-m K-band AO (0.05/0.07/0.1/0.2/0.5/1" -> 1/2/3/4.5/6/7 mag in Ks;
    generic assumption).
(4) FPP what-ifs with the fixed sign: speckle null, two-band transit depth ratio 1 +- 5% (+ pair A resolved), both,
    both + BD channels removed by RV.
(5) RV trend from the PTP companion: |dv/dt| of the M3V for the posterior-weighted companions (random orbital phase,
    e from r04, isotropic orientation) -> fraction with |dv/dt| > 10, 30 m/s/yr.
Output: out/s07_ptp_followup.json"""
import sys, os
sys.path.insert(0, 'tic214311908/analysis/referee3_bd/scripts')
from rlib3 import *
sys.path.insert(0, f'{B}/referee2_fpp/scripts')
cwd = os.getcwd(); os.chdir(f'{B}/referee2_fpp')
import r05_fpp as R
os.chdir(cwd)
rng0 = np.random.default_rng(7)
def edraw(kind, n, rng):
    if kind == 'kipping': return rng.beta(0.867, 3.03, n)
    if kind == 'rayleigh19': return rng.rayleigh(0.19, n)
    if kind == 'bdsp':
        a = np.clip(rng.normal(0.426, 0.092, n), 0.05, None); b = np.clip(rng.normal(2.354, 0.73, n), 0.1, None); return rng.beta(a, b)
def ecc_fixed(kind, n=2500):
    rng = np.random.default_rng(12345); nn = 10000
    e = np.clip(edraw(kind, nn, rng), 0, 0.9); w = rng.uniform(0, 2*np.pi, nn)
    g = (1 + e*np.sin(w))/np.sqrt(1 - e**2); geo = (1 + e*np.sin(w))/(1 - e**2)
    return 1/g, geo
orig_build = R.build
def build_patched(ds, kind, alpha):
    if '+' in kind:
        kp, kb = kind.split('+'); tp = orig_build(ds, kp, alpha); tb = orig_build(ds, kb, alpha)
        out = dict(tp); out['BD'] = tb['BD']; return out
    return orig_build(ds, kind, alpha)
R.ecc = ecc_fixed; R.build = build_patched
CP = R.CP
def ptp_weights(par):
    """re-computes the per-companion PTP weight exactly as r05.evaluate does (wl * p2 * ~twin)"""
    tab = R.build(par['ds'], par['ecc'], par['alpha'])
    Lc = np.ones(R.NC); Lc *= CP[par['ipd']]; Lc *= CP['L_ast']; Lc *= CP['L_dec']; Lc *= CP['L_wide']
    q = CP['q']; qmin = 0.075/CP['m1']; gq = par['gamma_q']; normq = (1 - qmin**(gq + 1))/(gq + 1)
    wq = ((1 - par['F_twin'])*q**gq/normq + par['F_twin']*(q >= 0.95)/0.05)*(1 - qmin)
    lg = CP['loga']
    wa = (np.exp(-0.5*((lg - par['mu_loga'])/par['sig_loga'])**2)/par['sig_loga'])/(np.exp(-0.5*((lg - 1.30)/1.16)**2)/1.16)
    Lcmd = R.cmd_like(par, R.DX_C, R.DR_C)
    wl = wq*wa*Lc*Lcmd
    Sa = np.where(CP['a'] <= 10, par['S_close'], np.where(CP['a'] >= 200, 1.0, par['S_close'] + (1 - par['S_close'])*np.log(CP['a']/10)/np.log(20)))
    occ = lambda cl, M: par['f_' + cl]/R.LNP*(np.asarray(M)/0.37)**par['beta_' + cl]
    Z2 = {cl: R.Zof(tab, cl, CP['Dc'], CP['rho2'], CP['R2']) for cl in ('sub', 'N', 'S', 'J', 'BD')}
    p2 = sum(occ(cl, CP['m2'])*Z2[cl] for cl in ('sub', 'N', 'S', 'J'))*Sa*CP['stab2']
    twin = CP['dT'] <= par['dT_twin']
    wcls = {cl: occ(cl, CP['m2'])*Z2[cl]*Sa*CP['stab2']*~twin*wl for cl in ('sub', 'N', 'S', 'J')}
    return wl*p2*~twin, wcls, wl*occ('BD', CP['m2'])*Z2['BD']*CP['stab2']
def wq_(x, w, ps=(5, 16, 50, 84, 95)):
    o = np.argsort(x); c = np.cumsum(w[o]); c /= c[-1]
    return [float(np.interp(p/100, c, x[o])) for p in ps]
CURVES = {'speckle_8m_832nm': ([0.03, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0], [2.0, 3.0, 4.0, 4.5, 5.0, 5.5, 6.0], 'dI'),
          'speckle_4m_I': ([0.05, 0.1, 0.15, 0.5, 1.0], [0.5, 1.5, 2.5, 4.0, 4.5], 'dI'),
          'AO_8m_Ks': ([0.05, 0.07, 0.1, 0.2, 0.5, 1.0], [1.0, 2.0, 3.0, 4.5, 6.0, 7.0], 'dKs')}
def detected(curve):
    hs, hv, band = curve
    lim = np.where(CP['sep'] >= hs[0], np.interp(CP['sep'], hs, hv), -99)
    return CP[band] < lim
out = {}
for ds in ('ref_beta', 'inv_beta'):
    for pe in ('kipping', 'rayleigh19'):
        par = dict(R.BASE); par.update(ds=ds, ecc=f'{pe}+bdsp', f_BD=1.3e-4)
        R.rng = np.random.default_rng(11)
        base = R.evaluate(par, detail=True)
        w, wcls, wbd = ptp_weights(par)
        tag = f'{ds}|{pe}'
        o = dict(FPP=base['FPP'], FPP_system=base['FPP_system'], PTP=base['post']['PTP_companion'], BD_total=base['post']['BD_on_M3V_single'] + base['post']['BD_on_primary_of_binary'] + base['post']['BDHEB_on_companion'],
                 NTP=base['post']['NTP_pairA'], WD=base['post']['EB_WD_on_M3V'])
        o['PTP_class_shares'] = {cl: float(v.sum()/w.sum()) for cl, v in wcls.items()}
        for key in ('sep', 'a', 'm2', 'dT', 'dI', 'dKs', 'Dc', 'e'):
            o[f'PTP_{key}_5_16_50_84_95'] = wq_(CP[key], w)
        o['PTP_frac_sep_lt'] = {f'{s}': float(np.sum(w*(CP['sep'] < s))/w.sum()) for s in (0.03, 0.05, 0.1, 0.2)}
        o['PTP_detect_frac'] = {k: float(np.sum(w*detected(c))/w.sum()) for k, c in CURVES.items()}
        o['BDcomp_detect_frac'] = {k: float(np.sum(wbd*detected(c))/wbd.sum()) for k, c in CURVES.items()}
        # RV trend of the M3V due to the companion (posterior-weighted)
        n = len(CP['a']); rng = np.random.default_rng(3)
        M_ = rng.uniform(0, 2*np.pi, n); e_ = np.clip(CP['e'], 0, 0.95)
        E_ = M_.copy()
        for _ in range(30): E_ = E_ - (E_ - e_*np.sin(E_) - M_)/(1 - e_*np.cos(E_))
        r_ = CP['a']*(1 - e_*np.cos(E_))*AU
        acc = G*CP['m2']*MSUN/r_**2          # cm/s^2 of the M3V
        proj = np.abs(rng.normal(size=n))/np.sqrt(3)  # crude random line-of-sight projection
        cosi = rng.uniform(-1, 1, n); proj = np.abs(np.sqrt(1 - cosi**2)*np.sin(rng.uniform(0, 2*np.pi, n)))
        dvdt = acc*proj*YR/100.      # m/s/yr
        o['PTP_rv_trend_m_s_yr_5_16_50_84_95'] = wq_(dvdt, w)
        o['PTP_frac_trend_gt'] = {f'{x}': float(np.sum(w*(dvdt > x))/w.sum()) for x in (3, 10, 30, 100)}
        # what-ifs with the fixed sign
        HRI = CURVES['speckle_8m_832nm'][:2]
        wi = {}
        for nm_, kw in (('speckle null 8m 832nm', dict(hri=HRI)), ('two-band 5% + pair resolved', dict(chroma=0.05, lr_pairA=1e-8)),
                        ('speckle + two-band', dict(hri=HRI, chroma=0.05, lr_pairA=1e-8)),
                        ('speckle + two-band + BD excluded by RV', dict(hri=HRI, chroma=0.05, lr_pairA=1e-8, f_BD=0.0))):
            p2_ = dict(par); p2_.update(kw); R.rng = np.random.default_rng(11)
            r_ = R.evaluate(p2_)
            wi[nm_] = dict(FPP=r_['FPP'], PTP=r_['post']['PTP_companion'], dominant=r_['dominant_FP'])
        o['whatif'] = wi
        out[tag] = o
        print(tag, json.dumps({k: v for k, v in o.items()}, default=float)[:2500], flush=True)
jdump(dict(method=__doc__, results=out), f'{R3}/out/s07_ptp_followup.json')
