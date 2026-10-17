"""s06: own marginal-likelihood (evidence) ratio of a BD on the M3V versus a Neptune on the M3V, from my own likelihood
grid of the CLEANED light curves (lc3; beta and white noise), independent of the r03/r05 machinery.
Grid: k (planet 0.095-0.175; BD 0.200-0.275; diluted-host bodies 0.15-0.45), b uniform in [0, 1+k) (64 steps), rho_c log
grid 0.5-80 g/cm3 (46); D = 1 for the M3V hosts and D in {0.12, 0.18, 0.25, 0.35} for a companion host; LD nominal;
dT0 fixed at the best planet fit; per-window linear baselines profiled.
Evidence for a body class on a host of density rho_h and radius R_h:  Z = sum_k p(k) sum_rho W(rho) sum_b L(k,b,rho) db
with W(rho bin) = E[(1+e sin w)/(1-e^2) 1(rho_h g^3 in bin)] (g = (1+e sin w)/sqrt(1-e^2); correct photo-eccentric sign),
and the geometric factor R_h/a common to classes on the same host. L = exp(-(chi2 - chi2_best_planet)/2).
Priors: Neptunes flat in ln R over 4-7 R_E; BDs R ~ N(0.88, 0.08) R_J (>= 0.72); M3V rho ~ N(9.42, 0.88).
e priors: Kipping13, Van Eylen singles, Sagear & Ballard 23 Rayleigh 0.19 (planets); Ferreira+26 short/long-period
Beta (BDs). Output ratio R = Z_BD/Z_TP, so P(BD on M3V)/P(planet on M3V) = (f_BD/f_N) R.
Also: a representative PTP companion (0.20 Msun, R 0.22 Rsun, rho 26.5, D 0.25) - evidence of a 4-10 R_E planet on it
relative to the M3V planet, per unit configuration prior, for each planet e prior (tests the e-dependence of PTP).
Output: out/s06_evidence.json"""
import sys; sys.path.insert(0, 'tic214311908/analysis/referee3_bd/scripts')
from lc3 import *
import json, time
from scipy.stats import norm
t00 = time.time()
KP = np.round(np.arange(0.095, 0.1751, 0.005), 4); KB = np.round(np.arange(0.200, 0.2751, 0.015), 4); KC = np.round(np.arange(0.15, 0.451, 0.025), 4)
NB = 64
RHO = np.exp(np.linspace(np.log(0.5), np.log(80.), 46)); dl = np.log(RHO[1]/RHO[0])
DC = np.array([0.12, 0.18, 0.25, 0.35])
fr = json.load(open(f'{R3}/out/s03_fit_real.json'))['fits']
def grid(lc, dT0):
    out = {}
    for name, ks, Ds in (('P', KP, [1.0]), ('BD', KB, [1.0]), ('C', KC, list(DC))):
        G2 = np.full((len(Ds), len(ks), NB, len(RHO)), np.nan)
        for ik, k in enumerate(ks):
            bs = (np.arange(NB) + 0.5)/NB*(1 + k)
            S = []
            for b in bs:
                for rho in RHO:
                    aR = float(aR_of_rho(rho))
                    S.append(lc.template(binned(k, b, aR), dT0) if b < aR else np.zeros(len(lc.t)))
            S = np.column_stack(S)
            Sp = lc.perp(S)
            for iD, D in enumerate(Ds):
                c2 = np.sum(lc.w[:, None]*(lc.yp[:, None] - D*Sp)**2, axis=0)
                G2[iD, ik] = c2.reshape(NB, len(RHO))
        out[name] = G2
    return out
rng = np.random.default_rng(5)
NE = 400000
def edr(kind):
    if kind == 'kipping': return rng.beta(0.867, 3.03, NE)
    if kind == 'vaneylen': return np.abs(rng.normal(0, 0.32, NE))
    if kind == 'rayleigh19': return rng.rayleigh(0.19, NE)
    if kind == 'circ': return np.zeros(NE)
    if kind == 'bdsp':
        a = np.clip(rng.normal(0.426, 0.092, NE), 0.05, None); b = np.clip(rng.normal(2.354, 0.73, NE), 0.1, None); return rng.beta(a, b)
    if kind == 'bdlp':
        a = np.clip(rng.normal(1.879, 0.92, NE), 0.05, None); b = np.clip(rng.normal(2.470, 1.25, NE), 0.1, None); return rng.beta(a, b)
def Wrho(kind, rho_h_draws):
    e = np.clip(edr(kind), 0, 0.9); w = rng.uniform(0, 2*np.pi, NE)
    g = (1 + e*np.sin(w))/np.sqrt(1 - e**2); geo = (1 + e*np.sin(w))/(1 - e**2)
    rc = rho_h_draws*g**3
    edges = np.exp(np.r_[np.log(RHO) - dl/2, np.log(RHO[-1]) + dl/2])
    W, _ = np.histogram(rc, bins=edges, weights=geo)
    return W/NE, float(np.mean(geo*((rc < edges[0]) | (rc > edges[-1]))))
res = {}
for ds in ('ref', 'inv'):
    for nm in ('beta', 'white'):
        t0 = time.time()
        lc = load_ref(beta=(nm == 'beta')) if ds == 'ref' else load_inv(beta=(nm == 'beta'))
        F = fr[f'{ds}_{nm}']
        dT0 = F['planet']['dT0_min']/1440.
        Gd = grid(lc, dT0)
        c2ref = float(np.nanmin(Gd['P']))
        L = {k_: np.exp(-(v - c2ref)/2) for k_, v in Gd.items()}
        db = lambda ks: ((1 + np.asarray(ks))/NB)[:, None]
        # sum over b -> [D, k, rho]
        Jb = {'P': np.sum(L['P'][0], axis=1)*db(KP), 'BD': np.sum(L['BD'][0], axis=1)*db(KB), 'C': np.stack([np.sum(L['C'][i], axis=1)*db(KC) for i in range(len(DC))])}
        # k priors
        R_P = KP*RSTAR*RSUN/REARTH; dlk = np.log(KP[1]/KP[0])
        pk_P = np.where((R_P >= 4) & (R_P <= 7), dlk/np.log(7/4.), 0.)
        RBD = KB*RSTAR*RSUN/RJUP
        pk_BD = norm.pdf(RBD, 0.88, 0.08)*(RBD >= 0.72); pk_BD = pk_BD/np.sum(pk_BD)
        rho_m3v = rng.normal(9.42, 0.88, NE)
        R = dict(chi2_best_planet_grid=c2ref, chi2_best_BD_grid=float(np.nanmin(Gd['BD'])), dchi2_BD_grid=float(np.nanmin(Gd['BD']) - c2ref))
        Z = {}
        for kind in ('kipping', 'vaneylen', 'rayleigh19', 'circ', 'bdsp', 'bdlp'):
            W, wout = Wrho(kind, rho_m3v)
            Z[f'TP_{kind}'] = float(np.sum(pk_P[:, None]*Jb['P']*W[None, :]))
            Z[f'BD_{kind}'] = float(np.sum(pk_BD[:, None]*Jb['BD']*W[None, :]))
            Z[f'W_out_of_grid_{kind}'] = wout
        R['Z'] = Z
        R['ratio_BD_over_TP'] = {f'BD e {be} / TP e {pe}': Z[f'BD_{be}']/Z[f'TP_{pe}'] for be in ('kipping', 'bdsp', 'bdlp') for pe in ('kipping', 'vaneylen', 'rayleigh19')}
        # where does the BD evidence come from? fraction of Z_BD from rho_c < 3
        W, _ = Wrho('bdsp', rho_m3v)
        zb = pk_BD[:, None]*Jb['BD']*W[None, :]
        R['BD_bdsp_frac_from_rho_lt_3'] = float(zb[:, RHO < 3].sum()/zb.sum())
        # representative PTP companion: 0.20 Msun, R 0.22 Rsun, D = 0.25; planets 4-10 R_E flat in ln R
        Rc, rhoc = 0.22, 0.20*MSUN/(4/3*np.pi*(0.22*RSUN)**3)
        R_C = KC*Rc*RSUN/REARTH; dlkc = np.log(KC[1]/KC[0])
        pk_C = np.where((R_C >= 4) & (R_C <= 10), dlkc/np.log(10/4.), 0.)
        rho_c_draws = rng.normal(rhoc, 0.1*rhoc, NE)
        iD = int(np.argmin(np.abs(DC - 0.25)))
        aR_m3v = float(aR_of_rho(9.42)); aR_c = float(aR_of_rho(rhoc))
        ptp = {}
        for kind in ('kipping', 'vaneylen', 'rayleigh19', 'circ'):
            W, _ = Wrho(kind, rho_c_draws)
            Zc = float(np.sum(pk_C[:, None]*Jb['C'][iD]*W[None, :]))
            ptp[kind] = dict(Z_C=Zc, ratio_vs_TP_same_e=Zc/Z[f'TP_{kind}']*(aR_m3v/aR_c))   # includes R/a geometric factor ratio
            for iDD, DD in enumerate(DC):
                ptp[kind][f'Z_C_D{DD}'] = float(np.sum(pk_C[:, None]*Jb['C'][iDD]*W[None, :]))
        R['PTP_representative'] = dict(rho_host=rhoc, D=float(DC[iD]), by_planet_e=ptp,
                                       note='ratio_vs_TP_same_e = [Z(planet 4-10 R_E on companion) (R_c/a_c)] / [Z(Neptune on M3V) (R*/a)], per unit configuration prior and equal per-star occurrence')
        R['chi2_min_C_by_D'] = {f'{DD}': float(np.nanmin(Gd['C'][i]) - c2ref) for i, DD in enumerate(DC)}
        res[f'{ds}_{nm}'] = R
        print(ds, nm, 'dchi2 BD grid', round(R['dchi2_BD_grid'], 2), 'ratios', {k_: f'{v:.2e}' for k_, v in R['ratio_BD_over_TP'].items()}, 'BD frac rho<3', round(R['BD_bdsp_frac_from_rho_lt_3'], 3), round(time.time() - t0), flush=True)
        print('   PTP rep ratio', {k_: f"{v['ratio_vs_TP_same_e']:.2e}" for k_, v in ptp.items()}, 'C dchi2 by D', {k_: round(v, 1) for k_, v in R['chi2_min_C_by_D'].items()}, flush=True)
        jdump(dict(method=__doc__, results=res, grids=dict(KP=KP, KB=KB, KC=KC, RHO=RHO, NB=NB, DC=DC)), f'{R3}/out/s06_evidence.json')
print('done', round(time.time() - t00))
