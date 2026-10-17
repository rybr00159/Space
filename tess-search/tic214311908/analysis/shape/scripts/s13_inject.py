"""s13: injection-recovery calibration of the shape analysis on the investigator's full light curve (INV).
Truth models are injected into the raw (not detrended) M3V flux at 10 fake ephemerides (phase offsets 0.14-0.86,
avoiding the real transits and phase 0.5), the series is re-detrended exactly as in ../phot/scripts/p06_lc.py
(wotan biweight 0.75 d, break tolerance 0.5 d, real and fake windows +-0.1 d masked, trend interpolated), and
windows of +-0.3 d around the fake epochs are fitted with the same models and noise weights as the real data.
Truths: 'planet' = the s03 INV 'single' max-likelihood planet (rho 9.4 g/cm3 host); 'wd' = the s06 INV 'single'
best white-dwarf-occultation model (if present); 'grazing' = the s05 best grazing fit with k = 0.27 on the M3V.
Fitted models per injection: planet free (rho free, LD prior), WD occultation by the M3V (s06 'single'
parameterisation, WD radius tied to its mass; light curve only), V-shape (planet model with flat part forced to zero, b = 1 - k... i.e. b >= 1 - k).
Output: out/inject.json, out/inject_rows.csv
usage: s13_inject.py [truths...]"""
import sys, json, numpy as np, pandas as pd
sys.path.insert(0, 'tic214311908/analysis/shape/scripts')
from shapelib import *
from wotan import flatten
from scipy.optimize import minimize
truths = sys.argv[1:] if len(sys.argv) > 1 else ['planet', 'wd']
z = dict(np.load(f'{SD}/data/inv_full.npz'))
ld = json.load(open(f'{SD}/out/ld_claret17.json'))['adopted_single']
fs = json.load(open(f'{SD}/out/fit_inv_single.json'))['best']
P = fs['P']; T0r = fs['T0']
def build_lc(tmodel, T0f):
    """inject tmodel (function of t returning flux in 979 units, centred on T0f) into raw flux; re-detrend; windows"""
    t = z['t']; raw = z['raw'].copy()
    raw = raw + (tmodel(t) - 1)*979/935
    f = np.full(len(t), np.nan)
    real = np.abs(((t - T0r)/P + 0.5) % 1 - 0.5)*P < 0.1
    fake = np.abs(((t - T0f)/P + 0.5) % 1 - 0.5)*P < 0.1
    for s in np.unique(z['sector']):
        m = z['sector'] == s
        use = m & ~real & ~fake & np.isfinite(raw)
        tr = flatten(t[use], raw[use], method='biweight', window_length=0.75, return_trend=True, break_tolerance=0.5, edge_cutoff=0)[1]
        ok = np.isfinite(tr)
        trend = np.interp(t[m], t[use][ok], tr[ok])
        f[m] = raw[m]/trend
    f = 1 + (f - 1)*935/979
    ep = np.round((t - T0f)/P); dt = t - (T0f + ep*P)
    w = (np.abs(dt) < 0.3) & np.isfinite(f) & ~(np.abs(((t - T0r)/P + 0.5) % 1 - 0.5)*P < 0.3)
    d = dict(t=t[w], f=f[w], sig=z['sig'][w], cad=z['cad'][w], sector=z['sector'][w], epoch=ep[w].astype(int))
    return LC(d)
def fit_planet(lc, T0f, rho_prior=None):
    def c2(x):
        T0, k, b, lr, u1, u2 = x
        if not (u1 + u2 < 1 and u1 + 2*u2 > 0): return 1e9
        aR = aR_from_rho(10**lr, P)
        v = lc.chi2(lc.transit(T0, P, k, aR, b, u1, u2)) + ((u1 - ld['u1'])/0.1)**2 + ((u2 - ld['u2'])/0.1)**2
        if rho_prior: v += ((10**lr - rho_prior[0])/rho_prior[1])**2
        return v
    bounds = [(T0f-0.01, T0f+0.01), (0.03, 0.45), (0, 1.4), (-1, 2.3), (-0.3, 0.8), (-0.2, 1.0)]
    best = None
    for b0 in (0.15, 0.55, 0.85):
        for lr0 in (0.8, 1.2):
            r = minimize(c2, [T0f, 0.12, b0, lr0, ld['u1'], ld['u2']], method='Powell', bounds=bounds, options=dict(xtol=1e-5, ftol=1e-6, maxfev=5000))
            if best is None or r.fun < best.fun: best = r
    T0, k, b, lr, u1, u2 = best.x
    aR = aR_from_rho(10**lr, P); T14, T23 = t14_t23(P, k, aR, b)
    return dict(chi2=float(best.fun), k=float(k), b=float(b), rho=float(10**lr), T14_h=float(T14*24), flat=float(T23/T14 if T14 > 0 else 0))
def fit_vshape(lc, T0f):
    """planet model forced to have no flat part: b >= 1 - k"""
    def c2(x):
        T0, k, bx, lr, u1, u2 = x
        b = (1 - k) + bx*(2*k)  # between 1-k (V) and 1+k
        if not (u1 + u2 < 1 and u1 + 2*u2 > 0): return 1e9
        aR = aR_from_rho(10**lr, P)
        return lc.chi2(lc.transit(T0, P, k, aR, b, u1, u2)) + ((u1 - ld['u1'])/0.1)**2 + ((u2 - ld['u2'])/0.1)**2
    bounds = [(T0f-0.01, T0f+0.01), (0.05, 0.9), (0, 0.6), (-1, 2.3), (-0.3, 0.8), (-0.2, 1.0)]
    best = None
    for k0 in (0.2, 0.4):
        for lr0 in (0.5, 1.0):
            r = minimize(c2, [T0f, k0, 0.05, lr0, ld['u1'], ld['u2']], method='Powell', bounds=bounds, options=dict(xtol=1e-5, ftol=1e-6, maxfev=5000))
            if best is None or r.fun < best.fun: best = r
    return dict(chi2=float(best.fun), k=float(best.x[1]))
def wd_model_factory(lc_or_t, T0, MA, RA, bB, fA, RB, MB, e=0., w=90.):
    a = (G*(MA+MB)*MSUN*(P*DAY)**2/(4*np.pi**2))**(1/3.)/RSUN
    k = RB/RA; aRA = a/RA
    return k, aRA, bB*k
def r_wd(MA, xr):
    mu = MA/1.44; rn = 0.0112*np.sqrt(mu**(-2/3.) - mu**(2/3.))
    return rn*(1 + xr*((2.5 if MA < 0.45 else 1.3) - 1))
def fit_wd(lc, T0f):
    def c2(x):
        T0, MA, xr, bB, fA, RB, MB = x
        k, aRA, bA = wd_model_factory(None, T0, MA, r_wd(MA, xr), bB, fA, RB, MB)
        m = 1 + fA*(lc.transit(T0, P, k, aRA, bA, 0.1, 0.2) - 1)
        return lc.chi2(m) + ((RB - 0.381)/0.012)**2 + ((MB - 0.370)/0.011)**2
    bounds = [(T0f-0.01, T0f+0.01), (0.15, 1.35), (0, 1), (0, 1.05), (0.005, 0.06), (0.37, 0.40), (0.34, 0.40)]
    best = None
    for MA0 in (0.25, 0.6):
        for bB0 in (0.1, 0.6):
            r = minimize(c2, [T0f, MA0, 0.3, bB0, 0.016, 0.381, 0.370], method='Powell', bounds=bounds, options=dict(xtol=1e-6, ftol=1e-6, maxfev=6000))
            if best is None or r.fun < best.fun: best = r
    return dict(chi2=float(best.fun), MA=float(best.x[1]), RA=float(r_wd(best.x[1], best.x[2])), bB=float(best.x[3]), fA=float(best.x[4]))
# truth models
def truth_fn(name, T0f):
    if name == 'planet':
        aR = aR_from_rho(10**fs['log10rho'], P)
        def fn(t):
            d = dict(t=t, f=np.ones(len(t)), sig=np.ones(len(t)), cad=z['cad'], sector=z['sector'], epoch=np.zeros(len(t), int))
            L = LC(d, use_offsets=False)
            return L.transit(T0f, P, fs['k'], aR, fs['b'], fs['u1'], fs['u2'])
        return fn
    if name == 'wd':
        wd = json.load(open(f'{SD}/out/compact_inv_single.json'))
        a = wd['a_Rsun']; k = wd['k']; aRA = wd['aR_A']; bA = wd['b_occ']*k
        def fn(t):
            d = dict(t=t, f=np.ones(len(t)), sig=np.ones(len(t)), cad=z['cad'], sector=z['sector'], epoch=np.zeros(len(t), int))
            L = LC(d, use_offsets=False)
            return 1 + wd['f_A']*(L.transit(T0f, P, k, aRA, bA, 0.1, 0.2) - 1)
        return fn
    if name == 'grazing':
        g = json.load(open(f'{SD}/out/dilution_profile_inv.json'))
        row = [r for r in g['k_profile'] if abs(r['k'] - 0.27) < 1e-6][0]['M3V_rho']
        def fn(t):
            d = dict(t=t, f=np.ones(len(t)), sig=np.ones(len(t)), cad=z['cad'], sector=z['sector'], epoch=np.zeros(len(t), int))
            L = LC(d, use_offsets=False)
            u1 = 2*np.sqrt(row['q1'])*row['q2']; u2 = np.sqrt(row['q1'])*(1 - 2*row['q2'])
            return L.transit(T0f, P, row['k'], row['aR'], row['b'], u1, u2)
        return fn
rows = []
phases = [0.15, 0.24, 0.33, 0.42, 0.58, 0.67, 0.76, 0.85]
for tr in truths:
    for phc in phases:
        T0f = T0r + phc*P
        lc = build_lc(truth_fn(tr, T0f), T0f)
        fp = fit_planet(lc, T0f); fw = fit_wd(lc, T0f); fv = fit_vshape(lc, T0f)
        row = dict(truth=tr, phase=phc, n=lc.n, nwin=int(lc.nw), chi2_flat=float(lc.chi2(np.ones(lc.n))), planet=fp, wd=fw, vshape=fv,
                   dchi2_wd_minus_planet=fw['chi2'] - fp['chi2'], dchi2_v_minus_planet=fv['chi2'] - fp['chi2'])
        rows.append(row)
        print(tr, phc, round(row['chi2_flat'] - fp['chi2'], 1), {k: round(v, 3) for k, v in fp.items()}, 'dWD', round(row['dchi2_wd_minus_planet'], 2), 'dV', round(row['dchi2_v_minus_planet'], 2), flush=True)
flat = []
for r in rows:
    flat.append(dict(truth=r['truth'], phase=r['phase'], n=r['n'], nwin=r['nwin'], dchi2_transit=r['chi2_flat'] - r['planet']['chi2'],
                     **{f'pl_{k}': v for k, v in r['planet'].items()}, **{f'wd_{k}': v for k, v in r['wd'].items()}, vs_chi2=r['vshape']['chi2'],
                     dchi2_wd_minus_planet=r['dchi2_wd_minus_planet'], dchi2_v_minus_planet=r['dchi2_v_minus_planet']))
D = pd.DataFrame(flat)
tag = '_'.join(truths)
D.to_csv(f'{SD}/out/inject_rows_{tag}.csv', index=False)
summ = {}
for tr in truths:
    q = D[D.truth == tr]
    summ[tr] = {c: dict(mean=float(q[c].mean()), std=float(q[c].std()), min=float(q[c].min()), max=float(q[c].max())) for c in
                ['pl_rho', 'pl_b', 'pl_k', 'pl_T14_h', 'pl_flat', 'dchi2_wd_minus_planet', 'dchi2_v_minus_planet', 'dchi2_transit']}
json.dump(summ, open(f'{SD}/out/inject_{tag}.json', 'w'), indent=1)
print(json.dumps(summ, indent=1))
