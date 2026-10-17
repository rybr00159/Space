"""s09: per-sector (and per-cadence) consistency of depth and duration.
For each data set (INV model B, INV model A = all-free deblending, REF) the global best-fit shape of the free fit
(s03, max-likelihood sample) is used; per sector we fit (a) the radius ratio k alone (depth), and (b) k plus a
duration scale s (a/R* -> a/R*/s, b fixed in units of R*, i.e. T14 -> ~s*T14). Errors from the Delta-chi2 = 1
profile with the beta-inflated errors (s01). Consistency chi2 about the weighted mean. Also odd/even depth.
Output: out/per_sector.json, out/per_sector.csv"""
import sys, json, numpy as np, pandas as pd
sys.path.insert(0, 'tic214311908/analysis/shape/scripts')
from shapelib import *
from scipy.optimize import minimize, minimize_scalar, brentq
def prof_err(f, x0, step):
    """1-sigma from chi2 profile of a 1-D function f (chi2) around its minimum x0"""
    c0 = f(x0)
    def g(x): return f(x) - c0 - 1
    try:
        hi = brentq(g, x0, x0 + step*20)
    except Exception: hi = np.nan
    try:
        lo = brentq(g, x0 - step*20, x0)
    except Exception: lo = np.nan
    return 0.5*((hi - x0) + (x0 - lo))
rows = []; summ = {}
for ds, path, extra in [('inv', 'data/inv_win.npz', None), ('invA', 'data/inv_win.npz', 'A'), ('ref', 'data/ref_win.npz', None)]:
    fit = json.load(open(f'{SD}/out/fit_{"inv" if ds == "invA" else ds}_free.json'))
    bst = fit['best']
    T0, P, k0, b0, rho0, u1, u2 = bst['T0'], bst['P'], bst['k'], bst['b'], 10**bst['log10rho'], bst['u1'], bst['u2']
    aR0 = aR_from_rho(rho0, P)
    z = np.load(f'{SD}/{path}')
    secs = sorted(np.unique(z['sector']))
    if extra == 'A':
        # model-A light curve: same windows, fA/eA; per-sector white scale x beta from OOT
        pass
    res = []
    for s in secs + ['all', 'odd', 'even', '30min', '10min', '200s']:
        if s == 'all': m = np.ones(len(z['t']), bool)
        elif s in ('odd', 'even'): m = (z['epoch'] % 2 == (1 if s == 'odd' else 0))
        elif s == '30min': m = z['cad'] == 1800
        elif s == '10min': m = z['cad'] == 600
        elif s == '200s': m = z['cad'] == 200
        else: m = z['sector'] == s
        lc = LC(f'{SD}/{path}', mask=m)
        if extra == 'A':
            fA = z['fA'][m]; eA = z['eA'][m]
            ok = np.isfinite(fA) & np.isfinite(eA)
            lc = LC(f'{SD}/{path}', mask=np.where(m)[0][ok] if False else (m & np.isfinite(z['fA']) & np.isfinite(z['eA'])))
            mm = m & np.isfinite(z['fA']) & np.isfinite(z['eA'])
            lc.f = z['fA'][mm]
            # error: OOT robust std per sector x beta of model B (same pixels/time sampling)
            sig = np.zeros(lc.n)
            nz = json.load(open(f'{SD}/out/noise.json'))['INV']
            for ss in np.unique(lc.sector):
                q = lc.sector == ss; oot = np.abs(z['dt'][mm][q]) > 0.1
                r = lc.f[q][oot] - 1
                sig[q] = 1.4826*np.median(np.abs(r - np.median(r)))*nz[str(int(ss))]['beta']
            lc.sig = sig; lc.w = 1/sig**2; lc.sw = np.bincount(lc.wi, lc.w, lc.nw)
        f1 = lambda kk: lc.chi2(lc.transit(T0, P, kk, aR0, b0, u1, u2))
        r1 = minimize_scalar(f1, bounds=(0.01, 0.4), method='bounded')
        ek = prof_err(f1, r1.x, 0.005)
        f2 = lambda x: lc.chi2(lc.transit(T0, P, x[0], aR0/x[1], b0, u1, u2))
        r2 = minimize(f2, [r1.x, 1.0], method='Nelder-Mead', options=dict(xatol=1e-5, fatol=1e-4))
        es = prof_err(lambda sc: min(minimize_scalar(lambda kk: f2([kk, sc]), bounds=(0.01, 0.4), method='bounded').fun for _ in [0]), r2.x[1], 0.02)
        T14, T23 = t14_t23(P, r2.x[0], aR0/r2.x[1], b0)
        dep = 1 - lc.transit(T0, P, r1.x, aR0, b0, u1, u2).min()
        dep_e = dep*2*ek/r1.x
        row = dict(dataset=ds, subset=str(s), n=int(lc.n), k=float(r1.x), k_err=float(ek), depth_ppt=float(dep*1e3), depth_err_ppt=float(dep_e*1e3),
                   dur_scale=float(r2.x[1]), dur_scale_err=float(es), T14_h=float(T14*24), T14_err_h=float(T14*24*es/r2.x[1]) if np.isfinite(es) else np.nan,
                   chi2_flat_minus_fit=float(lc.chi2(np.ones(lc.n)) - r1.fun))
        res.append(row); rows.append(row)
        print(row, flush=True)
    R = pd.DataFrame(res); S = R[R.subset.isin([str(s) for s in secs])]
    w = 1/S.depth_err_ppt**2; mu = np.sum(w*S.depth_ppt)/w.sum()
    wd = 1/S.dur_scale_err**2; mud = np.sum(wd*S.dur_scale)/wd.sum()
    oe = R.set_index('subset').loc[['odd', 'even']]
    summ[ds] = dict(depth_wmean_ppt=float(mu), depth_wmean_err=float(1/np.sqrt(w.sum())), depth_chi2=float(np.sum(w*(S.depth_ppt-mu)**2)), dof=len(S)-1,
                    dur_wmean=float(mud), dur_chi2=float(np.sum(wd*(S.dur_scale-mud)**2)),
                    odd_minus_even_ppt=float(oe.depth_ppt.iloc[0]-oe.depth_ppt.iloc[1]), odd_even_err=float(np.hypot(oe.depth_err_ppt.iloc[0], oe.depth_err_ppt.iloc[1])))
    print(ds, summ[ds])
pd.DataFrame(rows).to_csv(f'{SD}/out/per_sector.csv', index=False)
json.dump(summ, open(f'{SD}/out/per_sector.json', 'w'), indent=1)
