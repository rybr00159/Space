"""b05: (1) red-noise factors beta per sector from fake epochs (single-window depth-scale scatter of a fixed planet
template / formal error; REF: the referee's fake windows; INV: the 48 re-detrended fake ephemerides); (2) ML fits of
the real REF and INV-clean light curves for three shape classes (fitlib: planet / grazBD / compBD), with white and
beta noise; (3) profile of chi2(grazBD) - chi2(planet) in k and in rho_c (where must a grazing BD sit?).
Outputs: out/b05_beta.json, out/b05_fits.json"""
import sys; sys.path.insert(0, 'tic214311908/analysis/bd/scripts')
from blib import *
from lclib import *
import fitlib as F
from scipy.optimize import minimize
NP = {'ref': 2, 'inv': 3}
# ---------------- (1) beta
def beta_from_fake(ds):
    z = np.load(f'{BD}/data/{ds}_fake.npz')
    dd = Data(z['t'], z['f'], z['e'], z['sec'], z['cad'], z['win'], z['tc'], z['oot'], npoly=NP[ds])
    s = dd.template(cad_curves(0.12, 0.5, float(aR_of_rho(9.3))))
    # per-window amplitude: (proj-residual) least squares of y on s within each window
    out = {}
    amps, errs, secs = [], [], []
    for wv in range(dd.nw):
        m = dd.wi == wv
        if m.sum() < 20 or np.sum(s[m] != 0) < 3: continue
        X = np.vstack([dd.dt[m]**p for p in range(dd.npoly)] + [s[m]]).T
        W = dd.w[m]
        A = X.T @ (X*W[:, None]); c = np.linalg.solve(A, X.T @ (W*dd.y[m]))
        amps.append(c[-1]); errs.append(np.sqrt(np.linalg.inv(A)[-1, -1])); secs.append(dd.sec[m][0])
    amps, errs, secs = map(np.array, (amps, errs, secs))
    for sct in np.unique(secs):
        m = secs == sct
        z_ = amps[m]/errs[m]
        out[int(sct)] = dict(n=int(m.sum()), beta_std=float(np.std(z_)), beta_robust=float(rstd(z_)), beta=float(max(1.0, 0.5*(np.std(z_) + rstd(z_)))),
                             mean_amp=float(np.mean(amps[m])), white_scale=float(dd.scale[int(sct)]))
    return out
BETA = {ds: beta_from_fake(ds) for ds in ('ref', 'inv')}
jdump(BETA, f'{BD}/out/b05_beta.json'); print(json.dumps(BETA, indent=0), flush=True)
# ---------------- (2) fits
res = {}
for ds in ('ref', 'inv'):
    for nm in ('beta', 'white'):
        bt = {s: v['beta'] for s, v in BETA[ds].items()} if nm == 'beta' else None
        d = load_npz(f'{BD}/data/{ds}_real.npz', npoly=NP[ds], beta=bt)
        R = {}
        R['planet'] = F.fit(d, 'planet')
        dT = R['planet']['x'][0]
        R['grazBD'] = F.fit(d, 'grazBD', dT)
        R['compBD'] = F.fit(d, 'compBD', dT)
        R['Lambda_grazBD'] = R['grazBD']['chi2'] - R['planet']['chi2']
        R['Lambda_compBD'] = R['compBD']['chi2'] - R['planet']['chi2']
        R['chi2_null'] = float(d.Q(d.y, d.y)); R['N'] = len(d.t)
        res[f'{ds}_{nm}'] = R
        print(ds, nm, {k: (round(v['chi2'], 2), round(v['k'], 3), round(v['b'], 3), round(v['rho_c'], 2), round(v['D'], 3), round(v['T14_h'], 2), round(v['flat'], 2), round(v['dT0_min'], 2)) for k, v in R.items() if isinstance(v, dict)},
              'Lambda G/C', round(R['Lambda_grazBD'], 2), round(R['Lambda_compBD'], 2), flush=True)
        # ---------------- (3) profiles for the beta noise model
        if nm == 'beta':
            prof_k = []
            for kk in (0.19, 0.21, 0.23, 0.25, 0.27, 0.29, 0.31, 0.33):
                F.LIM['grazBD_fixk'] = dict(k=(kk*0.999, kk*1.001), D=(1.0, 1.0))
                st = [[dT, np.log(kk), bb, np.log(r), U1, U2] for bb, r in ((0.85, 2.), (0.88, 1.2), (0.9, 4.), (0.83, 0.8))]
                f_ = lambda x: F.chi2(d, x, 'grazBD_fixk')
                best = None
                for s0 in st:
                    r_ = minimize(f_, s0, method='Nelder-Mead', options=dict(maxiter=2500, maxfev=2500, xatol=1e-6, fatol=1e-4, adaptive=True))
                    if best is None or r_.fun < best.fun: best = r_
                _, k_, bb_, rho_, *_ = F.unpack(best.x, 'grazBD')
                prof_k.append(dict(k=kk, dchi2=float(best.fun - R['planet']['chi2']), b=float(bb_*(1 + k_)), rho_c=float(rho_)))
            prof_r = []
            for lr in np.log([0.5, 0.8, 1.2, 1.7, 2.5, 3.5, 5., 7., 9.3]):
                F.LIM['grazBD_fixr'] = dict(k=(F.KG_LO, F.KG_HI), D=(1.0, 1.0))
                f_ = lambda x: F.chi2(d, np.r_[x[0], x[1], x[2], lr, x[3], x[4]], 'grazBD_fixr')
                best = None
                for s0 in ([dT, np.log(0.22), 0.88, U1, U2], [dT, np.log(0.28), 0.86, U1, U2], [dT, np.log(0.32), 0.85, U1, U2]):
                    r_ = minimize(f_, s0, method='Nelder-Mead', options=dict(maxiter=2500, maxfev=2500, xatol=1e-6, fatol=1e-4, adaptive=True))
                    if best is None or r_.fun < best.fun: best = r_
                prof_r.append(dict(rho_c=float(np.exp(lr)), dchi2=float(best.fun - R['planet']['chi2']), k=float(np.exp(best.x[1])), b=float(best.x[2]*(1 + np.exp(best.x[1])))))
            R['profile_k_grazBD'] = prof_k; R['profile_rho_grazBD'] = prof_r
            print('  profile k', [(p_['k'], round(p_['dchi2'], 2), round(p_['rho_c'], 2)) for p_ in prof_k], flush=True)
            print('  profile rho', [(round(p_['rho_c'], 2), round(p_['dchi2'], 2)) for p_ in prof_r], flush=True)
jdump(dict(method=__doc__, beta=BETA, fits=res), f'{BD}/out/b05_fits.json')
