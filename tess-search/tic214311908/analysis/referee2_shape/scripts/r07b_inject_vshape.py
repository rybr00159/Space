"""r07b: calibration of dchi2(V-shape - planet) under a V-shaped (grazing) truth: the r03 'vshape' ML model of each
data set is injected into circularly shifted residuals (as r07 'shift'); planet_free and vshape are refitted.
usage: r07b_inject_vshape.py DATASET NTRIAL -> out/inject_{ds}_vshape.json/.csv"""
import sys, json, numpy as np, pandas as pd
sys.argv = [sys.argv[0], sys.argv[1], 'shift', '0'] + sys.argv[2:]
NT2 = int(sys.argv[4])
src = open('tic214311908/analysis/referee2_shape/scripts/r07_inject.py').read()
src = src[:src.index('rows = []')]
exec(src)
rng = np.random.default_rng(77 if ds == 'inv' else 78)
sup0 = Super(lc0)
_, res = lc0.chi2(planet_model(lc0, xp), ret=True)
truth = planet_model(lc0, xv, vshape=True)
rows = []
for j in range(NT2):
    sh = np.empty(lc0.n)
    for wv in range(lc0.nw):
        ii = np.where(lc0.wi == wv)[0]; sh[ii] = np.roll(res[ii], rng.integers(len(ii)))
    lc = LC(dict(D0, f=truth + sh))
    def c2p(x, vshape=False):
        T0, P, k, bx, lr, u1, u2 = x
        if not (0.01 < k < 0.9 and -1 < lr < 2.5): return 1e9
        if not (u1 + u2 < 1 and u1 > 0 and u1 + 2*u2 > 0): return 1e9
        if not vshape and abs(bx) >= 1 + k: return 1e9
        return lc.chi2(planet_model(lc, x, vshape)) + ((u1 - U1)/SU)**2 + ((u2 - U2)/SU)**2
    best = {}
    for lab, f, starts in (('p', c2p, [xp, xs, np.r_[xp[:3], 0.8, 0.6, U1, U2]]), ('v', lambda x: c2p(x, True), [xv, np.r_[xv[:3], 0.3, 0.9, U1, U2]])):
        bb = None
        for s in starts:
            r = minimize(f, s, method='Nelder-Mead', options=dict(maxiter=4000, maxfev=4000, xatol=1e-7, fatol=1e-4))
            r = minimize(f, r.x, method='Nelder-Mead', options=dict(maxiter=3000, maxfev=3000, xatol=1e-8, fatol=1e-5))
            if bb is None or r.fun < bb.fun: bb = r
        best[lab] = bb
    rows.append(dict(trial=j, dVP=best['v'].fun - best['p'].fun, b_planetfit=abs(best['p'].x[3]), k_planetfit=best['p'].x[2]))
    print(rows[-1], flush=True)
R = pd.DataFrame(rows); R.to_csv(f'{RD}/out/inject_{ds}_vshape.csv', index=False)
realv = S['vshape']['chi2'] - S['planet_free']['chi2']
o = dict(dataset=ds, real_dVP=realv, vshape_truth_dVP=dict(mean=float(R.dVP.mean()), std=float(R.dVP.std()), max=float(R.dVP.max()), frac_ge_real=float(np.mean(R.dVP >= realv)), n=len(R)))
print(o); json.dump(o, open(f'{RD}/out/inject_{ds}_vshape.json', 'w'), indent=1)
