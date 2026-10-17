"""r07: injection calibration of the light-curve shape statistic dchi2 = chi2(WD occultation) - chi2(planet) and
of dchi2(V-shape) - chi2(planet), for BOTH photometries.
Method 'shift' (INV and REF): residuals of the real data about the ML planet model (r03 planet_free) are
circularly shifted inside each window by a random number of cadences (keeps the noise and its correlations, breaks
the link with the real signal); a truth model is added: 'planet' = r03 planet_single ML, 'wd' = r03 wd_cool ML
(both of the same data set). Method 'fake' (INV only): truths injected into the raw flux at fake epochs and
re-detrended with the windows masked (as r01). Each fake light curve is fitted with planet_free (LD prior, rho
free), wd_cool (light curve + M-dwarf priors) and vshape, with a few starts around the truth and alternatives.
usage: r07_inject.py DATASET METHOD NTRIAL -> out/inject_{ds}_{method}.csv/.json"""
import sys, json, numpy as np, pandas as pd
sys.path.insert(0, 'tic214311908/analysis/referee2_shape/scripts')
from rlib import *
from rwd import *
from scipy.optimize import minimize
ds, method, NT = sys.argv[1], sys.argv[2], int(sys.argv[3])
rng = np.random.default_rng({'inv': 5, 'ref': 6}[ds] + (100 if method == 'fake' else 0))
S = json.load(open(f'{RD}/out/scen_{ds}_fake.json'))
D0 = load(ds, 'fake'); lc0 = LC(D0)
U1, U2, SU = 0.125, 0.533, 0.10
xp = np.array(S['planet_free']['x']); xs = np.array(S['planet_single']['x']); xw = np.array(S['wd_cool']['x']); xv = np.array(S['vshape']['x'])
def planet_model(lc, x, vshape=False):
    T0, P, k, bx, lr, u1, u2 = x
    b = (1 - k) + np.clip(bx, 0, 1)*2*k if vshape else abs(bx)
    return lc.transit(T0, P, k, aR_from_rho(10**lr, P), b, u1, u2)
def wd_model(lc, sup, y):
    T0, P, lMA, xf, bocc, fA, RB, MB = y
    MA = 10**lMA; fmax = 1.3
    RA = r_nauenberg(MA)*(1 + np.clip(xf, 0, 1)*(fmax - 1))
    a = kepler_a_rsun(MA + MB, P); cosi = abs(bocc)*RB/a
    if cosi >= 1: return None
    return wd_occultation(sup, T0, P, a, cosi, RA, RB, fA)
def fits(lc, sup, dT0=0.):
    def c2p(x, vshape=False):
        T0, P, k, bx, lr, u1, u2 = x
        if not (0.01 < k < 0.9 and -1 < lr < 2.5): return 1e9
        if not (u1 + u2 < 1 and u1 > 0 and u1 + 2*u2 > 0): return 1e9
        if not vshape and abs(bx) >= 1 + k: return 1e9
        return lc.chi2(planet_model(lc, x, vshape)) + ((u1 - U1)/SU)**2 + ((u2 - U2)/SU)**2
    def c2w(y):
        T0, P, lMA, xf, bocc, fA, RB, MB = y
        if not (np.log10(0.15) <= lMA <= np.log10(1.38) and 0 <= xf <= 1 and 0 <= fA < 0.2 and 0.05 < RB < 0.6 and 0.05 < MB < 0.6): return 1e9
        m = wd_model(lc, sup, y)
        if m is None: return 1e9
        return lc.chi2(m) + ((RB - 0.381)/0.012)**2 + ((MB - 0.370)/0.011)**2
    def ms(f, starts):
        best = None
        for s in starts:
            r = minimize(f, s, method='Nelder-Mead', options=dict(maxiter=4000, maxfev=4000, xatol=1e-7, fatol=1e-4))
            r = minimize(f, r.x, method='Nelder-Mead', options=dict(maxiter=3000, maxfev=3000, xatol=1e-8, fatol=1e-5))
            if best is None or r.fun < best.fun: best = r
        return best
    sh = np.array([dT0, 0, 0, 0, 0, 0, 0])
    p = ms(c2p, [xp + sh, xs + sh, np.r_[xp[:3] + sh[:3], 0.8, 0.6, U1, U2]])
    v = ms(lambda x: c2p(x, True), [xv + sh, np.r_[xv[:3] + sh[:3], 0.3, 0.9, U1, U2]])
    shw = np.r_[dT0, np.zeros(7)]
    w = ms(c2w, [xw + shw, np.r_[xw[:2] + shw[:2], np.log10(0.3), 0.0, 0.3, 0.0157, 0.381, 0.370], np.r_[xw[:2] + shw[:2], np.log10(0.8), 0.0, 0.1, 0.0157, 0.381, 0.370]])
    T0, P, k, bx, lr, u1, u2 = p.x
    aR = aR_from_rho(10**lr, P); T14, T23 = t14(P, k, aR, abs(bx))
    return dict(chi2_p=p.fun, chi2_w=w.fun, chi2_v=v.fun, k=k, b=abs(bx), rho=10**lr, T14_h=float(T14*24), flat=float(T23/T14) if T14 > 0 else 0., M_WD=10**w.x[2], b_occ=abs(w.x[4]))
rows = []
if method == 'shift':
    sup0 = Super(lc0)
    mp = planet_model(lc0, xp); _, res = lc0.chi2(mp, ret=True)
    truths = dict(planet=planet_model(lc0, xs), wd=wd_model(lc0, sup0, xw))
    for j in range(NT):
        sh = np.empty(lc0.n)
        for wv in range(lc0.nw):
            ii = np.where(lc0.wi == wv)[0]
            sh[ii] = np.roll(res[ii], rng.integers(len(ii)))
        for tr, tm in truths.items():
            lc = LC(dict(D0, f=tm + sh)); sup = Super(lc)
            r = fits(lc, sup); r.update(truth=tr, trial=j); rows.append(r)
            print(tr, j, round(r['chi2_w'] - r['chi2_p'], 2), round(r['chi2_v'] - r['chi2_p'], 2), round(r['rho'], 2), round(r['b'], 2), flush=True)
else:
    from wotan import flatten
    d = pd.read_csv(f'{B}/phot/lc_M3V_all.csv')
    d = d[(d.quality_ok == 1) & (d.flare_flag == 0) & np.isfinite(d.flux_raw_norm) & np.isfinite(d.flux_err)].copy()
    t = d.time_btjd_tdb.values; raw0 = d.flux_raw_norm.values; e = d.flux_err.values*935/979; sec = d.sector.values; cad = d.cadence_s.values.astype(int)
    real = np.abs(((t - T0REF)/PREF + 0.5) % 1 - 0.5)*PREF < 0.1
    bj = json.load(open(f'{RD}/out/beta_inv.json'))
    sig = e*np.array([bj[str(int(s))]['scale']*bj[str(int(s))]['beta'] for s in sec])
    phs = rng.permutation(np.r_[np.linspace(0.12, 0.40, 15), np.linspace(0.60, 0.88, 15)])[:NT]
    full = LC(dict(t=t, f=np.ones(len(t)), sig=np.ones(len(t)), sector=sec, cad=cad, epoch=np.zeros(len(t), int), dt=np.zeros(len(t)), wid=sec))
    fsup = Super(full)
    for j, ph in enumerate(phs):
        T0f = T0REF + ph*PREF; dT0 = ph*PREF
        dtf = ((t - T0f)/PREF + 0.5) % 1*PREF - 0.5*PREF
        fullj = LC(dict(t=t, f=np.ones(len(t)), sig=np.ones(len(t)), sector=sec, cad=cad, epoch=np.zeros(len(t), int), dt=dtf, wid=sec)); supj = Super(fullj)
        for tr in ('planet', 'wd'):
            if tr == 'planet':
                tm = planet_model(fullj, np.r_[xs[0] + dT0, xs[1:]])
            else:
                tm = wd_model(fullj, supj, np.r_[xw[0] + dT0, xw[1:]])
            raw = raw0 + (tm - 1)*979/935
            fake = np.abs(dtf) < 0.1
            f = np.full(len(t), np.nan)
            for s in np.unique(sec):
                m = sec == s; use = m & ~real & ~fake
                trd = flatten(t[use], raw[use], method='biweight', window_length=0.75, return_trend=True, break_tolerance=0.5, edge_cutoff=0)[1]
                ok = np.isfinite(trd)
                f[m] = raw[m]/np.interp(t[m], t[use][ok], trd[ok])
            f = 1 + (f - 1)*935/979
            ep = np.round((t - T0f)/PREF).astype(int)
            w = (np.abs(dtf) < 0.3) & np.isfinite(f) & ~(np.abs(((t - T0REF)/PREF + 0.5) % 1 - 0.5)*PREF < 0.3)
            D = dict(t=t[w], f=f[w], sig=sig[w], sector=sec[w], cad=cad[w], epoch=ep[w], dt=dtf[w], wid=sec[w]*1000 + ep[w])
            lc = LC(D); sup = Super(lc)
            r = fits(lc, sup, dT0); r.update(truth=tr, trial=j, phase=float(ph)); rows.append(r)
            print(tr, j, round(ph, 3), round(r['chi2_w'] - r['chi2_p'], 2), round(r['chi2_v'] - r['chi2_p'], 2), round(r['rho'], 2), round(r['b'], 2), flush=True)
R = pd.DataFrame(rows); R['dWP'] = R.chi2_w - R.chi2_p; R['dVP'] = R.chi2_v - R.chi2_p
R.to_csv(f'{RD}/out/inject_{ds}_{method}.csv', index=False)
real = dict(dWP=S['wd_cool']['chi2'] - S['planet_free']['chi2'], dVP=S['vshape']['chi2'] - S['planet_free']['chi2'])
out = dict(dataset=ds, method=method, real=real)
for tr in ('planet', 'wd'):
    q = R[R.truth == tr]
    out[tr] = {c: dict(mean=float(q[c].mean()), std=float(q[c].std()), p16=float(q[c].quantile(0.16)), p84=float(q[c].quantile(0.84)), n=int(len(q))) for c in ['dWP', 'dVP', 'rho', 'b', 'T14_h', 'k', 'M_WD']}
# likelihood ratio of the real dWP under the two truths (Gaussian approximation to each injection distribution, and KDE)
from scipy.stats import norm, gaussian_kde
for stat in ('dWP',):
    a = R[R.truth == 'planet'][stat].values; b = R[R.truth == 'wd'][stat].values; x = real[stat]
    out[f'LR_planet_over_wd_{stat}_gauss'] = float(norm.pdf(x, a.mean(), a.std())/norm.pdf(x, b.mean(), b.std()))
    try: out[f'LR_planet_over_wd_{stat}_kde'] = float(gaussian_kde(a)(x)[0]/gaussian_kde(b)(x)[0])
    except Exception as ex: out['kde_error'] = str(ex)
    out[f'frac_wd_truth_ge_real_{stat}'] = float(np.mean(b >= x)); out[f'frac_planet_truth_ge_real_{stat}'] = float(np.mean(a >= x))
print(json.dumps(out, indent=1))
json.dump(out, open(f'{RD}/out/inject_{ds}_{method}.json', 'w'), indent=1)
