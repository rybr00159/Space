"""s03: nested-sampling (dynesty) transit fits of the M3V light curves.
usage: s03_planet_fit.py DATASET MODEL [nlive]
DATASET: inv (../phot/lc_M3V_all.csv windows +-0.3 d, rescaled to 979 e-/s) | ref (../referee_phot windows)
         | invA (investigator model-A photometry) -- all with per-sector white-scale x beta errors (s01) and
         per-window constant offsets profiled.
MODEL:
  free      : P, T0, k, b, log10 rho_circ (U[-1, 2.3], g/cm3), u1, u2 ~ N(Claret17, 0.1); circular
  freeLD    : as free with Kipping q1, q2 ~ U(0,1)
  single    : rho ~ N(9.41, 0.88) (single-star value), circular
  twin      : planet around one component of an equal twin: rho ~ N(15.5, 1.45), the light is diluted by 2
              (observed = 1 + (transit - 1)/2), k = undiluted radius ratio, LD of the 3237 K twin; circular
  single_ecc / twin_ecc : as single / twin with e, w free (sqrt(e) cos w, sqrt(e) sin w ~ U, e < 0.9)
  wd        : (scenario test, see s06) not here
Output: out/fit_{DATASET}_{MODEL}.npz (equal-weight posterior samples + logZ), out/fit_{DATASET}_{MODEL}.json"""
import sys, json, numpy as np
sys.path.insert(0, 'tic214311908/analysis/shape/scripts')
from shapelib import *
import dynesty
from dynesty import utils as dyfunc
from scipy.stats import norm
ds, model = sys.argv[1], sys.argv[2]
nlive = int(sys.argv[3]) if len(sys.argv) > 3 else 500
path = {'inv': f'{SD}/data/inv_win.npz', 'ref': f'{SD}/data/ref_win.npz'}[ds]
lc = LC(path)
ld = json.load(open(f'{SD}/out/ld_claret17.json'))
twin = model.startswith('twin'); ecc = model.endswith('_ecc')
U1, U2 = (ld['adopted_twin']['u1'], ld['adopted_twin']['u2']) if twin else (ld['adopted_single']['u1'], ld['adopted_single']['u2'])
names = ['T0', 'P', 'k', 'b', 'log10rho', 'u1', 'u2'] + (['secosw', 'sesinw'] if ecc else [])
nd = len(names)
def ptform(x):
    th = np.empty(nd)
    th[0] = T0REF + (x[0]-0.5)*0.03
    th[1] = PREF + (x[1]-0.5)*0.0003
    th[2] = 0.02 + x[2]*(0.6 if twin else 0.4)
    th[3] = x[3]*(1 + th[2])
    if model in ('free', 'freeLD'):
        th[4] = -1 + x[4]*3.3
    elif model.startswith('single'):
        th[4] = np.log10(max(0.5, norm.ppf(x[4], 9.41, 0.88)))
    else:
        th[4] = np.log10(max(0.5, norm.ppf(x[4], 15.5, 1.45)))
    if model == 'freeLD':
        q1, q2 = x[5], x[6]; th[5] = 2*np.sqrt(q1)*q2; th[6] = np.sqrt(q1)*(1-2*q2)
    else:
        th[5] = norm.ppf(x[5], U1, 0.1); th[6] = norm.ppf(x[6], U2, 0.1)
    if ecc:
        th[7] = -1 + 2*x[7]; th[8] = -1 + 2*x[8]
    return th
def unpack(th):
    T0, P, k, b, lr, u1, u2 = th[:7]
    if ecc:
        e = th[7]**2 + th[8]**2; w = np.degrees(np.arctan2(th[8], th[7]))
    else:
        e, w = 0., 90.
    return T0, P, k, b, 10**lr, u1, u2, e, w
def loglike(th):
    T0, P, k, b, rho, u1, u2, e, w = unpack(th)
    if e >= 0.9: return -1e300
    if not (u1 + u2 < 1 and u1 + 2*u2 > 0 and u1 > -0.5): return -1e300
    aR = aR_from_rho(rho, P)
    # eccentric: rho here is the TRUE stellar density; aR the true semi-major axis
    if aR*(1-e) < 1.5: return -1e300
    m = lc.transit(T0, P, k, aR, b, u1, u2, e, w)
    if twin: m = 1 + (m - 1)/2
    return -0.5*lc.chi2(m)
rng = np.random.default_rng(42)
sampler = dynesty.NestedSampler(loglike, ptform, nd, nlive=nlive, bound='multi', sample='rwalk', walks=40, rstate=rng)
sampler.run_nested(dlogz=0.05, print_progress=False)
res = sampler.results
import pickle
pickle.dump(dict(samples=res.samples, logwt=res.logwt, logl=res.logl, logz=res.logz, logzerr=res.logzerr, ncall=res.ncall), open(f'{SD}/out/raw_{ds}_{model}.pkl', 'wb'))
w = np.exp(res.logwt - res.logz[-1])
samp = dyfunc.resample_equal(res.samples, w/w.sum(), rstate=rng)
T0, P, k, b, rho, u1, u2, e, wv = np.array([unpack(s) for s in samp]).T
aR = aR_from_rho(rho, P)
fac = (1 + e*np.sin(np.radians(wv)))/(1 - e**2)
inc = np.degrees(np.arccos(np.clip(b/aR*fac, 0, 1)))
T14, T23 = t14_t23(P, k, aR, b, e, wv)
# circular-equivalent density (what a circular fit would infer): rho_circ = rho * g^3, g = (1+e sin w)/sqrt(1-e^2)
g = (1 + e*np.sin(np.radians(wv)))/np.sqrt(1 - e**2)
rho_circ = rho*g**3
lmax = res.logl.max()
best = res.samples[np.argmax(res.logl)]
# chi2 of best
chi2best = -2*lmax
# observed (diluted) depth of the model at mid-transit: compute from best-fit-like samples quickly (k^2 proxy for posterior)
dep = k**2*(1 - U1*0 )  # proxy; precise model depth below for the max-likelihood sample
mb = lc.transit(*[unpack(best)[i] for i in (0, 1, 2)], aR_from_rho(unpack(best)[4], unpack(best)[1]), unpack(best)[3], unpack(best)[5], unpack(best)[6], unpack(best)[7], unpack(best)[8])
depth_best = 1 - lc.transit(unpack(best)[0], unpack(best)[1], unpack(best)[2], aR_from_rho(unpack(best)[4], unpack(best)[1]), unpack(best)[3], unpack(best)[5], unpack(best)[6], unpack(best)[7], unpack(best)[8]).min()
np.savez(f'{SD}/out/fit_{ds}_{model}.npz', T0=T0, P=P, k=k, b=b, rho=rho, rho_circ=rho_circ, u1=u1, u2=u2, e=e, w=wv, aR=aR, inc=inc, T14=T14, T23=T23,
         logz=res.logz[-1], logzerr=res.logzerr[-1], best=best, names=names, lmax=lmax)
def q(x):
    p = np.percentile(x, [2.275, 15.865, 50, 84.135, 97.725, 0.135, 99.865]); return dict(med=p[2], lo1=p[1], hi1=p[3], lo2=p[0], hi2=p[4], lo3=p[5], hi3=p[6])
summ = dict(dataset=ds, model=model, nlive=nlive, n=lc.n, logz=float(res.logz[-1]), logzerr=float(res.logzerr[-1]), chi2_best=float(chi2best), chi2_flat=float(lc.chi2(np.ones(lc.n))),
            ncall=int(np.sum(res.ncall)), best=dict(zip(names, map(float, best))), depth_best_model=float(depth_best),
            post={nm: q(v) for nm, v in dict(T0=T0, P=P, k=k, b=b, rho=rho, rho_circ=rho_circ, u1=u1, u2=u2, e=e, aR=aR, inc=inc, T14_h=T14*24, T23_h=T23*24, flat=T23/T14).items()})
json.dump(summ, open(f'{SD}/out/fit_{ds}_{model}.json', 'w'), indent=1, default=float)
print(json.dumps({k: (v if k != 'post' else {kk: [round(vv['med'], 5), round(vv['lo1'], 5), round(vv['hi1'], 5)] for kk, vv in v.items()}) for k, v in summ.items()}, indent=0, default=float))
