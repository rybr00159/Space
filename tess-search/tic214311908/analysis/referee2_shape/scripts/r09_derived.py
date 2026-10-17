"""r09: derived quantities from the r02 emcee posteriors.
(a) Vanderburg+2019 / Seager & Mallen-Ornelas maximum dilution: k <= T12/T13 for any transit geometry, so the
    undiluted depth is <= (T12/T13)^2 (uniform source; LD makes the true bound slightly looser) and
    dT_max = 2.5 log10((T12/T13)^2 / depth_obs). Posterior percentiles of dT_max (free fits).
(b) Photo-eccentric host comparison: rho_circ posterior (uniform prior in log rho) -> marginal-likelihood ratio of
    the single host N(9.41, 0.88) vs the twin host N(15.5, 1.45) (rho only, circular orbits).
(c) Adopted planet parameters from the 'single' and 'twin' posteriors (R* 0.381 +- 0.012 / 0.289 +- 0.009,
    Teff 3350 +- 55 / 3237, M* 0.370 / 0.264), mean of INV and REF.
Output: out/derived.json"""
import sys, json, numpy as np
sys.path.insert(0, 'tic214311908/analysis/referee2_shape/scripts')
from rlib import *
import batman
from scipy.stats import gaussian_kde
rng = np.random.default_rng(4)
out = {}
def obs_depth(T0, P, k, aR, b, u1, u2, dil=1.0):
    p = batman.TransitParams(); p.t0 = 0; p.per = P; p.rp = k; p.a = aR; p.inc = np.degrees(np.arccos(min(1, b/aR))); p.ecc = 0; p.w = 90
    p.limb_dark = 'quadratic'; p.u = [u1, u2]
    return (1 - batman.TransitModel(p, np.array([0.0])).light_curve(p)[0])/dil
def pct(x, q=(0.135, 2.275, 15.865, 50, 84.135, 97.725, 99.865)):
    return [float(v) for v in np.percentile(x, q)]
for ds in ('inv', 'ref'):
    z = np.load(f'{RD}/out/mcmc_{ds}_free_fake.npz')
    ch = z['chain']; T14_, T23_ = z['T14'], z['T23']
    T0, P, k, b, lr, u1, u2 = ch.T
    T12 = (T14_ - T23_)/2; T13 = T14_ - T12
    ratio = T12/T13
    idx = rng.choice(len(ch), 3000, replace=False)
    dep = np.array([obs_depth(T0[i], P[i], k[i], aR_from_rho(10**lr[i], P[i]), b[i], u1[i], u2[i]) for i in idx])
    dTmax = 2.5*np.log10(ratio[idx]**2/dep)
    rho = 10**lr
    kd = gaussian_kde(lr)
    def ml(mu, sd):
        s = rng.normal(mu, sd, 20000); s = s[s > 0.5]
        return float(np.mean(kd(np.log10(s))))
    Ls, Lt = ml(9.41, 0.88), ml(15.5, 1.45)
    out[ds] = dict(depth_obs_ppt=pct(dep*1e3), T12_over_T13=pct(ratio), dTmax=pct(dTmax), frac_dTmax_gt_3p5=float(np.mean(dTmax > 3.5)), frac_dTmax_gt_3=float(np.mean(dTmax > 3.0)),
                   rho_circ=pct(rho), B_single_over_twin_rho_only=Ls/Lt, lnB=float(np.log(Ls/Lt)),
                   pct_single_rho=float(np.mean(rho < 9.41)*100), pct_twin_rho=float(np.mean(rho < 15.5)*100))
    print(ds, json.dumps(out[ds]))
# adopted parameters
ad = {}
for host, R, sR, Teff, M in (('single', 0.381, 0.012, 3350., 0.370), ('twin', 0.289, 0.009, 3237., 0.264)):
    vals = {}
    for ds in ('inv', 'ref'):
        z = np.load(f'{RD}/out/mcmc_{ds}_{host}_fake.npz'); ch = z['chain']
        T0, P, k, b, lr, u1, u2 = ch.T
        n = len(k); Rs = rng.normal(R, sR, n); Ts = rng.normal(Teff, 55, n)
        aR = z['aR']
        vals[ds] = dict(k=k, Rp=k*Rs*109.08, b=b, aR=aR, inc=z['inc'], T14=z['T14']*24, Teq=Ts*np.sqrt(1/(2*aR)), dep=None)
    res = {}
    for q in ('k', 'Rp', 'b', 'aR', 'inc', 'T14', 'Teq'):
        m_ = [np.median(vals[ds][q]) for ds in ('inv', 'ref')]
        s_ = [0.5*(np.percentile(vals[ds][q], 84.135) - np.percentile(vals[ds][q], 15.865)) for ds in ('inv', 'ref')]
        res[q] = dict(inv=float(m_[0]), ref=float(m_[1]), mean=float(np.mean(m_)), err=float(np.hypot(max(s_), 0.5*abs(m_[0] - m_[1]))))
    ad[host] = res
    print(host, {q: (round(v['mean'], 4), round(v['err'], 4)) for q, v in res.items()})
out['adopted'] = ad
json.dump(out, open(f'{RD}/out/derived.json', 'w'), indent=1)
