"""Combine the six sectors: sky offset of the eclipsing source with a lateral-colour (chromatic) nuisance model.

Data: per sector s and eclipse type t (deep, shallow) the measured offset o_st = (dE, dN) of the difference-image
source from the G star (fast PRF fit, real_fast.json), bias-corrected with the empirical-image injection at the
target (inj_s*.json, delta = 0), with the injection covariance (scatter over 500 fake-epoch trials) as its error.

Chromatic model: a source whose light is redder than the G star appears displaced by
   chi_st = c_t * [k_rad,s * u_rad,s + k_tan,s * u_tan,s] / k_ref
where u_rad,s = unit sky vector from the star TOWARDS the camera boresight in sector s, k_rad,s and k_tan,s are the
field-star lateral-colour coefficients of that sector (x11, mpx per mag of BP-RP, mean of two FFIs), k_ref = mean k_rad,
and c_t is a free 'colour lever' per eclipse type (arcsec of radial shift at k_ref; negative = away from boresight).
Models:
   M0: o = delta                       (2 params)
   M1: o = chi(c_d, c_s)               (2 params)
   M2: o = delta + chi(c_d, c_s)       (4 params)
   MT: tangential components only: o.u_tan = delta.u_tan   (2 params, colour-free apart from small k_tan)
Empirical calibration: 4000 simulated data sets, each taking for every sector/type one random injection trial's
recovered-minus-injected offset (target host, empirical image) as noise around the M2 best fit; refit -> delta
distribution (calibrated error ellipse) and chi2 distributions (goodness of fit).
"""
import sys, os, json, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import least_squares
from loclib import R, SECTORS, to_json

TYPES = ['deep', 'shallow']
real = json.load(open(f'{R}/results/real_fast.json'))
bs = json.load(open(f'{R}/results/x02_boresight.json'))
ct = json.load(open(f'{R}/results/x11_colour_terms.json'))
krad = {}; ktan = {}
for s in SECTORS:
    rr = [c for c in ct if c['sec'] == s]
    krad[s] = np.mean([c['k_rad'] for c in rr]); ktan[s] = np.mean([c['k_tan'] for c in rr])
kref = np.mean(list(krad.values()))
obs, cov, U, Tn, KR, KT, lab, noise = [], [], [], [], [], [], [], []
for s in SECTORS:
    inj = json.load(open(f'{R}/results/inj_s{s}.json'))
    pa = np.radians(bs[str(s)]['pa_to_boresight'])
    u = np.array([np.sin(pa), np.cos(pa)])          # towards boresight (E, N)
    tn = np.array([u[1], -u[0]])
    for t in TYPES:
        r = real[str(s)][t]
        s0 = inj[t]['summary'][0]
        assert s0['sep'] == 0.0
        b = np.array(s0['emp']['bias']); C = np.array(s0['emp']['cov'])
        obs.append(np.array([r['dE'], r['dN']]) - b); cov.append(C)
        U.append(u); Tn.append(tn); KR.append(krad[s] / kref); KT.append(ktan[s] / kref); lab.append((s, t))
        z = np.load(f'{R}/results/inj_s{s}_{t}.npz')
        noise.append(z['emp'][0] - np.nanmean(z['emp'][0], 0))   # (NT, 2) zero-mean noise draws
obs = np.array(obs); cov = np.array(cov); U = np.array(U); Tn = np.array(Tn); KR = np.array(KR); KT = np.array(KT)
W = np.array([np.linalg.cholesky(np.linalg.inv(c)) for c in cov])   # whitening: r^T C^-1 r = |L^T r|^2
N = len(obs)
tidx = np.array([TYPES.index(t) for _, t in lab])


def model(p, kind):
    d = np.zeros((N, 2))
    if kind in ('M0', 'M2', 'MT'):
        d += p[:2]
    if kind in ('M1', 'M2'):
        c = np.array(p[-2:])[tidx]
        d += (c * KR)[:, None] * U + (c * KT)[:, None] * Tn
    return d


def resid(p, kind, o):
    r = o - model(p, kind)
    if kind == 'MT':
        rt = (r * Tn).sum(1)
        st = np.sqrt(np.array([t @ c @ t for t, c in zip(Tn, cov)]))
        return rt / st
    return np.concatenate([W[i].T @ r[i] for i in range(N)])


def fit(kind, o):
    p0 = {'M0': [0, 0], 'M1': [-2, -4], 'M2': [0, 0, -2, -4], 'MT': [0, 0]}[kind]
    f = least_squares(resid, p0, args=(kind, o))
    J = f.jac
    try:
        C = np.linalg.inv(J.T @ J)
    except np.linalg.LinAlgError:
        C = np.full((len(p0), len(p0)), np.nan)
    return f.x, float((f.fun ** 2).sum()), C


res = {}
for kind in ['M0', 'M1', 'M2', 'MT']:
    p, c2, C = fit(kind, obs)
    dof = (N * 2 if kind != 'MT' else N) - len(p)
    res[kind] = dict(p=p.tolist(), chi2=c2, dof=dof, formal_cov=C.tolist())
    print(f'{kind}: p = {np.round(p, 3)}  chi2 = {c2:.1f} / {dof}  formal sd {np.round(np.sqrt(np.diag(C)), 3)}')

# empirical calibration around the M2 best fit
rng = np.random.default_rng(7)
pM2 = np.array(res['M2']['p'])
truth = model(pM2, 'M2')
NS = int(os.environ.get('COMB_NS', 4000))
sims = {k: [] for k in ['M0', 'M1', 'M2', 'MT']}
chi = {k: [] for k in ['M0', 'M1', 'M2', 'MT']}
for _ in range(NS):
    o = truth + np.array([nz[rng.integers(len(nz))] for nz in noise])
    for kind in ['M2', 'MT', 'M1']:
        p, c2, _ = fit(kind, o)
        sims[kind].append(p); chi[kind].append(c2)
for k in ['M2', 'MT', 'M1']:
    sims[k] = np.array(sims[k]); chi[k] = np.array(chi[k])
# null for M1 (delta = 0): simulate from M1 best fit to get chi2 distribution of M1 when true
pM1 = np.array(res['M1']['p'])
truth1 = model(pM1, 'M1')
chi_M1_null = []; dM2_null = []
for _ in range(NS):
    o = truth1 + np.array([nz[rng.integers(len(nz))] for nz in noise])
    p, c2, _ = fit('M1', o); chi_M1_null.append(c2)
    p2, c22, _ = fit('M2', o); dM2_null.append(p2[:2])
chi_M1_null = np.array(chi_M1_null); dM2_null = np.array(dM2_null)
d = sims['M2'][:, :2]
Cd = np.cov(d.T)
out = dict(kref_mpx_per_mag=kref, krad=krad, ktan=ktan, labels=[f'S{s}_{t}' for s, t in lab], obs=obs.tolist(),
           cov=cov.tolist(), models=res)
out['M2_delta'] = pM2[:2].tolist()
out['M2_delta_emp_cov'] = Cd.tolist()
out['M2_delta_emp_sd'] = np.sqrt(np.diag(Cd)).tolist()
out['M2_c_emp_sd'] = sims['M2'][:, 2:].std(0).tolist()
out['M2_chi2_emp_pct'] = [float(np.percentile(chi['M2'], q)) for q in (16, 50, 84, 95, 99)]
out['MT_delta'] = res['MT']['p']
out['MT_delta_emp_sd'] = sims['MT'].std(0).tolist()
out['MT_delta_emp_cov'] = np.cov(sims['MT'].T).tolist()
out['M1_chi2_null_pct'] = [float(np.percentile(chi_M1_null, q)) for q in (16, 50, 84, 95, 99, 99.9)]
out['M1_chi2_obs_p_emp'] = float((chi_M1_null >= res['M1']['chi2']).mean())
# delta = 0 test under M2: Mahalanobis distance with empirical covariance, and from the M1-truth simulations
md = float(np.sqrt(pM2[:2] @ np.linalg.inv(Cd) @ pM2[:2]))
r_null = np.sqrt(np.einsum('ij,jk,ik->i', dM2_null, np.linalg.inv(np.cov(dM2_null.T)), dM2_null))
out['M2_delta0_mahalanobis'] = md
out['M2_delta0_p_emp'] = float((r_null >= md).mean())
sep = np.hypot(*d.T); sepT = np.hypot(*pM2[:2])
out['M2_sep'] = float(sepT); out['M2_sep_pct'] = [float(np.percentile(sep, q)) for q in (50, 68.3, 95, 99.7)]
# separation upper limits from the delta=0 truth simulations (what radius contains 95/99.7% if source on target)
sep0 = np.hypot(*dM2_null.T)
out['sep_if_on_target_pct'] = [float(np.percentile(sep0, q)) for q in (68.3, 95, 99.7)]
print('\nM2 delta =', np.round(pM2[:2], 2), 'emp sd', np.round(np.sqrt(np.diag(Cd)), 2), 'corr', round(Cd[0, 1] / np.sqrt(Cd[0, 0] * Cd[1, 1]), 2))
print('   colour levers c_deep, c_shallow =', np.round(pM2[2:], 2), '+-', np.round(sims['M2'][:, 2:].std(0), 2), '(arcsec radial at k_ref)')
print('   chi2 M2 obs', round(res['M2']['chi2'], 1), 'emp pct 16/50/84/95/99', np.round(out['M2_chi2_emp_pct'], 1))
print('   delta=0: Mahalanobis', round(md, 2), ' empirical p', out['M2_delta0_p_emp'])
print('   |delta| =', round(sepT, 2), '"; if on target, |delta_hat| 68/95/99.7% =', np.round(out['sep_if_on_target_pct'], 2))
print('MT delta =', np.round(res['MT']['p'], 2), 'emp sd', np.round(sims['MT'].std(0), 2))
print('M1 chi2 obs', round(res['M1']['chi2'], 1), ' null pct 16/50/84/95/99/99.9', np.round(out['M1_chi2_null_pct'], 1), ' p_emp', out['M1_chi2_obs_p_emp'])
print('M0 chi2 obs', round(res['M0']['chi2'], 1), '/', res['M0']['dof'])
np.savez_compressed(f'{R}/results/s05_sims.npz', M2=sims['M2'], MT=sims['MT'], M1=sims['M1'], dM2_null=dM2_null, chiM2=chi['M2'], chiM1null=chi_M1_null)
to_json(out, f'{R}/results/s05_combine.json')
