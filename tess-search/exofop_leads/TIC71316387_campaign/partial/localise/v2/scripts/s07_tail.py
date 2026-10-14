"""Empirical tail calibration of the separation exclusion.

For a hypothetical true source position delta_true (sep r at the PA that minimises the observed profile dchi2),
simulate data = delta_true + chromatic(c_best) + noise, where noise for each sector/type is a randomly chosen
injection trial (recovered - injected, empirical-image injection at the target, 500 trials per sector/type).
Statistic: profile dchi2(delta_true) = chi2(delta fixed at delta_true) - chi2_min.  p = fraction of simulations with
dchi2 >= observed dchi2(delta_true).  Also reports the noise kurtosis per sector/type.
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from s06_variants import build, fit, SECTORS, TYPES, R
from loclib import to_json
from scipy.stats import kurtosis

D = build()
pb, c2b, _ = fit(D)
z = np.load(f'{R}/results/s06_chi2map.npz'); g = z['g']; M = z['M']
EE, NN = np.meshgrid(g, g); rr = np.hypot(EE, NN)
noise = []
for s in SECTORS:
    for t in TYPES:
        e = np.load(f'{R}/results/inj_s{s}_{t}.npz')['emp'][0]
        noise.append(e - np.nanmean(e, 0))
kur = [float(np.mean(kurtosis(n, axis=0))) for n in noise]
frac3 = []
for n, C in zip(noise, D['cov']):
    m = np.einsum('ij,jk,ik->i', n, np.linalg.inv(C), n)
    frac3.append(float((m > 11.83).mean()))
print('noise excess kurtosis per sector/type:', np.round(kur, 2))
print('fraction of single-sector noise draws beyond the 2-D 3-sigma contour (Gaussian 0.0027):', np.round(frac3, 4))
rng = np.random.default_rng(11)
NS = int(os.environ.get('TAIL_NS', 2000))
out = dict(kurtosis=kur, frac_beyond_3sigma=frac3, tests=[])
for r0 in [0.0, 1.0, 1.25, 1.5, 2.0, 3.0]:
    if r0 == 0:
        dt = np.zeros(2)
    else:
        s = np.abs(rr - r0) < 0.03
        k = np.argmin(np.where(s, M, np.inf)); dt = np.array([EE.flat[k], NN.flat[k]])
    _, c2fix, _ = fit(D, fix_delta=dt)
    obs_d = c2fix - c2b
    c_best = pb[2:]
    tmap = D['ti']
    chrom = (c_best[tmap] * D['KR'])[:, None] * D['U'] + (c_best[tmap] * D['KT'])[:, None] * D['Tn']
    truth = dt + chrom
    stats = []
    for _ in range(NS):
        o = truth + np.array([n[rng.integers(len(n))] for n in noise])
        Ds = dict(D); Ds['obs'] = o
        p, c2m, _ = fit(Ds)
        _, c2f, _ = fit(Ds, fix_delta=dt)
        stats.append(c2f - c2m)
    stats = np.array(stats)
    p_emp = float((stats >= obs_d).mean())
    out['tests'].append(dict(sep=r0, delta_true=dt.tolist(), pa=float(np.degrees(np.arctan2(dt[0], dt[1])) % 360) if r0 else None,
                             obs_dchi2=float(obs_d), p_emp=p_emp, n_sims=NS, sim_dchi2_pct=[float(np.percentile(stats, q)) for q in (50, 95, 99.7)]))
    print(f'sep {r0:.2f}" at {dt.round(2)}: observed dchi2 {obs_d:6.1f}; sims median {np.median(stats):.2f}, 95% {np.percentile(stats,95):.2f}, '
          f'99.7% {np.percentile(stats,99.7):.2f}; p_emp = {p_emp:.4f} ({(stats >= obs_d).sum()}/{NS})', flush=True)
to_json(out, f'{R}/results/s07_tail.json')
