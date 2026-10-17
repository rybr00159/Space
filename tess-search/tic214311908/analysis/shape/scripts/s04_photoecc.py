"""s04: stellar-density (photo-eccentric) test and host-model comparison.
Inputs: s03 posteriors. For each data set:
1. rho_circ posterior (free fit, prior uniform in log10 rho on [-1, 2.3]) vs the single-star value N(9.41, 0.88)
   and the twin value N(15.5, 1.45): marginal-likelihood ratios E_post[p_host(rho)/prior(rho)] (circular orbits,
   ignoring the shape change of the twin dilution); the full nested-sampling evidences of the s03 'single', 'twin',
   'single_ecc', 'twin_ecc' models (which include the dilution of the twin case) are also tabulated.
2. Photo-eccentric posterior on e (Dawson & Johnson 2012): (e, w) ~ uniform (e < 0.9), rho_true ~ host, weight =
   smoothed-histogram density (0.03 dex) of the rho_circ posterior (in log rho) at rho_true g^3 / prior, g = (1 + e sin w)/sqrt(1 - e^2).
   Minimum eccentricity needed: g = (rho_circ/rho_host)^(1/3) -> e_min = |1 - g^2|/(1 + g^2).
Output: out/photoecc.json"""
import sys, json, numpy as np
sys.path.insert(0, 'tic214311908/analysis/shape/scripts')
from shapelib import *
from scipy.stats import gaussian_kde, norm
rng = np.random.default_rng(11)
HOST = dict(single=(9.41, 0.88), twin=(15.5, 1.45), single_feh=(9.16, 0.85))
out = {}
for ds in ['inv', 'ref']:
    z = np.load(f'{SD}/out/fit_{ds}_free.npz')
    lr = np.log10(z['rho_circ'])
    # density of log10 rho_circ: Gaussian-smoothed histogram (fast; bandwidth 0.03 dex)
    edges = np.linspace(-1, 2.3, 661); hd, _ = np.histogram(lr, edges, density=True)
    from scipy.ndimage import gaussian_filter1d
    hd = gaussian_filter1d(hd, 0.03/(edges[1]-edges[0]))
    cen = 0.5*(edges[1:] + edges[:-1])
    kde = lambda x, _c=cen, _h=hd: np.interp(x, _c, _h, left=0, right=0)
    prior_dens = 1/3.3  # per unit log10 rho
    res = dict(rho_circ_median=float(np.median(z['rho_circ'])), rho_circ_16_84=list(map(float, np.percentile(z['rho_circ'], [16, 84]))),
               rho_circ_2p5_97p5=list(map(float, np.percentile(z['rho_circ'], [2.5, 97.5]))), rho_circ_0p15_99p85=list(map(float, np.percentile(z['rho_circ'], [0.135, 99.865]))))
    for h, (mu, sd) in HOST.items():
        rr = rng.normal(mu, sd, 200000); rr = rr[rr > 0.5]
        # circular: marginal likelihood ratio relative to the free fit
        L = np.mean(kde(np.log10(rr))/prior_dens)
        res[f'{h}_circ_marglike_ratio'] = float(L)
        res[f'{h}_frac_post_within'] = float(np.mean(np.abs(z['rho_circ'] - mu) < sd))
        # PIT / tension: fraction of posterior beyond host value
        res[f'{h}_post_frac_below'] = float(np.mean(z['rho_circ'] < mu))
        # eccentric
        N = 400000
        e = rng.uniform(0, 0.9, N); w = rng.uniform(0, 2*np.pi, N); rt = rng.normal(mu, sd, N)
        ok = rt > 0.5; e, w, rt = e[ok], w[ok], rt[ok]
        g = (1 + e*np.sin(w))/np.sqrt(1 - e**2)
        rc = rt*g**3
        inside = (np.log10(rc) > -1) & (np.log10(rc) < 2.3)
        wt = np.where(inside, kde(np.log10(np.clip(rc, 0.1, 199)))/prior_dens, 0)
        wt /= wt.sum()
        idx = rng.choice(len(e), 20000, p=wt)
        es = e[idx]
        res[f'{h}_ecc_marglike_ratio'] = float(np.mean(np.where(inside, kde(np.log10(np.clip(rc, 0.1, 199)))/prior_dens, 0)))
        res[f'{h}_e_post'] = dict(median=float(np.median(es)), p16=float(np.percentile(es, 16)), p84=float(np.percentile(es, 84)), p95=float(np.percentile(es, 95)),
                                  P_e_lt_0p1=float(np.mean(es < 0.1)), P_e_lt_0p2=float(np.mean(es < 0.2)))
        gneed = (np.median(z['rho_circ'])/mu)**(1/3.)
        res[f'{h}_g_needed_median'] = float(gneed)
        res[f'{h}_e_min_median'] = float(abs(1 - gneed**2)/(1 + gneed**2))
        for q in [2.5, 97.5]:
            gq = (np.percentile(z['rho_circ'], q)/mu)**(1/3.)
            res[f'{h}_e_min_at_rho_p{q}'] = float(abs(1 - gq**2)/(1 + gq**2))
    res['bayes_single_over_twin_circ_from_rho_only'] = res['single_circ_marglike_ratio']/res['twin_circ_marglike_ratio']
    # nested-sampling evidences
    ev = {}
    for m in ['free', 'freeLD', 'single', 'twin', 'single_ecc', 'twin_ecc']:
        try:
            j = json.load(open(f'{SD}/out/fit_{ds}_{m}.json')); ev[m] = dict(logz=j['logz'], logzerr=j['logzerr'], chi2_best=j['chi2_best'],
                                                                               rho=j['post']['rho']['med'], rho_circ=j['post']['rho_circ']['med'], b=j['post']['b']['med'],
                                                                               k=j['post']['k']['med'], T14_h=j['post']['T14_h']['med'], e=j['post']['e']['med'])
        except FileNotFoundError:
            pass
    res['evidence'] = ev
    if 'single' in ev and 'twin' in ev:
        res['lnB_single_twin_circ'] = ev['single']['logz'] - ev['twin']['logz']
    if 'single_ecc' in ev and 'twin_ecc' in ev:
        res['lnB_single_twin_ecc'] = ev['single_ecc']['logz'] - ev['twin_ecc']['logz']
    if 'single' in ev and 'single_ecc' in ev:
        res['lnB_single_circ_over_ecc'] = ev['single']['logz'] - ev['single_ecc']['logz']
    out[ds] = res
    print(ds, json.dumps(res, indent=1))
json.dump(out, open(f'{SD}/out/photoecc.json', 'w'), indent=1)
