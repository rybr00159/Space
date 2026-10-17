"""r01: pixel-localisation likelihood ratios LR_X = p(s_obs | X dims) / p(s_obs | M3V dims), s = chi2(X) - chi2(M3V),
re-derived from the referee_loc injection sets (1000 injections per source and analysis variant) with several density
models, because s_obs lies far outside every injection on the alternative sources (0/1000), so the LR is an
extrapolation and its size is model-dependent:
  gauss  : Gaussian fit to each injection cloud (analyst's 'Gaussian')
  t      : Student-t ML fit (analyst's adopted 't')
  gpd    : generalised-Pareto fit to the upper 10% tail of the X-truth cloud (extreme-value extrapolation)
  sysinfl: Gaussian with the X-truth width inflated in quadrature by the variant-to-variant scatter of the real s_obs
           (a crude allowance for analysis systematics, which injections do not sample)
For each available variant (inj_*.npz with a matching fit_*.json) the LR is computed with that variant's own
injections and s_obs. Output: out/r01_loc_lr.json"""
import numpy as np, glob, os, sys, json
from scipy import stats
sys.path.insert(0, 'tic214311908/analysis/referee2_fpp/scripts')
from rlib import *
SRC = ['M3V', 'pairA', 'pairB', 'Mira', 'G175', 'G177', 'G170']
summ = jload('referee_loc/out/summary.json')
files = sorted(glob.glob(f'{B}/referee_loc/out/inj_*.npz'))
res = {}
real_by_var = {X: [] for X in SRC[1:]}
for f in files:
    tag = os.path.basename(f)[4:-4]
    key = 'base' if tag == 'base' else tag
    z = np.load(f)
    if 'real_comb_dchi2' not in z.files: continue
    real = z['real_comb_dchi2']
    out = {}
    for j, X in enumerate(SRC):
        if X == 'M3V' or f'{X}__comb_c' not in z.files: continue
        sX = (z[f'{X}__comb_c'][:, j] - z[f'{X}__comb_c'][:, 0]).astype(float)
        sM = (z['M3V__comb_c'][:, j] - z['M3V__comb_c'][:, 0]).astype(float)
        so = float(real[j]); real_by_var[X].append(so)
        lM = stats.norm(*stats.norm.fit(sM)).logpdf(so)
        lX_g = stats.norm(*stats.norm.fit(sX)).logpdf(so)
        tX = stats.t.fit(sX); tM = stats.t.fit(sM)
        lX_t = stats.t(*tX).logpdf(so); lM_t = stats.t(*tM).logpdf(so)
        # GPD on the upper 10% of the X cloud: density at so = 0.1 * gpd.pdf(so - u)
        u = np.percentile(sX, 90); exc = sX[sX > u] - u
        c, loc, sc = stats.genpareto.fit(exc, floc=0)
        lX_gpd = np.log(0.1) + stats.genpareto(c, 0, sc).logpdf(so - u) if (c >= 0 or so - u < -sc/c) else -np.inf
        out[X] = dict(s_obs=so, X_mean=float(sX.mean()), X_sd=float(sX.std()), z_gauss=float((so - sX.mean())/sX.std()),
                      M_pct=float(np.mean(sM <= so)), nu_t_X=float(tX[0]), gpd_shape=float(c), gpd_scale=float(sc),
                      log10LR_gauss=float((lX_g - lM)/np.log(10)), log10LR_t=float((lX_t - lM_t)/np.log(10)),
                      log10LR_gpd=float((lX_gpd - lM)/np.log(10)))
    res[tag] = out
# systematics-inflated Gaussian, base injections + variant scatter of s_obs
z = np.load(f'{B}/referee_loc/out/inj_base.npz')
sysinfl = {}
for j, X in enumerate(SRC):
    if X == 'M3V': continue
    sX = (z[f'{X}__comb_c'][:, j] - z[f'{X}__comb_c'][:, 0]).astype(float)
    sM = (z['M3V__comb_c'][:, j] - z['M3V__comb_c'][:, 0]).astype(float)
    vals = [v['real'][X] for k, v in summ.items() if isinstance(v, dict) and 'real' in v and X in v['real']]
    sd_sys = float(np.std(vals)); so_min = float(np.min(vals))
    s_tot = np.hypot(sX.std(), sd_sys)
    lX = stats.norm(sX.mean(), s_tot).logpdf(z['real_comb_dchi2'][j])
    lM = stats.norm(sM.mean(), np.hypot(sM.std(), sd_sys)).logpdf(z['real_comb_dchi2'][j])
    lX_min = stats.norm(sX.mean(), s_tot).logpdf(so_min)
    lM_min = stats.norm(sM.mean(), np.hypot(sM.std(), sd_sys)).logpdf(so_min)
    sysinfl[X] = dict(n_variants=len(vals), s_obs_range=[float(np.min(vals)), float(np.max(vals))], sd_between_variants=sd_sys,
                      log10LR_sysinfl=float((lX - lM)/np.log(10)), log10LR_sysinfl_at_min_variant=float((lX_min - lM_min)/np.log(10)))
summary = {}
for X in SRC[1:]:
    rows = [r[X] for r in res.values() if X in r]
    summary[X] = {k: dict(base=res['base'][X][k], min=float(np.min([r[k] for r in rows])), median=float(np.median([r[k] for r in rows])),
                          max=float(np.max([r[k] for r in rows])))
                  for k in ('s_obs', 'z_gauss', 'log10LR_gauss', 'log10LR_t', 'log10LR_gpd')}
    summary[X]['n_variants_with_injections'] = len(rows)
    summary[X]['sysinfl'] = sysinfl[X]
    print(X, json.dumps(summary[X], default=float))
jdump(dict(method=__doc__.split('\n')[0], per_variant=res, summary=summary), 'out/r01_loc_lr.json')
