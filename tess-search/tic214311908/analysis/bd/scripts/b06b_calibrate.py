"""b06b: calibrated likelihood ratios of the shape statistics from the fake-epoch injections (b06).
For each data set (REF, INV-clean) and alternative truth X (grazBD: Lambda_G; compBD40: Lambda_C):
  observed Lambda_obs (b05, beta noise; Lambda_C at k = 0.40 from b05b, the companion-hosted case that the injections
  represent, and the free-k Lambda_C);
  sampling distributions of Lambda under the planet truth and under X;
  LR_cal(X:planet) = p(Lambda_obs | X)/p(Lambda_obs | planet), with three density estimates: Gaussian KDE, Gaussian fit,
  and a conservative tail version [ (n_X(Lambda >= Lambda_obs) + 1)/(N_X + 1) ] / [ n_P(Lambda >= Lambda_obs)/N_P ];
  the adopted LR is the LARGEST of the three (most favourable to the false positive).
  LR_formal = exp(-Lambda_obs/2) (point hypotheses at the ML solutions). C_factor = LR_adopted/LR_formal.
Output: out/b06_calibration.json"""
import sys; sys.path.insert(0, 'tic214311908/analysis/bd/scripts')
from blib import *
import pandas as pd, os
from scipy.stats import gaussian_kde, norm
FJ = jload(f'{BD}/out/b05_fits.json'); K40 = jload(f'{BD}/out/b05_comp_k40.json')
out = {}
for ds in ('ref', 'inv'):
    out[ds] = {}
    fp = f'{BD}/out/b06_inject_{ds}_planet.csv'
    if not os.path.exists(fp): continue
    Pp = pd.read_csv(fp)
    for truth, stat, lobs in (('grazBD', 'Lambda_G', FJ['fits'][f'{ds}_beta']['Lambda_grazBD']),
                              ('compBD40', 'Lambda_C', FJ['fits'][f'{ds}_beta']['Lambda_compBD'])):
        fx = f'{BD}/out/b06_inject_{ds}_{truth}.csv'
        if not os.path.exists(fx): continue
        X = pd.read_csv(fx)
        xp, xx = Pp[stat].values, X[stat].values
        r = dict(N_planet=len(xp), N_X=len(xx), Lambda_obs=float(lobs),
                 planet_truth=dict(mean=float(np.mean(xp)), sd=float(np.std(xp)), p16=float(np.percentile(xp, 16)), p50=float(np.median(xp)), p84=float(np.percentile(xp, 84)),
                                   frac_ge_obs=float(np.mean(xp >= lobs))),
                 X_truth=dict(mean=float(np.mean(xx)), sd=float(np.std(xx)), p16=float(np.percentile(xx, 16)), p50=float(np.median(xx)), p84=float(np.percentile(xx, 84)),
                              max=float(np.max(xx)), frac_ge_obs=float(np.mean(xx >= lobs)), n_ge_obs=int(np.sum(xx >= lobs))))
        kp, kx = gaussian_kde(xp), gaussian_kde(xx)
        lr_kde = float(kx(lobs)[0]/max(kp(lobs)[0], 1e-300))
        lr_gau = float(norm.pdf(lobs, np.mean(xx), np.std(xx))/norm.pdf(lobs, np.mean(xp), np.std(xp)))
        tp = max(np.mean(xp >= lobs), 1/(len(xp) + 1))
        lr_tail = float(((np.sum(xx >= lobs) + 1)/(len(xx) + 1))/tp)
        lr_formal = float(np.exp(-0.5*lobs))
        adopted = max(lr_kde, lr_gau, lr_tail)
        r.update(LR_kde=lr_kde, LR_gauss=lr_gau, LR_tail_conservative=lr_tail, LR_formal=lr_formal, LR_adopted=adopted, C_factor=adopted/lr_formal,
                 p_value_X=float((np.sum(xx >= lobs) + 1)/(len(xx) + 1)))
        if truth == 'compBD40':
            lk40 = float(K40[ds]['dchi2'])
            r['Lambda_obs_k040'] = lk40
            r['note'] = 'Lambda_C is the free-k (k >= 0.25) statistic; Lambda at fixed k = 0.40 on the real data is Lambda_obs_k040'
        out[ds][truth] = r
        print(ds, truth, json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if not isinstance(v, dict)}), flush=True)
        print('   planet', r['planet_truth'], '\n   X', r['X_truth'], flush=True)
jdump(dict(method=__doc__, **out), f'{BD}/out/b06_calibration.json')
