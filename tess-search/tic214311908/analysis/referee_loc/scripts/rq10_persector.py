"""Per-sector consistency of the real free-fit positions with M3V injections (baseline run): Mahalanobis distance
of the real per-sector free fit from the M3V-injection cloud (robust covariance), the fraction of injections farther
out, and per-sector chi2(pairA)-chi2(M3V) / chi2(G175)-chi2(M3V) percentiles under M3V and under the alternative.
Output: out/persector.json"""
import sys; sys.path.insert(0, 'tic214311908/analysis/referee_loc/scripts')
from rl_lib import *
d = json.load(open(f'{ROOT}/out/fit_base.json')); Z = np.load(f'{ROOT}/out/inj_base.npz')
MAIN = d['main']; out = {}
for si, s in enumerate(d['secs']):
    r = d['real']['per_sector'][str(s)]
    fr = np.array(r['free_pos'])
    F = Z['M3V__free'][:, si, :].astype(float)
    med = np.median(F, 0); C = np.cov((F - med)[np.all(np.abs(F - med) < 15, 1)].T)
    md = lambda x: float(np.sqrt((x - med) @ np.linalg.inv(C) @ (x - med)))
    mds = np.array([md(x) for x in F])
    o = dict(free_real=fr.tolist(), inj_median=med.tolist(), inj_sd=np.sqrt(np.diag(C)).tolist(), mahal_real=md(fr),
             frac_inj_mahal_ge=float(np.mean(mds >= md(fr))))
    for h in ('pairA', 'Mira', 'G175'):
        ii = MAIN.index(h)
        x3 = Z['M3V__chi'][:, si, ii] - Z['M3V__chi'][:, si, 0]; xh = Z[f'{h}__chi'][:, si, ii] - Z[f'{h}__chi'][:, si, 0]
        o[h] = dict(real=r['dchi2'][h], frac_M3Vinj_le=float(np.mean(x3 <= r['dchi2'][h])), frac_trueinj_ge=float(np.mean(xh >= r['dchi2'][h])))
    out[s] = o
    print(s, f"free {np.round(fr,1)} M3V-inj median {np.round(med,1)} sd {np.round(np.sqrt(np.diag(C)),1)}  P(farther) {o['frac_inj_mahal_ge']:.3f}  " +
          '  '.join(f"{h}: {o[h]['real']:+.1f} (M3V-inj<= {o[h]['frac_M3Vinj_le']:.2f}, true-inj>= {o[h]['frac_trueinj_ge']:.3f})" for h in ('pairA', 'Mira', 'G175')))
jdump(out, f'{ROOT}/out/persector.json')
