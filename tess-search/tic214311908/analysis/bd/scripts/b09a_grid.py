"""b09a: light-curve likelihood grid for 'dark body + dilution' scenarios on the REAL data (REF and INV-clean; beta and
white noise). For each template (k, beta_b = b/(1+k), rho_c) the transit curve (unit dilution, M3V LD (U1, U2),
dT0 fixed at the planet ML value of each data set from b05) gives s; chi2(D) = Qyy - 2 D Qys + D^2 Qss with per-window
baselines profiled (lclib.Data.Q). Stored: Qys, Qss (float64) for every template and data/noise combination, plus Qyy and
the planet ML chi2 (b05) as reference.
Grid: ln k 0.05..1.2 (160), beta_b 0..0.6 step 0.02 + 0.6..0.998 step 0.004 (131), ln rho_c 0.4..100 (62).
Output: data/b09a_grid.npz"""
import sys; sys.path.insert(0, 'tic214311908/analysis/bd/scripts')
from blib import *
from lclib import *
import time
FJ = jload(f'{BD}/out/b05_fits.json')
NP = {'ref': 2, 'inv': 3}
D = {}
for ds in ('ref', 'inv'):
    bt = {int(s): v['beta'] for s, v in FJ['beta'][ds].items()}
    D[f'{ds}_beta'] = load_npz(f'{BD}/data/{ds}_real.npz', npoly=NP[ds], beta=bt)
    D[f'{ds}_white'] = load_npz(f'{BD}/data/{ds}_real.npz', npoly=NP[ds])
dT0 = {ds: FJ['fits'][f'{ds}_beta']['planet']['x'][0] for ds in ('ref', 'inv')}
ks = np.exp(np.linspace(np.log(0.05), np.log(1.2), 160))
bbs = np.r_[np.arange(0, 0.6, 0.02), np.arange(0.6, 0.9981, 0.004)]
rhos = np.exp(np.linspace(np.log(0.4), np.log(100.), 62))
shape = (len(ks), len(bbs), len(rhos))
QYS = {k: np.zeros(shape) for k in D}; QSS = {k: np.zeros(shape) for k in D}
QYY = {k: d.Q(d.y, d.y) for k, d in D.items()}
UY = {k: np.einsum('wij,wj->wi', d.Ainv, d.proj(d.y)) for k, d in D.items()}
t1 = time.time()
for i, k in enumerate(ks):
    for j, bb in enumerate(bbs):
        b = bb*(1 + k)
        for l, rho in enumerate(rhos):
            aR = float(aR_of_rho(rho))
            if b >= aR: continue
            cv = cad_curves(k, b, aR)
            for ds in ('ref', 'inv'):
                s = D[f'{ds}_beta'].template(cv, dT0[ds])
                if not np.any(s): continue
                for nm in ('beta', 'white'):
                    d = D[f'{ds}_{nm}']
                    pa = d.proj(s)
                    QYS[f'{ds}_{nm}'][i, j, l] = float(np.sum(d.w*d.y*s) - np.sum(UY[f'{ds}_{nm}']*pa))
                    QSS[f'{ds}_{nm}'][i, j, l] = float(np.sum(d.w*s*s) - np.einsum('wi,wij,wj->', pa, d.Ainv, pa))
    print(i, round(k, 4), round(time.time() - t1), 's', flush=True)
np.savez_compressed(f'{BD}/data/b09a_grid.npz', ks=ks, bbs=bbs, rhos=rhos, dT0_ref=dT0['ref'], dT0_inv=dT0['inv'],
                    **{f'QYS_{k}': v for k, v in QYS.items()}, **{f'QSS_{k}': v for k, v in QSS.items()},
                    **{f'QYY_{k}': np.array(v) for k, v in QYY.items()},
                    **{f'chi2ML_{k}': np.array(FJ['fits'][k]['planet']['chi2']) for k in QYY})
print('done', round(time.time() - t1), 's')
