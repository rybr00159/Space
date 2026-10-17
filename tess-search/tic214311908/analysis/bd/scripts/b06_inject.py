"""b06: injection calibration of the shape statistics at FAKE EPOCHS.
Lambda_G = chi2_min(grazing BD on the M3V; D = 1, k in [0.19, 0.33]) - chi2_min(planet; D = 1, k free)
Lambda_C = chi2_min(BD/giant on a companion; k in [0.25, 1.2], D free) - chi2_min(planet)
Truths (ML fits of the real data, out/b05_fits.json, beta noise): 'planet', 'grazBD', and 'compBD40' (best companion fit
with k fixed at 0.40, i.e. a BD of ~0.9 R_J on a ~0.2 Msun companion).
REF: each trial draws, per sector, as many of the referee's fake windows as there are real windows (3,4,3,4,4,3),
     without replacement inside a trial. INV: each trial is one of the 48 re-detrended fake ephemerides (b04).
The truth model (exposure-integrated, same code as the fits) multiplies the fake-window flux. Each fake light curve is
fitted exactly as the real one (same noise model: per-sector OOT scale x b05 beta).
usage: b06_inject.py DS TRUTH NTRIAL SEED -> out/b06_inject_{ds}_{truth}.csv"""
import sys; sys.path.insert(0, 'tic214311908/analysis/bd/scripts')
from blib import *
from lclib import *
import fitlib as F
import pandas as pd
ds, truth, NT, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
rng = np.random.default_rng(seed)
NP = {'ref': 2, 'inv': 3}
FJ = jload(f'{BD}/out/b05_fits.json')
R = FJ['fits'][f'{ds}_beta']
bt = {int(s): v['beta'] for s, v in FJ['beta'][ds].items()}
if truth == 'compBD40':
    xt = jload(f'{BD}/out/b05_comp_k40.json')[ds]['x']; scen_t = 'compBD'
else:
    xt = R[truth]['x']; scen_t = truth
z = np.load(f'{BD}/data/{ds}_fake.npz')
NREAL = {12: 3, 39: 4, 66: 3, 93: 4, 103: 4, 104: 3}
rows = []
out = f'{BD}/out/b06_inject_{ds}_{truth}.csv'
for it in range(NT):
    if ds == 'ref':
        sel = np.zeros(len(z['t']), bool)
        for s, n in NREAL.items():
            wins = np.unique(z['win'][z['sec'] == s])
            sel |= np.isin(z['win'], rng.choice(wins, n, replace=False))
    else:
        j = it % int(z['fset'].max() + 1)
        sel = z['fset'] == j
    d = Data(z['t'][sel], z['f'][sel], z['e'][sel], z['sec'][sel], z['cad'][sel], z['win'][sel], z['tc'][sel], z['oot'][sel], npoly=NP[ds], beta=bt)
    s_t = F.model_s(d, np.array(xt), scen_t)
    d.f = d.f*(1 + s_t); d.y = d.f - 1
    rp = F.fit(d, 'planet')
    dT = rp['x'][0]
    rg = F.fit(d, 'grazBD', dT)
    rc = F.fit(d, 'compBD', dT)
    rows.append(dict(trial=it, truth=truth, chi2_p=rp['chi2'], chi2_g=rg['chi2'], chi2_c=rc['chi2'], Lambda_G=rg['chi2'] - rp['chi2'],
                     Lambda_C=rc['chi2'] - rp['chi2'], k_p=rp['k'], b_p=rp['b'], rho_p=rp['rho_c'], T14_p=rp['T14_h'], flat_p=rp['flat'],
                     k_g=rg['k'], rho_g=rg['rho_c'], k_c=rc['k'], D_c=rc['D'], rho_c=rc['rho_c']))
    print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(out, index=False)
