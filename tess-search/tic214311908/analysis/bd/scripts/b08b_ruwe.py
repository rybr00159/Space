"""b08b: RUWE likelihood ratio for a BD on the M3V with a CLEAN single-star comparison sample (the b08 sample of
../referee_star/data/mdw_cmd.csv includes unresolved binaries, median RUWE 1.15, which makes the target's 1.076 look
'low' and penalises any extra wobble). Clean sample: |G - 14.28| < 0.6, ipd_frac_multi_peak = 0, non_single_star = 0,
RUWE < 1.4. Same wobble model as b08. Updates out/b08_gaia.json: 'ruwe' -> clean sample (adopted), 'ruwe_allsample'
-> the b08 version."""
import sys; sys.path.insert(0, 'tic214311908/analysis/bd/scripts')
from blib import *
import pandas as pd
from scipy.stats import gaussian_kde
J = jload(f'{BD}/out/b08_gaia.json')
if 'ruwe_allsample' not in J: J['ruwe_allsample'] = J['ruwe']
m = pd.read_csv(f'{B}/referee_star/data/mdw_cmd.csv')
rs = m[(np.abs(m.phot_g_mean_mag - 14.281) < 0.6) & (m.ipd_frac_multi_peak == 0) & (m.non_single_star == 0) & (m.ruwe < 1.4)].ruwe.values
kr = gaussian_kde(np.log(rs)); pr = lambda x: kr(np.log(x))/x
RU, plx = 1.076204, 15.497
tab = {}
for sccd in (0.19, 0.22, 0.25):
    row = []
    for Mb in (13, 20, 30, 40, 50, 60, 70, 78.5):
        q = Mb/1047.57; a1 = 0.0448*q/(0.370 + q)*plx
        lrs = []
        for cth in (0.15, 0.25, 0.35):
            x0 = RU**2 - a1**2/2*cth*2/sccd**2
            lrs.append(float(pr(np.sqrt(x0))[0]*RU/np.sqrt(x0)/pr(RU)[0]))
        row.append(dict(M_bd=float(Mb), alpha_mas=float(a1), LR_by_projection=lrs, LR_mean=float(np.mean(lrs))))
    tab[f'sigma_ccd_{sccd}'] = row
    print(sccd, [(r['M_bd'], round(r['LR_mean'], 3)) for r in row])
J['ruwe'] = dict(sample='clean single-star comparison (b08b)', n_comparison=int(len(rs)), ruwe_pct={p: float(np.percentile(rs, p)) for p in (16, 50, 84, 95)},
                 target=RU, target_percentile=float(np.mean(rs < RU)*100), table=tab)
jdump(J, f'{BD}/out/b08_gaia.json')
