"""b12: diagnostic figures. (a) injection distributions of Lambda_G (grazing BD minus planet chi2) under planet and
grazing-BD truths with the observed value, REF and INV-clean; (b) profile of Lambda_G versus rho_circ (where a grazing
BD would have to sit) with the M3V density and the e needed; (c) same for Lambda_C under planet / compBD40 truths.
Output: figs/b12_shape_calibration.png"""
import sys; sys.path.insert(0, 'tic214311908/analysis/bd/scripts')
from blib import *
import pandas as pd, os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
C1, C2, C3 = '#2a78d6', '#eb6834', '#1baf7a'
INK, INK2, GRID, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.facecolor': SURF, 'figure.facecolor': SURF, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6,
                     'axes.spines.top': False, 'axes.spines.right': False})
FJ = jload(f'{BD}/out/b05_fits.json')
fig, ax = plt.subplots(2, 3, figsize=(11, 6.4))
for i, ds in enumerate(('ref', 'inv')):
    a = ax[i, 0]
    for truth, col, lab in (('planet', C1, 'planet truth'), ('grazBD', C2, 'grazing-BD truth')):
        fp = f'{BD}/out/b06_inject_{ds}_{truth}.csv'
        if os.path.exists(fp):
            x = pd.read_csv(fp).Lambda_G.values
            a.hist(x, bins=np.linspace(-6, 40, 47), color=col, alpha=0.75, label=f'{lab} (n={len(x)})', edgecolor=SURF, linewidth=1)
    lo = FJ['fits'][f'{ds}_beta']['Lambda_grazBD']
    a.axvline(lo, color=INK, lw=2); a.text(lo + 0.5, a.get_ylim()[1]*0.9, f'observed {lo:.1f}', color=INK, fontsize=8)
    a.set_xlabel('Lambda_G = chi2(grazing BD) - chi2(planet)'); a.set_ylabel('injections'); a.set_title(f'{ds.upper()}: grazing BD vs planet', color=INK, fontsize=9, loc='left')
    a.legend(frameon=False, fontsize=7)
    a = ax[i, 1]
    pr = FJ['fits'][f'{ds}_beta']['profile_rho_grazBD']
    r = np.array([p['rho_c'] for p in pr]); d = np.array([p['dchi2'] for p in pr])
    a.plot(r, d, color=C2, lw=2, marker='o', ms=4)
    a.axvspan(RHO_STAR - RHO_ERR, RHO_STAR + RHO_ERR, color=C1, alpha=0.15, lw=0)
    a.text(RHO_STAR*0.97, d.max()*0.85, 'M3V density', color=INK2, fontsize=8, ha='right')
    a.set_xscale('log'); a.set_xlabel('rho_circ of the grazing-BD fit (g/cm3)'); a.set_ylabel('Lambda_G (profile)')
    a.set_title(f'{ds.upper()}: grazing BD needs rho_circ ~ 1.7', color=INK, fontsize=9, loc='left')
    a = ax[i, 2]
    for truth, col, lab in (('planet', C1, 'planet truth'), ('compBD40', C3, 'k = 0.4 companion-host truth')):
        fp = f'{BD}/out/b06_inject_{ds}_{truth}.csv'
        if os.path.exists(fp):
            x = pd.read_csv(fp).Lambda_C.values
            a.hist(x, bins=np.linspace(-4, 12, 33), color=col, alpha=0.75, label=f'{lab} (n={len(x)})', edgecolor=SURF, linewidth=1)
    lo = FJ['fits'][f'{ds}_beta']['Lambda_compBD']
    a.axvline(lo, color=INK, lw=2); a.text(lo + 0.3, a.get_ylim()[1]*0.9, f'observed {lo:.1f}', color=INK, fontsize=8)
    a.set_xlabel('Lambda_C = chi2(k>=0.25 on companion) - chi2(planet)'); a.set_ylabel('injections')
    a.set_title(f'{ds.upper()}: companion-hosted large body vs planet', color=INK, fontsize=9, loc='left'); a.legend(frameon=False, fontsize=7)
fig.tight_layout()
os.makedirs(f'{BD}/figs', exist_ok=True)
fig.savefig(f'{BD}/figs/b12_shape_calibration.png', dpi=130)
print('saved')
