"""Robustness table over analysis variants (all diff_*.json files).
For each variant: per-sector s_par (bootstrap sd), f_N, dchi2(T-N), and an inverse-variance
(bootstrap) weighted combination of s_par with z-scores against target (0) and neighbour (1).
Outputs: variants_table.json, fig_variants.png
"""
import sys, os, glob, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import HERE, SECTORS, to_json

desc = {
    'prf_main': 'TESS PRF, r_fit=3 px, registration local r=4 (baseline)',
    'gauss_main': 'Gaussian PSF (fitted width), registration with N/T flux ratio fixed',
    'prf_pointing': 'PRF + per-pixel baseline decorrelated against measured pointing (dx,dy)',
    'gauss_pointing': 'Gaussian + pointing decorrelation',
    'prf_r2': 'PRF, fit radius 2 px', 'prf_r2.5': 'PRF, fit radius 2.5 px', 'prf_r4': 'PRF, fit radius 4 px',
    'prf_regext': 'PRF, registration from other stars only (target region masked)',
    'prf_regratio': 'PRF, registration with N/T flux ratio fixed to Gaia',
    'prf_regr2.5': 'PRF, registration local r=2.5 px',
    'prf_constbase': 'PRF, constant (not linear) OOT baseline',
    'prf_ootfar': 'PRF, OOT windows 1.5-3.0 T14 from mid-transit',
    'prf_in50': 'PRF, in-transit = central 50% of T14',
    'prf_plane': 'PRF, planar background in all fits',
    'prf_clean': 'PRF, windows with scattered-light change > 2 e/s/px dropped',
    'prf_clean_plane': 'PRF, clean windows + planar background',
}
rows = {}
for f in sorted(glob.glob(f'{HERE}/diff_*.json')):
    tag = os.path.basename(f)[5:-5]
    r = json.load(open(f))
    per = {}
    w = []; x = []
    for s in SECTORS:
        v = r[str(s)]
        st = v['stacked']; b = v['bootstrap']
        per[s] = dict(n=v['n_transits'], s_par=st['s_par'], s_sd=b['s_par']['std'], offT_E=st['offT_E'], offT_N=st['offT_N'],
                      f_N=st['f_N'], f_N_sd=b['f_N']['std'], dchi2=st['dchi2_TminusN'], perp=st['perp_arcsec'])
        sd = max(b['s_par']['std'], 0.03)
        w.append(1 / sd ** 2); x.append(st['s_par'])
    w = np.array(w); x = np.array(x)
    comb = float((w * x).sum() / w.sum()); err = float(1 / np.sqrt(w.sum()))
    chi2_consist = float((w * (x - comb) ** 2).sum())
    rows[tag] = dict(desc=desc.get(tag, tag), sectors=per, s_par_comb=comb, s_par_comb_err=err,
                     z_vs_target=comb / err, z_vs_neigh=(comb - 1) / err, chi2_between_sectors=chi2_consist,
                     dchi2_sum=float(sum(per[s]['dchi2'] for s in SECTORS)))
    print(f"{tag:16s} s_par_comb={comb:.3f}+-{err:.3f}  z_T={comb/err:5.1f} z_N={(comb-1)/err:+5.1f}  "
          f"chi2_sect={chi2_consist:5.1f}/4  sum dchi2={rows[tag]['dchi2_sum']:6.1f}  per-sector: " +
          ' '.join(f"{per[s]['s_par']:.2f}" for s in SECTORS))
to_json(rows, f'{HERE}/variants_table.json')
fig, ax = plt.subplots(figsize=(11, 0.45 * len(rows) + 1.5))
tags = list(rows)
for i, t in enumerate(tags):
    for k, s in enumerate(SECTORS):
        p = rows[t]['sectors'][s]
        ax.errorbar(p['s_par'], i + (k - 2) * 0.1, xerr=p['s_sd'], fmt='.', color=f'C{k}', alpha=0.6, label=f'S{s}' if i == 0 else None)
    ax.errorbar(rows[t]['s_par_comb'], i, xerr=rows[t]['s_par_comb_err'], fmt='D', color='k', ms=6, label='combined' if i == 0 else None)
ax.axvline(0, color='r', label='target'); ax.axvline(1, color='m', label='neighbour')
ax.set_yticks(range(len(tags))); ax.set_yticklabels([f'{t}' for t in tags], fontsize=8)
ax.set_xlabel('s_par (projection of difference-image position on target->neighbour axis)')
ax.legend(fontsize=7, ncol=4); ax.set_xlim(-0.5, 2.0)
ax.set_title('Robustness of the free-position localisation across analysis variants (bootstrap errors)')
plt.tight_layout(); plt.savefig(f'{HERE}/fig_variants.png', dpi=110); plt.close()
