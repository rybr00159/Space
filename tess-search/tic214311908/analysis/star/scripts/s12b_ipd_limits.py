"""Which resolved companions do the Gaia DR3 image diagnostics of TIC 214311908 exclude?
Empirical calibration on high-resolution imaging samples:
  (1) the TIC 71316387 campaign sample (SOAR speckle Tokovinin 2018-2024, Ziegler 2020/2021 TOIs, Furlan 2017 KOIs,
      Ziegler 2017 Robo-AO KOIs) -> tic713/localise/v3/speckle/pairs_gaia.csv, singles_gaia.csv;
  (2) the Robo-AO M-dwarf survey (Lamman+20) -> data/mdw_imaging_sample.csv (s12a).
Target data D: ipd_frac_multi_peak = 0, ipd_gof_harmonic_amplitude = 0.0132, RUWE = 1.076, no other Gaia source
within 3" (nearest 3.13"), BP blended in 2/38 transits; G = 14.28.
'Clean like the target' = fmp <= 1 and gofha <= 0.025 and RUWE <= 1.2 and n_gaia3 == 0. Calibrators restricted to
G 11.5-16.5 (gofha and fmp sensitivity depend on brightness). For each (separation, contrast) cell:
P(clean | pair in cell) with a Jeffreys 68% interval, compared with P(clean | single) from imaging non-detections.
Cells with P(clean|pair) / P(clean|single) < 0.1 (and the 84% upper bound < 0.2) are 'excluded'; < 0.33 'disfavoured'.
Contrasts are in I / i' / LP600 / 692-880 nm (close to TESS); Delta G is larger for red companions (see README).
Output: out/ipd_limits.json, figs/ipd_limits.png"""
import json, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from scipy.stats import beta
B = 'tic214311908/analysis/star'
V3 = '<local>/tic713/localise/v3/speckle'
P = pd.read_csv(f'{V3}/pairs_gaia.csv'); S = pd.read_csv(f'{V3}/singles_gaia.csv')
for d in (P, S):
    d.rename(columns={'Gmag': 'G', 'IPDfmp': 'fmp', 'IPDgofha': 'gofha', 'RUWE': 'ruwe'}, inplace=True)
    d['src'] = 'speckle_compilation'
M = pd.read_csv(f'{B}/data/mdw_imaging_sample.csv').rename(columns={'phot_g_mean_mag': 'G', 'ipd_frac_multi_peak': 'fmp', 'ipd_gof_harmonic_amplitude': 'gofha'})
M['src'] = 'Lamman20'
MP = M[np.isfinite(M.comp_sep)].rename(columns={'comp_sep': 'sep', 'comp_dm': 'dm'})
MS = M[~M.resolved.astype(bool)]
PP = pd.concat([P[['G', 'sep', 'dm', 'fmp', 'gofha', 'ruwe', 'n_gaia3', 'bprp', 'src']], MP[['G', 'sep', 'dm', 'fmp', 'gofha', 'ruwe', 'n_gaia3', 'bp_rp', 'src']].rename(columns={'bp_rp': 'bprp'})], ignore_index=True)
SS = pd.concat([S[['G', 'fmp', 'gofha', 'ruwe', 'n_gaia3', 'bprp', 'src']], MS[['G', 'fmp', 'gofha', 'ruwe', 'n_gaia3', 'bp_rp', 'src']].rename(columns={'bp_rp': 'bprp'})], ignore_index=True)
gsel = lambda d: (d.G > 11.5) & (d.G < 16.5)
PP = PP[gsel(PP) & np.isfinite(PP.sep) & np.isfinite(PP.dm)]; SS = SS[gsel(SS)]
def clean(d, fmpmax=1, gmax=0.025, rmax=1.2):
    return (d.fmp <= fmpmax) & (d.gofha <= gmax) & (d.ruwe <= rmax) & (d.n_gaia3 == 0)
def jeff(k, n):
    return [float(beta.ppf(0.16, k + 0.5, n - k + 0.5)), float(beta.ppf(0.84, k + 0.5, n - k + 0.5))] if n else [np.nan, np.nan]
out = dict(n_pairs=int(len(PP)), n_singles=int(len(SS)), by_src_pairs=PP.src.value_counts().to_dict(), by_src_singles=SS.src.value_counts().to_dict())
ks, ns = int(clean(SS).sum()), len(SS)
Ls = (ks + 0.5) / (ns + 1)
out['single'] = dict(k=ks, n=ns, L=Ls, ci68=jeff(ks, ns))
# also single-star rate for red singles and for G 13.5-15.5
for lab, m in [('red_bprp>2', SS.bprp > 2), ('G13.5-15.5', (SS.G > 13.5) & (SS.G < 15.5))]:
    k, n = int(clean(SS[m]).sum()), int(m.sum()); out[f'single_{lab}'] = dict(k=k, n=n, L=(k + 0.5) / (n + 1), ci68=jeff(k, n))
print('P(clean|single):', out['single'], out['single_red_bprp>2'], out['single_G13.5-15.5'])
SEPS = [0.08, 0.15, 0.25, 0.4, 0.6, 0.9, 1.3, 2.0, 2.5]
DMS = [0, 1, 2, 3, 4, 5, 7]
cells = []
for s0, s1 in zip(SEPS[:-1], SEPS[1:]):
    line = f'{s0:.2f}-{s1:.2f}"'
    for d0, d1 in zip(DMS[:-1], DMS[1:]):
        m = (PP.sep >= s0) & (PP.sep < s1) & (PP.dm >= d0) & (PP.dm < d1)
        n = int(m.sum()); k = int(clean(PP[m]).sum())
        L = (k + 0.5) / (n + 1) if n else np.nan; ci = jeff(k, n)
        LR = L / Ls if n else np.nan; LRhi = ci[1] / Ls if n else np.nan
        status = 'no data' if n < 4 else ('excluded' if (LR < 0.1 and LRhi < 0.2) else ('disfavoured' if LR < 0.33 else 'allowed'))
        cells.append(dict(sep=[s0, s1], dm=[d0, d1], n=n, k=k, L=L, ci68=ci, LR=LR, LR_hi84=LRhi, status=status,
                          fmp_med=float(PP.fmp[m].median()) if n else None, gofha_med=float(PP.gofha[m].median()) if n else None,
                          ruwe_med=float(PP.ruwe[m].median()) if n else None, frac_gaia_resolved=float((PP.n_gaia3[m] > 0).mean()) if n else None))
        line += f'  {k:3d}/{n:<4d}{status[:4]:>5s}'
    print(line)
out['cells'] = cells
# separation in au at the target distance
out['dist_pc'] = 64.31
# figure
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(7, 4.5))
col = {'excluded': '#c0392b', 'disfavoured': '#e59866', 'allowed': '#7fb3d5', 'no data': '#dddddd'}
for c in cells:
    ax.add_patch(plt.Rectangle((c['sep'][0], c['dm'][0]), c['sep'][1] - c['sep'][0], c['dm'][1] - c['dm'][0], color=col[c['status']], ec='w'))
    if c['n']:
        ax.text(np.sqrt(c['sep'][0] * c['sep'][1]), np.mean(c['dm']), f"{c['k']}/{c['n']}", ha='center', va='center', fontsize=6)
ax.set_xscale('log'); ax.set_xlim(0.08, 2.5); ax.set_ylim(7, 0); ax.set_xlabel('separation (arcsec)'); ax.set_ylabel('contrast (I / i\' / 600-900 nm, mag)')
sec = ax.secondary_xaxis('top', functions=(lambda x: x * 64.31, lambda x: x / 64.31)); sec.set_xlabel('projected separation (au)')
ax.set_title('Gaia DR3 IPD/RUWE diagnostics: pairs as clean as TIC 214311908 (k/n)', fontsize=9)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=v, label=k) for k, v in col.items()], fontsize=7, loc='lower left')
fig.tight_layout(); fig.savefig(f'{B}/figs/ipd_limits.png', dpi=120)
# sensitivity: different thresholds
for lab, kw in [('strict', dict(fmpmax=0, gmax=0.02, rmax=1.15)), ('loose', dict(fmpmax=3, gmax=0.04, rmax=1.4))]:
    ks2, ns2 = int(clean(SS, **kw).sum()), len(SS); Ls2 = (ks2 + 0.5) / (ns2 + 1)
    ex = []
    for s0, s1 in zip(SEPS[:-1], SEPS[1:]):
        for d0, d1 in zip(DMS[:-1], DMS[1:]):
            m = (PP.sep >= s0) & (PP.sep < s1) & (PP.dm >= d0) & (PP.dm < d1); n = int(m.sum()); k = int(clean(PP[m], **kw).sum())
            if n >= 4:
                L = (k + 0.5) / (n + 1); ex.append(dict(sep=[s0, s1], dm=[d0, d1], k=k, n=n, LR=L / Ls2))
    out[f'variant_{lab}'] = dict(single=[ks2, ns2], cells=ex)
json.dump(out, open(f'{B}/out/ipd_limits.json', 'w'), indent=1, default=float)
# --- variant restricted to red calibrators (BP-RP > 1.8 for singles; pairs: combined BP-RP > 1.8), coarser grid
SEPR = [0.08, 0.2, 0.4, 0.8, 1.3, 2.5]; DMR = [0, 2, 3.5, 5, 7]
red_p = PP[PP.bprp > 1.8]; red_s = SS[SS.bprp > 1.8]
ksr, nsr = int(clean(red_s).sum()), len(red_s); Lsr = (ksr + 0.5) / (nsr + 1)
print(f'\nred calibrators: pairs {len(red_p)}, singles {nsr}, P(clean|single, red) = {Lsr:.3f}')
cr = []
for s0, s1 in zip(SEPR[:-1], SEPR[1:]):
    line = f'{s0:.2f}-{s1:.2f}"'
    for d0, d1 in zip(DMR[:-1], DMR[1:]):
        m = (red_p.sep >= s0) & (red_p.sep < s1) & (red_p.dm >= d0) & (red_p.dm < d1); n = int(m.sum()); k = int(clean(red_p[m]).sum())
        L = (k + 0.5) / (n + 1) if n else np.nan; ci = jeff(k, n)
        st = 'no data' if n < 4 else ('excluded' if (L / Lsr < 0.1 and ci[1] / Lsr < 0.2) else ('disfavoured' if L / Lsr < 0.33 else 'allowed'))
        cr.append(dict(sep=[s0, s1], dm=[d0, d1], n=n, k=k, L=L, LR=L / Lsr if n else np.nan, LR_hi84=ci[1] / Lsr if n else np.nan, status=st))
        line += f'  {k:3d}/{n:<4d}{st[:4]:>5s}'
    print(line)
out['red_variant'] = dict(single=dict(k=ksr, n=nsr, L=Lsr), cells=cr)
json.dump(out, open(f'{B}/out/ipd_limits.json', 'w'), indent=1, default=float)
