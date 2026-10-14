"""Empirical calibration of Gaia DR3 IPD diagnostics with high-resolution imaging (v01 sample).

Target: G = 11.78, RUWE 1.026, ipd_frac_multi_peak (fmp) = 41, ipd_gof_harmonic_amplitude (gofha) = 0.029,
phase 18.2 deg, ipd_frac_odd_win 0, no Gaia source within 18".

Primary = brightest Gaia DR3 source within 4" of the catalogue position (WDS positions can be ~1-3" off).
Pairs: companions measured at 0.03-6"; per primary, the 'IPD-relevant' companion = smallest delta-mag among those at
0.08-2.5" (median over repeat measurements, I-band for SOAR; i/880/692/LP600 for KOIs).
'Gaia-unresolved' = no other Gaia source within 3" of the primary (as for the target, whose nearest is 18.2").
Singles = targets with no companion in any catalogue (speckle/AO non-detections).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from v3lib import V3, to_json
from scipy.stats import beta

O = V3 + '/speckle'
R = pd.read_csv(f'{O}/rows.csv')
X = pd.read_csv(f'{O}/xm.csv')
X = X.rename(columns={'BP-RP': 'bprp', 'E(BP/RP)': 'bpexc'})
# primary per position
Xc = X[X.angDist < 4.0].sort_values(['pos', 'Gmag'])
prim = Xc.groupby('pos').first().reset_index()
# Gaia neighbours of the primary within 3" (from the 6" xmatch list around the same position)
nb = []
for pos, grp in X.groupby('pos'):
    p = prim[prim.pos == pos]
    if len(p) == 0:
        continue
    p = p.iloc[0]
    d = np.hypot((grp.RAdeg - p.RAdeg) * np.cos(np.radians(p.DEdeg)), grp.DEdeg - p.DEdeg) * 3600
    o = (grp.Source != p.Source)
    nb.append((pos, int((o & (d < 3.0)).sum()), float(d[o].min()) if o.any() else 99.0,
               float((grp.Gmag[o & (d < 3.0)]).min() - p.Gmag) if (o & (d < 3.0)).any() else np.nan))
nb = pd.DataFrame(nb, columns=['pos', 'n_gaia3', 'd_gaia_min', 'dG_gaia3'])
prim = prim.merge(nb, on='pos', how='left')
R = R.merge(prim, on='pos', how='inner')
print('rows with Gaia primary:', len(R), 'unique primaries:', R.Source.nunique())
# band choice
soar = R.src.str.contains('155/235|158/48|160/7|162/41|164/58|166/139|168/28|159/19|162/192')
R['use_dm'] = R.dm
R.loc[soar & ~R.filt.fillna('').str.strip().isin(['I', '']), 'use_dm'] = np.nan     # SOAR: I band only (Ziegler tables are I)
res = R[(R.kind == 'res') & np.isfinite(R.sep) & (R.sep > 0.03) & (R.sep < 6)]
# per primary: IPD-relevant companion
pairs = []
for sid, g in res.groupby('Source'):
    g2 = g[(g.sep >= 0.08) & (g.sep <= 2.5) & np.isfinite(g.use_dm)]
    if len(g2) == 0:
        continue
    # group repeat measurements of the same companion by separation (within 30%)
    g2 = g2.sort_values('use_dm')
    best = g2.iloc[0]
    same = g2[np.abs(g2.sep - best.sep) < 0.3 * best.sep + 0.05]
    pa = np.radians(same.pa.dropna().values * 2)
    pairs.append(dict(Source=sid, sep=float(same.sep.median()), dm=float(same.use_dm.median()), pa=float(np.degrees(np.arctan2(np.sin(pa).mean(), np.cos(pa).mean()) / 2) % 180) if len(pa) else np.nan,
                      nmeas=len(same), src=';'.join(sorted(set(same.src))), ncomp=int(len(g.sep.round(1).unique()))))
pairs = pd.DataFrame(pairs).merge(prim.drop_duplicates('Source'), on='Source', how='left')
resolved_ids = set(res.Source)
sing = R[(R.kind == 'unres') & ~R.Source.isin(resolved_ids)].drop_duplicates('Source')
print('pairs (companion 0.08-2.5"):', len(pairs), ' singles:', len(sing))


def rate(k, n):
    lo = beta.ppf(0.16, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(0.84, k + 1, n - k) if k < n else 1.0
    return dict(k=int(k), n=int(n), frac=(k / n if n else None), ci68=[float(lo), float(hi)])


out = {}
Gsel = lambda d: (d.Gmag > 10) & (d.Gmag < 13)
ru = lambda d: d.RUWE < 1.4
fmpT = lambda d: (d.IPDfmp >= 30) & (d.IPDfmp <= 50)
unres = lambda d: d.n_gaia3 == 0
# (1) singles
for lab, m in [('G10-13', Gsel(sing)), ('G10-13_ruwe<1.4', Gsel(sing) & ru(sing)), ('allG', np.ones(len(sing), bool))]:
    s = sing[m]
    out[f'singles_{lab}'] = dict(fmp_ge30=rate(((s.IPDfmp >= 30)).sum(), len(s)), fmp_30_50=rate(fmpT(s).sum(), len(s)),
                                fmp_ge10=rate((s.IPDfmp >= 10).sum(), len(s)),
                                gofha_q=[float(np.nanpercentile(s.IPDgofha, q)) for q in (50, 84, 95)],
                                fmp_q=[float(np.nanpercentile(s.IPDfmp, q)) for q in (50, 95, 99)],
                                blendBP_frac_med=float(np.nanmedian(s.NBPblend / s.o_BPmag)),
                                blendBP_frac_q=[float(np.nanpercentile(s.NBPblend / s.o_BPmag, q)) for q in (5, 16, 50, 84)])
    print(lab, 'singles', out[f'singles_{lab}'])
# (2) pairs like the target
sel_sets = {'G10-13,fmp30-50,ruwe<1.4': Gsel(pairs) & fmpT(pairs) & ru(pairs),
            'G10-13,fmp30-50,ruwe<1.4,gaia-unresolved': Gsel(pairs) & fmpT(pairs) & ru(pairs) & unres(pairs),
            'G10-13,fmp30-50,ruwe<1.4,gaia-unresolved,gofha<0.06': Gsel(pairs) & fmpT(pairs) & ru(pairs) & unres(pairs) & (pairs.IPDgofha < 0.06),
            'G9-14,fmp25-60,ruwe<1.4,gaia-unresolved': (pairs.Gmag > 9) & (pairs.Gmag < 14) & (pairs.IPDfmp >= 25) & (pairs.IPDfmp <= 60) & ru(pairs) & unres(pairs),
            'allG,fmp30-50,ruwe<1.4,gaia-unresolved': fmpT(pairs) & ru(pairs) & unres(pairs)}
for lab, m in sel_sets.items():
    p = pairs[m]
    out[f'pairs_{lab}'] = dict(n=int(len(p)), sep_q=[float(np.percentile(p.sep, q)) for q in (5, 16, 50, 84, 95)] if len(p) else None,
                               dm_q=[float(np.percentile(p.dm, q)) for q in (5, 16, 50, 84, 95)] if len(p) else None,
                               G_med=float(p.Gmag.median()) if len(p) else None,
                               list=p[['Source', 'Gmag', 'sep', 'dm', 'pa', 'IPDfmp', 'RUWE', 'IPDgofha', 'IPDgofhp', 'src']].round(3).to_dict('records') if len(p) <= 60 else None)
    print(f'\n{lab}: n={len(p)} sep q5/16/50/84/95 = {np.round(out[f"pairs_{lab}"]["sep_q"], 2) if len(p) else None}  dm = {np.round(out[f"pairs_{lab}"]["dm_q"], 2) if len(p) else None}')
# (3) likelihood grid: P(fmp in 30-50 & RUWE<1.4 & Gaia-unresolved | sep, dm), G 9.5-13.5
grid = []
seps = [0.08, 0.15, 0.25, 0.4, 0.6, 0.85, 1.2, 1.7, 2.5]
dms = [0, 1.5, 2.5, 3.5, 4.5, 6.5]
gsel = (pairs.Gmag > 9.5) & (pairs.Gmag < 13.5)
print('\nP(fmp30-50 & RUWE<1.4 & Gaia-unresolved | sep, dm) for G 9.5-13.5 pairs [k/n]:')
hdr = '        ' + ''.join(f'  dm{a:.1f}-{b:.1f}' for a, b in zip(dms[:-1], dms[1:]))
print(hdr)
for s0, s1 in zip(seps[:-1], seps[1:]):
    line = f'{s0:.2f}-{s1:.2f}'
    for d0, d1 in zip(dms[:-1], dms[1:]):
        m = gsel & (pairs.sep >= s0) & (pairs.sep < s1) & (pairs.dm >= d0) & (pairs.dm < d1)
        k = int((m & fmpT(pairs) & ru(pairs) & unres(pairs)).sum()); n = int(m.sum())
        kf = int((m & (pairs.IPDfmp >= 30)).sum())
        r = rate(k, n)
        grid.append(dict(sep=[s0, s1], dm=[d0, d1], n=n, k=k, frac=r['frac'], ci68=r['ci68'], k_fmp_ge30=kf,
                         fmp_med=float(pairs.IPDfmp[m].median()) if n else None, ruwe_med=float(pairs.RUWE[m].median()) if n else None,
                         frac_gaia_resolved=float((pairs.n_gaia3[m] > 0).mean()) if n else None))
        line += f'   {k:3d}/{n:<4d}  '
    print(line)
out['likelihood_grid_G9.5-13.5'] = grid
# (4) harmonic phase vs companion PA for Gaia-unresolved pairs with significant multi-peak / elongation
hp = pairs[(pairs.sep > 0.1) & (pairs.sep < 1.5) & unres(pairs) & np.isfinite(pairs.pa) & np.isfinite(pairs.IPDgofhp)]
for lab, m in [('gofha>0.05', hp.IPDgofha > 0.05), ('gofha<0.05&fmp>=20', (hp.IPDgofha < 0.05) & (hp.IPDfmp >= 20)),
               ('gofha<0.04', hp.IPDgofha < 0.04)]:
    h = hp[m]
    d = (h.pa - h.IPDgofhp) % 180
    # circular statistics on doubled angles
    ang = np.radians(2 * d)
    Rbar = np.hypot(np.cos(ang).mean(), np.sin(ang).mean()) if len(h) else np.nan
    mean = (np.degrees(np.arctan2(np.sin(ang).mean(), np.cos(ang).mean())) / 2) % 180 if len(h) else np.nan
    hist = np.histogram(d, bins=np.arange(0, 181, 22.5))[0].tolist()
    # Rayleigh p-value for doubled angles
    n = len(h); z = n * Rbar ** 2 if n else np.nan
    p = float(np.exp(-z) * (1 + (2 * z - z ** 2) / (4 * n) - (24 * z - 132 * z ** 2 + 76 * z ** 3 - 9 * z ** 4) / (288 * n ** 2))) if n > 5 else None
    # fraction within +-22.5 deg of the mean
    within = float(np.mean(np.abs(((d - mean + 90) % 180) - 90) < 22.5)) if n else None
    out[f'phase_vs_pa_{lab}'] = dict(n=n, mean_offset_deg=mean, Rbar=Rbar, rayleigh_p=p, hist_22p5deg=hist, frac_within_22p5=within)
    print(f'phase test {lab}: n={n} mean (PA - phase) mod 180 = {mean:.1f} deg, Rbar={Rbar:.2f}, Rayleigh p={p}, hist={hist}')
pairs.to_csv(f'{O}/pairs_gaia.csv', index=False); sing.to_csv(f'{O}/singles_gaia.csv', index=False)
to_json(out, f'{V3}/results/v02_calib.json')
