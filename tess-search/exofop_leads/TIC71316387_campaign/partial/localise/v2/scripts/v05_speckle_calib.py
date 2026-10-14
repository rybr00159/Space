"""Empirical calibration of Gaia DR3 IPD diagnostics against SOAR speckle results.

Per Gaia primary (brightest Gaia DR3 source within 2.5" of the speckle position):
  * binaries: closest resolved speckle pair with 0.03" < sep < 3" (measurement epoch closest to 2016; I-band preferred);
    'gaia_resolved' if another Gaia source lies within max(0.3", 0.5*sep) of the predicted secondary position;
  * single controls: only 'unresolved' speckle rows (no resolved pair in any table) with Delta-m(1") >= 3.5.
Statistics in the target's regime (G 10-13.5, RUWE < 1.4).
"""
import numpy as np, pandas as pd, json, os
O = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v2/speckle'
R = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v2'
sp = pd.read_csv(f'{O}/speckle_rows.csv')
x = pd.read_csv(f'{O}/xmatch_gaia.csv')
x = x.rename(columns={'ra': 'ra_in', 'dec': 'dec_in'})
x['key'] = list(zip(x.ra_in.round(4), x.dec_in.round(4)))
sp['key'] = list(zip(sp.ra.round(4), sp.dec.round(4)))
# primary per key
prim = {}
allg = {}
for k, g in x.groupby('key'):
    g = g.sort_values('Gmag')
    allg[k] = g
    near = g[g.angDist < 2.5]
    if len(near):
        prim[k] = near.iloc[0]
sp['src'] = [int(prim[k].Source) if k in prim else -1 for k in sp.key]
sp = sp[sp.src > 0]
gcols = ['Source', 'RAdeg', 'DEdeg', 'Gmag', 'BP-RP', 'RUWE', 'IPDfmp', 'IPDfow', 'IPDgofha', 'IPDgofhp', 'epsi', 'sepsi', 'Plx', 'pmRA', 'pmDE',
         'E(BP/RP)', 'NBPblend', 'o_BPmag', 'NRPblend', 'o_RPmag', 'Dup', 'Nper', 'Solved']
G = pd.concat([prim[k][gcols].to_frame().T for k in prim]).drop_duplicates('Source').set_index('Source')
for c in gcols[1:]:
    G[c] = pd.to_numeric(G[c], errors='coerce')


def gaia_secondary(src, sep, pa):
    """Is there a Gaia source near the predicted secondary position?"""
    p = G.loc[src]
    for k, g in allg.items():
        pass
    return None


# index of all gaia sources by primary for the secondary check
by_src = {}
for k, g in allg.items():
    if k in prim:
        by_src.setdefault(int(prim[k].Source), []).append(g)
out = []
for src, g in sp.groupby('src'):
    res = g[g.resolved & (g.sep > 0.03) & (g.sep < 3.0)]
    unr = g[~g.resolved]
    rec = dict(src=src)
    if len(res):
        smin = res.sep.min()
        cl = res[res.sep < 1.3 * smin + 0.02].copy()
        cl['de'] = np.abs(cl.epoch - 2016.0)
        cl['pri'] = (cl.filt != 'I').astype(int)
        cl = cl.sort_values(['pri', 'de'])
        r0 = cl.iloc[0]
        dmI = cl[cl.filt == 'I'].dmag.median() if (cl.filt == 'I').any() else np.nan
        dmy = cl[cl.filt == 'y'].dmag.median() if (cl.filt == 'y').any() else np.nan
        rec.update(kind='binary', sep=r0.sep, pa=r0.pa, dmag=r0.dmag, filt=r0.filt, epoch=r0.epoch, dmI=dmI, dmy=dmy,
                   n_pairs=int(res.comp.nunique()), cat=r0['cat'])
        # Gaia-resolved?
        p = G.loc[src]
        gg = pd.concat(by_src.get(src, []))
        gg = gg[gg.Source != src].drop_duplicates('Source')
        if len(gg):
            dE = (gg.RAdeg - p.RAdeg) * np.cos(np.radians(p.DEdeg)) * 3600; dN = (gg.DEdeg - p.DEdeg) * 3600
            d = np.hypot(dE, dN)
            ok = np.abs(d - r0.sep) < max(0.3, 0.5 * r0.sep)
            rec['gaia_resolved'] = bool(ok.any())
            rec['nearest_gaia'] = float(d.min())
        else:
            rec['gaia_resolved'] = False; rec['nearest_gaia'] = np.nan
    elif len(unr):
        rec.update(kind='single', dm1=unr.dm1.max(), dm015=unr.dm015.max(), minsep=unr.minsep.min(), cat=unr['cat'].iloc[0])
        gg = pd.concat(by_src.get(src, []))
        gg = gg[gg.Source != src].drop_duplicates('Source')
        rec['nearest_gaia'] = float(gg.angDist.min()) if len(gg) else np.nan
    out.append(rec)
S = pd.DataFrame(out).set_index('src').join(G, how='left')
S.to_csv(f'{O}/calib_sample.csv')
print(S.kind.value_counts())

summary = {}
reg = (S.Gmag > 10) & (S.Gmag < 13.5) & (S.RUWE < 1.4)
B = S[(S.kind == 'binary') & reg & ~S.gaia_resolved.astype(bool)]
C = S[(S.kind == 'single') & reg & (S.dm1 >= 3.5)]
print(f'regime G 10-13.5, RUWE<1.4: binaries (Gaia-unresolved) {len(B)}, single controls {len(C)}')
for lo, hi in [(0, 1), (1, 10), (10, 20), (20, 30), (30, 50), (50, 101)]:
    nb = ((B.IPDfmp >= lo) & (B.IPDfmp < hi)).sum(); nc = ((C.IPDfmp >= lo) & (C.IPDfmp < hi)).sum()
    print(f'  IPDfmp {lo:3d}-{hi:3d}: binaries {nb:4d}  singles {nc:4d}')
summary['n_binaries_regime'] = int(len(B)); summary['n_singles_regime'] = int(len(C))
summary['singles_frac_fmp_ge30'] = float((C.IPDfmp >= 30).mean()); summary['singles_n_fmp_ge30'] = int((C.IPDfmp >= 30).sum())
summary['singles_frac_fmp_ge20'] = float((C.IPDfmp >= 20).mean())
summary['singles_frac_fmp_ge10'] = float((C.IPDfmp >= 10).mean())
summary['binaries_frac_fmp_ge30'] = float((B.IPDfmp >= 30).mean())
sel = B[(B.IPDfmp >= 30) & (B.IPDfmp <= 50)]
print(f'\nBinaries with IPDfmp 30-50 (G 10-13.5, RUWE<1.4, Gaia-unresolved): n={len(sel)}')
print(sel[['sep', 'dmag', 'filt', 'dmI', 'Gmag', 'RUWE', 'IPDfmp', 'IPDgofha', 'IPDgofhp', 'pa', 'NBPblend', 'o_BPmag', 'cat']].sort_values('sep').to_string())
summary['fmp30_50'] = dict(n=int(len(sel)), sep_q=[float(q) for q in np.nanpercentile(sel.sep, [5, 16, 50, 84, 95])] if len(sel) else None,
                           dmag_q=[float(q) for q in np.nanpercentile(sel.dmag, [5, 16, 50, 84, 95])] if len(sel) else None)
# broader: all G, RUWE<1.4
sel2 = S[(S.kind == 'binary') & (S.RUWE < 1.4) & (S.IPDfmp >= 30) & (S.IPDfmp <= 50) & ~S.gaia_resolved.astype(bool)]
summary['fmp30_50_allG'] = dict(n=int(len(sel2)), sep_q=[float(q) for q in np.nanpercentile(sel2.sep, [5, 16, 50, 84, 95])],
                                dmag_q=[float(q) for q in np.nanpercentile(sel2.dmag, [5, 16, 50, 84, 95])])
print('\nall G, RUWE<1.4, IPDfmp 30-50, Gaia-unresolved: n', len(sel2), 'sep pct', summary['fmp30_50_allG']['sep_q'], 'dmag pct', summary['fmp30_50_allG']['dmag_q'])
json.dump(summary, open(f'{R}/results/v05_summary.json', 'w'), indent=1)
