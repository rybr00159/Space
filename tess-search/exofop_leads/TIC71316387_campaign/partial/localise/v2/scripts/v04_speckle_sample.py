"""Unified SOAR speckle sample (resolved pairs + unresolved stars) cross-matched to Gaia DR3 (CDS XMatch, I/355/gaiadr3).
Output: speckle/speckle_rows.csv (all speckle rows normalised) and speckle/xmatch_gaia.csv (all Gaia sources within 6")."""
import glob, os
import numpy as np
import pandas as pd
from astropy.table import Table
import astropy.units as u
from astroquery.xmatch import XMatch
import warnings; warnings.filterwarnings('ignore')

O = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v2/speckle'
rows = []


def col(t, names):
    for n in names:
        if n in t.colnames:
            return np.array(t[n].filled(np.nan) if hasattr(t[n], 'filled') else t[n], dtype=float)
    return np.full(len(t), np.nan)


def scol(t, names):
    for n in names:
        if n in t.colnames:
            return np.array([str(x) for x in t[n]])
    return np.array([''] * len(t))


spec = {
    # file: (kind, sep cols, pa cols, dmag cols, epoch cols, filter col, minsep cols, dm015, dm1, id cols)
    'J_AJ_159_19_table6': ('mixed', ['rho'], ['theta'], ['Dmag'], ['Year'], 'Filt', ['minrho'], ['Dm0.15'], ['Dm1.0'], ['TIC']),
    'J_AJ_162_192_table1': ('mixed', ['Sep'], ['PA'], ['Con'], ['Obs'], 'Filt', ['MinSep'], ['Con0.15'], ['Con1'], ['TIC']),
    'J_AJ_162_192_table2': ('mixed', ['Sep'], ['PA'], ['Con'], ['Obs'], None, ['MinSep'], ['Con0.15'], ['Con1'], ['TIC']),
    'J_AJ_165_180_table1': ('mixed', ['Sep'], ['PA'], ['DmagI'], ['Date'], None, ['rho'], ['dI0.15'], ['dI1'], ['MSC']),
    'J_AJ_168_28_table1': ('res', ['Sep'], ['PA'], ['Delm'], ['Epoch'], 'Filt', [], [], [], ['WDS']),
    'J_AJ_168_28_table2': ('unres', [], [], [], ['Epoch'], 'Filt', ['MinSep'], ['Delm0.15'], ['Delm1'], ['WDS']),
    'J_AJ_166_139_table2': ('res', ['Sep'], ['PA'], ['Delm'], ['Epoch.Y'], 'Filt', [], [], [], ['WDS']),
    'J_AJ_166_139_table3': ('unres', [], [], [], ['Epoch'], 'Filt', ['MinSep'], ['Delm0.15'], ['Delm1'], ['WDS']),
    'J_AJ_164_58_table2': ('res', ['Sep'], ['PA'], ['Dmag'], ['Epoch'], 'Filt', [], [], [], ['WDS']),
    'J_AJ_164_58_table3': ('unres', [], [], [], ['Epoch'], 'Filt', ['AR'], ['Dmag0.15'], ['Dmag1'], ['WDS']),
    'J_AJ_162_41_table2': ('res', ['Sep'], ['PA'], ['dmag'], ['Tobs'], 'Filt', [], [], [], ['WDS']),
    'J_AJ_162_41_table3': ('unres', [], [], [], ['Tobs'], 'Filt', ['Res'], ['dmag0.15'], ['dmag1'], ['WDS']),
    'J_AJ_160_7_table2': ('res', ['Sep'], ['PA'], ['Delm'], ['Epoch.Y'], 'Filt', [], [], [], ['WDS']),
    'J_AJ_160_7_table3': ('unres', [], [], [], ['Epoch.Y'], 'Filt', ['Sepmin'], ['Dm0.15'], ['Dm1'], ['WDS']),
    'J_AJ_158_48_table2': ('res', ['Sep'], ['PA'], ['Dm'], ['Epoch'], 'Filt', [], [], [], ['WDS']),
    'J_AJ_158_48_table3': ('unres', [], [], [], ['Epoch'], 'Filt', ['Sepmin'], ['Dm0.15'], ['Dm1'], ['WDS']),
    'J_AJ_155_235_table3': ('res', ['Sep'], ['PA'], ['Dmag'], ['Epoch'], 'Filt', [], [], [], ['WDS']),
}
for f in sorted(glob.glob(f'{O}/J_AJ_*.ecsv')):
    key = os.path.basename(f).replace('.ecsv', '')
    if key not in spec:
        continue
    kind, sc, pc, dc, ec, fc, ms, d15, d1, idc = spec[key]
    t = Table.read(f)
    ra = col(t, ['RAJ2000', '_RA']); de = col(t, ['DEJ2000', '_DE'])
    sep = col(t, sc) if sc else np.full(len(t), np.nan)
    pa = col(t, pc) if pc else np.full(len(t), np.nan)
    dm = col(t, dc) if dc else np.full(len(t), np.nan)
    ep = col(t, ec)
    ep = np.where(ep < 100, ep + 2000, ep)            # Ziegler 'Obs' given as 19.2063 etc.
    filt = scol(t, [fc]) if fc else np.array(['I'] * len(t))
    msep = col(t, ms) if ms else np.full(len(t), np.nan)
    m15 = col(t, d15) if d15 else np.full(len(t), np.nan)
    m1 = col(t, d1) if d1 else np.full(len(t), np.nan)
    ids = scol(t, idc)
    comp = scol(t, ['Comp', 'Disc', 'Design', 'DesLet', 'Dis'])
    for i in range(len(t)):
        resolved = np.isfinite(sep[i]) and sep[i] > 0
        rows.append(dict(cat=key, id=ids[i], comp=comp[i], ra=ra[i], dec=de[i], epoch=ep[i], filt=filt[i].strip(),
                         resolved=bool(resolved), sep=sep[i], pa=pa[i], dmag=dm[i], minsep=msep[i], dm015=m15[i], dm1=m1[i]))
df = pd.DataFrame(rows)
df = df[np.isfinite(df.ra) & np.isfinite(df.dec)].reset_index(drop=True)
df['row'] = np.arange(len(df))
df.to_csv(f'{O}/speckle_rows.csv', index=False)
print(len(df), 'rows;', df.resolved.sum(), 'resolved;', (~df.resolved).sum(), 'unresolved/limits')
# unique positions for xmatch
pos = df.groupby([df.ra.round(4), df.dec.round(4)]).row.first().reset_index()
pos.columns = ['ra', 'dec', 'row']
tab = Table.from_pandas(pos)
print(len(tab), 'unique positions')
res = []
for s in range(0, len(tab), 3000):
    r = XMatch.query(cat1=tab[s:s + 3000], cat2='vizier:I/355/gaiadr3', max_distance=6 * u.arcsec, colRA1='ra', colDec1='dec')
    res.append(r.to_pandas())
    print('xmatch chunk', s, len(r), flush=True)
x = pd.concat(res)
x.to_csv(f'{O}/xmatch_gaia.csv', index=False)
print(len(x), 'gaia matches')
