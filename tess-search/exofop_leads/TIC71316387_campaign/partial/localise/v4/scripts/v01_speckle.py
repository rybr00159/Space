"""v4: high-resolution-imaging calibration sample (own parser) + own CDS XMatch to Gaia DR3 (I/355/gaiadr3, 5").

Input: VizieR table copies downloaded by v3 (v3/speckle/*.ecsv; catalogue data only, no derived products reused):
  SOAR speckle (Tokovinin et al. 2018-2024): J/AJ/155/235, 158/48, 160/7, 162/41, 164/58, 166/139, 168/28
     resolved-pair tables (I and y filters; only I kept for contrasts) and non-resolution tables (Dm at 0.15" and 1")
  SOAR speckle of TESS TOIs (Ziegler et al. 2020, 2021): J/AJ/159/19 t6, J/AJ/162/192 t1, t2   (I band)
  Kepler KOI imaging compilation (Furlan et al. 2017): J/AJ/153/71 t8 (companions; i/880nm/692nm/LP600 contrasts preferred)
     + t2 (all KOIs; Nc = 0 -> no companion found)
  Robo-AO KOIs (Ziegler et al. 2017): J/AJ/153/66 t4, t6 (LP600 contrasts)
Output v4/speckle/rows.csv (one row per measurement or non-detection) and xm.csv (all Gaia DR3 sources within 5").
"""
import os, time
import numpy as np, pandas as pd
from astropy.table import Table
import astropy.units as u
from astroquery.xmatch import XMatch
import warnings; warnings.filterwarnings('ignore')

SRC = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v3/speckle'
O = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v4/speckle'


def T(name):
    return Table.read(f'{SRC}/{name}.ecsv')


def col(t, *names):
    for n in names:
        if n in t.colnames:
            x = t[n]
            if hasattr(x, 'filled'):
                try:
                    return np.array(x.filled(np.nan), float)
                except Exception:
                    x = x.filled('')
            out = []
            for v in x:
                try:
                    out.append(float(v))
                except Exception:
                    out.append(np.nan)
            return np.array(out)
    return np.full(len(t), np.nan)


def scol(t, *names):
    for n in names:
        if n in t.colnames:
            x = t[n]
            if hasattr(x, 'filled'):
                x = x.filled('')
            return np.array([str(v).strip() for v in x])
    return np.array([''] * len(t))


def pos(t):
    ra = col(t, 'RAJ2000', '_RAJ2000'); de = col(t, 'DEJ2000', '_DEJ2000')
    # some tables store sexagesimal RAJ2000 strings -> fall back to VizieR-computed _RAJ2000
    bad = ~np.isfinite(ra) | ~np.isfinite(de)
    if bad.any():
        ra2 = col(t, '_RAJ2000'); de2 = col(t, '_DEJ2000')
        ra[bad] = ra2[bad]; de[bad] = de2[bad]
    return ra, de


rows = []
soar_res = {'J_AJ_155_235_table3': 'Dmag', 'J_AJ_158_48_table2': 'Dm', 'J_AJ_160_7_table2': 'Delm', 'J_AJ_162_41_table2': 'dmag',
            'J_AJ_164_58_table2': 'Dmag', 'J_AJ_166_139_table2': 'Delm', 'J_AJ_168_28_table1': 'Delm'}
soar_unres = {'J_AJ_158_48_table3': ('Dm0.15', 'Dm1'), 'J_AJ_160_7_table3': ('Dm0.15', 'Dm1'), 'J_AJ_162_41_table3': ('dmag0.15', 'dmag1'),
              'J_AJ_164_58_table3': ('Dmag0.15', 'Dmag1'), 'J_AJ_166_139_table3': ('Delm0.15', 'Delm1'), 'J_AJ_168_28_table2': ('Delm0.15', 'Delm1')}
for name, dmc in soar_res.items():
    t = T(name); ra, de = pos(t)
    rows.append(pd.DataFrame(dict(cat=name, kind='res', ra=ra, dec=de, band=scol(t, 'Filt'), sep=col(t, 'Sep'), pa=col(t, 'PA'),
                                  dm=col(t, dmc), lim015=np.nan, lim1=np.nan, name=scol(t, 'WDS'))))
for name, (l1, l2) in soar_unres.items():
    t = T(name); ra, de = pos(t)
    rows.append(pd.DataFrame(dict(cat=name, kind='unres', ra=ra, dec=de, band=scol(t, 'Filt'), sep=np.nan, pa=np.nan, dm=np.nan,
                                  lim015=col(t, l1), lim1=col(t, l2), name=scol(t, 'WDS'))))
# Ziegler et al. SOAR TOI tables: rows with a separation are companions, others non-detections
for name, sc, pc, dc, l1, l2 in [('J_AJ_159_19_table6', 'rho', 'theta', 'Dmag', 'Dm0.15', 'Dm1.0'),
                                 ('J_AJ_162_192_table1', 'Sep', 'PA', 'Con', 'Con0.15', 'Con1'),
                                 ('J_AJ_162_192_table2', 'Sep', 'PA', 'Con', 'Con0.15', 'Con1')]:
    t = T(name); ra, de = pos(t)
    sep = col(t, sc)
    rows.append(pd.DataFrame(dict(cat=name, kind=np.where(np.isfinite(sep) & (sep > 0), 'res', 'unres'), ra=ra, dec=de,
                                  band=np.where(scol(t, 'Filt') == '', 'I', scol(t, 'Filt')), sep=sep, pa=col(t, pc), dm=col(t, dc),
                                  lim015=col(t, l1), lim1=col(t, l2), name='TIC' + scol(t, 'TIC'))))
# Furlan et al. 2017 KOIs
t2 = T('J_AJ_153_71_table2'); t8 = T('J_AJ_153_71_table8')
ra2, de2 = pos(t2)
koi2 = pd.DataFrame(dict(KOI=col(t2, 'KOI'), ra=ra2, dec=de2, Nc=col(t2, 'Nc')))
dm = np.full(len(t8), np.nan); band = np.array([''] * len(t8), dtype=object)
for c, lab in (('Dimag', 'i'), ('D880nm', '880nm'), ('D692nm', '692nm'), ('DLP600', 'LP600'), ('DF775W', 'F775W')):
    x = col(t8, c); m = np.isfinite(x) & ~np.isfinite(dm)
    dm[m] = x[m]; band[m] = lab
k8 = pd.DataFrame(dict(KOI=col(t8, 'KOI'), sep=col(t8, 'Sep'), pa=col(t8, 'PA'), dm=dm, band=band, dK=col(t8, 'DKmag')))
k8 = k8.merge(koi2[['KOI', 'ra', 'dec']], on='KOI', how='left')
rows.append(pd.DataFrame(dict(cat='J_AJ_153_71_table8', kind='res', ra=k8.ra, dec=k8.dec, band=k8.band, sep=k8.sep, pa=k8.pa, dm=k8.dm,
                              lim015=np.nan, lim1=np.nan, name='KOI' + k8.KOI.astype('Int64').astype(str))))
u0 = koi2[koi2.Nc == 0]
rows.append(pd.DataFrame(dict(cat='J_AJ_153_71_table2', kind='unres', ra=u0.ra, dec=u0.dec, band='', sep=np.nan, pa=np.nan, dm=np.nan,
                              lim015=np.nan, lim1=np.nan, name='KOI' + u0.KOI.astype('Int64').astype(str))))
for name in ('J_AJ_153_66_table4', 'J_AJ_153_66_table6'):
    t = T(name); ra, de = pos(t)
    rows.append(pd.DataFrame(dict(cat=name, kind='res', ra=ra, dec=de, band=scol(t, 'Filter'), sep=col(t, 'Sep'), pa=col(t, 'PA'),
                                  dm=col(t, 'Dmag'), lim015=np.nan, lim1=np.nan, name='KOI' + scol(t, 'KOI'))))
R = pd.concat(rows, ignore_index=True)
R = R[np.isfinite(R.ra) & np.isfinite(R.dec)].reset_index(drop=True)
R['ra'] = R.ra.round(6); R['dec'] = R.dec.round(6)
U = R[['ra', 'dec']].drop_duplicates().reset_index(drop=True); U['pos'] = np.arange(len(U))
R = R.merge(U, on=['ra', 'dec'], how='left')
R.to_csv(f'{O}/rows.csv', index=False)
print('rows', len(R), 'positions', len(U)); print(R.groupby(['cat', 'kind']).size().to_string())
outs = []
for i0 in range(0, len(U), 5000):
    sub = Table.from_pandas(U.iloc[i0:i0 + 5000])
    for k in range(4):
        try:
            r = XMatch.query(cat1=sub, cat2='vizier:I/355/gaiadr3', max_distance=5 * u.arcsec, colRA1='ra', colDec1='dec')
            break
        except Exception as e:
            print('retry', e); time.sleep(10)
    outs.append(r.to_pandas())
    print('xmatched', i0, len(r), flush=True)
X = pd.concat(outs, ignore_index=True)
X.to_csv(f'{O}/xm.csv', index=False)
print('xm rows', len(X), 'columns', list(X.columns))
