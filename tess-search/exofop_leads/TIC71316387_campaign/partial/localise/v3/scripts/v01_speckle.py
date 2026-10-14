"""Build an empirical calibration sample of high-resolution-imaging targets (resolved pairs and non-detections) and
cross-match them to Gaia DR3 (CDS XMatch, VizieR I/355/gaiadr3, all sources within 6").

Sources (VizieR):
  SOAR speckle, Tokovinin et al.:  J/AJ/155/235 (2018), J/AJ/158/48 (2019), J/AJ/160/7 (2020), J/AJ/162/41 (2021),
                                   J/AJ/164/58 (2022), J/AJ/166/139 (2023), J/AJ/168/28 (2024)   [resolved + unresolved]
  SOAR speckle of TESS TOIs, Ziegler et al.: J/AJ/159/19 (2020), J/AJ/162/192 (2021)       [mixed rows]
  Kepler KOI imaging compilation, Furlan et al. 2017: J/AJ/153/71 table8 (companions), table2 (all KOIs, Nc)
  Robo-AO KOI survey, Ziegler et al. 2017: J/AJ/153/66 table4/6 (companions <2.5"), table9 (all observed)
Output: v3/speckle/rows.csv (normalised rows) and v3/speckle/xm.csv (Gaia sources within 6" of each row).
"""
import os, glob, time
import numpy as np
import pandas as pd
from astropy.table import Table
import astropy.units as u
from astroquery.vizier import Vizier
from astroquery.xmatch import XMatch
import warnings; warnings.filterwarnings('ignore')

O = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v3/speckle'
v = Vizier(row_limit=-1, columns=['**', '_RAJ2000', '_DEJ2000'], timeout=180)


def get(t):
    fn = f'{O}/{t.replace("/", "_")}.ecsv'
    if os.path.exists(fn):
        return Table.read(fn)
    for k in range(3):
        try:
            r = v.get_catalogs(t)
            tab = r[0]
            tab.write(fn, overwrite=True, format='ascii.ecsv')
            print('downloaded', t, len(tab), flush=True)
            return tab
        except Exception as e:
            print('retry', t, e); time.sleep(5)
    raise RuntimeError(t)


def c(t, names, default=np.nan):
    for n in names:
        if n in t.colnames:
            x = t[n]
            try:
                return np.array(x.filled(np.nan) if hasattr(x, 'filled') else x, dtype=float)
            except Exception:
                return np.array([float(s) if str(s).strip() not in ('', '--') else np.nan for s in x])
    return np.full(len(t), default)


def s(t, names):
    for n in names:
        if n in t.colnames:
            return np.array([str(x) for x in t[n]])
    return np.array([''] * len(t))


def radec(t):
    ra = c(t, ['RAJ2000', '_RAJ2000', '_RA']); de = c(t, ['DEJ2000', '_DEJ2000', '_DE'])
    return ra, de


rows = []
# --- Tokovinin SOAR: resolved tables and unresolved tables
tok = {'J/AJ/155/235/table3': 'res', 'J/AJ/158/48/table2': 'res', 'J/AJ/158/48/table3': 'unres',
       'J/AJ/160/7/table2': 'res', 'J/AJ/160/7/table3': 'unres', 'J/AJ/162/41/table2': 'res', 'J/AJ/162/41/table3': 'unres',
       'J/AJ/164/58/table2': 'res', 'J/AJ/164/58/table3': 'unres', 'J/AJ/166/139/table2': 'res', 'J/AJ/166/139/table3': 'unres',
       'J/AJ/168/28/table1': 'res', 'J/AJ/168/28/table2': 'unres'}
for name, kind in tok.items():
    t = get(name)
    ra, de = radec(t)
    filt = s(t, ['Filt'])
    ep = c(t, ['Epoch', 'Epoch.Y', 'Tobs'])
    df = pd.DataFrame(dict(src=name, kind=kind, ra=ra, dec=de, filt=filt, epoch=ep,
                           sep=c(t, ['Sep']) if kind == 'res' else np.nan, pa=c(t, ['PA']) if kind == 'res' else np.nan,
                           dm=c(t, ['Dmag', 'Dm', 'Delm', 'dmag']) if kind == 'res' else np.nan,
                           lim015=c(t, ['Dm0.15', 'dmag0.15', 'Dmag0.15', 'Delm0.15']), lim1=c(t, ['Dm1', 'dmag1', 'Dmag1', 'Delm1']),
                           ident=s(t, ['WDS']) + ' ' + s(t, ['Comp', 'DesLet', 'Disc', 'Dis', 'Design'])))
    rows.append(df)
# --- Ziegler SOAR TESS (rows with Sep = companion; rows without Sep = non-detection)
for name, sc, pc, dc, fc in [('J/AJ/159/19/table6', 'rho', 'theta', 'Dmag', 'Filt'),
                             ('J/AJ/162/192/table1', 'Sep', 'PA', 'Con', 'Filt'),
                             ('J/AJ/162/192/table2', 'Sep', 'PA', 'Con', None)]:
    t = get(name)
    ra, de = radec(t)
    sep = c(t, [sc]); pa = c(t, [pc]); dm = c(t, [dc])
    kind = np.where(np.isfinite(sep), 'res', 'unres')
    df = pd.DataFrame(dict(src=name, kind=kind, ra=ra, dec=de, filt=s(t, [fc]) if fc else 'I', epoch=np.nan, sep=sep, pa=pa, dm=dm,
                           lim015=c(t, ['Dm0.15', 'Con0.15']), lim1=c(t, ['Dm1.0', 'Con1']), ident='TIC ' + s(t, ['TIC'])))
    rows.append(df)
# --- Furlan et al. 2017 (KOIs): companions table8 joined to coordinates of table2
t2 = get('J/AJ/153/71/table2'); t8 = get('J/AJ/153/71/table8')
k2 = pd.DataFrame(dict(KOI=c(t2, ['KOI']), ra=radec(t2)[0], dec=radec(t2)[1], Nc=c(t2, ['Nc'])))
k8 = pd.DataFrame(dict(KOI=c(t8, ['KOI']), sep=c(t8, ['Sep']), pa=c(t8, ['PA'])))
# delta-mag: prefer i, 880nm, 692nm, LP600, F775W (closest to Gaia G/RP for red companions)
dmi = np.full(len(t8), np.nan); fsel = np.array([''] * len(t8), dtype=object)
for col, lab in [('Dimag', 'i'), ('D880nm', '880'), ('DF775W', 'F775W'), ('D692nm', '692'), ('DLP600', 'LP600')]:
    x = c(t8, [col]); m = np.isfinite(x) & ~np.isfinite(dmi)
    dmi[m] = x[m]; fsel[m] = lab
k8['dm'] = dmi; k8['filt'] = fsel; k8['dK'] = c(t8, ['DKmag'])
k8 = k8.merge(k2[['KOI', 'ra', 'dec']], on='KOI', how='left')
rows.append(pd.DataFrame(dict(src='J/AJ/153/71/table8', kind='res', ra=k8.ra, dec=k8.dec, filt=k8.filt, epoch=np.nan, sep=k8.sep,
                              pa=k8.pa, dm=k8.dm, lim015=np.nan, lim1=np.nan, ident='KOI ' + k8.KOI.astype(int).astype(str))))
un = k2[k2.Nc == 0]
rows.append(pd.DataFrame(dict(src='J/AJ/153/71/table2', kind='unres', ra=un.ra, dec=un.dec, filt='', epoch=np.nan, sep=np.nan,
                              pa=np.nan, dm=np.nan, lim015=np.nan, lim1=np.nan, ident='KOI ' + un.KOI.astype(int).astype(str))))
# --- Robo-AO KOIs (Ziegler et al. 2017): secure + likely detections within 2.5"
for name in ['J/AJ/153/66/table4', 'J/AJ/153/66/table6']:
    t = get(name)
    ra, de = radec(t)
    rows.append(pd.DataFrame(dict(src=name, kind='res', ra=ra, dec=de, filt=s(t, ['Filter']), epoch=np.nan, sep=c(t, ['Sep']),
                                  pa=c(t, ['PA']), dm=c(t, ['Dmag']), lim015=np.nan, lim1=np.nan, ident='KOI ' + s(t, ['KOI']))))
R = pd.concat(rows, ignore_index=True)
R = R[np.isfinite(R.ra) & np.isfinite(R.dec)].reset_index(drop=True)
R['row'] = np.arange(len(R))
R.to_csv(f'{O}/rows.csv', index=False)
print('rows', len(R), R.groupby(['src', 'kind']).size().to_string(), flush=True)
# unique positions -> XMatch
U = R.drop_duplicates(subset=['ra', 'dec'])[['ra', 'dec']].reset_index(drop=True)
U['pos'] = np.arange(len(U))
R = R.merge(U, on=['ra', 'dec'], how='left'); R.to_csv(f'{O}/rows.csv', index=False)
outs = []
for i0 in range(0, len(U), 4000):
    sub = Table.from_pandas(U.iloc[i0:i0 + 4000])
    for k in range(3):
        try:
            r = XMatch.query(cat1=sub, cat2='vizier:I/355/gaiadr3', max_distance=6 * u.arcsec, colRA1='ra', colDec1='dec')
            break
        except Exception as e:
            print('xmatch retry', e); time.sleep(10)
    keep = ['pos', 'angDist', 'Source', 'RAdeg', 'DEdeg', 'Plx', 'pmRA', 'pmDE', 'Gmag', 'BPmag', 'RPmag', 'BP-RP', 'RUWE', 'IPDfmp',
            'IPDgofha', 'IPDgofhp', 'IPDfow', 'epsi', 'sepsi', 'Solved', 'Nper', 'E(BP/RP)', 'o_BPmag', 'NBPblend', 'o_RPmag',
            'NRPblend', 'Dup', 'RV', 'RVamp', 'RVgof', 'NAL', 'NgAL']
    outs.append(r[[k for k in keep if k in r.colnames]].to_pandas())
    print('xmatched', i0, len(r), flush=True)
X = pd.concat(outs, ignore_index=True)
X.to_csv(f'{O}/xm.csv', index=False)
print('xm rows', len(X))
