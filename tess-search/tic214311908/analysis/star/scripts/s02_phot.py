"""Archival photometry around TIC 214311908 from VizieR (all sources within 20" of the M3V's 2016.0 position).
For each catalogue, the M3V's position is propagated to the catalogue's mean epoch with Gaia DR3 PM, and every
source's offset from that position is listed, so mismatches with the background star 3.1" NE (G 17.5) or with the
M-dwarf pair (11.7") are visible.
Output: data/vizier_<cat>.csv for each catalogue; data/phot_raw.json (summary of matched rows).
"""
import os, json, time, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord
D = 'tic214311908/analysis/star/data'
ra0, de0, pmra, pmde = 260.969674, -51.323821, -113.671152, 85.54355   # Gaia DR3, epoch 2016.0
c0 = SkyCoord(ra0, de0, unit='deg')
cats = {  # name: (vizier id, mean epoch for propagation)
    '2mass': ('II/246/out', 1999.5), 'allwise': ('II/328/allwise', 2010.5), 'catwise': ('II/365/catwise', 2015.4),
    'unwise': ('II/363/unwise', 2013.0), 'skymapper4': ('II/379/smssdr4', 2016.0), 'skymapper2': ('II/358/smss', 2016.0),
    'apass9': ('II/336/apass9', 2012.0), 'vhs5': ('II/367/vhs_dr5', 2013.0), 'vhs6': ('II/359/vhs_dr4', 2013.0),
    'denis': ('B/denis/denis', 1998.0), 'ucac4': ('I/322A/out', 2001.0), 'galex_guv': ('II/335/gufuv', 2007.0),
    'galex_ais': ('II/312/ais', 2007.0), 'gsc242': ('I/353/gsc242', 2016.0), 'nomad': ('I/297/out', 2000.0),
    'usnob1': ('I/284/out', 1980.0), 'decaps2': ('II/373/decaps2', 2018.0), 'vvvx': ('II/376/vvvx', 2016.0),
    'erass1': ('J/A+A/682/A34/erass1-m', 2020.0), 'rass2': ('J/A+A/588/A103/cat2rxs', 1990.8),
    'xmm4': ('IX/69/xmm4d13s', 2010.0), 'sxps2': ('IX/58/2sxps', 2012.0), 'tic82': ('IV/39/tic82', 2000.0),
    'tess_ctl': ('IV/39/tic82', 2000.0), 'akari': ('II/297/irc', 2006.5), 'gaia_dr3_vari': ('I/358/vclassre', 2016.0),
    'andrae23': ('J/ApJS/267/8/table1', 2016.0), 'ucd_mann': ('J/ApJ/804/64/table5', 2016.0),
    'ElBadry21': ('J/MNRAS/506/2269/table', 2016.0), 'wds': ('B/wds/wds', 2000.0), 'vsx': ('B/vsx/vsx', 2000.0),
    'asassn_v': ('II/366/catalog', 2016.0), 'atlas_var': ('J/AJ/156/241/table4', 2017.0), 'lepine_mdw': ('J/AJ/142/138/Mdwarfs', 2000.0),
    'cns5': ('J/A+A/670/A19/cns5', 2016.0), 'gaiaucd': ('J/A+A/669/A139/table1', 2016.0),
}
v = Vizier(columns=['**', '_RAJ2000', '_DEJ2000'], row_limit=-1, timeout=120)
summ = {}
for name, (cid, ep) in cats.items():
    fn = f'{D}/vizier_{name}.csv'
    if os.path.exists(fn):
        t = pd.read_csv(fn)
    else:
        t = None
        for k in range(2):
            try:
                r = v.query_region(c0, radius=20 * u.arcsec, catalog=cid)
                t = r[0].to_pandas() if len(r) else pd.DataFrame()
                break
            except Exception as e:
                print(name, 'retry', e); time.sleep(3)
        if t is None:
            summ[name] = 'query failed'; continue
        if len(t):
            dt = ep - 2016.0
            rap = ra0 + pmra * dt / 3.6e6 / np.cos(np.radians(de0)); dep = de0 + pmde * dt / 3.6e6
            ra = t['_RAJ2000'].astype(float); de = t['_DEJ2000'].astype(float)
            t['dE_arcsec'] = (ra - rap) * np.cos(np.radians(de0)) * 3600; t['dN_arcsec'] = (de - dep) * 3600
            t['sep_m3v_epoch'] = np.hypot(t.dE_arcsec, t.dN_arcsec)
            t = t.sort_values('sep_m3v_epoch')
        t.to_csv(fn, index=False)
    n = len(t)
    summ[name] = dict(cat=cid, n_within_20=n, nearest=(t.iloc[0].dropna().to_dict() if n else None))
    print(f'{name:12s} {cid:28s} n={n:3d}' + (f'  nearest sep {t.sep_m3v_epoch.iloc[0]:.2f}"' if n else ''), flush=True)
json.dump(summ, open(f'{D}/phot_raw.json', 'w'), indent=1, default=str)
