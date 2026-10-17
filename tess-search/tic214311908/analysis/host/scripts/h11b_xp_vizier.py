"""h11b: Gaia DR3 XP *sampled* mean spectra (official, 336-1020 nm, 2 nm; VizieR I/355/xpsample) for the sources in
data/xp/xp_ids.csv plus the M3V, the A(+B) entry and the Mira, fetched with multi-position VizieR cone queries
(0.5", Gaia 2016 positions) in chunks of 100. (The Gaia-archive DataLink route of h11 delivered ~50 spectra per 4 min
and was abandoned.) Resumable. Output: data/xp/xps_chunk_XXX.npz (source_id, flux, flux_error on the 343-point grid)"""
import os, numpy as np, pandas as pd, warnings, time
warnings.filterwarnings('ignore')
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import astropy.units as u
B = 'tic214311908/analysis/host'
ids = pd.read_csv(f'{B}/data/xp/xp_ids.csv')
C = pd.read_csv(f'{B}/data/calibrators_gaia.csv').drop_duplicates('source_id').set_index('source_id')[['gra', 'gdec']]
pos = C.rename(columns={'gra': 'ra', 'gdec': 'dec'})
if os.path.exists(f'{B}/data/field_mdwarfs.csv'):
    F = pd.read_csv(f'{B}/data/field_mdwarfs.csv').set_index('source_id')[['ra', 'dec']]
    pos = pd.concat([pos, F[~F.index.isin(pos.index)]])
extra = pd.DataFrame(dict(ra=[260.969674, 260.974873, 260.9776003047465], dec=[-51.323821, -51.324055, -51.3237201328113]),
                     index=[5925373272850784768, 5925373272850782208, 5925373268528684288])
pos = pd.concat([extra, pos[~pos.index.isin(extra.index)]])
want = list(extra.index) + [s for s in ids.source_id.values if s in pos.index]
wl = np.arange(336, 1021, 2.0)
v = Vizier(columns=['Source', 'lambda', 'Flux', 'e_Flux'], row_limit=-1, timeout=600)
CH = 100
for k in range(0, len(want), CH):
    fn = f'{B}/data/xp/xps_chunk_{k // CH:03d}.npz'
    sub = want[k:k + CH]
    if os.path.exists(fn):
        z = np.load(fn)
        if set(z['requested']) == set(sub): continue
    c = SkyCoord(pos.loc[sub, 'ra'].values * u.deg, pos.loc[sub, 'dec'].values * u.deg)
    for att in range(5):
        try:
            r = v.query_region(c, radius=0.5 * u.arcsec, catalog='I/355/xpsample'); break
        except Exception as e:
            print('retry', k, str(e)[:100], flush=True); time.sleep(20); r = None
    if r is None or len(r) == 0:
        print('failed', k, flush=True); continue
    t = r[0].to_pandas(); t = t[t.Source.isin(sub)]
    sids, fl, fe = [], [], []
    for s_, g in t.groupby('Source'):
        g = g.sort_values('lambda')
        if len(g) != len(wl): continue
        sids.append(s_); fl.append(g.Flux.values); fe.append(g.e_Flux.values)
    np.savez_compressed(fn, source_id=np.array(sids, dtype=np.int64), flux=np.array(fl, dtype=np.float64), flux_error=np.array(fe, dtype=np.float64),
                        requested=np.array(sub, dtype=np.int64))
    print('chunk', k // CH, len(sids), 'of', len(sub), flush=True)
print('done')
