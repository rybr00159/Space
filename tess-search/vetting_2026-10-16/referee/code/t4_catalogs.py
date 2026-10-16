"""TIC 55829034 (and 214311906): search variable-star catalogues near the target for
periodic variables matching P (or 2P, P/2)."""
import sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/five_referee/code')
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.vizier import Vizier
from astroquery.gaia import Gaia
from reflib import save_json, BASE

TARGETS = {55829034: (272.03163, -35.69140, 0.616146), 214311906: (260.97765, -51.32369, 5.69149)}
CATS = {'VSX': 'B/vsx/vsx', 'GaiaDR3_EB': 'I/358/veb', 'GaiaDR3_vclassre': 'I/358/vclassre', 'OGLE_BLG_ECL': 'J/AcA/66/405',
        'ATLAS_var': 'J/AJ/156/241', 'ASASSN_var': 'II/366', 'ZTF_per': 'J/ApJS/249/18', 'Chen_WISE': 'J/ApJS/237/28'}
OUT = {}
v = Vizier(columns=['**', '+_r'], row_limit=-1)
for tic, (ra, dec, P) in TARGETS.items():
    c = SkyCoord(ra, dec, unit='deg')
    print(f'=== TIC {tic}  (l,b)=({c.galactic.l.deg:.2f},{c.galactic.b.deg:.2f})  P={P}')
    OUT[tic] = {}
    for name, cat in CATS.items():
        try:
            r = v.query_region(c, radius=3 * u.arcmin, catalog=cat)
        except Exception as e:
            print(f'  {name}: query failed {e}')
            continue
        if len(r) == 0:
            print(f'  {name}: none within 3\'')
            OUT[tic][name] = []
            continue
        rows = []
        for tab in r:
            pcols = [cn for cn in tab.colnames if cn.lower() in ('period', 'per', 'p', 'pf', 'frequency', 'freq', 'period_ls', 'perlc')]
            for row in tab:
                d = {cn: (str(row[cn]) if not isinstance(row[cn], (int, float, np.floating, np.integer)) else float(row[cn])) for cn in tab.colnames[:25]}
                rows.append(d)
                per = None
                for pc in pcols:
                    try:
                        per = float(row[pc])
                        if pc.lower().startswith('freq'):
                            per = 1 / per
                        break
                    except Exception:
                        pass
                match = ''
                if per:
                    for k, lab in [(1, 'P'), (2, '2P'), (0.5, 'P/2')]:
                        if abs(per / (k * P) - 1) < 0.01:
                            match = f'  <== matches {lab}'
                print(f"  {name}: r={float(row['_r'])*60 if float(row['_r'])<1 else float(row['_r']):.1f} {'arcsec' if float(row['_r'])<1 else 'arcmin?'} "
                      f"P={per} " + ' '.join(f'{cn}={row[cn]}' for cn in tab.colnames[1:7]) + match)
        OUT[tic][name] = rows
    # Gaia DR3 vari_eclipsing_binary / vari_summary within 3'
    q = (f"SELECT v.source_id, v.frequency, 1/v.frequency AS period, v.global_ranking, g.ra, g.dec, g.phot_g_mean_mag, "
         f"DISTANCE(POINT('ICRS', g.ra, g.dec), POINT('ICRS', {ra}, {dec}))*3600 AS sep "
         f"FROM gaiadr3.vari_eclipsing_binary AS v JOIN gaiadr3.gaia_source AS g ON v.source_id=g.source_id "
         f"WHERE 1=CONTAINS(POINT('ICRS', g.ra, g.dec), CIRCLE('ICRS', {ra}, {dec}, 0.05))")
    try:
        t = Gaia.launch_job(q).get_results()
        print(f'  Gaia DR3 vari_eclipsing_binary within 3\': {len(t)}')
        for row in t:
            per = float(row['period'])
            m = ''.join(f' <== matches {lab}' for k, lab in [(1, 'P'), (2, '2P'), (0.5, 'P/2')] if abs(per / (k * P) - 1) < 0.01)
            print(f"    {row['source_id']} sep={float(row['sep']):.1f}\" G={float(row['phot_g_mean_mag']):.2f} P={per:.6f}{m}")
        OUT[tic]['gaia_veb'] = [dict(source_id=int(r_['source_id']), sep=float(r_['sep']), P=float(r_['period']), G=float(r_['phot_g_mean_mag'])) for r_ in t]
    except Exception as e:
        print('  Gaia veb query failed', e)
save_json(OUT, f'{BASE}/results/t4_catalogs.json')
