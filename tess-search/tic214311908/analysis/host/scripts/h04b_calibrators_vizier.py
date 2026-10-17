"""h04b: builds data/calibrators_gaia.csv without the Gaia archive (whose upload/IN-list jobs failed with HTTP 500/503
on 2026-10-17): Gaia DR3 columns from the cached CDS XMatch (data/calibrators_xmatch.csv, h04; VizieR I/355/gaiadr3),
match accepted when the Gaia position propagated to J2000 is < 2" from the literature position; then CDS XMatch of the
Gaia positions (propagated to 1999.5 / 2010.5) with 2MASS PSC (II/246/out) and AllWISE (II/328/allwise), radius 2",
nearest. Columns renamed to the Gaia-archive names used by the other scripts. The large XMatch cache is then
reduced to the matched rows only (disk).
Output: data/calibrators_gaia.csv"""
import numpy as np, pandas as pd, warnings, time
warnings.filterwarnings('ignore')
from astroquery.xmatch import XMatch
from astropy.table import Table
import astropy.units as u
B = 'tic214311908/analysis/host'
C = pd.read_csv(f'{B}/data/calibrators_raw.csv')
X = pd.read_csv(f'{B}/data/calibrators_xmatch.csv')
if 'source_id' not in X:   # full cache from h04 -> reduce
  X = X[np.isfinite(X.pmRA) & (X.Gmag < 18)]
  cosd = np.cos(np.radians(X.dec))
  X['ra00'] = X.RAdeg - X.pmRA * 16 / 3.6e6 / cosd; X['de00'] = X.DEdeg - X.pmDE * 16 / 3.6e6
  X['d00'] = np.hypot((X.ra00 - X.ra) * cosd, X.de00 - X.dec) * 3600
  X = X[X.d00 < 2.0].sort_values(['row', 'd00']).drop_duplicates('row')
ren = {'Source': 'source_id', 'RAdeg': 'ra', 'DEdeg': 'dec', 'Plx': 'parallax', 'e_Plx': 'parallax_error', 'pmRA': 'pmra', 'pmDE': 'pmdec',
       'Gmag': 'phot_g_mean_mag', 'BPmag': 'phot_bp_mean_mag', 'RPmag': 'phot_rp_mean_mag', 'BP-RP': 'bp_rp', 'E(BP/RP)': 'phot_bp_rp_excess_factor',
       'RUWE': 'ruwe', 'IPDfmp': 'ipd_frac_multi_peak', 'IPDgofha': 'ipd_gof_harmonic_amplitude', 'RV': 'radial_velocity', 'e_RV': 'radial_velocity_error',
       'GRVSmag': 'grvs_mag', 'NBPblend': 'phot_bp_n_blended_transits', 'epsi': 'astrometric_excess_noise', 'Nper': 'visibility_periods_used'}
if 'source_id' not in X:
  keep = ['row', 'd00'] + list(ren.keys())
  X = X[keep].rename(columns=ren)
  X.to_csv(f'{B}/data/calibrators_xmatch.csv', index=False)   # shrink cache
X = X.rename(columns={'ra': 'gra', 'dec': 'gdec'})
M = C.merge(X, on='row', how='inner')
print('matched literature rows', len(M), M.src.value_counts().to_dict())
G = M.drop_duplicates('source_id')[['source_id', 'gra', 'gdec', 'pmra', 'pmdec']].copy()
cd = np.cos(np.radians(G.gdec))
for lab, ep, cat, rcol in (('tm', 1999.5, 'vizier:II/246/out', None), ('aw', 2010.5, 'vizier:II/328/allwise', None)):
    G[f'ra_{lab}'] = G.gra + G.pmra * (ep - 2016) / 3.6e6 / cd; G[f'de_{lab}'] = G.gdec + G.pmdec * (ep - 2016) / 3.6e6
    out = []
    for k in range(0, len(G), 20000):
        sub = G.iloc[k:k + 20000][['source_id', f'ra_{lab}', f'de_{lab}']]
        for att in range(4):
            try:
                out.append(XMatch.query(cat1=Table.from_pandas(sub), cat2=cat, max_distance=2 * u.arcsec, colRA1=f'ra_{lab}', colDec1=f'de_{lab}').to_pandas()); break
            except Exception as e:
                print('retry', lab, str(e)[:100]); time.sleep(15)
    o = pd.concat(out).sort_values(['source_id', 'angDist']).drop_duplicates('source_id')
    if lab == 'tm':
        o = o.rename(columns={'Jmag': 'j_m', 'e_Jmag': 'j_msigcom', 'Hmag': 'h_m', 'e_Hmag': 'h_msigcom', 'Kmag': 'ks_m', 'e_Kmag': 'ks_msigcom', 'Qfl': 'ph_qual'})
        o = o[['source_id', 'j_m', 'j_msigcom', 'h_m', 'h_msigcom', 'ks_m', 'ks_msigcom', 'ph_qual']]
    else:
        o = o.rename(columns={'W1mag': 'w1mpro', 'e_W1mag': 'w1mpro_error', 'W2mag': 'w2mpro', 'e_W2mag': 'w2mpro_error', 'W3mag': 'w3mpro', 'ccf': 'wise_cc', 'ex': 'wise_ext', 'qph': 'wise_ph_qual'})
        o = o[['source_id', 'w1mpro', 'w1mpro_error', 'w2mpro', 'w2mpro_error', 'w3mpro', 'wise_cc', 'wise_ext', 'wise_ph_qual']]
    M = M.merge(o, on='source_id', how='left')
    print(lab, 'matched', o.source_id.nunique())
M['has_xp_continuous'] = M.phot_g_mean_mag < 17.65
M.to_csv(f'{B}/data/calibrators_gaia.csv', index=False)
print('saved', len(M))
