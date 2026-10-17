"""M-dwarf high-resolution imaging calibration sample for Gaia DR3 image/astrometric diagnostics.
Sources:
  Lamman et al. 2020 (AJ 159, 139; Robo-AO M-dwarf survey): J/AJ/159/139 table5 (LSPM -> Gaia DR2, observed flag),
     table3 (multiple systems) and table4 (doubles): Sep, Dmag (i' / LP600).
  Vrijmoet et al. 2022 (AJ 163, 178; SOAR speckle of 25-pc M dwarfs): table1 targets, table2 results (resolved or
     not, Sep, Dmag in I, detection limits).
Gaia DR3 counterparts: Lamman via dr2_neighbourhood (DR2 id -> DR3); Vrijmoet via CDS XMatch (15" around the
J2000 position, nearest source with parallax within 30% of the catalogue parallax and G < 17).
For every primary: RUWE, ipd_frac_multi_peak, ipd_gof_harmonic_amplitude, BP-blend fraction, G, BP-RP, parallax,
and the number of other Gaia DR3 sources within 3" (n_gaia3).
Output: data/mdw_imaging_sample.csv (one row per primary; columns comp_sep, comp_dm for the brightest companion at
0.08-2.5", resolved flag; NaN for non-detections)"""
import time, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from astroquery.vizier import Vizier
from astroquery.gaia import Gaia
from astroquery.xmatch import XMatch
from astropy.table import Table
import astropy.units as u
Gaia.ROW_LIMIT = -1
B = 'tic214311908/analysis/star'
v = Vizier(row_limit=-1, columns=['**', '_RAJ2000', '_DEJ2000'], timeout=300)
L5 = v.get_catalogs('J/AJ/159/139/table5')[0].to_pandas()
L3 = v.get_catalogs('J/AJ/159/139/table3')[0].to_pandas(); L4 = v.get_catalogs('J/AJ/159/139/table4')[0].to_pandas()
print('Lamman table5', len(L5), L5.Obs.value_counts().to_dict() if 'Obs' in L5 else '')
comp = pd.concat([L3[['LSPM', 'Sep', 'Dmag']], L4[['LSPM', 'Sep', 'Dmag']]], ignore_index=True)
comp['LSPM'] = comp.LSPM.astype(str).str.strip()
L5['LSPM'] = L5.LSPM.astype(str).str.strip()
obs = L5[L5.Obs.astype(str).str.strip().isin(['1', 'Y', 'y', 'True', 'yes'])] if 'Obs' in L5 else L5
if len(obs) < 100: obs = L5
print('Lamman observed targets used:', len(obs))
# DR2 -> DR3
ids = obs.Gaia.dropna().astype('int64').unique()
rows = []
for i in range(0, len(ids), 800):
    chunk = ','.join(str(x) for x in ids[i:i + 800])
    q = (f"select n.dr2_source_id, n.dr3_source_id, n.angular_distance, n.magnitude_difference from gaiadr3.dr2_neighbourhood as n "
         f"where n.dr2_source_id in ({chunk})")
    for k in range(3):
        try:
            rows.append(Gaia.launch_job(q).get_results().to_pandas()); break
        except Exception as e:
            print('retry', e); time.sleep(10)
nb = pd.concat(rows).sort_values(['dr2_source_id', 'angular_distance'])
nb['absdm'] = nb.magnitude_difference.abs()
nb = nb.sort_values(['dr2_source_id', 'absdm']).drop_duplicates('dr2_source_id')
obs = obs.merge(nb[['dr2_source_id', 'dr3_source_id']], left_on=obs.Gaia.astype('Int64'), right_on='dr2_source_id', how='left')
lam = []
for _, r in obs.iterrows():
    c = comp[comp.LSPM == r.LSPM]
    c = c[(c.Sep >= 0.08) & (c.Sep <= 2.5)]
    if len(c):
        b = c.sort_values('Dmag').iloc[0]
        lam.append(dict(survey='Lamman20', name=r.LSPM, dr3=r.dr3_source_id, resolved=True, comp_sep=b.Sep, comp_dm=b.Dmag, dm_band="i'/LP600"))
    else:
        anyc = (comp.LSPM == r.LSPM).any()
        lam.append(dict(survey='Lamman20', name=r.LSPM, dr3=r.dr3_source_id, resolved=bool(anyc), comp_sep=np.nan, comp_dm=np.nan, dm_band="i'/LP600",
                        only_wide=bool(anyc)))
lam = pd.DataFrame(lam)
# Vrijmoet 2022
V1 = v.get_catalogs('J/AJ/163/178/table1')[0].to_pandas(); V2 = v.get_catalogs('J/AJ/163/178/table2')[0].to_pandas()
V1['WDS'] = V1.WDS.astype(str).str.strip(); V2['WDS'] = V2.WDS.astype(str).str.strip()
print('Vrijmoet targets', len(V1), 'results', len(V2), V2.Res.value_counts().to_dict())
tab = Table.from_pandas(V1[['recno', '_RAJ2000', '_DEJ2000', 'plx']].rename(columns={'_RAJ2000': 'ra', '_DEJ2000': 'dec'}).astype({'ra': float, 'dec': float}))
xm = XMatch.query(cat1=tab, cat2='vizier:I/355/gaiadr3', max_distance=20 * u.arcsec, colRA1='ra', colDec1='dec').to_pandas()
xm = xm[(np.abs(xm.Plx - xm.plx) < 0.3 * xm.plx) & (xm.Gmag < 17)].sort_values(['recno', 'Gmag']).drop_duplicates('recno')
V1 = V1.merge(xm[['recno', 'Source']], on='recno', how='left')
vr = []
for _, r in V1.iterrows():
    res = V2[V2.WDS == r.WDS]
    rr = res[(res.Sep >= 0.08) & (res.Sep <= 2.5) & np.isfinite(res.Dmag)]
    resolved_any = res.Res.astype(str).str.contains('R|Y|1', regex=True).any() if len(res) else False
    if len(rr):
        b = rr.sort_values('Dmag').iloc[0]
        vr.append(dict(survey='Vrijmoet22', name=r.Name, dr3=r.Source, resolved=True, comp_sep=b.Sep, comp_dm=b.Dmag, dm_band='I'))
    elif len(res):
        vr.append(dict(survey='Vrijmoet22', name=r.Name, dr3=r.Source, resolved=bool(np.isfinite(res.Sep).any()), comp_sep=np.nan, comp_dm=np.nan, dm_band='I',
                       only_wide=bool(np.isfinite(res.Sep).any()), lim015=float(np.nanmax(res['Dmag0.15'])) if 'Dmag0.15' in res else np.nan))
vr = pd.DataFrame(vr)
S = pd.concat([lam, vr], ignore_index=True)
S = S[S.dr3.notna()].copy(); S['dr3'] = S.dr3.astype('int64')
# Gaia DR3 diagnostics + neighbours within 3"
ids = S.dr3.unique(); G = []
for i in range(0, len(ids), 800):
    chunk = ','.join(str(x) for x in ids[i:i + 800])
    q = (f"select source_id, ra, dec, phot_g_mean_mag, bp_rp, parallax, parallax_error, ruwe, ipd_frac_multi_peak, ipd_gof_harmonic_amplitude, "
         f"ipd_frac_odd_win, phot_bp_n_blended_transits, phot_bp_n_obs, phot_bp_rp_excess_factor, astrometric_params_solved, visibility_periods_used "
         f"from gaiadr3.gaia_source where source_id in ({chunk})")
    G.append(Gaia.launch_job(q).get_results().to_pandas())
G = pd.concat(G)
S = S.merge(G, left_on='dr3', right_on='source_id', how='inner')
tab = Table.from_pandas(S[['source_id', 'ra', 'dec']].drop_duplicates('source_id'))
xn = XMatch.query(cat1=tab, cat2='vizier:I/355/gaiadr3', max_distance=3.0 * u.arcsec, colRA1='ra', colDec1='dec').to_pandas()
xn = xn[xn.Source != xn.source_id]
cnt = xn.groupby('source_id').size()
S['n_gaia3'] = S.source_id.map(cnt).fillna(0).astype(int)
S.to_csv(f'{B}/data/mdw_imaging_sample.csv', index=False)
print(S.groupby(['survey', 'resolved']).size())
print('with close companion (0.08-2.5"):', int(np.isfinite(S.comp_sep).sum()))
