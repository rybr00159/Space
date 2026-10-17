"""h04: M-dwarf metallicity calibrators from the literature (spectroscopic [Fe/H], mostly on the Mann+13 scale):
  Mann+15 (J/ApJ/804/64), Mann+13 wide-binary M companions (J/AJ/145/52 table1, [Fe/H] of the FGK primary),
  Newton+14 (J/AJ/147/20 table6 + cpm table with primary [Fe/H]), Rojas-Ayala+12 (J/ApJ/748/93), Terrien+15
  (J/ApJS/220/16, mean of Mann+13 H- and K-band [Fe/H]), Gaidos+14 (J/MNRAS/443/2561), Birky+20 APOGEE/Cannon
  (J/ApJ/892/31 table6), Ding+22 LAMOST M dwarfs (J/ApJS/260/45, separate scale; only 3100<Teff<3650 K, SNR>30, random 30000).
Cross-match to Gaia DR3 (CDS XMatch, 30", Gaia position propagated back to J2000 with its PM, accept < 2").
Then Gaia archive: same columns as h03 (incl. 2MASS, AllWISE best neighbours) by source_id.
Output: data/calibrators_raw.csv (literature rows), data/calibrators_gaia.csv (merged)"""
import numpy as np, pandas as pd, warnings, time
warnings.filterwarnings('ignore')
from astroquery.vizier import Vizier
from astroquery.xmatch import XMatch
from astroquery.gaia import Gaia
from astropy.table import Table
import astropy.units as u
B = 'tic214311908/analysis/host'
v = Vizier(row_limit=-1, columns=['**', '_RAJ2000', '_DEJ2000'], timeout=600)
def g(cat):
    return v.get_catalogs(cat)[0].to_pandas()
L = []
def add(src, name, ra, dec, feh, efeh, scale='M13'):
    L.append(pd.DataFrame(dict(src=src, name=np.asarray(name).astype(str), ra=np.asarray(ra, float), dec=np.asarray(dec, float),
                               feh=np.asarray(feh, float), efeh=np.asarray(efeh, float), scale=scale)))
t = g('J/ApJ/804/64/stars'); add('Mann15', t.Name, t._RAJ2000, t._DEJ2000, t['[Fe/H]'], t['e_[Fe/H]'])
t = g('J/AJ/145/52/table1'); add('Mann13wb', t.Name, t._RAJ2000, t._DEJ2000, t['[Fe/H]'], t['e_[Fe/H]'])
t = g('J/AJ/147/20/table6'); add('Newton14', t.Name, t._RAJ2000, t._DEJ2000, t['[Fe/H]'], t['e_[Fe/H]'])
t = g('J/AJ/147/20/cpm'); add('Newton14cpm', t.Name2, t._RAJ2000, t._DEJ2000, t['[Fe/H]1'], t['e_[Fe/H]1'])
t = g('J/ApJ/748/93/stars'); add('RojasAyala12', t.Name, t['_RA'], t['_DE'], t['[Fe/H]'], t['e_[Fe/H]'])
t = g('J/ApJS/220/16/table4')
fe = np.nanmean(np.vstack([t.MFeHH.astype(float), t.MFeHK.astype(float)]), axis=0)
add('Terrien15', t['2MASS'], t._RAJ2000, t._DEJ2000, fe, np.full(len(t), 0.12))
t = g('J/MNRAS/443/2561/stars'); add('Gaidos14', t.PM, t._RAJ2000, t._DEJ2000, t['[Fe/H]'], t['e_[Fe/H]'])
t = g('J/ApJ/892/31/table6'); add('Birky20', t['2MASS'], t._RAJ2000, t._DEJ2000, t['[Fe/H]'], np.full(len(t), 0.12))
try:
    t = g('J/ApJS/260/45/tableb1'); print('Ding22 rows', len(t))
    t = t[(t.Teff > 3100) & (t.Teff < 3650) & (t.SNR > 30)]
    t = t.sample(n=min(len(t), 30000), random_state=1); print('Ding22 kept (3100<Teff<3650, SNR>30, random <=30000)', len(t))
    add('Ding22', t.LAMOST, t._RAJ2000, t._DEJ2000, t['[Fe/H]'], np.full(len(t), 0.15), scale='LAMOST')
except Exception as e:
    print('Ding22 failed', e)
C = pd.concat(L, ignore_index=True)
C = C[np.isfinite(C.feh) & np.isfinite(C.ra) & np.isfinite(C.dec)].reset_index(drop=True)
C['row'] = np.arange(len(C))
print(C.src.value_counts().to_dict())
C.to_csv(f'{B}/data/calibrators_raw.csv', index=False)
# XMatch in chunks (cached)
import os
XF = f'{B}/data/calibrators_xmatch.csv'
X = []
for k in (range(0, len(C), 20000) if not os.path.exists(XF) else []):
    sub = C.iloc[k:k + 20000]
    for attempt in range(3):
        try:
            x = XMatch.query(cat1=Table.from_pandas(sub[['row', 'ra', 'dec']]), cat2='vizier:I/355/gaiadr3', max_distance=30 * u.arcsec,
                             colRA1='ra', colDec1='dec').to_pandas()
            X.append(x); break
        except Exception as e:
            print('retry', e); time.sleep(10)
    print('xmatch', k, flush=True)
if os.path.exists(XF):
    X = pd.read_csv(XF)
else:
    X = pd.concat(X, ignore_index=True); X.to_csv(XF, index=False)
X = X[np.isfinite(X.pmRA) & (X.Gmag < 18)]
cosd = np.cos(np.radians(X.dec))
X['ra00'] = X.RAdeg - X.pmRA * 16 / 3.6e6 / cosd; X['de00'] = X.DEdeg - X.pmDE * 16 / 3.6e6
X['d00'] = np.hypot((X.ra00 - X.ra) * cosd, X.de00 - X.dec) * 3600
X = X[X.d00 < 2.0].sort_values(['row', 'd00']).drop_duplicates('row')
C = C.merge(X[['row', 'Source', 'd00']], on='row', how='inner').rename(columns={'Source': 'source_id'})
C['source_id'] = C.source_id.astype(np.int64)
print('matched', len(C), C.src.value_counts().to_dict())
ids = C.source_id.unique()
cols = """g.source_id, g.ra, g.dec, g.l, g.b, g.parallax, g.parallax_error, g.pmra, g.pmdec, g.phot_g_mean_mag, g.phot_bp_mean_mag,
 g.phot_rp_mean_mag, g.bp_rp, g.phot_bp_rp_excess_factor, g.phot_bp_n_blended_transits, g.phot_rp_n_blended_transits, g.phot_bp_n_obs,
 g.ruwe, g.ipd_frac_multi_peak, g.ipd_gof_harmonic_amplitude, g.astrometric_excess_noise, g.non_single_star, g.visibility_periods_used,
 g.radial_velocity, g.radial_velocity_error, g.rv_method_used, g.rv_nb_transits, g.rv_nb_deblended_transits, g.rv_expected_sig_to_noise,
 g.rv_template_teff, g.rv_chisq_pvalue, g.rv_renormalised_gof, g.rv_amplitude_robust, g.grvs_mag, g.has_xp_continuous,
 t.j_m, t.j_msigcom, t.h_m, t.h_msigcom, t.ks_m, t.ks_msigcom, t.ph_qual,
 w.w1mpro, w.w1mpro_error, w.w2mpro, w.w2mpro_error, w.w3mpro, w.cc_flags as wise_cc, w.ext_flag as wise_ext, w.ph_qual as wise_ph_qual"""
G = []
from astropy.table import Table as _T
for k in range(0, len(ids), 5000):
    up = _T({'sid': np.array(ids[k:k + 5000], dtype=np.int64)})
    up.write('tic214311908/analysis/host/data/_upload_ids.xml', format='votable', overwrite=True)
    q = f"""select {cols} from tap_upload.ids as u join gaiadr3.gaia_source as g on g.source_id = u.sid
left outer join gaiadr3.tmass_psc_xsc_best_neighbour as xt on xt.source_id = g.source_id
left outer join gaiadr1.tmass_original_valid as t on t.designation = xt.original_ext_source_id
left outer join gaiadr3.allwise_best_neighbour as xw on xw.source_id = g.source_id
left outer join gaiadr1.allwise_original_valid as w on w.designation = xw.original_ext_source_id"""
    for attempt in range(4):
        try:
            G.append(Gaia.launch_job_async(q, upload_resource='tic214311908/analysis/host/data/_upload_ids.xml', upload_table_name='ids').get_results().to_pandas()); break
        except Exception as e:
            print('retry', e); time.sleep(10)
    print('gaia', k, flush=True)
G = pd.concat(G, ignore_index=True).drop_duplicates('source_id')
M = C.merge(G, on='source_id', how='left')
M.to_csv(f'{B}/data/calibrators_gaia.csv', index=False)
print('saved', len(M))
