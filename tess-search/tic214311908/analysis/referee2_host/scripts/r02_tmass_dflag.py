"""r02: Is 2MASS photometry flagged 'd' (diffraction-spike contamination, cc_flg) biased bright at the 0.03-mag level that
the colour-residual test r(G-Ks) needs? The M3V has Cflg = 'ddd' (Mira, Ks 6.68, 17.8" E).
Method: clean field M dwarfs (host/data/field_mdwarfs.csv raw Gaia rows; RUWE<1.1, IPD clean, |C*|<3, 2MASS AAA,
BP-RP 2.45-3.10) are cross-matched (CDS XMatch, positions propagated to epoch 2000.0, 2") with the 2MASS PSC to get
Cflg and the 2MASS neighbour info. Colour residuals r(G-J), r(G-H), r(G-Ks) at fixed BP-RP (running-median ridge)
are compared between 'd'-flagged and clean ('000') stars, and for 'd' stars with a neighbour 3-5 mag brighter in Ks
within 10-30" (the M3V's situation) - matched on dMKs (|dMKs| < 0.3 normal stars; the d-flag is independent of
metallicity).
Output: out/r02_tmass_dflag.json"""
import json, numpy as np, pandas as pd, warnings, io
warnings.filterwarnings('ignore')
from astropy.table import Table
from astroquery.xmatch import XMatch
import astropy.units as u
H = 'tic214311908/analysis/host'
B = 'tic214311908/analysis/referee2_host'
F = pd.read_csv(f'{H}/data/field_mdwarfs.csv')
x = F.bp_rp; fx = 1.162004 + 0.011464 * x + 0.049255 * x ** 2 - 0.005879 * x ** 3
F['Cs'] = (F.phot_bp_rp_excess_factor - fx) / (0.0059898 + 8.817481e-12 * F.phot_g_mean_mag ** 7.618399)
cl = (F.ruwe < 1.1) & (F.ipd_frac_multi_peak == 0) & (F.parallax / F.parallax_error > 20) & (F.ph_qual == 'AAA') & (F.Cs.abs() < 3) & F.bp_rp.between(2.45, 3.10)
D = F[cl].copy().reset_index(drop=True)
dt = -16.0
D['ra0'] = D.ra + D.pmra * dt / 3.6e6 / np.cos(np.radians(D.dec)); D['de0'] = D.dec + D.pmdec * dt / 3.6e6
up = Table.from_pandas(D[['source_id', 'ra0', 'de0']])
import os
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
cache = f'{B}/data/r02_vizier_2mass_cflg.csv'
if os.path.exists(cache): X = pd.read_csv(cache)
else:
    V = Vizier(columns=['2MASS', 'RAJ2000', 'DEJ2000', 'Jmag', 'Hmag', 'Kmag', 'Qflg', 'Cflg', 'prox', 'pxPA'], row_limit=-1, timeout=600)
    res = []
    for k in range(0, len(D), 2000):
        sub = D.iloc[k:k + 2000]
        c = SkyCoord(sub.ra0.values, sub.de0.values, unit='deg')
        for att in range(3):
            try:
                r = V.query_region(c, radius=2 * u.arcsec, catalog='II/246/out')[0].to_pandas(); break
            except Exception as e:
                print('retry', k, str(e)[:80], flush=True)
        r['source_id'] = sub.source_id.values[r['_q'].values - 1]
        cr = SkyCoord(r.RAJ2000.values, r.DEJ2000.values, unit='deg'); cq = SkyCoord(sub.ra0.values[r['_q'].values - 1], sub.de0.values[r['_q'].values - 1], unit='deg')
        r['angDist'] = cr.separation(cq).arcsec
        res.append(r); print('vizier chunk', k, len(r), flush=True)
    X = pd.concat(res, ignore_index=True); X.to_csv(cache, index=False)
print(list(X.columns))
X = X.sort_values('angDist').drop_duplicates('source_id')
cols = [c for c in X.columns if c in ('source_id', 'angDist', 'Jmag', 'Hmag', 'Kmag', 'Qflg', 'Cflg', 'prox', 'Bflg', 'Rflg')]
X = X[cols]
D = D.merge(X, on='source_id', how='inner')
print('matched', len(D), D.Cflg.value_counts().head(10).to_dict())
D['MK'] = D.ks_m + 5 * np.log10(D.parallax / 100)
D['GK'] = D.phot_g_mean_mag - D.ks_m; D['GJ'] = D.phot_g_mean_mag - D.j_m; D['GH'] = D.phot_g_mean_mag - D.h_m
def ridge(xv, yv):
    g = np.arange(2.45, 3.101, 0.02); m = [np.median(yv[np.abs(xv - g0) < 0.03]) for g0 in g]; return lambda z: np.interp(z, g, m)
base = D.Cflg.astype(str).eq('000')
for q in ('MK', 'GK', 'GJ', 'GH'):
    f = ridge(D.bp_rp[base].values, D[q][base].values); D['r' + q] = D[q] - f(D.bp_rp)
out = dict(n_matched=int(len(D)), cflg_counts=D.Cflg.astype(str).value_counts().head(15).to_dict())
norm_ = D.rMK.abs() < 0.3
groups = {'000': D.Cflg.astype(str).eq('000'), 'K_d (3rd char d)': D.Cflg.astype(str).str[2].eq('d'), 'ddd': D.Cflg.astype(str).eq('ddd'),
          'K_c': D.Cflg.astype(str).str[2].eq('c'), 'K_h': D.Cflg.astype(str).str[2].eq('h')}
for g, sel in groups.items():
    k = sel & norm_
    if k.sum() < 5: out[g] = dict(n=int(k.sum())); continue
    e = {}
    for q in ('rGK', 'rGJ', 'rGH', 'rMK'):
        v = D[q][k].values; e[q] = dict(median=float(np.median(v)), mean_clipped=float(np.mean(v[np.abs(v - np.median(v)) < 0.15])),
                                          err_median=float(1.2533 * 1.4826 * np.median(np.abs(v - np.median(v))) / np.sqrt(len(v))), robust_sd=float(1.4826 * np.median(np.abs(v - np.median(v)))))
    out[g] = dict(n=int(k.sum()), **e)
    print(g, k.sum(), {q: (round(e[q]['median'], 4), round(e[q]['err_median'], 4)) for q in e})
# bootstrap difference d-flag minus 000 for rGK
rng = np.random.default_rng(3)
a = D.rGK[groups['K_d (3rd char d)'] & norm_].values; b = D.rGK[groups['000'] & norm_].values
if len(a) > 5:
    bs = [np.median(rng.choice(a, len(a))) - np.median(rng.choice(b, len(b))) for _ in range(2000)]
    out['dflag_minus_clean_rGK'] = dict(median=float(np.median(a) - np.median(b)), boot_16_84=[float(np.percentile(bs, 16)), float(np.percentile(bs, 84))],
                                        boot_95=float(np.percentile(bs, 95)), frac_a_gt_clean_p84=float(np.mean(a > np.percentile(b, 84))))
    print('d minus 000 rGK', out['dflag_minus_clean_rGK'])
D[['source_id', 'bp_rp', 'rMK', 'rGK', 'rGJ', 'rGH', 'Cflg', 'prox', 'ks_m']].to_csv(f'{B}/data/r02_field_2mass_cflg.csv', index=False)
json.dump(out, open(f'{B}/out/r02_tmass_dflag.json', 'w'), indent=1)
