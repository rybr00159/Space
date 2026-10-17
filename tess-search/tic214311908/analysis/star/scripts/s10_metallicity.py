"""Photometric metallicity / over-luminosity calibration for TIC 214311908.
Calibration stars with spectroscopic [Fe/H]: Mann+15 (J/ApJ/804/64, 183 stars), Newton+14 (J/AJ/147/20 table6,
MEarth M dwarfs, NIR [Fe/H]), Rojas-Ayala+12 (J/ApJ/748/93). Gaia DR3 counterparts: CDS XMatch within 40" of the
J2000 position, choosing the source whose position propagated back to J2000 with its Gaia PM is within 1.5".
For calibrators with BP-RP 2.2-3.3, parallax_over_error > 10, RUWE < 1.4: offset from the s07 empirical ridge
(dMKs at fixed BP-RP; same ridge construction, recomputed here on data/mdwarf_sample.csv over BP-RP 2.45-3.05, with a
wider-range ridge for 2.2-3.3) vs [Fe/H]; linear fit with intrinsic scatter. Inverted for the target
(dMKs = -0.72) -> photometric [Fe/H] if single. Also: the wide companion TIC 214311910 A (G resolved, BP/RP blended)
under assumed colours of B.
Output: out/metallicity.json, figs/metallicity.png"""
import json, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from astroquery.vizier import Vizier
from astroquery.xmatch import XMatch
from astropy.table import Table
import astropy.units as u
B = 'tic214311908/analysis/star'
v = Vizier(row_limit=-1, columns=['**', '_RAJ2000', '_DEJ2000'], timeout=300)
cal = []
t = v.get_catalogs('J/ApJ/804/64/stars')[0].to_pandas(); cal.append(pd.DataFrame(dict(src='Mann15', name=t.Name, ra=t._RAJ2000, dec=t._DEJ2000, feh=t['[Fe/H]'], efeh=t['e_[Fe/H]'])))
t = v.get_catalogs('J/AJ/147/20/table6')[0].to_pandas(); cal.append(pd.DataFrame(dict(src='Newton14', name=t.Name, ra=t._RAJ2000, dec=t._DEJ2000, feh=t['[Fe/H]'], efeh=t['e_[Fe/H]'])))
t = v.get_catalogs('J/ApJ/748/93/stars')[0].to_pandas(); cal.append(pd.DataFrame(dict(src='RojasAyala12', name=t.Name, ra=t._RAJ2000 if '_RAJ2000' in t else t._RA, dec=t._DEJ2000 if '_DEJ2000' in t else t._DE, feh=t['[Fe/H]'], efeh=t['e_[Fe/H]'])))
C = pd.concat(cal, ignore_index=True); C = C[np.isfinite(C.feh) & np.isfinite(C.ra)].reset_index(drop=True); C['row'] = np.arange(len(C))
print(C.src.value_counts().to_dict())
X = XMatch.query(cat1=Table.from_pandas(C[['row', 'ra', 'dec']].astype({'ra': float, 'dec': float})), cat2='vizier:I/355/gaiadr3', max_distance=40 * u.arcsec, colRA1='ra', colDec1='dec').to_pandas()
X = X[np.isfinite(X.pmRA) & (X.Gmag < 17)]
cosd = np.cos(np.radians(X.dec))
X['ra00'] = X.RAdeg - X.pmRA * 16 / 3.6e6 / cosd; X['de00'] = X.DEdeg - X.pmDE * 16 / 3.6e6
X['d00'] = np.hypot((X.ra00 - X.ra) * cosd, X.de00 - X.dec) * 3600
X = X[X.d00 < 1.5].sort_values(['row', 'd00']).drop_duplicates('row')
C = C.merge(X[['row', 'Source', 'Plx', 'e_Plx', 'Gmag', 'BP-RP', 'RUWE']], on='row', how='inner')
# 2MASS Ks via Gaia archive best neighbour
from astroquery.gaia import Gaia
ids = ','.join(str(int(s)) for s in C.Source.unique())
k = Gaia.launch_job(f"select x.source_id, t.ks_m, t.ks_msigcom, t.ph_qual from gaiadr3.tmass_psc_xsc_best_neighbour as x join gaiadr1.tmass_original_valid as t on x.original_ext_source_id = t.designation where x.source_id in ({ids})").get_results().to_pandas()
C = C.merge(k, left_on='Source', right_on='source_id', how='left')
C['MK'] = C.ks_m + 5 * np.log10(C.Plx / 100)
C.to_csv(f'{B}/data/feh_calibrators.csv', index=False)
# ridge over a wider colour range from a d<50 pc M-dwarf sample: reuse s04 sample (2.45-3.05) and extend with a fit
m = pd.read_csv(f'{B}/data/mdwarf_sample.csv')
m['MK'] = m.ks_m + 5 * np.log10(m.parallax / 100)
good = (m.ruwe < 1.4) & (m.ipd_frac_multi_peak <= 5) & m.ph_qual.astype(str).str[2].eq('A') & (m.ks_msigcom < 0.05)
s = m[good]
x, y = s.bp_rp.values, s.MK.values
msk = np.isfinite(x) & np.isfinite(y)
x, y = x[msk], y[msk]; keep = np.ones(len(x), bool)
for it in range(10):
    c = np.polyfit(x[keep], y[keep], 2); r = y - np.polyval(c, x); sd = 1.4826 * np.median(np.abs(r[keep] - np.median(r[keep]))); keep = np.abs(r - np.median(r[keep])) < 2.5 * sd
h, e = np.histogram(r[np.abs(r) < 1.5], bins=np.arange(-1.5, 1.51, 0.05)); mode = 0.5 * (e[np.argmax(h)] + e[np.argmax(h) + 1])
ridge = lambda xx: np.polyval(c, xx) + mode
cal_sel = C[(C['BP-RP'] > 2.45) & (C['BP-RP'] < 3.05) & (C.Plx / C.e_Plx > 10) & (C.RUWE < 1.4) & C.ph_qual.astype(str).str[2].eq('A')].copy()
cal_sel['dMK'] = cal_sel.MK - ridge(cal_sel['BP-RP'])
cal_sel = cal_sel.drop_duplicates('Source')
print('calibrators in colour range:', len(cal_sel), cal_sel.src.value_counts().to_dict())
# fit [Fe/H] = a + b dMK, with intrinsic scatter; and the reverse dMK = c + d [Fe/H]
A = np.column_stack([np.ones(len(cal_sel)), cal_sel.dMK]); cf, *_ = np.linalg.lstsq(A, cal_sel.feh, rcond=None)
res = cal_sel.feh - A @ cf; sdr = np.std(res)
Ar = np.column_stack([np.ones(len(cal_sel)), cal_sel.feh]); cr, *_ = np.linalg.lstsq(Ar, cal_sel.dMK, rcond=None); sdm = np.std(cal_sel.dMK - Ar @ cr)
dMK_t = 10.482 + 5 * np.log10(15.5494 / 100) - ridge(2.756149)
feh_t = cf[0] + cf[1] * dMK_t
# bootstrap
bs = []
rng = np.random.default_rng(5)
for kk in range(2000):
    ii = rng.integers(0, len(cal_sel), len(cal_sel)); Ab = A[ii]; cb, *_ = np.linalg.lstsq(Ab, cal_sel.feh.values[ii], rcond=None); bs.append(cb[0] + cb[1] * dMK_t)
bs = np.array(bs)
print(f'[Fe/H] = {cf[0]:.3f} + {cf[1]:.3f} dMK  (scatter {sdr:.3f} dex);  dMK = {cr[0]:.3f} + {cr[1]:.3f} [Fe/H] (scatter {sdm:.3f} mag)')
print(f'target dMK {dMK_t:+.3f} -> [Fe/H]_phot = {feh_t:+.2f} +- {np.hypot(sdr, np.std(bs)):.2f} (if single)')
# fraction of calibrators as over-luminous and their [Fe/H]
ov = cal_sel[cal_sel.dMK < -0.5]
print('calibrators with dMK < -0.5:', len(ov), ' their [Fe/H] median', ov.feh.median(), ov[['name', 'src', 'feh', 'dMK', 'RUWE']].round(2).to_dict('records'))
# metallicity needed to explain dMK = -0.72 entirely: from reverse relation
feh_needed = (dMK_t - cr[0]) / cr[1]
# what dMK does the most metal-rich calibrators show
mr = cal_sel[cal_sel.feh > 0.3]
out = dict(n_cal=int(len(cal_sel)), fit_feh_vs_dMK=cf.tolist(), scatter_feh=sdr, fit_dMK_vs_feh=cr.tolist(), scatter_dMK=sdm, target_dMK=dMK_t,
           feh_phot_if_single=[float(feh_t), float(np.hypot(sdr, np.std(bs)))], feh_needed_from_reverse=float(feh_needed),
           metal_rich_cal_dMK=[float(mr.dMK.median()) if len(mr) else None, int(len(mr))],
           frac_cal_dMK_below_target=float(np.mean(cal_sel.dMK <= dMK_t)), overluminous_cal=ov[['name', 'src', 'feh', 'dMK', 'RUWE']].round(3).to_dict('records'))
# wide companion TIC 214311910 A: deblend BP/RP assuming B colours
pairA = {}
for bprpB, grpB in [(2.9, 1.18), (3.06, 1.25), (3.2, 1.30), (3.4, 1.36)]:
    GB = 15.740896; RPB = GB - grpB; BPB = RPB + bprpB
    RPA = -2.5 * np.log10(10 ** (-0.4 * 13.099223) - 10 ** (-0.4 * RPB)); BPA = -2.5 * np.log10(10 ** (-0.4 * 15.977035) - 10 ** (-0.4 * BPB))
    cA = BPA - RPA; MGA = 14.654037 + 5 * np.log10(15.5494 / 100)
    pairA[f'B_bprp_{bprpB}'] = dict(bprp_A=float(cA), MG_A=float(MGA))
out['pairA_deblended'] = pairA
print('pair A deblended colours:', pairA)
json.dump(out, open(f'{B}/out/metallicity.json', 'w'), indent=1, default=float)
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(5.5, 4))
for sname, mk in [('Mann15', 'o'), ('Newton14', 's'), ('RojasAyala12', '^')]:
    q = cal_sel[cal_sel.src == sname]; ax.plot(q.feh, q.dMK, mk, ms=4, alpha=0.7, label=sname)
xx = np.linspace(-0.8, 0.6, 10); ax.plot(xx, cr[0] + cr[1] * xx, 'k-'); ax.axhline(dMK_t, color='C3', ls='--', label='TIC 214311908')
ax.axhline(dMK_t + 0.753, color='C3', ls=':', label='target if equal-mass binary')
ax.invert_yaxis(); ax.set_xlabel('[Fe/H] (spectroscopic)'); ax.set_ylabel('dMKs from ridge at fixed BP-RP (mag)'); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(f'{B}/figs/metallicity.png', dpi=120)
