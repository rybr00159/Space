"""f02b: magnitudes of the six stars that matter for ground follow-up: the M3V, the Mira, pair A, pair B, the 3.1" star
(Gaia DR3 5925373272827260544) and the 5.5" star (5925373272827262208), plus the 12.0" G=17.0 star.
Sources (all read-only):
* Gaia DR3 G/BP/RP and XP synthetic photometry (gaiadr3.synthetic_photometry_gspc; data/gaia_25arcsec_gspc.csv from f02a);
  Riello+2021 (EDR3 doc. 5.5.1 Table 5.9) G->V, Ic polynomials where no XP photometry exists;
* SkyMapper DR4, 2MASS, VHS DR6 rows already downloaded by ../star (star/data/vizier_*.csv, 20" cone), DECaPS2 object
  catalogue downloaded by ../star (star/data/decaps2_90.csv);
* positions propagated with Gaia proper motions to each catalogue's epoch before matching (match radius 1.0").
Output: out/mags_table.csv (star x band, every source), out/mags_adopted.json (adopted AB/Vega magnitude per band)"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flib import *
import numpy as np, pandas as pd
SD = f'{B}/star/data'
g = pd.read_csv(f'{F}/data/gaia_25arcsec_gspc.csv', dtype={'source_id': np.int64})
STARS = {'M3V': 5925373272850784768, 'pairA': 5925373272850782208, 'pairB': 5925373272827253632, 'Mira': 5925373268528684288,
         'star3.1': 5925373272827260544, 'star5.5': 5925373272827262208, 'star12.0_G17.0': 5925373272827256576}
gs = g.set_index('source_id')
# pair B has a 2-parameter solution: give it pair A's proper motion (bound, 0.83")
def pos_at(name, jyear):
    r = gs.loc[STARS[name]]
    pmra, pmdec = (r.pmra, r.pmdec)
    if not np.isfinite(pmra):
        ra_ = gs.loc[STARS['pairA']]
        pmra, pmdec = (ra_.pmra, ra_.pmdec) if name == 'pairB' else (0.0, 0.0)
    dt = jyear - 2016.0
    dec = r.dec + pmdec*dt/3.6e6
    ra = r.ra + pmra*dt/3.6e6/np.cos(np.radians(r.dec))
    return ra, dec
def sep(ra1, de1, ra2, de2):
    return 3600*np.hypot((ra1 - ra2)*np.cos(np.radians(de1)), de1 - de2)
mjd2jy = lambda mjd: 2000.0 + (mjd - 51544.5)/365.25
rows = []
def add(star, band, mag, err, source, note=''):
    if mag is None or not np.isfinite(mag): return
    rows.append(dict(star=star, band=band, mag=float(mag), err=(float(err) if err is not None and np.isfinite(err) else np.nan), source=source, note=note))
for name, sid in STARS.items():
    r = gs.loc[sid]
    add(name, 'G', r.phot_g_mean_mag, 0.003, 'Gaia DR3', 'Vega')
    add(name, 'BP', r.phot_bp_mean_mag, 0.01, 'Gaia DR3', 'Vega')
    add(name, 'RP', r.phot_rp_mean_mag, 0.01, 'Gaia DR3', 'Vega')
    for b, c in (('V', 'v_jkc_mag'), ('Rc', 'r_jkc_mag'), ('Ic', 'i_jkc_mag'), ('g', 'g_sdss_mag'), ('r', 'r_sdss_mag'), ('i', 'i_sdss_mag'), ('z', 'z_sdss_mag'), ('y', 'y_ps1_mag')):
        add(name, b, r[c], 0.02, 'Gaia DR3 XP synthetic (gspc)', 'JKC Vega' if b in ('V', 'Rc', 'Ic') else 'SDSS/PS1 AB')
    x = r.bp_rp
    if np.isfinite(x):
        add(name, 'V', r.phot_g_mean_mag - (-0.02704 + 0.01424*x - 0.2156*x**2 + 0.01426*x**3), 0.05, 'Gaia G,BP-RP -> V (Riello+21)', 'sd 0.046 of the fit')
        add(name, 'Ic', r.phot_g_mean_mag - (0.01753 + 0.76*x - 0.0991*x**2), 0.05, 'Gaia G,BP-RP -> Ic (Riello+21)', 'sd 0.037 of the fit')
# SkyMapper DR4
s = pd.read_csv(f'{SD}/vizier_skymapper4.csv')
for row in s.itertuples():
    jy = mjd2jy(row.EpMean)
    for name in STARS:
        ra, de = pos_at(name, jy)
        if sep(ra, de, row.RAICRS, row.DEICRS) < 1.0:
            for b in 'griz':
                add(name, b, getattr(row, f'{b}PSF'), getattr(row, f'e_{b}PSF'), 'SkyMapper DR4 PSF', f'AB; FWHM {row.FWHM}"; {b}mmvar {getattr(row, b + "mmvar")}')
# DECaPS2
d = pd.read_csv(f'{SD}/decaps2_90.csv')
for row in d.itertuples():
    jy = mjd2jy(row.epochmean)
    for name in STARS:
        ra, de = pos_at(name, jy)
        if sep(ra, de, row.ra, row.dec) < 1.0:
            for b in 'grizy':
                m = getattr(row, f'mean_mag_{b}')
                if np.isfinite(m) and m < 30:
                    add(name, b if b != 'y' else 'Y', m, getattr(row, f'stdev_{b}', np.nan), 'DECaPS2 mean', f'AB; nmag_ok_{b}={getattr(row, "nmag_ok_" + b)}; stdev_{b} is flux stdev')
# 2MASS
t = pd.read_csv(f'{SD}/vizier_2mass.csv')
for row in t.itertuples():
    jy = 2000.0 + (row.JD - 2451545.0)/365.25
    for name in STARS:
        ra, de = pos_at(name, jy)
        if sep(ra, de, row.RAJ2000, row.DEJ2000) < 1.5:
            for b, c in (('J', 'Jmag'), ('H', 'Hmag'), ('Ks', 'Kmag')):
                add(name, b, getattr(row, c), getattr(row, 'e_' + c), '2MASS PSC', f'Vega; Qflg {row.Qflg}, Cflg {row.Cflg}')
# VHS DR6
v = pd.read_csv(f'{SD}/vizier_vhs6.csv')
for row in v.itertuples():
    for name in STARS:
        jy = mjd2jy(row.JMJD) if np.isfinite(row.JMJD) else 2014.6
        ra, de = pos_at(name, jy)
        if sep(ra, de, row.RAJ2000, row.DEJ2000) < 1.0:
            add(name, 'J', row.Jpmag, row.e_Jpmag, 'VHS DR6 Jpmag', 'Vega; saturates near J~12 (M3V/Mira/pair values unreliable)')
            ks = getattr(row, 'Kspmag', np.nan)
            add(name, 'Ks', ks, getattr(row, 'e_Kspmag', np.nan), 'VHS DR6 Kspmag', 'Vega; saturates near Ks~11')
T = pd.DataFrame(rows)
T.to_csv(f'{OUT}/mags_table.csv', index=False)
# ---------- adopted values ----------
pref = {'g': ['Gaia DR3 XP synthetic (gspc)', 'DECaPS2 mean', 'SkyMapper DR4 PSF'], 'r': ['Gaia DR3 XP synthetic (gspc)', 'DECaPS2 mean', 'SkyMapper DR4 PSF'],
        'i': ['Gaia DR3 XP synthetic (gspc)', 'DECaPS2 mean', 'SkyMapper DR4 PSF'], 'z': ['Gaia DR3 XP synthetic (gspc)', 'DECaPS2 mean', 'SkyMapper DR4 PSF'],
        'V': ['Gaia DR3 XP synthetic (gspc)', 'Gaia G,BP-RP -> V (Riello+21)'], 'Ic': ['Gaia DR3 XP synthetic (gspc)', 'Gaia G,BP-RP -> Ic (Riello+21)'],
        'J': ['2MASS PSC', 'VHS DR6 Jpmag'], 'Ks': ['2MASS PSC', 'VHS DR6 Kspmag'], 'G': ['Gaia DR3'], 'BP': ['Gaia DR3'], 'RP': ['Gaia DR3']}
bad2mass = {('pairA', '2MASS PSC'), ('pairB', '2MASS PSC')}       # blended (Cflg c, Qflg F/D/U)
badvhs = {'M3V', 'Mira', 'pairA'}
ad = {}
for name in STARS:
    ad[name] = {}
    for b, srcs in pref.items():
        for src in srcs:
            if (name, src) in bad2mass: continue
            if src.startswith('VHS') and name in badvhs: continue
            sel = T[(T.star == name) & (T.band == b) & (T.source == src)]
            if len(sel):
                ad[name][b] = dict(mag=float(sel.mag.iloc[0]), source=src); break
# 3.1"/5.5" stars: no XP photometry; DECaPS2 griz if matched, else colour transform from G and BP-RP (flagged)
for name in ('star3.1', 'star5.5'):
    for b in 'griz':
        if b not in ad[name]:
            ad[name][b] = dict(mag=None, source='not matched')
jdump(dict(adopted=ad, star_ids={k: str(v) for k, v in STARS.items()},
           note='g r i z AB (SDSS-like); V Ic Vega (JKC); J Ks Vega (2MASS); the Mira values are single-epoch or mean values of a '
                'large-amplitude variable (see f03_mira for its range)'), f'{OUT}/mags_adopted.json')
pd.set_option('display.width', 250)
print(T.pivot_table(index=['star', 'source'], columns='band', values='mag', aggfunc='first').round(2).to_string())
for k, v in ad.items(): print(k, {b: (round(x['mag'], 2) if x['mag'] is not None else None) for b, x in v.items()})
