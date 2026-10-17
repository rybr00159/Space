"""h16: calibrate published XP-based metallicities on M dwarfs and read off the M3V's value.
 - Gaia Net (Huson+25, J/ApJ/984/58/table1): XP coefficients only (no photometry, no astrometry -> luminosity-free;
   M3V: [Fe/H] = -0.03, logTeff 3.521, logg 4.80);
 - SHBoost (Khalatyan+24, V/160/shboost): XP + photometry + parallax (NOT luminosity-free; M3V: [M/H] = -0.16,
   mass 0.50, logg 4.45);
 - StarHorse DR2 (Anders+19, I/349): photometry+parallax only (M3V: met50 = +0.29) - uses the over-luminosity.
Calibrators: Mann-scale spectroscopic [Fe/H] (data/calibrators_gaia.csv, h04b; 2.45 < BP-RP < 3.10, RUWE < 1.4,
plx/err > 20), matched by Gaia source_id via CDS XMatch on the Gaia DR3 positions (1").
Also the field stars at the M3V's colour (data/field_mdwarfs.csv): 'normal' (|dMKs| < 0.2) vs 'over-luminous'
(-1.0 < dMKs < -0.55), RUWE < 1.2 - empirical twin null for Gaia Net.
Output: out/h16_lit_xp_feh.json, figs/h16_lit_xp_feh.png"""
import json, numpy as np, pandas as pd, warnings, time
warnings.filterwarnings('ignore')
from astroquery.xmatch import XMatch
from astropy.table import Table
import astropy.units as u
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
B = 'tic214311908/analysis/host'
C = pd.read_csv(f'{B}/data/calibrators_gaia.csv'); C = C[(C.scale != 'LAMOST') & np.isfinite(C.feh)]
C['w'] = 1 / np.clip(C.efeh.fillna(0.12), 0.05, 0.5) ** 2
g = C.groupby('source_id').apply(lambda d: pd.Series(dict(feh=np.sum(d.w * d.feh) / np.sum(d.w), src='+'.join(sorted(set(d.src)))))).reset_index()
C = g.merge(C.drop_duplicates('source_id')[['source_id', 'gra', 'gdec', 'bp_rp', 'ruwe', 'parallax', 'parallax_error', 'ks_m', 'phot_g_mean_mag']], on='source_id')
C['MKx'] = C.ks_m + 5 * np.log10(C.parallax / 100)
_m = C.bp_rp.between(2.45, 3.10) & np.isfinite(C.MKx) & np.isfinite(C.bp_rp)
_rr = np.polyfit(C.bp_rp[_m], C.MKx[_m], 2)
C = C[C.bp_rp.between(2.45, 3.10) & (C.ruwe < 1.4) & (C.parallax / C.parallax_error > 20) & (np.abs(C.MKx - np.polyval(_rr, C.bp_rp)) < 1.5)]
F = pd.read_csv(f'{B}/data/field_mdwarfs.csv')
F['MK'] = F.ks_m + 5 * np.log10(F.parallax / 100)
s0 = F[(F.ruwe < 1.2) & (F.ph_qual == 'AAA') & (np.abs(F.bp_rp - 2.756) < 0.10)]
c1 = np.polyfit(s0.bp_rp, s0.MK, 1); r = s0.MK - np.polyval(c1, s0.bp_rp)
h, e = np.histogram(r, bins=np.arange(-1.5, 1.5, 0.05)); mode = 0.5 * (e[np.argmax(h)] + e[np.argmax(h) + 1])
F['dMK'] = F.MK - np.polyval(c1, F.bp_rp) - mode
Fs = F[(np.abs(F.bp_rp - 2.756) < 0.10) & (F.ruwe < 1.2) & (F.ph_qual == 'AAA') & ((F.dMK.abs() < 0.2) | F.dMK.between(-1.0, -0.55))].copy()
Fs['grp'] = np.where(Fs.dMK.abs() < 0.2, 'normal', 'overlum')
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
def xm(df, ra, de, cat):
    # multi-position VizieR cone queries (0.5", Gaia 2016 positions); match back by nearest position
    cat = cat.replace('vizier:', '')
    v = Vizier(columns=['**', '_q'], row_limit=-1, timeout=600)
    out = []
    for k in range(0, len(df), 300):
        sub = df.iloc[k:k + 300]
        c = SkyCoord(sub[ra].values * u.deg, sub[de].values * u.deg)
        for att in range(4):
            try:
                r = v.query_region(c, radius=0.5 * u.arcsec, catalog=cat); break
            except Exception as e:
                print('retry', str(e)[:100]); time.sleep(15); r = None
        if r is None or len(r) == 0: continue
        t = r[0].to_pandas(); t['source_id'] = sub.source_id.values[t['_q'].values - 1]
        out.append(t)
    o = pd.concat(out)
    o['angDist'] = 0.0
    return o.drop_duplicates('source_id')
out = {}
for lab, cat in (('gaianet', 'vizier:J/ApJ/984/58/table1'), ('shboost', 'vizier:V/160/shboost')):
    xc = xm(C, 'gra', 'gdec', cat); xf = xm(Fs.rename(columns={'ra': 'fra', 'dec': 'fdec'}), 'fra', 'fdec', cat)
    print(lab, 'columns', list(xc.columns)[:30])
    col = '[Fe/H]' if lab == 'gaianet' else 'Met'
    mc = C.merge(xc[['source_id', col]], on='source_id'); mc = mc[np.isfinite(mc[col])]
    d = mc[col] - mc.feh
    A = np.polyfit(mc.feh, mc[col], 1)
    loc = np.abs(mc.bp_rp - 2.756) < 0.10
    res = dict(n=int(len(mc)), median_offset=float(np.median(d)), robust_sd=float(1.4826 * np.median(np.abs(d - np.median(d)))),
               corr=float(np.corrcoef(mc.feh, mc[col])[0, 1]), slope_pred_vs_true=float(A[0]), intercept=float(A[1]),
               n_local=int(loc.sum()), local_median_offset=float(np.median(d[loc])) if loc.sum() else None,
               local_robust_sd=float(1.4826 * np.median(np.abs(d[loc] - np.median(d[loc])))) if loc.sum() > 3 else None,
               local_corr=float(np.corrcoef(mc.feh[loc], mc[col][loc])[0, 1]) if loc.sum() > 3 else None)
    # true [Fe/H] of calibrators whose catalogue value is close to the M3V's
    tv = -0.03 if lab == 'gaianet' else -0.161819
    nb = loc & (np.abs(mc[col] - tv) < 0.07)
    res['target_value'] = tv
    res['true_feh_given_similar_value'] = dict(n=int(nb.sum()), median=float(np.median(mc.feh[nb])) if nb.sum() else None,
                                               p16=float(np.percentile(mc.feh[nb], 16)) if nb.sum() > 3 else None, p84=float(np.percentile(mc.feh[nb], 84)) if nb.sum() > 3 else None,
                                               frac_gt_0p25=float(np.mean(mc.feh[nb] > 0.25)) if nb.sum() else None)
    # inverse calibration: [Fe/H]_true = a + b * value (regression of true on catalogue value, local colour)
    if loc.sum() > 10:
        Bc = np.polyfit(mc[col][loc], mc.feh[loc], 1); rr = mc.feh[loc] - np.polyval(Bc, mc[col][loc])
        res['inverse_fit_local'] = dict(a=float(Bc[1]), b=float(Bc[0]), scatter=float(np.std(rr)), target_true_est=float(np.polyval(Bc, tv)))
    # metal-rich calibrators: what does the catalogue give for true [Fe/H] > +0.25 at the target colour?
    mr = loc & (mc.feh > 0.25)
    res['catalogue_value_for_true_gt_0p25'] = dict(n=int(mr.sum()), median=float(np.median(mc[col][mr])) if mr.sum() else None,
                                                   frac_le_target=float(np.mean(mc[col][mr] <= tv)) if mr.sum() else None)
    sol = loc & (mc.feh.abs() < 0.1)
    res['catalogue_value_for_true_solar'] = dict(n=int(sol.sum()), median=float(np.median(mc[col][sol])) if sol.sum() else None,
                                                 frac_le_target=float(np.mean(mc[col][sol] <= tv)) if sol.sum() else None)
    mf = Fs.merge(xf[['source_id', col]], on='source_id'); mf = mf[np.isfinite(mf[col])]
    for gname in ('normal', 'overlum'):
        v = mf[col][mf.grp == gname]
        res[f'field_{gname}'] = dict(n=int(len(v)), median=float(np.median(v)), p16=float(np.percentile(v, 16)), p84=float(np.percentile(v, 84)),
                                     frac_le_target=float(np.mean(v <= tv)))
    out[lab] = res
    print(lab, json.dumps(res, indent=0)[:2500])
    if lab == 'gaianet':
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
        ax[0].scatter(mc.feh, mc[col], s=4, alpha=0.4, label='calibrators'); ax[0].scatter(mc.feh[loc], mc[col][loc], s=6, c='orange', label='|dBP-RP|<0.1')
        ax[0].plot([-1, 0.6], [-1, 0.6], 'k-'); ax[0].axhline(tv, c='r', ls='--', label='M3V'); ax[0].set_xlabel('[Fe/H] spectroscopic'); ax[0].set_ylabel('Gaia Net [Fe/H] (XP only)'); ax[0].legend()
        for gname, cl in (('normal', 'k'), ('overlum', 'm')):
            ax[1].hist(mf[col][mf.grp == gname], bins=np.arange(-0.8, 0.6, 0.04), histtype='step', density=True, color=cl, label=f'field {gname}')
        ax[1].axvline(tv, c='r', ls='--'); ax[1].set_xlabel('Gaia Net [Fe/H]'); ax[1].legend()
        plt.tight_layout(); plt.savefig(f'{B}/figs/h16_lit_xp_feh.png', dpi=110)
json.dump(out, open(f'{B}/out/h16_lit_xp_feh.json', 'w'), indent=1)
