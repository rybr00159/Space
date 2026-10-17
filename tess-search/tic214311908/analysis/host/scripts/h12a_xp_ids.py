"""h12a: choose the Gaia sources whose XP spectra are downloaded (h11):
 (a) metallicity calibrators (data/calibrators_gaia.csv; 2.40 < BP-RP < 3.15, has_xp_continuous, plx/err > 20,
     RUWE < 1.4, |C*| < 3): all Mann-scale non-Birky stars, Birky+20-only stars (random <= 1500), LAMOST Ding+22
     (random <= 700, separate scale);
 (b) field M dwarfs at the M3V colour (|BP-RP - 2.756| < 0.10, 12.0 < G < 15.0 [sampled XP exist only for G < ~15], RUWE < 1.2, has_xp, 2MASS AAA):
     'normal' CMD position (|dMKs| < 0.2; random 300) and 'over-luminous' (-1.0 < dMKs < -0.55; all, <= 300).
Usage: 'h12a_xp_ids.py cal' (stage 1, calibrators) then 'h12a_xp_ids.py' (appends field stars).
Output: data/xp/xp_ids.csv (source_id, set)"""
import numpy as np, pandas as pd
B = 'tic214311908/analysis/host'
rng = np.random.default_rng(12)
def Cs(d):
    x = d.bp_rp; f = 1.162004 + 0.011464 * x + 0.049255 * x ** 2 - 0.005879 * x ** 3
    return (d.phot_bp_rp_excess_factor - f) / (0.0059898 + 8.817481e-12 * d.phot_g_mean_mag ** 7.618399)
C = pd.read_csv(f'{B}/data/calibrators_gaia.csv')
C['Cs'] = Cs(C)
q = C.bp_rp.between(2.40, 3.15) & (C.has_xp_continuous.astype(str).str.lower() == 'true') & (C.parallax / C.parallax_error > 20) & (C.ruwe < 1.4) & (C.Cs.abs() < 3)
C = C[q]
srcs = C.groupby('source_id').src.apply(lambda s: set(s))
m13 = [s for s, v in srcs.items() if (v - {'Birky20', 'Ding22'})]
birky_only = [s for s, v in srcs.items() if v == {'Birky20'}]
lam_only = [s for s, v in srcs.items() if v == {'Ding22'}]
sel = [(s, 'cal_m13') for s in m13] + [(s, 'cal_birky') for s in rng.choice(birky_only, min(1500, len(birky_only)), replace=False)] + \
      [(s, 'cal_lamost') for s in rng.choice(lam_only, min(700, len(lam_only)), replace=False)]
import os, sys
if len(sys.argv) > 1 and sys.argv[1] == 'cal':   # stage 1: calibrators only
    out = pd.DataFrame(sel, columns=['source_id', 'set']).drop_duplicates('source_id')
    out.to_csv(f'{B}/data/xp/xp_ids.csv', index=False); print(out.set.value_counts().to_dict()); sys.exit()
prev = pd.read_csv(f'{B}/data/xp/xp_ids.csv') if os.path.exists(f'{B}/data/xp/xp_ids.csv') else None
F = pd.read_csv(f'{B}/data/field_mdwarfs.csv')
F['MK'] = F.ks_m + 5 * np.log10(F.parallax / 100)
s0 = F[(F.ruwe < 1.2) & (F.ph_qual == 'AAA') & (np.abs(F.bp_rp - 2.756) < 0.10)]
c = np.polyfit(s0.bp_rp, s0.MK, 1)
# mode of residual
r = s0.MK - np.polyval(c, s0.bp_rp); h, e = np.histogram(r, bins=np.arange(-1.5, 1.5, 0.05)); mode = 0.5 * (e[np.argmax(h)] + e[np.argmax(h) + 1])
F['dMK'] = F.MK - np.polyval(c, F.bp_rp) - mode
base = (np.abs(F.bp_rp - 2.756) < 0.10) & F.phot_g_mean_mag.between(12.0, 15.0) & (F.ruwe < 1.2) & (F.has_xp_continuous.astype(str).str.lower() == 'true') & (F.ph_qual == 'AAA')
nor = F[base & (F.dMK.abs() < 0.2)].source_id.values; ovl = F[base & F.dMK.between(-1.0, -0.55)].source_id.values
sel += [(s, 'field_normal') for s in rng.choice(nor, min(300, len(nor)), replace=False)] + [(s, 'field_overlum') for s in rng.choice(ovl, min(300, len(ovl)), replace=False)]
out = pd.DataFrame(sel, columns=['source_id', 'set']).drop_duplicates('source_id')
if prev is not None:   # keep the stage-1 order so that already-downloaded chunks stay valid; append the field stars
    out = pd.concat([prev, out[~out.source_id.isin(prev.source_id)]], ignore_index=True)
out.to_csv(f'{B}/data/xp/xp_ids.csv', index=False)
print(out.set.value_counts().to_dict(), 'field normal avail', len(nor), 'overlum avail', len(ovl))
