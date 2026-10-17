"""h12d: is the M3V's XP spectrum better described by a two-star composite than by a single star, compared with
other stars? Templates: colour-binned (0.04 mag) medians of the XP shapes of stars with a normal CMD position
(|dMKs| < 0.25; field + calibrators, data/xp/xps_chunk_*.npz). Single model: template at a free colour. Composite:
primary colour c1 (2.45-2.80) + secondary colour c2 >= c1, each weighted by its absolute G flux from the field
M_G(BP-RP) ridge (so the companion's brightness follows its colour). chi2 uses the star-to-star scatter of
colour-matched normal stars as noise. Delta chi2 = chi2_single - chi2_composite for the M3V and, as calibration, for
normal field stars, over-luminous field stars and calibrators with [Fe/H] > +0.2 or < -0.2 near the M3V's colour
(templates rebuilt without the star being tested).
Output: out/h12d_xp_composite.json"""
import glob, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
B = 'tic214311908/analysis/host'
Z = [np.load(f) for f in sorted(glob.glob(f'{B}/data/xp/xps_chunk_*.npz'))]
sid = np.concatenate([z['source_id'] for z in Z]); FL = np.vstack([z['flux'] for z in Z if len(z['source_id'])])
_, ui = np.unique(sid, return_index=True); ui = np.sort(ui); sid, FL = sid[ui], FL[ui]
wl = np.arange(336, 1021, 2.0); win = (wl >= 420) & (wl <= 1000); wn = (wl >= 620) & (wl <= 1000)
TID = 5925373272850784768
C = pd.read_csv(f'{B}/data/calibrators_gaia.csv').drop_duplicates('source_id').set_index('source_id')
Cfe = pd.read_csv(f'{B}/data/calibrators_gaia.csv'); Cfe = Cfe[Cfe.scale != 'LAMOST'].groupby('source_id').feh.mean()
F = pd.read_csv(f'{B}/data/field_mdwarfs.csv').set_index('source_id')
ids = pd.read_csv(f'{B}/data/xp/xp_ids.csv').set_index('source_id')
def get(s, col_c, col_f):
    if s in F.index: return F.loc[s, col_f]
    if s in C.index: return C.loc[s, col_c]
    return np.nan
bprp = np.array([get(s, 'bp_rp', 'bp_rp') for s in sid]); plx = np.array([get(s, 'parallax', 'parallax') for s in sid])
ks = np.array([get(s, 'ks_m', 'ks_m') for s in sid]); G = np.array([get(s, 'phot_g_mean_mag', 'phot_g_mean_mag') for s in sid])
MK = ks + 5 * np.log10(plx / 100); MG = G + 5 * np.log10(plx / 100)
ok = np.isfinite(bprp) & np.isfinite(MK) & (sid != TID)
rK = np.polyfit(bprp[ok & (bprp > 2.3) & (bprp < 3.3)], MK[ok & (bprp > 2.3) & (bprp < 3.3)], 3)
dMK = MK - np.polyval(rK, bprp)
Ff = F[(F.ruwe < 1.2) & (F.parallax / F.parallax_error > 20)]
rMG = np.polyfit(Ff.bp_rp, Ff.phot_g_mean_mag + 5 * np.log10(Ff.parallax / 100), 3)
LS = np.log10(np.clip(FL, 1e-30, None)); LS = LS - np.log10(FL[:, wn].mean(1))[:, None]
it = np.where(sid == TID)[0][0]
nrm = ok & (np.abs(dMK) < 0.25)
bins = np.arange(2.30, 3.42, 0.04)
def build(exclude):
    tc, ts = [], []
    for b0 in bins:
        k = nrm & (np.abs(bprp - b0) < 0.03) & (np.arange(len(sid)) != exclude)
        if k.sum() >= 6: tc.append(b0); ts.append(np.median(LS[k][:, win], 0))
    return np.array(tc), np.array(ts)
def tmpl(tc, ts, c):
    j = np.clip(np.searchsorted(tc, c) - 1, 0, len(tc) - 2); f = (c - tc[j]) / (tc[j + 1] - tc[j]); return ts[j] * (1 - f) + ts[j + 1] * f
cm = nrm & (np.abs(bprp - 2.756) < 0.03)
med = np.median(LS[cm][:, win], 0); sd = 1.4826 * np.median(np.abs(LS[cm][:, win] - med), 0)
def dchi(i):
    tc, ts = build(i); St = LS[i, win]
    c0 = bprp[i] if i != it else 2.756
    sing = min(np.sum(((St - tmpl(tc, ts, c)) / sd) ** 2) for c in np.linspace(max(2.32, c0 - 0.2), min(3.38, c0 + 0.2), 41))
    best = np.inf
    for c1 in np.linspace(max(2.32, c0 - 0.3), min(3.3, c0 + 0.05), 15):
        for c2 in np.linspace(c1, 3.38, 20):
            w1 = 10 ** (-0.4 * np.polyval(rMG, c1)); w2 = 10 ** (-0.4 * np.polyval(rMG, c2))
            f = w1 * 10 ** tmpl(tc, ts, c1) + w2 * 10 ** tmpl(tc, ts, c2)
            best = min(best, np.sum(((St - np.log10(f / f[wn[win]].mean())) / sd) ** 2))
    return sing, best
res = {}
s_, c_ = dchi(it); res['M3V'] = dict(chi2_single=float(s_), chi2_comp=float(c_), dchi2=float(s_ - c_))
print('M3V', res['M3V'])
rng = np.random.default_rng(1)
groups = {'field_normal': [i for i in np.where(ok)[0] if ids.set.get(sid[i], '') == 'field_normal' and abs(bprp[i] - 2.756) < 0.06],
          'field_overlum': [i for i in np.where(ok)[0] if ids.set.get(sid[i], '') == 'field_overlum' and abs(bprp[i] - 2.756) < 0.06 and -1.0 < dMK[i] < -0.55],
          'cal_feh_gt_0.2': [i for i in np.where(ok)[0] if sid[i] in Cfe.index and Cfe[sid[i]] > 0.2 and abs(bprp[i] - 2.756) < 0.12],
          'cal_feh_lt_-0.2': [i for i in np.where(ok)[0] if sid[i] in Cfe.index and Cfe[sid[i]] < -0.2 and abs(bprp[i] - 2.756) < 0.12]}
for g, lst in groups.items():
    lst = list(rng.choice(lst, min(60, len(lst)), replace=False)) if len(lst) else []
    d = np.array([np.subtract(*dchi(i)) for i in lst])
    res[g] = dict(n=int(len(d)), median_dchi2=float(np.median(d)) if len(d) else None, p84=float(np.percentile(d, 84)) if len(d) else None,
                  frac_ge_M3V=float(np.mean(d >= res['M3V']['dchi2'])) if len(d) else None)
    print(g, res[g], flush=True)
json.dump(res, open(f'{B}/out/h12d_xp_composite.json', 'w'), indent=1)
