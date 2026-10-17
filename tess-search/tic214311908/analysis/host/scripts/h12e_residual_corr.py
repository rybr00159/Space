"""h12e: correlation of the M3V's XP shape residual (relative to the median of colour-matched normal field stars) with
the mean residual patterns of metal-rich ([Fe/H] > +0.25) and metal-poor ([Fe/H] < -0.25) calibrators within 0.05 mag
of its BP-RP, and with the mean pattern of over-luminous field stars; significance from the same correlation computed
for normal field stars (null distribution). Output: out/h12e_residual_corr.json"""
import glob, json, numpy as np, pandas as pd
B = 'tic214311908/analysis/host'
Z = [np.load(f) for f in sorted(glob.glob(f'{B}/data/xp/xps_chunk_*.npz'))]
sid = np.concatenate([z['source_id'] for z in Z]); FL = np.vstack([z['flux'] for z in Z if len(z['source_id'])])
_, ui = np.unique(sid, return_index=True); ui = np.sort(ui); sid, FL = sid[ui], FL[ui]
wl = np.arange(336, 1021, 2.0); win = (wl >= 420) & (wl <= 1000); wn = (wl >= 620) & (wl <= 1000)
S = np.log10(np.clip(FL[:, win], 1e-30, None) / FL[:, wn].mean(1)[:, None])
P = pd.read_csv(f'{B}/out/h12b_xp_pred_field.csv').set_index('source_id')
C = pd.read_csv(f'{B}/data/calibrators_gaia.csv'); C = C[C.scale != 'LAMOST']
fe = C.groupby('source_id').feh.mean(); col = C.drop_duplicates('source_id').set_index('source_id').bp_rp
T = S[sid == 5925373272850784768][0]
idx = {s: i for i, s in enumerate(sid)}
nor = [s for s in P.index if P.set[s] == 'field_normal' and abs(P.bp_rp[s] - 2.756) < 0.03 and s in idx]
ovl = [s for s in P.index if P.set[s] == 'field_overlum' and abs(P.bp_rp[s] - 2.756) < 0.05 and -1.0 < P.dMK_local[s] < -0.55 and s in idx]
med = np.median(S[[idx[s] for s in nor]], 0); sd = 1.4826 * np.median(np.abs(S[[idx[s] for s in nor]] - med), 0)
res = lambda v: (v - med) / sd
mr = [s for s in fe.index if fe[s] > 0.25 and s in idx and abs(col[s] - 2.756) < 0.05]
mp = [s for s in fe.index if fe[s] < -0.25 and s in idx and abs(col[s] - 2.756) < 0.05]
pat = {'metal_rich': np.mean([res(S[idx[s]]) for s in mr], 0), 'metal_poor': np.mean([res(S[idx[s]]) for s in mp], 0),
       'overluminous_field': np.mean([res(S[idx[s]]) for s in ovl], 0)}
pat['metal_rich_minus_poor'] = pat['metal_rich'] - pat['metal_poor']
out = dict(n_normal=len(nor), n_metal_rich=len(mr), n_metal_poor=len(mp), n_overlum=len(ovl))
rt = res(T)
for k, v in pat.items():
    c = float(np.corrcoef(rt, v)[0, 1])
    null = np.array([np.corrcoef(res(S[idx[s]]), v)[0, 1] for s in nor])
    out[k] = dict(corr_M3V=c, null_median=float(np.median(null)), null_sd=float(np.std(null)), pct_of_normal_stars_below=float(np.mean(null < c)))
    print(k, out[k])
json.dump(out, open(f'{B}/out/h12e_residual_corr.json', 'w'), indent=1)
