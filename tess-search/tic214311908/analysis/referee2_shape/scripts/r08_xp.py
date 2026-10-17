"""r08: Gaia DR3 XP blue-excess test of the 'cool white dwarf occulted by the M3V' scenario (missed test).
A WD that supplies 1.57% of the TESS-band light (the eclipse depth) at the Teff required by its radius (4500-7500 K
for 0.2-1.2 Msun) adds a much larger fraction of the light at 400-550 nm. The target's XP sampled spectrum is
compared with comparison M dwarfs of the same colour (../referee_star/data/mdw_cmd.csv: BP-RP 2.66-2.86,
G 13.0-15.5, RUWE < 1.3, ipd_frac_multi_peak <= 2), all with XP spectra from the Gaia archive (datalink XP_SAMPLED).
Statistic: blue index X = ln[F(420-540 nm)/F(740-860 nm)] regressed on BP-RP (and on the red index ln F(640-700)/F(740-860));
target residual vs the robust scatter of the comparison residuals. Prediction for each WD mass: the WD flux
(blackbody of radius R_WD, Teff solving the 1.57% TESS share against the observed XP spectrum) added to an M3V
spectrum => predicted change of X.
Output: out/xp_test.json, data/xp_*.csv (downloaded XP spectra; data only)"""
import sys, json, os, numpy as np, pandas as pd
from astroquery.gaia import Gaia
B = 'tic214311908/analysis'
RD = f'{B}/referee2_shape'
TID = 5925373272850784768
m = pd.read_csv(f'{B}/referee_star/data/mdw_cmd.csv')
c = m[(m.bp_rp > 2.66) & (m.bp_rp < 2.86) & (m.phot_g_mean_mag > 13.0) & (m.phot_g_mean_mag < 15.5) & (m.ruwe < 1.3) & (m.ipd_frac_multi_peak <= 2) & (m.source_id != TID)]
c = c.sample(min(len(c), 400), random_state=2)
ids = [TID] + c.source_id.astype(np.int64).tolist()
cache = f'{RD}/data/xp_sampled.csv'
if not os.path.exists(cache):
    rows = []
    for i in range(0, len(ids), 45):
        chunk = ids[i:i+45]
        res = Gaia.load_data(ids=chunk, data_release='Gaia DR3', retrieval_type='XP_SAMPLED', data_structure='DATAMODEL_GAIA', format='csv')
        for k, tabs in res.items():
            for tb in tabs:
                for row in tb:
                    fl = np.array(str(row['flux']).strip('()[]{}').split(','), float)
                    fe = np.array(str(row['flux_error']).strip('()[]{}').split(','), float)
                    rows.append(dict(source_id=int(row['source_id']), flux=json.dumps(fl.tolist()), flux_error=json.dumps(fe.tolist())))
        print('downloaded', len(rows), flush=True)
    pd.DataFrame(rows).to_csv(cache, index=False)
X = pd.read_csv(cache)
print(X.shape)
# ---------------- analysis
wl = np.arange(336, 1021, 2.0)            # XP sampled grid (nm)
F = {int(r.source_id): (np.array(json.loads(r.flux)), np.array(json.loads(r.flux_error))) for r in X.itertuples()}
assert TID in F, 'target XP missing'
def band(f, lo, hi):
    s = (wl >= lo) & (wl < hi); return np.sum(f[s])
def idx(f):
    red = band(f, 740, 860)
    return dict(blue=np.log(band(f, 420, 540)/red), green=np.log(band(f, 540, 640)/red), mid=np.log(band(f, 640, 700)/red), far=np.log(band(f, 900, 1000)/red))
cc = c.set_index('source_id')
rows = []
for sid, (f, e) in F.items():
    if sid == TID or sid not in cc.index: continue
    if np.any(~np.isfinite(f)): continue
    d = idx(f); d.update(source_id=sid, bp_rp=float(cc.loc[sid, 'bp_rp']), G=float(cc.loc[sid, 'phot_g_mean_mag']), ruwe=float(cc.loc[sid, 'ruwe'])); rows.append(d)
Cmp = pd.DataFrame(rows)
ft, et = F[TID]; T = idx(ft); T['bp_rp'] = 2.756149
print('n comparison', len(Cmp), 'target', T)
def rstd(x): return 1.4826*np.median(np.abs(x - np.median(x)))
out = dict(n_comp=int(len(Cmp)), target=T)
# regressions: blue ~ a + b*(bp_rp) + c*mid  (mid captures TiO depth / Teff at fixed colour)
for name, cols in (('bprp', ['bp_rp']), ('bprp_mid', ['bp_rp', 'mid']), ('bprp_green', ['bp_rp', 'green'])):
    A = np.column_stack([np.ones(len(Cmp))] + [Cmp[k].values for k in cols])
    for it in range(3):
        coef = np.linalg.lstsq(A, Cmp.blue.values, rcond=None)[0]
        res = Cmp.blue.values - A @ coef
        keep = np.abs(res - np.median(res)) < 4*rstd(res)
        A2 = A[keep]; coef = np.linalg.lstsq(A2, Cmp.blue.values[keep], rcond=None)[0]
    res = Cmp.blue.values - A @ coef
    tres = T['blue'] - (np.r_[1, [T[k] for k in cols]] @ coef)
    out[f'fit_{name}'] = dict(coef=coef.tolist(), robust_sd=float(rstd(res)), std=float(np.std(res)), target_resid=float(tres),
                             target_percentile=float(np.mean(res < tres)*100), z_robust=float(tres/rstd(res)))
    print(name, out[f'fit_{name}'])
# predicted WD effect on the blue index: WD blackbody with 1.57% of the TESS-band photons of the observed spectrum
filt = np.loadtxt(f'{B}/shape/data/filt_TESS.Red.dat')   # (Angstrom, transmission); data file read only
tw = np.interp(wl*10, filt[:, 0], filt[:, 1], left=0, right=0)
def planck(T):
    lam = wl*1e-9; return 1/(lam**5*(np.exp(1.4388e-2/(lam*T)) - 1))
fA = 0.0157
pred = {}
for TWD in [4500, 5000, 5500, 6000, 6500, 7000, 8000]:
    Bw = planck(TWD)
    # photon-weighted TESS share (XP flux is energy per nm; photons ~ f*lambda); the XP grid ends at 1020 nm (TESS to 1050)
    s = fA/(1 - fA)*np.sum(ft*wl*tw)/np.sum(Bw*wl*tw)
    fw = s*Bw
    Ti = idx(ft); Tw = idx(ft + fw)
    pred[str(TWD)] = dict(d_blue=float(Tw['blue'] - Ti['blue']), d_mid=float(Tw['mid'] - Ti['mid']), share_420_540=float(band(fw, 420, 540)/band(ft + fw, 420, 540)),
                          share_G_approx=float(np.sum(fw*wl)/np.sum((ft + fw)*wl)))
    print(TWD, pred[str(TWD)])
out['wd_prediction'] = pred
out['note'] = ('d_blue is the shift of the target blue index the WD would cause (target treated as M3V-only + WD; '
               'the WD Teff needed for M_WD = 0.3/0.45/0.6/0.8 is ~5100/5400-5700/5900-6300/6300-6800 K, ../shape s06b/s16)')
json.dump(out, open(f'{RD}/out/xp_test.json', 'w'), indent=1)
# ---------------- consistent prediction of the locus residual (the WD also changes the regressors)
def synth_bprp(f):
    # crude BP/RP proxies from the XP grid (photon sums 330-680 and 630-1020 nm)
    bp = np.sum((f*wl)[(wl >= 330) & (wl < 680)]); rp = np.sum((f*wl)[(wl >= 630) & (wl <= 1020)])
    return -2.5*np.log10(bp/rp)
regs = {}
for name, cols in (('bprp', ['bp_rp']), ('red_only', ['mid', 'far'])):
    A = np.column_stack([np.ones(len(Cmp))] + [Cmp[k].values for k in cols])
    keep = np.ones(len(Cmp), bool)
    for it in range(3):
        coef = np.linalg.lstsq(A[keep], Cmp.blue.values[keep], rcond=None)[0]
        res = Cmp.blue.values - A @ coef
        keep = np.abs(res - np.median(res)) < 4*rstd(res)
    regs[name] = (cols, coef, rstd(res), res)
def resid(f, name, bprp_obs=2.756149, f_ref=None):
    cols, coef, sd, _ = regs[name]
    d = idx(f)
    if 'bp_rp' in cols:
        d['bp_rp'] = bprp_obs + (synth_bprp(f) - synth_bprp(f_ref) if f_ref is not None else 0.)
    return d['blue'] - np.r_[1, [d[k] for k in cols]] @ coef
out['consistent'] = {}
for name in regs:
    cols, coef, sd, res = regs[name]
    r_obs = resid(ft, name)
    row = dict(regressors=cols, robust_sd=float(sd), target_resid=float(r_obs), target_percentile=float(np.mean(res < r_obs)*100), pred={})
    for TWD in [4500, 5000, 5500, 6000, 6500, 7000, 8000]:
        Bw = planck(TWD); s = fA/(1 - fA)*np.sum(ft*wl*tw)/np.sum(Bw*wl*tw); fw = s*Bw
        m3 = ft - fw                     # the M3V alone under the WD hypothesis
        dres = r_obs - resid(m3, name, f_ref=ft)     # residual the WD adds (observed minus M3V-only)
        zW = (r_obs - dres)/sd           # residual of the inferred M3V-only spectrum vs locus
        zP = r_obs/sd
        row['pred'][str(TWD)] = dict(wd_adds_resid=float(dres), z_m3v_only_if_wd=float(zW), dchi2_wd_minus_planet=float(zW**2 - zP**2),
                                     pct_m3v_only_if_wd=float(np.mean(res < r_obs - dres)*100))
    out['consistent'][name] = row
    print(name, json.dumps(row))
json.dump(out, open(f'{RD}/out/xp_test.json', 'w'), indent=1)
