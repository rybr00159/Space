"""s01: analysis data sets and noise characterisation.
INV  = ../phot/lc_M3V_all.csv, quality_ok==1, flare_flag==0, finite. Fractional flux rescaled from the 935 e-/s
       (catalogue) to the 979 e-/s (PRF-unit) normalisation: f = 1 + (flux_norm - 1) * 935/979, err * 935/979.
       Transit-window subset |dt| < 0.30 d (dt from the reference ephemeris), and the full series (for phase 0.5 /
       phase-curve work and injections).
REF  = ../referee_phot/lc_M3V_referee_windows.csv (already 1 + dflux/979; per-pixel quadratic baselines fitted with
       +-0.07 d masked).
Noise: per sector, the white scale s = robust std(OOT residual)/median(err) and the time-averaging red-noise factor
beta (Pont+06, Winn+08) on bins of 20, 30, 45, 60, 90 min of the OOT data (INV: full series outside +-0.1 d of
transits and of phase 0.5). For REF the in-window beta (bins inside windows) is biased low by the per-pixel
baselines, so beta = mean(err_fake, err_fake_robust)/err_white from the referee's fake-window calibration
(../referee_phot/out/per_sector.csv) is adopted (1.0-2.0). Adopted per-point sigma =
err * s * beta, with beta = median over the bin sizes (floored at 1).
Output: data/inv_win.npz, data/inv_full.npz, data/ref_win.npz, out/noise.json"""
import numpy as np, pandas as pd, json
B = 'tic214311908/analysis'
T0, P = 3840.91414, 5.6915087
SC = 935.0/979.0
def rstd(x):
    x = x[np.isfinite(x)]; return 1.4826*np.median(np.abs(x-np.median(x)))
def beta_factor(t, r, cad_s, groups, binmins=(20,30,45,60,90)):
    out = {}
    s1 = np.std(r)
    for bm in binmins:
        N = max(1, int(round(bm*60/cad_s)))
        if N < 2: continue
        means = []
        for g in np.unique(groups):
            rg = r[groups==g]; tg = t[groups==g]
            # contiguous chunks of N points without gaps > 1.5 cadences
            br = np.where(np.diff(tg) > 1.5*cad_s/86400)[0]+1
            for seg in np.split(np.arange(len(rg)), br):
                for i in range(0, len(seg)-N+1, N):
                    means.append(rg[seg[i:i+N]].mean())
        means = np.array(means); M = len(means)
        if M < 8: continue
        exp = s1/np.sqrt(N)*np.sqrt(M/(M-1.))
        out[bm] = float(np.std(means)/exp)
    return out
noise = {}
# ---------------- INV
d = pd.read_csv(f'{B}/phot/lc_M3V_all.csv')
q = d[(d.quality_ok==1) & np.isfinite(d.flux_norm) & (d.flare_flag==0)].copy()
q['f'] = 1 + (q.flux_norm - 1)*SC; q['e'] = q.flux_err*SC
ep = np.round((q.time_btjd_tdb - T0)/P); q['epoch'] = ep.astype(int); q['dt'] = q.time_btjd_tdb - (T0 + ep*P)
ph = ((q.time_btjd_tdb - T0)/P) % 1
q['phase'] = ph
oot = (np.abs(q.dt) > 0.10) & (np.abs(ph-0.5) > 0.02)
noise['INV'] = {}
sig = np.zeros(len(q))
for s in sorted(q.sector.unique()):
    m = (q.sector==s).values
    r = (q.f-1).values[m & oot.values]; e = q.e.values[m & oot.values]
    sc = rstd(r)/np.median(e)
    bf = beta_factor(q.time_btjd_tdb.values[m & oot.values], r - np.median(r), q.cadence_s.values[m][0], np.zeros(len(r)))
    beta = max(1.0, float(np.median(list(bf.values()))))
    noise['INV'][int(s)] = dict(scale=float(sc), beta_bins=bf, beta=beta, rstd_ppt=float(rstd(r)*1e3), median_err_ppt=float(np.median(e)*1e3), n_oot=int(len(r)))
    sig[m] = q.e.values[m]*sc*beta
    print('INV', s, noise['INV'][int(s)])
q['sig'] = sig
win = np.abs(q.dt) < 0.30
cols = dict(t=q.time_btjd_tdb.values, f=q.f.values, e=q.e.values, sig=q.sig.values, sector=q.sector.values, cad=q.cadence_s.values, epoch=q.epoch.values, dt=q.dt.values, phase=q.phase.values,
            fA=1 + (q.flux_norm_modelA.values - 1)*SC, eA=q.flux_err_modelA.values*SC, raw=q.flux_raw_norm.values, trend=q.trend.values, dflux=q.dflux_e_s.values)
np.savez('data/inv_full.npz', **cols)
np.savez('data/inv_win.npz', **{k: v[win.values] for k, v in cols.items()})
print('INV windows', win.sum(), 'full', len(q))
# ---------------- REF
r = pd.read_csv(f'{B}/referee_phot/lc_M3V_referee_windows.csv')
r = r[np.isfinite(r.flux_norm)].copy()
ep = np.round((r.time_btjd_tdb - T0)/P); r['dt'] = r.time_btjd_tdb - (T0 + ep*P)
noise['REF'] = {}
sig = np.zeros(len(r)); 
for s in sorted(r.sector.unique()):
    m = (r.sector==s).values; o = m & (r.oot.values==1)
    res = r.flux_norm.values[o]-1; e = r.flux_err.values[o]
    sc = rstd(res)/np.median(e)
    bf = beta_factor(r.time_btjd_tdb.values[o], res - np.median(res), r.cadence_s.values[m][0], r.window_id.values[o]*1000+r.epoch.values[o])
    beta_win = max(1.0, float(np.median(list(bf.values())))) if len(bf) else 1.0
    # in-window beta is biased low (quadratic per-pixel baselines absorb slow noise); adopt the referee's
    # fake-window calibration instead: beta = mean(err_fake, err_fake_robust)/err_white (../referee_phot/out/per_sector.csv)
    ps = pd.read_csv(f'{B}/referee_phot/out/per_sector.csv').set_index('sector')
    beta = max(1.0, float(0.5*(ps.loc[s,'err_fake'] + ps.loc[s,'err_fake_robust'])/ps.loc[s,'err_white']))
    noise['REF'][int(s)] = dict(beta_in_window=beta_win, scale=float(sc), beta_bins=bf, beta=beta, rstd_ppt=float(rstd(res)*1e3), median_err_ppt=float(np.median(e)*1e3), n_oot=int(o.sum()))
    sig[m] = r.flux_err.values[m]*sc*beta
    print('REF', s, noise['REF'][int(s)])
np.savez('data/ref_win.npz', t=r.time_btjd_tdb.values, f=r.flux_norm.values, e=r.flux_err.values, sig=sig, sector=r.sector.values, cad=r.cadence_s.values,
         epoch=r.epoch.values, dt=r.dt.values, win=r.window_id.values, oot=r.oot.values)
print('REF', len(r))
json.dump(noise, open('out/noise.json','w'), indent=1)
