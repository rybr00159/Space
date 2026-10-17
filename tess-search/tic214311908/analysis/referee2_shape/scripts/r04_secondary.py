"""r04: secondary-eclipse depth at phase 0.5 on the INV full series, re-detrended with the tested window masked.
For the test phase (0.5) and calibration phases (0.12-0.42, 0.58-0.88, step 0.01), the raw model-B flux is
re-detrended (wotan biweight 0.75 d, real transits and the tested windows +-0.1 d masked, as phot/p06), windows of
+-0.3 d are cut and a template is fitted for its depth with per-window offsets (white errors x per-sector scale).
Templates: (a) the planet shape (ML INV: k 0.12, b 0.5, rho 9.5, LD; T14 1.72 h), (b) the WD-transit shape (a small
body crossing the M3V on the same chord as a 1.35-h box primary: k 0.04, b 0.0, a/R from M_WD 0.45, LD).
Injection test: a 2.0-ppt template injected at phase 0.5 in the raw flux is recovered after re-detrending.
Output: out/secondary.json, out/secondary_rows.csv"""
import sys, json, numpy as np, pandas as pd
sys.path.insert(0, 'tic214311908/analysis/referee2_shape/scripts')
from rlib import *
from wotan import flatten
d = pd.read_csv(f'{B}/phot/lc_M3V_all.csv')
d = d[(d.quality_ok == 1) & (d.flare_flag == 0) & np.isfinite(d.flux_raw_norm) & np.isfinite(d.flux_err)].copy()
t = d.time_btjd_tdb.values; raw0 = d.flux_raw_norm.values; e = d.flux_err.values*935/979; sec = d.sector.values; cad = d.cadence_s.values.astype(int)
real = np.abs(((t - T0REF)/PREF + 0.5) % 1 - 0.5)*PREF < 0.1
fl = 1 + (d.flux_norm.values - 1)*935/979
scale = {int(s): rstd(fl[(sec == s) & ~real & np.isfinite(fl)] - 1)/np.median(e[(sec == s) & ~real & np.isfinite(fl)]) for s in np.unique(sec)}
sigw = e*np.array([scale[int(s)] for s in sec])
bj = json.load(open(f'{RD}/out/beta_inv.json'))
def measure(ph, raw):
    T0f = T0REF + ph*PREF
    fake = np.abs(((t - T0f)/PREF + 0.5) % 1 - 0.5)*PREF < 0.1
    f = np.full(len(t), np.nan)
    for s in np.unique(sec):
        m = sec == s; use = m & ~real & ~fake
        tr = flatten(t[use], raw[use], method='biweight', window_length=0.75, return_trend=True, break_tolerance=0.5, edge_cutoff=0)[1]
        ok = np.isfinite(tr)
        f[m] = raw[m]/np.interp(t[m], t[use][ok], tr[ok])
    f = 1 + (f - 1)*935/979
    ep = np.round((t - T0f)/PREF).astype(int); dt = t - (T0f + ep*PREF)
    w = (np.abs(dt) < 0.3) & np.isfinite(f) & ~(np.abs(((t - T0REF)/PREF + 0.5) % 1 - 0.5)*PREF < 0.3)
    D = dict(t=t[w], f=f[w], sig=sigw[w], sector=sec[w], cad=cad[w], epoch=ep[w], dt=dt[w], wid=sec[w]*1000 + ep[w])
    lc = LC(D)
    res = {}
    for name, args in (('planet', (0.12, aR_from_rho(9.5, PREF), 0.5, 0.125, 0.533)), ('wdtransit', (0.04, 34.0, 0.0, 0.125, 0.533))):
        tm = lc.transit(T0f, PREF, *args) - 1; tm = tm/(-tm.min())
        for wmode in ('white', 'beta'):
            ww = lc.w.copy()
            if wmode == 'beta':
                ww = ww/np.array([bj[str(int(s))]['beta'] for s in lc.sector])**2
            r = lc.f - 1; x = tm.copy(); wi = lc.wi
            sw = np.bincount(wi, ww, lc.nw)
            r = r - (np.bincount(wi, ww*r, lc.nw)/sw)[wi]; x = x - (np.bincount(wi, ww*x, lc.nw)/sw)[wi]
            den = np.sum(ww*x*x)
            res[f'{name}_{wmode}'] = (float(np.sum(ww*x*r)/den), float(1/np.sqrt(den)))   # depth (positive = dip): x is the template normalised to -1 at mid-event
    return res
rows = []
r05 = measure(0.5, raw0); print('phase 0.5', r05, flush=True)
# injection at phase 0.5: 2 ppt of each template shape
inj = {}
for name, args in (('planet', (0.12, aR_from_rho(9.5, PREF), 0.5, 0.125, 0.533)), ('wdtransit', (0.04, 34.0, 0.0, 0.125, 0.533))):
    D = dict(t=t, f=np.ones(len(t)), sig=np.ones(len(t)), sector=sec, cad=cad, epoch=np.zeros(len(t), int), dt=((t - T0REF - 0.5*PREF)/PREF + 0.5) % 1*PREF - 0.5*PREF, wid=sec)
    L = LC(D); tm = L.transit(T0REF + 0.5*PREF, PREF, *args) - 1; tm = tm/(-tm.min())
    raw = raw0 + 0.002*tm*979/935
    inj[name] = measure(0.5, raw)[f'{name}_white'][0] - r05[f'{name}_white'][0]
print('injected 2 ppt -> recovered', inj, flush=True)
phs = np.r_[np.arange(0.12, 0.425, 0.01), np.arange(0.58, 0.885, 0.01)]
for ph in phs:
    rr = measure(ph, raw0); rows.append(dict(phase=ph, **{f'{k}_d': v[0] for k, v in rr.items()}, **{f'{k}_e': v[1] for k, v in rr.items()}))
R = pd.DataFrame(rows); R.to_csv(f'{RD}/out/secondary_rows.csv', index=False)
out = dict(phase05={k: dict(depth_ppt=v[0]*1e3, formal_ppt=v[1]*1e3) for k, v in r05.items()}, injection_recovery_of_2ppt=inj, n_cal=len(R))
for k in r05:
    dd = R[f'{k}_d'].values
    emp = float(np.std(dd)); empr = float(rstd(dd))
    out['phase05'][k].update(empirical_std_ppt=emp*1e3, empirical_robust_ppt=empr*1e3, z=float(r05[k][0]/max(emp, empr)),
                             upper3_ppt=float((r05[k][0] + 3*max(emp, empr))*1e3), cal_mean_ppt=float(np.mean(dd)*1e3))
print(json.dumps(out, indent=1))
json.dump(out, open(f'{RD}/out/secondary.json', 'w'), indent=1)
