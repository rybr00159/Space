"""r01: per-sector red-noise factor of the INV light curve from fake epochs (independent of ../shape s01).
For each of NF fake ephemerides (phase offsets 0.10-0.40 and 0.60-0.90, real transits and phase 0.5 avoided), the RAW
model-B flux (flux_raw_norm) is re-detrended with wotan biweight 0.75 d with BOTH the real (+-0.1 d) and fake (+-0.1 d)
windows masked, exactly as phot/p06; windows of +-0.30 d around the fake epochs are cut, and a fixed planet template
(k 0.119, b 0.5, rho 9.5, Claret LD) is fitted for a per-sector depth scale with per-window offsets and WHITE errors
(err x robust scale). beta_sector = std(fitted depth over fake epochs)/median(formal error) (also a robust version);
also the same for single windows (per-transit beta).
Output: out/beta_inv.json, out/beta_inv_rows.csv"""
import sys, json, numpy as np, pandas as pd
sys.path.insert(0, 'tic214311908/analysis/referee2_shape/scripts')
from rlib import *
from wotan import flatten
d = pd.read_csv(f'{B}/phot/lc_M3V_all.csv')
d = d[(d.quality_ok == 1) & (d.flare_flag == 0) & np.isfinite(d.flux_raw_norm) & np.isfinite(d.flux_err)].copy()
t = d.time_btjd_tdb.values; raw = d.flux_raw_norm.values; e = d.flux_err.values*935/979; sec = d.sector.values; cad = d.cadence_s.values.astype(int)
real = np.abs(((t - T0REF)/PREF + 0.5) % 1 - 0.5)*PREF < 0.1
phs = np.r_[np.linspace(0.10, 0.40, 31), np.linspace(0.60, 0.90, 31)]
rows = []
# white scale per sector from the delivered detrended flux (OOT)
fl = 1 + (d.flux_norm.values - 1)*935/979
scale = {}
for s in np.unique(sec):
    k = (sec == s) & ~real & np.isfinite(fl)
    scale[int(s)] = rstd(fl[k] - 1)/np.median(e[k])
sigw = e*np.array([scale[int(s)] for s in sec])
for ph in phs:
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
    tmpl = lc.transit(T0f, PREF, 0.119, aR_from_rho(9.5, PREF), 0.5, 0.125, 0.533) - 1
    # per-window offset profiled linear fit of depth scale A: f = 1 + o_w + A*tmpl
    wts = lc.w
    def fitA(sel):
        # remove weighted window means from data and template
        r = lc.f[sel] - 1; x = tmpl[sel]; wi = lc.wi[sel]; ww = wts[sel]
        sw = np.bincount(wi, ww, lc.nw); sw[sw == 0] = 1
        r = r - (np.bincount(wi, ww*r, lc.nw)/sw)[wi]; x = x - (np.bincount(wi, ww*x, lc.nw)/sw)[wi]
        den = np.sum(ww*x*x)
        if den <= 0: return np.nan, np.nan
        return np.sum(ww*x*r)/den, 1/np.sqrt(den)
    for s in np.unique(lc.sector):
        A, sA = fitA(lc.sector == s)
        rows.append(dict(phase=ph, sector=int(s), win=-1, A=A, sA=sA))
        for wv in np.unique(lc.wi[lc.sector == s]):
            sel = lc.wi == wv
            if np.sum(np.abs(tmpl[sel]) > 1e-4) < 2: continue
            A1, sA1 = fitA(sel)
            rows.append(dict(phase=ph, sector=int(s), win=int(wv), A=A1, sA=sA1))
    print(round(ph, 3), flush=True)
R = pd.DataFrame(rows); R.to_csv(f'{RD}/out/beta_inv_rows.csv', index=False)
out = {}
for s in np.unique(R.sector):
    q = R[(R.sector == s) & (R.win == -1)]; q1 = R[(R.sector == s) & (R.win >= 0)]
    b = float(np.std(q.A)/np.median(q.sA)); br = float(rstd(q.A.values)/np.median(q.sA))
    b1 = float(np.std(q1.A/q1.sA)); b1r = float(rstd((q1.A/q1.sA).values))
    out[int(s)] = dict(scale=float(scale[int(s)]), beta_sector=b, beta_sector_robust=br, beta_single=b1, beta_single_robust=b1r, n_fake=int(len(q)), n_single=int(len(q1)),
                       beta=max(1., 0.5*(b + br)))
    print(s, out[int(s)])
json.dump(out, open(f'{RD}/out/beta_inv.json', 'w'), indent=1)
