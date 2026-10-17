"""b07: occultation (secondary-eclipse) depth scan for ECCENTRIC orbits. An eccentric BD occults at phase
0.5 + ~(2/pi) e cos w with a duration that differs from the transit's, so the phase-0.5 limit of ../referee2_shape is not
enough. Data: INV-clean full light curve (delivered model-B detrended flux minus the referee2_search RRab model,
../referee2_search/data/m3v_clean_r2.npz; flares and real transit windows +-0.1 d removed). For phase phi = 0.20-0.80
(step 0.0025) and box duration 0.4, 0.6, 0.85, 1.2, 1.7, 2.4 h: per-epoch box depth with a local linear baseline in
+-0.3 d, combined over epochs (inverse-variance), formal error; empirical calibration: the robust scatter of
depth/error over all phases (a phase scan is its own null, apart from a real occultation) -> sigma_emp.
Output: out/b07_secondary_scan.npz, out/b07_secondary_scan.json"""
import sys; sys.path.insert(0, 'tic214311908/analysis/bd/scripts')
from blib import *
import pandas as pd
z = np.load(f'{B}/referee2_search/data/m3v_clean_r2.npz')
d = pd.read_csv(f'{B}/phot/lc_M3V_all.csv', usecols=['time_btjd_tdb', 'quality_ok', 'flux_err', 'flare_flag'])
d = d[d.quality_ok == 1]; d['tk'] = np.round(d.time_btjd_tdb.values, 7)
e_map = pd.Series(d.flux_err.values*935.97/979, index=d.tk); fl_map = pd.Series(d.flare_flag.values, index=d.tk)
tk = np.round(z['t'], 7)
ok = np.isin(tk, e_map.index)
t = z['t'][ok]; f = z['clean'][ok]; s = z['s'][ok]; e = e_map.loc[tk[ok]].values; fl = fl_map.loc[tk[ok]].values.astype(bool)
ph0 = ((t - T0)/P) % 1
m = ~fl & np.isfinite(f) & (np.minimum(ph0, 1 - ph0)*P > 0.1)
t, f, s, e = t[m], f[m], s[m], e[m]
# per-sector white scale from robust scatter
sc = {int(x): rstd(f[s == x] - 1)/np.median(e[s == x]) for x in np.unique(s)}
sig = e*np.array([sc[int(x)] for x in s])
phis = np.arange(0.20, 0.8001, 0.0025)
durs = np.array([0.4, 0.6, 0.85, 1.2, 1.7, 2.4])/24
DEP = np.full((len(phis), len(durs)), np.nan); ERR = np.full_like(DEP, np.nan)
ep_all = np.floor((t - T0)/P).astype(int)
for i, ph in enumerate(phis):
    tc_all = T0 + (ep_all + ph)*P
    dt = t - tc_all
    for j, du in enumerate(durs):
        num = 0.; den = 0.
        for ep in np.unique(ep_all):
            w = (ep_all == ep) & (np.abs(dt) < 0.3)
            if w.sum() < 30: continue
            box = (np.abs(dt[w]) < du/2).astype(float)
            if box.sum() < 3 or box.sum() > w.sum() - 10: continue
            X = np.vstack([np.ones(w.sum()), dt[w], -box]).T
            W = 1/sig[w]**2
            A = X.T @ (X*W[:, None])
            try: c = np.linalg.solve(A, X.T @ (W*(f[w] - 1))); var = np.linalg.inv(A)[2, 2]
            except np.linalg.LinAlgError: continue
            num += c[2]/var; den += 1/var
        if den > 0: DEP[i, j] = num/den; ERR[i, j] = 1/np.sqrt(den)
    if i % 40 == 0: print(round(ph, 3), np.round(DEP[i]*1e3, 3), np.round(ERR[i]*1e3, 3), flush=True)
zs = DEP/ERR
beta = np.array([rstd(zs[:, j]) for j in range(len(durs))])
std_ = np.array([np.nanstd(zs[:, j]) for j in range(len(durs))])
emp = ERR*np.maximum(beta, std_)[None, :]
i5 = np.argmin(np.abs(phis - 0.5))
out = dict(phases=phis, durations_h=durs*24, beta_robust=beta, beta_std=std_, scale=sc,
           depth_ppt_at_0p5=DEP[i5]*1e3, err_emp_ppt_at_0p5=emp[i5]*1e3,
           max_abs_z=float(np.nanmax(np.abs(DEP/emp))), max_depth_z=float(np.nanmax(DEP/emp)),
           limit3sig_ppt_median_over_phase=[float(np.nanmedian(DEP[:, j] + 3*emp[:, j])*1e3) for j in range(len(durs))])
np.savez(f'{BD}/out/b07_secondary_scan.npz', phases=phis, durs_d=durs, depth=DEP, err=ERR, err_emp=emp)
jdump(dict(method=__doc__, **out), f'{BD}/out/b07_secondary_scan.json')
print(json.dumps({k: v for k, v in out.items() if k not in ('phases',)}, default=lambda o: o.tolist() if hasattr(o, 'tolist') else str(o), indent=0))
