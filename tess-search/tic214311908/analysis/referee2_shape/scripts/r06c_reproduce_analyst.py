"""r06c: reproduce ../shape s16's RUWE likelihood with ITS assumptions (random inclination, uniform scan angles,
41 transits, sigma 0.35 / 0.22 mas, threshold statistic P(RUWE<=1.10)/P0(RUWE<=1.10), no 1e-3 floor; 20000 orbits),
then switch the assumptions one at a time (edge-on; NSL scan angles) to show which drives the difference.
Output: out/ruwe_reproduce.json"""
import json, numpy as np, pandas as pd
rng = np.random.default_rng(21)
B = 'tic214311908/analysis'
m = pd.read_csv(f'{B}/referee_star/data/mdw_cmd.csv')
r0 = m[(m.phot_g_mean_mag > 13.5) & (m.phot_g_mean_mag < 15.5) & (m.ruwe < 1.25) & (m.ipd_frac_multi_peak <= 2)].ruwe.values
pR0 = np.mean(r0 <= 1.10)
P = 5.6915087; M1 = 0.370; plx = 15.55; NTR = 41; NCCD = 9
# NSL scan angles of this source (best nu0) from r06
import sys
src = open(f'{B}/referee2_shape/scripts/r06_ruwe.py').read(); src = src[:src.index('obs = {k:')]
exec(src)
times_nsl, psi_nsl, pf_nsl = simulate(np.radians(60.), 1)
def sim(M2, sig, incl='random', scan='uniform', eps=0.02, ns=20000):
    q = M2/M1; beta = q/(1+q) - eps/(1+eps); a = (M1+M2)**(1/3.)*(P/365.25)**(2/3.); alpha = beta*a*plx
    res = np.empty(ns)
    for k in range(ns):
        if scan == 'uniform':
            t = np.sort(rng.uniform(2014.6, 2017.4, NTR)); psi = rng.uniform(0, 2*np.pi, NTR)
        else:
            sel = np.sort(rng.choice(len(times_nsl), NTR if NTR <= len(times_nsl) else len(times_nsl), replace=False))
            t = 2014.56 + (times_nsl[sel] - times_nsl.min())/365.25; psi = psi_nsl[sel]
        lam = 2*np.pi*(t - 2015.0); pf = np.sin(lam)*np.sin(psi) + 0.6*np.cos(lam)*np.cos(psi)
        A = np.column_stack([np.sin(psi), np.cos(psi), np.sin(psi)*(t-2016), np.cos(psi)*(t-2016), pf])
        inc = np.arccos(rng.uniform(-1, 1)) if incl == 'random' else np.pi/2
        Om = rng.uniform(0, 2*np.pi); ph = rng.uniform(0, 2*np.pi)
        u = 2*np.pi*t/(P/365.25) + ph
        x, y = np.cos(u), np.sin(u)*np.cos(inc)
        dE = alpha*(x*np.cos(Om) - y*np.sin(Om)); dN = alpha*(x*np.sin(Om) + y*np.cos(Om))
        al = dE*np.sin(psi) + dN*np.cos(psi)
        rr = al - A @ np.linalg.lstsq(A, al, rcond=None)[0]
        res[k] = np.sqrt(rng.choice(r0)**2 + NCCD*np.sum(rr**2)/sig**2/(len(t)*NCCD - 5))
    p = np.mean(res <= 1.10)
    return dict(P_le_1p10=float(p), m2lnL=float(-2*np.log(max(p, 0.5/ns)/pR0)), median=float(np.median(res)))
out = {}
for M2 in (0.2, 0.3, 0.45, 0.6):
    for sig in (0.22, 0.35):
        for incl, scan in (('random', 'uniform'), ('edgeon', 'uniform'), ('edgeon', 'nsl')):
            key = f'M{M2}_s{sig}_{incl}_{scan}'; out[key] = sim(M2, sig, incl, scan); print(key, out[key], flush=True)
json.dump(out, open(f'{B}/referee2_shape/out/ruwe_reproduce.json', 'w'), indent=1)
