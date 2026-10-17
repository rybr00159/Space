"""r12: per-sector, odd/even and early/late depth consistency with a fixed transit shape (r03 planet_single ML
of each data set), linear depth scale per group with per-window offsets profiled; errors = beta-inflated
(per-sector fake-epoch calibration, rlib 'fake'). Also a per-sector duration scale (stretch of the time axis
about mid-transit) by a 1-D chi2 profile.
Output: out/per_sector.json"""
import sys, json, numpy as np
sys.path.insert(0, 'tic214311908/analysis/referee2_shape/scripts')
from rlib import *
out = {}
for ds in ('inv', 'ref'):
    S = json.load(open(f'{RD}/out/scen_{ds}_fake.json'))
    T0, P, k, b, lr, u1, u2 = S['planet_single']['x']
    D = load(ds, 'fake'); lc = LC(D)
    tm = lc.transit(T0, P, k, aR_from_rho(10**lr, P), abs(b), u1, u2) - 1
    d0 = -tm.min()
    def fit(sel):
        r = lc.f[sel] - 1; x = tm[sel]/d0; wi = lc.wi[sel]; ww = lc.w[sel]
        sw = np.bincount(wi, ww, lc.nw); sw[sw == 0] = 1
        r = r - (np.bincount(wi, ww*r, lc.nw)/sw)[wi]; x = x - (np.bincount(wi, ww*x, lc.nw)/sw)[wi]
        den = np.sum(ww*x*x); return float(-np.sum(ww*x*r)/den*1e3*-1), float(1/np.sqrt(den)*1e3)
    res = {}
    for s in np.unique(lc.sector):
        res[f'S{int(s)}'] = fit(lc.sector == s)
    ep = lc.epoch
    res['odd'] = fit(ep % 2 == 1); res['even'] = fit(ep % 2 == 0)
    res['early_S12_S66'] = fit(lc.sector <= 66); res['late_S93_S104'] = fit(lc.sector >= 93)
    res['all'] = fit(np.ones(lc.n, bool))
    sec = [res[f'S{int(s)}'] for s in np.unique(lc.sector)]
    d = np.array([v[0] for v in sec]); e = np.array([v[1] for v in sec])
    mu = np.sum(d/e**2)/np.sum(1/e**2)
    # duration scale per sector
    dur = {}
    for s in np.unique(lc.sector):
        sel = lc.sector == s
        grid = np.linspace(0.6, 1.4, 41); cs = []
        for g in grid:
            # stretch: same shape with a/R scaled by 1/g (duration x g), depth fixed to the global value
            m = lc.transit(T0, P, k, aR_from_rho(10**lr, P)/g, abs(b), u1, u2)
            r = lc.f - m; o = np.bincount(lc.wi, r*lc.w, lc.nw)/lc.sw; r = r - o[lc.wi]
            cs.append(np.sum((r*r*lc.w)[sel]))
        cs = np.array(cs); i = np.argmin(cs)
        ok = cs < cs[i] + 1
        dur[f'S{int(s)}'] = (float(grid[i]), float(grid[ok].min()), float(grid[ok].max()))
    out[ds] = dict(depth_ppt={kk: dict(depth=v[0], err=v[1]) for kk, v in res.items()}, sector_chi2=float(np.sum(((d - mu)/e)**2)), dof=len(d) - 1,
                   odd_minus_even=float(res['odd'][0] - res['even'][0]), odd_minus_even_err=float(np.hypot(res['odd'][1], res['even'][1])),
                   late_minus_early=float(res['late_S93_S104'][0] - res['early_S12_S66'][0]), late_minus_early_err=float(np.hypot(res['late_S93_S104'][1], res['early_S12_S66'][1])),
                   duration_scale=dur)
    print(ds, json.dumps(out[ds]))
json.dump(out, open(f'{RD}/out/per_sector.json', 'w'), indent=1)
