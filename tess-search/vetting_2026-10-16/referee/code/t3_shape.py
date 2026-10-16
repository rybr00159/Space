"""TIC 214311906/214311908: stacked transit shape from TESScut photometry centred
on the M dwarf TIC 214311908 (r=1 px): locally detrended, folded, binned;
grid fit of trapezoid (depth, T14, T12/T14) to test flat vs V, plus odd/even."""
import sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/five_referee/code')
import numpy as np
from reflib import *

P, T0 = 5.69149, 3840.913
dts, fls, ns = [], [], []
for sec in (12, 39, 66, 93, 103, 104):
    cut = load_cut(214311906, sec)
    g = gaia_table(214311906, btjd_to_jyear(np.median(cut['time'])))
    a = g[g.source_id == 5925373272850784768].iloc[0]
    t, f, e, ex = aperture_lc(214311906, sec, r=1.0, center=(a['ra_ep'], a['dec_ep']))
    if np.median(np.diff(t)) > 0.015:   # skip 30-min S12 for shape (keep for depth only)
        continue
    n0, n1 = int(np.floor((t.min() - T0) / P)), int(np.ceil((t.max() - T0) / P))
    for n in range(n0, n1 + 1):
        c = T0 + n * P
        s = np.abs(t - c) < 0.3
        oot = s & (np.abs(t - c) > 1.3 / 24)
        if oot.sum() < 20 or np.sum(s & (np.abs(t - c) < 0.5 / 24)) < 3:
            continue
        co = np.polyfit(t[oot] - c, f[oot], 2)
        dts.append(t[s] - c); fls.append(f[s] / np.polyval(co, t[s] - c)); ns.append(np.full(s.sum(), n))
dt = np.concatenate(dts); fl = np.concatenate(fls); nn = np.concatenate(ns)
bins = np.arange(-0.3, 0.3001, 10 / 1440.)
idx = np.digitize(dt, bins)
bc, bf, be = [], [], []
for i in range(1, len(bins)):
    m = idx == i
    if m.sum() > 5:
        bc.append(dt[m].mean()); bf.append(fl[m].mean()); be.append(np.std(fl[m]) / np.sqrt(m.sum()))
bc, bf, be = map(np.array, (bc, bf, be))
# rescale errors from out-of-transit scatter
oot = np.abs(bc) > 1.3 / 24
scale = np.std(bf[oot] - 1) / np.median(be[oot])
be = be * max(1, scale)
grid = []
for T14 in np.arange(1.0, 2.41, 0.05) / 24:
    for fr in np.arange(0.05, 0.501, 0.025):
        T12 = fr * T14
        sh = 1 - supersample(trapezoid, bc, 10 / 1440., t0=0.0, depth=1.0, T14=T14, T12=T12)
        # analytic best depth and offset
        X = np.vstack([np.ones_like(bc), -sh]).T
        b, *_ = np.linalg.lstsq(X / be[:, None], bf / be, rcond=None)
        chi = np.sum(((bf - X @ b) / be) ** 2)
        grid.append((T14 * 24, fr, b[1], chi))
grid = np.array(grid)
best = grid[np.argmin(grid[:, 3])]
print(f'stacked binned profile (S39-S104, {len(np.unique(nn))} events): best T14={best[0]:.2f} h, T12/T14={best[1]:.3f}, depth={best[2]*1e6:.0f} ppm, chi2={best[3]:.1f} for {len(bc)-4} dof (err scale {scale:.2f})')
# marginal min chi2 vs T12/T14
out = {}
for fr in (0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5):
    m = np.isclose(grid[:, 1], fr)
    out[fr] = float(grid[m, 3].min() - best[3])
print('delta-chi2 (min over T14) vs ingress fraction T12/T14: ' + ', '.join(f'{k}: {v:.1f}' for k, v in out.items()))
# odd/even on the stacked data
for par in (0, 1):
    m = (nn % 2 == par) & (np.abs(dt) < 0.5 / 24)
    mo = (nn % 2 == par) & (np.abs(dt) > 1.3 / 24)
    print(f'  parity {par}: mean in-transit (|dt|<0.5h) deficit {(1-np.mean(fl[m]))*1e6:.0f} +- {np.std(fl[mo])/np.sqrt(m.sum())*1e6:.0f} ppm (n_in={m.sum()})')
save_json(dict(best=best.tolist(), dchi2_vs_ingress_frac=out, binned=dict(t=bc.tolist(), f=bf.tolist(), e=be.tolist())), f'{BASE}/results/t3_shape.json')
