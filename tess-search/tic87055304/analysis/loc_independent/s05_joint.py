"""Joint all-sector localisation in the sky plane.

One point source at a common sky offset (dE, dN) from the target (registered frame of each sector),
free amplitude and constant per sector; chi2 summed over sectors using the stacked difference images
from s02.  Also a chi2 map on a sky grid (with each sector's chi2 rescaled by its own best-fit
chi2/dof so that badly-modelled sectors do not dominate), Delta chi2 at the Gaia stars, and a
bootstrap over transits (resampling within each sector).

usage: python s05_joint.py <real_tag> [prf|gauss] [reg_variant]
Outputs: joint_<tag>.json, fig_joint_chi2map_<tag>.png
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import least_squares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pipeline import *

tag = sys.argv[1] if len(sys.argv) > 1 else 'prf_main'
psfkind = sys.argv[2] if len(sys.argv) > 2 else 'prf'
regv = sys.argv[3] if len(sys.argv) > 3 else 'local_r4'
NB = 300
rng = np.random.default_rng(5)
arr = np.load(f'{HERE}/diff_{tag}.npz')
ctxs = {s: SectorContext(s, psfkind=psfkind, reg_variant=regv, rfit=3.0) for s in SECTORS}


def sector_chi2(ctx, D, sig, EN):
    p = ctx.pT + ctx.jinv @ np.asarray(EN)
    m = ctx.mask
    w = 1 / sig[m]
    A = np.column_stack([ctx.psf(p[0], p[1])[m], np.ones(m.sum())]) * w[:, None]
    b = D[m] * w
    c, *_ = np.linalg.lstsq(A, b, rcond=None)
    return A @ c - b, c


def joint_fit(data, scale=None, x0=(5.0, 5.0)):
    def resid(EN):
        rr = []
        for s in SECTORS:
            r, _ = sector_chi2(ctxs[s], data[s][0], data[s][1], EN)
            rr.append(r / (np.sqrt(scale[s]) if scale is not None else 1.0))
        return np.concatenate(rr)
    best = None
    for e in np.linspace(-15, 20, 8):
        for n in np.linspace(-15, 20, 8):
            c2 = (resid([e, n]) ** 2).sum()
            if best is None or c2 < best[0]:
                best = (c2, [e, n])
    r = least_squares(resid, best[1], diff_step=1e-3, x_scale=1.0)
    cov = np.linalg.inv(r.jac.T @ r.jac)
    return r.x, cov, float((r.fun ** 2).sum()), len(r.fun)


data = {s: (arr[f's{s}_D'], arr[f's{s}_sig']) for s in SECTORS}
# per-sector chi2/dof at its own best single-source fit, used to rescale
real = json.load(open(f'{HERE}/diff_{tag}.json'))
scale = {s: max(real[str(s)]['stacked']['free_chi2_per_dof'], 1e-3) for s in SECTORS}
xy, cov, c2, n = joint_fit(data, scale)
g = ctxs[66].S.gaia
neighEN = np.array([g.de[ctxs[66].inb], g.dn[ctxs[66].inb]])
c3EN = np.array([g.de[ctxs[66].ic3], g.dn[ctxs[66].ic3]])
# bootstrap over transits within each sector
boots = []
for b in range(NB):
    db = {}
    for s in SECTORS:
        Dk = arr[f's{s}_Dk']; Sk = arr[f's{s}_sigk']
        k = rng.integers(0, len(Dk), len(Dk))
        db[s] = stack(Dk[k], Sk[k], mask=ctxs[s].mask)
    xb, _, _, _ = joint_fit(db, scale, x0=xy)
    boots.append(xb)
boots = np.array(boots)
bcov = np.cov(boots.T)
vT = neighEN
s_par = float(xy @ vT / (vT @ vT))
s_par_b = boots @ vT / (vT @ vT)
u = vT / np.linalg.norm(vT)
perp = float(-xy[0] * u[1] + xy[1] * u[0])
perp_b = -boots[:, 0] * u[1] + boots[:, 1] * u[0]
# chi2 map
E = np.linspace(-12, 22, 69); Nn = np.linspace(-12, 22, 69)
cmap = np.zeros((len(Nn), len(E)))
for i, nn in enumerate(Nn):
    for j, ee in enumerate(E):
        cmap[i, j] = sum((sector_chi2(ctxs[s], *data[s], [ee, nn])[0] ** 2).sum() / scale[s] for s in SECTORS)
cmin = cmap.min()
def dchi_at(EN):
    return float(sum((sector_chi2(ctxs[s], *data[s], EN)[0] ** 2).sum() / scale[s] for s in SECTORS) - c2)
cand = {}
for i in np.where((g.sep_t < 60) & (g.tmag < 17.4))[0]:
    cand[int(g.source_id[i])] = dict(dE=float(g.de[i]), dN=float(g.dn[i]), tmag=float(g.tmag[i]),
                                     dchi2_vs_best=dchi_at([g.de[i], g.dn[i]]))
res = dict(tag=tag, psf=psfkind, best_EN=xy.tolist(), formal_cov=cov.tolist(), boot_cov=bcov.tolist(),
           boot_sd_EN=np.sqrt(np.diag(bcov)).tolist(), chi2_scaled=c2, n=n, chi2_scale_per_sector=scale,
           offset_from_target=float(np.hypot(*xy)), offset_from_neigh=float(np.hypot(*(xy - neighEN))),
           neigh_EN=neighEN.tolist(), c3_EN=c3EN.tolist(), s_par=s_par, s_par_boot_sd=float(np.std(s_par_b)),
           s_par_boot_p16_p84=np.percentile(s_par_b, [16, 84]).tolist(), perp=perp, perp_boot_sd=float(np.std(perp_b)),
           dchi2_target=dchi_at([0, 0]), dchi2_neigh=dchi_at(neighEN), dchi2_c3=dchi_at(c3EN), candidates=cand)
# significance of target exclusion from bootstrap (Mahalanobis distance)
d = xy - 0
res['mahal_target_boot'] = float(np.sqrt(d @ np.linalg.solve(bcov, d)))
d = xy - neighEN
res['mahal_neigh_boot'] = float(np.sqrt(d @ np.linalg.solve(bcov, d)))
d = xy - c3EN
res['mahal_c3_boot'] = float(np.sqrt(d @ np.linalg.solve(bcov, d)))
to_json(res, f'{HERE}/joint_{tag}.json')
print(json.dumps({k: v for k, v in res.items() if k != 'candidates'}, indent=1))
print('candidates dchi2:', {str(k)[-6:]: round(v['dchi2_vs_best'], 1) for k, v in cand.items()})

fig, ax = plt.subplots(figsize=(7.5, 6.5))
lv = [1, 2.3, 6.2, 11.8, 25, 50, 100, 200, 400]
im = ax.contourf(E, Nn, np.clip(cmap - cmin, 0, 600), levels=np.r_[0, lv, 600], cmap='viridis_r')
ax.contour(E, Nn, cmap - cmin, levels=[2.3, 6.2, 11.8], colors='w', linewidths=0.8)
plt.colorbar(im, label=r'$\Delta\chi^2$ (joint, each sector rescaled by its $\chi^2$/dof)')
sel = (g.sep_t < 40) & (g.phot_g_mean_mag < 20.5)
ax.scatter(g.de[sel], g.dn[sel], s=np.clip(4 * (21 - g.phot_g_mean_mag[sel]) ** 1.6, 4, 400), facecolor='none', edgecolor='r', lw=0.8)
ax.plot(0, 0, 'r+', ms=16, mew=2, label='target TIC 87055304')
ax.plot(*neighEN, 'mx', ms=13, mew=2, label='TIC 87055315 (neighbour)')
ax.plot(*c3EN, 'cD', ms=7, mfc='none', mew=1.5, label='Gaia ...830336 (G=18.0)')
ax.plot(boots[:, 0], boots[:, 1], 'k.', ms=1.5, alpha=0.4, label='bootstrap over transits')
ax.errorbar(xy[0], xy[1], xerr=np.sqrt(bcov[0, 0]), yerr=np.sqrt(bcov[1, 1]), fmt='o', color='orange', mec='k', ms=7, label='joint best fit')
ax.set_xlim(22, -12); ax.set_ylim(-12, 22)
ax.set_xlabel('offset East (arcsec)'); ax.set_ylabel('offset North (arcsec)')
ax.set_title(f'Joint 5-sector difference-image localisation ({psfkind})\ncircle size ~ Gaia G brightness; contours 1,2,3 sigma (2 dof)')
ax.legend(fontsize=7, loc='lower left'); ax.set_aspect('equal')
plt.tight_layout(); plt.savefig(f'{HERE}/fig_joint_chi2map_{tag}.png', dpi=120); plt.close()
