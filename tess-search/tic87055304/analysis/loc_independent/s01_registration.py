"""Step 2(a): registration of the mean out-of-transit direct image to Gaia DR3.

For each sector:
  * mean image of all good (QUALITY==0) cadences outside +-T14 of every transit;
  * model = sum_i F_i PSF(x_i+dx, y_i+dy) + b, with x_i,y_i = Gaia DR3 positions propagated to the
    sector epoch and projected with the TESScut WCS (APERTURE extension);
      - FREE fluxes (>=0): every star with T<17 within 60" of the target, plus T<15 stars whose
        position is within (R_fit + 1.5) px of the target;
      - all other T<17.5 stars within the stamp(+4 px): catalogue flux x one free common scale;
  * (dx,dy) free global offset solved by variable projection over a circular fit region.
Variants (to expose registration systematics):
  * R_fit = 2.5, 3, 4, 5 px ('local' registration, dominated by the target+neighbour blend);
  * 'external': target region (r<2.5 px) masked, annulus 2.5-8 px (registration from other stars only);
  * neighbour/target flux ratio fixed to the Gaia-predicted value;
  * PSF = TESS PRF (Bell & Higgins, MAST PRF files) or pixel-integrated Gaussian with fitted width.
  * WCS check: brightest stars (T<12.6) and target/neighbour fitted individually for (x,y,F)
    with all other stars held at the global solution.
Outputs: registration.json, fig_registration_sXX.png
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import least_squares, lsq_linear, minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import *


def mean_image(S):
    tcs, _ = S.transit_epochs()
    m = S.good.copy()
    for tc in tcs:
        m &= np.abs(S.t - tc) > T14
    img = S.flux[m].mean(0)
    err = S.flux[m].std(0) / np.sqrt(m.sum())
    return img, err, m.sum()


def star_sets(S, xt, yt, rfit):
    g = S.gaia
    d_t = np.hypot(g.x_wcs - xt, g.y_wcs - yt)
    instamp = (g.x_wcs > -4) & (g.x_wcs < S.nx + 3) & (g.y_wcs > -4) & (g.y_wcs < S.ny + 3)
    free = ((g.tmag < 17) & (g.sep_t < 60)) | ((g.tmag < 15) & (d_t < rfit + 1.5))
    free &= instamp
    fixed = instamp & (g.tmag < 17.5) & ~free & (d_t < rfit + 7)
    return np.where(free)[0], np.where(fixed)[0]


def build_design(S, psf, off, ifree, ifixed, mask, ratio_fix=None):
    g = S.gaia
    cols = []
    names = []
    free_list = list(ifree)
    if ratio_fix is not None:
        # target and neighbour tied: F_N = ratio * F_T
        free_list = [i for i in free_list if i not in (S.itarget, S.ineigh)]
        c = psf(g.x_wcs[S.itarget] + off[0], g.y_wcs[S.itarget] + off[1]) + \
            ratio_fix * psf(g.x_wcs[S.ineigh] + off[0], g.y_wcs[S.ineigh] + off[1])
        cols.append(c[mask]); names.append('T+N tied')
    for i in free_list:
        cols.append(psf(g.x_wcs[i] + off[0], g.y_wcs[i] + off[1])[mask]); names.append(int(i))
    fx = np.zeros(mask.shape)
    for i in ifixed:
        fx += 10 ** (-0.4 * (g.tmag[i] - 10.0)) * psf(g.x_wcs[i] + off[0], g.y_wcs[i] + off[1])
    cols.append(fx[mask]); names.append('faint_scale')
    cols.append(np.ones(mask.sum())); names.append('bkg')
    return np.column_stack(cols), names


def register(S, psf, img, sig, mask, ifree, ifixed, ratio_fix=None, off0=(0., 0.)):
    w = 1 / sig[mask]
    b = img[mask] * w

    def solve(off):
        A, names = build_design(S, psf, off, ifree, ifixed, mask, ratio_fix)
        Aw = A * w[:, None]
        lb = np.zeros(A.shape[1]); lb[-1] = -np.inf
        r = lsq_linear(Aw, b, bounds=(lb, np.inf), method='bvls')
        return Aw, r.x, names

    def resid(off):
        Aw, c, _ = solve(off)
        return Aw @ c - b

    r = least_squares(resid, np.array(off0, float), diff_step=1e-3, x_scale=0.05)
    Aw, coef, names = solve(r.x)
    J = r.jac
    cov = np.linalg.pinv(J.T @ J) * (r.fun ** 2).sum() / max(len(b) - len(coef) - 2, 1)
    model = (Aw @ coef) / w
    return dict(dx=float(r.x[0]), dy=float(r.x[1]), cov=cov, coef=coef, names=names,
                chi2=float((r.fun ** 2).sum()), npix=int(mask.sum()), model_masked=model)


def full_model(S, psf, off, ifree, ifixed, coef, names, ratio_fix=None):
    mask = np.ones((S.ny, S.nx), bool)
    A, _ = build_design(S, psf, off, ifree, ifixed, mask, ratio_fix)
    return (A @ coef).reshape(S.ny, S.nx)


def fit_gauss_width(S, img, sig, mask, ifree, ifixed):
    """Fit Gaussian sigma (x,y separately) together with offset (nested)."""
    best = {}

    def obj(p):
        sx, sy = p
        if sx < 0.3 or sy < 0.3 or sx > 2.5 or sy > 2.5:
            return 1e30
        psf = GaussPSF(S, sx, sy)
        r = register(S, psf, img, sig, mask, ifree, ifixed, off0=best.get('off', (0, 0)))
        best['off'] = (r['dx'], r['dy'])
        return r['chi2']

    res = minimize(obj, [0.8, 0.8], method='Nelder-Mead', options=dict(xatol=2e-3, fatol=1e-2, maxiter=80))
    return res.x


def single_star_fit(S, psf, img, sig, reg, ifree, ifixed, istar, rad=2.0):
    """Free (x,y,F) for one star; all other stars fixed at the global solution."""
    g = S.gaia
    off = (reg['dx'], reg['dy'])
    mfull = full_model(S, psf, off, ifree, ifixed, reg['coef'], reg['names'])
    # model of the star alone at the global solution
    names = reg['names']
    if istar in names:
        fstar = reg['coef'][names.index(istar)]
    else:
        fstar = reg['coef'][names.index('faint_scale')] * 10 ** (-0.4 * (g.tmag[istar] - 10))
    x0, y0 = g.x_wcs[istar] + off[0], g.y_wcs[istar] + off[1]
    others = mfull - fstar * psf(x0, y0)
    mask = circle_mask(img.shape, x0, y0, rad)
    res = fit_free_position(img - others, sig, psf, x0, y0, mask, bkg=True, box=1.0)
    dxy = np.array([res['x'] - x0, res['y'] - y0])
    return dict(source_id=int(g.source_id[istar]), tmag=float(g.tmag[istar]), x_pred=x0, y_pred=y0,
                dx=float(dxy[0]), dy=float(dxy[1]), dE_dN=(S.jac @ dxy).tolist(), flux=res['amp'])


def main():
    out = {}
    for sec in SECTORS:
        S = Sector(sec)
        img, err, nimg = mean_image(S)
        sig = np.sqrt(err ** 2 + (0.01 * np.abs(img)) ** 2 + 1.0)
        g = S.gaia
        xt, yt = g.x_wcs[S.itarget], g.y_wcs[S.itarget]
        tp = TessPRF(S)
        res = dict(sector=sec, camera=S.camera, ccd=S.ccd, col0=S.col0, row0=S.row0, epoch=S.epoch,
                   n_cad_mean=int(nimg), jac=S.jac.tolist(),
                   target_wcs_xy=[float(xt), float(yt)],
                   neigh_wcs_xy=[float(g.x_wcs[S.ineigh]), float(g.y_wcs[S.ineigh])],
                   neigh_offset_EN_arcsec=[float(g.de[S.ineigh]), float(g.dn[S.ineigh])],
                   variants={})
        # --- Gaussian width (fit once with R=4 local)
        ifree, ifixed = star_sets(S, xt, yt, 4.0)
        mask4 = circle_mask(img.shape, xt, yt, 4.0)
        sxy = fit_gauss_width(S, img, sig, mask4, ifree, ifixed)
        res['gauss_sigma_xy'] = [float(sxy[0]), float(sxy[1])]
        gp = GaussPSF(S, sxy[0], sxy[1])
        rT_exp = 10 ** (-0.4 * (g.tmag[S.ineigh] - g.tmag[S.itarget]))
        res['neigh_target_ratio_catalog'] = float(rT_exp)
        for psfname, psf in [('prf', tp), ('gauss', gp)]:
            for vname, rad, ext, rfix in [('local_r2.5', 2.5, False, None), ('local_r3', 3.0, False, None),
                                         ('local_r4', 4.0, False, None), ('local_r5', 5.0, False, None),
                                         ('local_r4_ratiofix', 4.0, False, rT_exp),
                                         ('external_2.5-8', 8.0, True, None)]:
                ifree, ifixed = star_sets(S, xt, yt, rad)
                mask = circle_mask(img.shape, xt, yt, rad)
                if ext:
                    mask &= ~circle_mask(img.shape, xt, yt, 2.5)
                    # do not let target/neighbour fluxes float on pixels that barely see them
                r = register(S, psf, img, sig, mask, ifree, ifixed, ratio_fix=rfix)
                names = r['names']
                def flux_of(i):
                    if rfix is not None and i in (S.itarget, S.ineigh):
                        ft = r['coef'][names.index('T+N tied')]
                        return ft if i == S.itarget else ft * rfix
                    return r['coef'][names.index(i)] if i in names else np.nan
                fT, fN = flux_of(S.itarget), flux_of(S.ineigh)
                dE, dN = S.jac @ np.array([r['dx'], r['dy']])
                rr = dict(dx=r['dx'], dy=r['dy'], sdx=float(np.sqrt(r['cov'][0, 0])), sdy=float(np.sqrt(r['cov'][1, 1])),
                          dE=float(dE), dN=float(dN), chi2=r['chi2'], npix=r['npix'], n_free=len(ifree),
                          n_fixed=len(ifixed), F_target=float(fT), F_neigh=float(fN),
                          ratio_NT=float(fN / fT) if fT > 0 else np.nan,
                          faint_scale=float(r['coef'][names.index('faint_scale')]), bkg=float(r['coef'][-1]))
                res['variants'][f'{psfname}:{vname}'] = rr
                print(sec, psfname, vname, 'dx=%.3f dy=%.3f  (dE,dN)=(%.2f,%.2f)"  FT=%.0f FN=%.0f N/T=%.3f (cat %.3f) chi2/n=%.1f' % (
                    r['dx'], r['dy'], dE, dN, fT, fN, fN / fT if fT > 0 else np.nan, rT_exp, r['chi2'] / r['npix']), flush=True)
                if vname == 'local_r4':
                    off = (r['dx'], r['dy'])
                    mod = full_model(S, psf, off, ifree, ifixed, r['coef'], names)
                    res[f'{psfname}_resid_rms_frac_r2'] = float(np.sqrt(np.mean(((img - mod) / img)[circle_mask(img.shape, xt, yt, 2.0)] ** 2)))
                    if psfname == 'prf':
                        model_prf, reg_prf, sets_prf = mod, r, (ifree, ifixed)
                    else:
                        model_g = mod
        # --- WCS check: individual fits of bright stars with PRF (global solution from external fit, r=8)
        ifree, ifixed = star_sets(S, xt, yt, 8.0)
        mask = circle_mask(img.shape, xt, yt, 8.0)
        regw = register(S, tp, img, sig, mask, ifree, ifixed)
        bright = np.where((g.tmag < 13.0) & (g.x_wcs > 1.5) & (g.x_wcs < S.nx - 2.5) & (g.y_wcs > 1.5) & (g.y_wcs < S.ny - 2.5))[0]
        cands = list(bright) + [S.ineigh]
        checks = []
        for i in cands:
            if i not in regw['names']:
                continue
            c = single_star_fit(S, tp, img, sig, regw, ifree, ifixed, i)
            c['is_target'] = bool(i == S.itarget); c['is_neigh'] = bool(i == S.ineigh)
            checks.append(c)
            print('   star', c['source_id'], 'T=%.2f' % c['tmag'], 'resid (dE,dN)=(%.2f,%.2f)"' % tuple(c['dE_dN']))
        # --- photocentre of the target+neighbour blend: mean image minus model of all OTHER stars,
        #     single PSF + constant at a free position, r<2.5 px  (registration-free reference point)
        for psfname, psf, vkey in [('prf', tp, 'prf:local_r4'), ('gauss', gp, 'gauss:local_r4_ratiofix')]:
            v = res['variants'][vkey]
            ifree, ifixed = star_sets(S, xt, yt, 4.0)
            mask = circle_mask(img.shape, xt, yt, 4.0)
            rfix = rT_exp if 'ratiofix' in vkey else None
            r = register(S, psf, img, sig, mask, ifree, ifixed, ratio_fix=rfix)
            off = (r['dx'], r['dy'])
            names = r['names']
            coef_o = r['coef'].copy()
            if rfix is None:
                for i in (S.itarget, S.ineigh):
                    coef_o[names.index(i)] = 0.0
            else:
                coef_o[names.index('T+N tied')] = 0.0
            others = full_model(S, psf, off, ifree, ifixed, coef_o, names, ratio_fix=rfix)
            m25 = circle_mask(img.shape, xt + off[0], yt + off[1], 2.5)
            fb = fit_free_position(img - others, sig, psf, xt + off[0], yt + off[1], m25, bkg=True, box=1.5)
            pT = np.array([xt + off[0], yt + off[1]]); pN = np.array([g.x_wcs[S.ineigh] + off[0], g.y_wcs[S.ineigh] + off[1]])
            vT = S.jac @ (pN - pT); oT = S.jac @ (np.array([fb['x'], fb['y']]) - pT)
            res[f'{psfname}_blend_photocentre'] = dict(x=fb['x'], y=fb['y'], offT_EN=oT.tolist(),
                                                       s_par=float(oT @ vT / (vT @ vT)),
                                                       s_par_expected_catalog=float(rT_exp / (1 + rT_exp)))
            print('   blend photocentre', psfname, 'offT (E,N)=(%.2f,%.2f)" s_par=%.3f (catalog flux ratio predicts %.3f)' % (
                oT[0], oT[1], oT @ vT / (vT @ vT), rT_exp / (1 + rT_exp)))
        res['wcs_check_global'] = dict(dx=regw['dx'], dy=regw['dy'], dE_dN=(S.jac @ np.array([regw['dx'], regw['dy']])).tolist())
        res['wcs_check_stars'] = checks
        nonT = [c for c in checks if not c['is_target'] and not c['is_neigh']]
        if len(nonT) > 1:
            arr = np.array([c['dE_dN'] for c in nonT])
            res['wcs_check_rms_EN_arcsec'] = arr.std(0, ddof=1).tolist()
            res['wcs_check_mean_EN_arcsec'] = arr.mean(0).tolist()
        out[str(sec)] = res
        # ---- figure
        fig, ax = plt.subplots(1, 4, figsize=(17, 4.4))
        ext = (-0.5, S.nx - 0.5, -0.5, S.ny - 0.5)
        v = np.percentile(img, [5, 99.5])
        ax[0].imshow(np.log10(np.clip(img, 1, None)), origin='lower', extent=ext, cmap='viridis')
        sel = g.tmag < 16
        off = (reg_prf['dx'], reg_prf['dy'])
        ax[0].scatter(g.x_wcs[sel] + off[0], g.y_wcs[sel] + off[1], s=np.clip(3 * (17 - g.tmag[sel]) ** 2, 2, 200),
                      facecolor='none', edgecolor='w', lw=0.6)
        ax[0].plot(xt + off[0], yt + off[1], 'r+', ms=12); ax[0].plot(g.x_wcs[S.ineigh] + off[0], g.y_wcs[S.ineigh] + off[1], 'mx', ms=10)
        for c in checks:
            ax[0].annotate('', xy=(c['x_pred'] + 20 * c['dx'], c['y_pred'] + 20 * c['dy']), xytext=(c['x_pred'], c['y_pred']),
                           arrowprops=dict(arrowstyle='->', color='orange'))
        ax[0].set_title(f'S{sec} mean image (log), Gaia T<16 (registered)\norange: per-star WCS residual x20')
        for k, (m, nm) in enumerate([(model_prf, 'TESS PRF'), (model_g, 'Gaussian')]):
            r = (img - m) / img
            ax[1 + k].imshow(r, origin='lower', extent=ext, cmap='RdBu_r', vmin=-0.1, vmax=0.1)
            ax[1 + k].plot(xt + off[0], yt + off[1], 'k+', ms=12)
            ax[1 + k].set_title(f'(data-model)/data, {nm}, R=4px fit')
            ax[1 + k].add_patch(plt.Circle((xt, yt), 4, fill=False, color='k', ls='--'))
        names_v = [k for k in res['variants']]
        dEs = [res['variants'][k]['dE'] for k in names_v]
        dNs = [res['variants'][k]['dN'] for k in names_v]
        for k, nm in enumerate(names_v):
            ax[3].plot(dEs[k], dNs[k], 'o' if nm.startswith('prf') else 's', label=nm, ms=5)
        ax[3].plot(g.de[S.ineigh] * 0 + 0, 0, 'r+', ms=14)
        ax[3].set_xlabel('registration offset dE (")'); ax[3].set_ylabel('dN (")')
        ax[3].legend(fontsize=5.5, ncol=2); ax[3].set_title('global WCS offset, all variants')
        ax[3].invert_xaxis()
        plt.tight_layout(); plt.savefig(f'{HERE}/fig_registration_s{sec:03d}.png', dpi=110); plt.close()
    to_json(out, f'{HERE}/registration.json')


if __name__ == '__main__':
    main()
