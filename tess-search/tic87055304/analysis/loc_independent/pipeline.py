"""Difference-image machinery shared by the real analysis (s02), injections (s03) and robustness runs.

A SectorContext fixes, for one sector: the cube, the PSF model (TESS PRF or Gaussian of fitted width),
the registered pixel positions of all Gaia stars (WCS + Gaia propagated to the sector epoch + global
registration offset from s01), the fit region, and pre-computed PSF images at the fixed star positions.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from common import *

T12_INJ = 0.045  # ingress duration of the injected trapezoid (d); b~0.85, k~0.067 -> T12/T14 ~ 0.24


class SectorContext:
    def __init__(self, sec, psfkind='prf', reg_variant='local_r4', rfit=3.0, gauss_sigma=None, S=None):
        self.S = S if S is not None else Sector(sec)
        S = self.S
        reg = json.load(open(f'{HERE}/registration.json'))[str(sec)]
        self.reg = reg
        regkey = f"{'prf' if psfkind == 'prf' else 'gauss'}:{reg_variant}"
        self.dx = reg['variants'][regkey]['dx']
        self.dy = reg['variants'][regkey]['dy']
        if psfkind == 'prf':
            self.psf = TessPRF(S)
        else:
            sx, sy = gauss_sigma if gauss_sigma is not None else reg['gauss_sigma_xy']
            self.psf = GaussPSF(S, sx, sy)
        self.psfkind = psfkind
        g = S.gaia
        self.x = g.x_wcs.values + self.dx
        self.y = g.y_wcs.values + self.dy
        self.it, self.inb = S.itarget, S.ineigh
        self.pT = np.array([self.x[self.it], self.y[self.it]])
        self.pN = np.array([self.x[self.inb], self.y[self.inb]])
        self.rfit = rfit
        self.mask = circle_mask((S.ny, S.nx), self.pT[0], self.pT[1], rfit)
        self.jac = S.jac
        self.jinv = np.linalg.inv(S.jac)
        # candidate hosts: Gaia stars within 30" (physical: T<17.4 can give the depth; 'all': G<20)
        sep = g.sep_t.values
        self.i_phys30 = [i for i in np.where((sep < 30) & (g.tmag.values < 17.4))[0]]
        self.i_all30 = [i for i in np.where((sep < 30) & (g.phot_g_mean_mag.values < 20.0))[0]]
        self.i_phys60 = [i for i in np.where((sep < 60) & (g.tmag.values < 17.4))[0]]
        for lst in (self.i_phys30, self.i_all30, self.i_phys60):
            for i in (self.it, self.inb):
                if i not in lst:
                    lst.insert(0, i)
        self.ic3 = int(np.where(g.source_id.values == C3_ID)[0][0])
        self._cache = {}
        self.plane = os.environ.get('LOC_PLANE', '0') == '1'
        yy, xx = np.mgrid[0:S.ny, 0:S.nx]
        self.bkg_cols = [np.ones(int(self.mask.sum()))]
        if self.plane:
            self.bkg_cols += [(xx - self.pT[0])[self.mask], (yy - self.pT[1])[self.mask]]
        self.vecTN = (self.pN - self.pT)

    def psf_at(self, i):
        if i not in self._cache:
            self._cache[i] = self.psf(self.x[i], self.y[i])
        return self._cache[i]

    # ------------------------------------------------------------ windows
    def real_windows(self, **kw):
        tcs, ns = self.S.transit_epochs()
        out = []
        for tc, n in zip(tcs, ns):
            w = self.S.window(tc, **kw)
            w['epoch'] = int(n)
            if w['ok']:
                out.append(w)
        return out

    def fake_windows(self, shifts=np.arange(0.25, 0.751, 0.05), **kw):
        """Windows at tc + phi*P; never overlapping a real transit window (+-OOT_OUT)."""
        tcs, ns = self.S.transit_epochs()
        out = []
        for phi in shifts:
            for tc, n in zip(tcs, ns):
                tf = tc + phi * P
                if np.min(np.abs(tf - tcs)) < 2 * OOT_OUT + 0.02:
                    continue
                w = self.S.window(tf, **kw)
                w['epoch'] = int(n); w['phi'] = float(phi)
                if w['ok']:
                    out.append(w)
        return out

    def diff(self, win, regressors=None):
        return difference_image(self.S.t, self.S.flux, win, regressors=regressors)

    def trapezoid(self, t, tc):
        a = np.abs(t - tc)
        h = T14 / 2
        return np.clip((h - a) / T12_INJ, 0, 1)

    def inject_diff(self, win, host_xy, amp, cadence_integrate=True):
        """Difference image of the injected signal alone (exact, since D is linear in the cube):
        cube_inj(t) = -amp * shape(t) * PSF(host). shape is integrated over the exposure."""
        idx = np.concatenate([win['lft'], win['inn'], win['rgt']])
        t = self.S.t
        if cadence_integrate:
            sub = np.linspace(-0.5, 0.5, 11) * self.S.exptime
            shape = np.mean([self.trapezoid(t + s, win['tc']) for s in sub], axis=0)
        else:
            shape = self.trapezoid(t, win['tc'])
        prf = self.psf(host_xy[0], host_xy[1]) if not isinstance(host_xy, np.ndarray) or host_xy.ndim == 1 else host_xy
        # only the window cadences matter: build a small cube
        small = -amp * shape[idx][:, None, None] * prf[None]
        tt = t[idx]
        loc = dict(tc=win['tc'], inn=np.searchsorted(idx, win['inn']), lft=np.searchsorted(idx, win['lft']),
                   rgt=np.searchsorted(idx, win['rgt']))
        D, _, _ = difference_image(tt, small, loc)
        return D

    def in_window_shape_mean(self, win):
        sub = np.linspace(-0.5, 0.5, 11) * self.S.exptime
        t = self.S.t[win['inn']]
        return float(np.mean([self.trapezoid(t + s, win['tc']) for s in sub]))

    # ------------------------------------------------------------ analysis of one difference image
    def analyse(self, D, sig, free=True, extra=True):
        m = self.mask
        out = {}
        if free:
            r = fit_free_position(D, sig, self.psf, self.pT[0], self.pT[1], m, bkg=True, box=2.5, plane=self.plane)
            p = np.array([r['x'], r['y']])
            offT = self.jac @ (p - self.pT)
            offN = self.jac @ (p - self.pN)
            vT = self.jac @ self.vecTN
            s_par = float(offT @ vT / (vT @ vT))
            u = vT / np.linalg.norm(vT)
            perp = float(offT[0] * (-u[1]) + offT[1] * u[0])
            covsky = self.jac @ r['cov'] @ self.jac.T
            out.update(free_x=r['x'], free_y=r['y'], free_amp=r['amp'], free_chi2=r['chi2'],
                       offT_E=float(offT[0]), offT_N=float(offT[1]), offT=float(np.hypot(*offT)),
                       offN_E=float(offN[0]), offN_N=float(offN[1]), offN=float(np.hypot(*offN)),
                       s_par=s_par, perp_arcsec=perp, free_at_bound=r['at_bound'],
                       free_cov_sky=covsky.tolist())
        # two-source fixed amplitudes
        w = 1 / sig[m]
        b = D[m] * w
        def lin(idx):
            A = np.column_stack([self.psf_at(i)[m] for i in idx] + self.bkg_cols) * w[:, None]
            c, *_ = np.linalg.lstsq(A, b, rcond=None)
            cov = np.linalg.pinv(A.T @ A)
            chi2 = float(((A @ c - b) ** 2).sum())
            c = np.r_[c[:len(idx)], c[len(idx)]]  # keep [amps..., const] layout (plane terms dropped)
            return c, cov, chi2
        c2, cov2, chi2_2 = lin([self.it, self.inb])
        aT, aN = c2[0], c2[1]
        out.update(a_T=float(aT), a_N=float(aN), sa_T=float(np.sqrt(cov2[0, 0])), sa_N=float(np.sqrt(cov2[1, 1])),
                   f_N=float(aN / (aT + aN)) if (aT + aN) != 0 else np.nan, chi2_TN=chi2_2)
        cT, _, chi2_T = lin([self.it])
        cN, _, chi2_N = lin([self.inb])
        out.update(chi2_T=chi2_T, chi2_N=chi2_N, dchi2_TminusN=chi2_T - chi2_N, amp_T_only=float(cT[0]),
                   amp_N_only=float(cN[0]))
        if extra:
            from scipy.optimize import lsq_linear
            def nnls(idx):
                A = np.column_stack([self.psf_at(i)[m] for i in idx] + self.bkg_cols) * w[:, None]
                lb = np.zeros(A.shape[1]); lb[len(idx):] = -np.inf
                r = lsq_linear(A, b, bounds=(lb, np.inf), method='bvls')
                return r.x[:len(idx) + 1], float(((A @ r.x - b) ** 2).sum())
            ids = lambda idx: [int(self.S.gaia.source_id[i]) for i in idx]
            c3, cov3, chi3 = lin([self.it, self.inb, self.ic3])
            out['TN3_ids'] = ids([self.it, self.inb, self.ic3])
            out['TN3_amp'] = c3[:-1].tolist(); out['TN3_amp_err'] = np.sqrt(np.diag(cov3))[:-1].tolist()
            out['TN3_frac'] = (c3[:-1] / c3[:-1].sum()).tolist()
            for name, idx in [('phys30', self.i_phys30), ('all30', self.i_all30)]:
                cl, _, chil = lin(idx)
                cn, chin = nnls(idx)
                out[f'{name}_ids'] = ids(idx)
                out[f'{name}_amp_lin'] = cl[:-1].tolist()
                out[f'{name}_frac_lin'] = (cl[:-1] / cl[:-1].sum()).tolist()
                out[f'{name}_amp_nnls'] = cn[:-1].tolist()
                out[f'{name}_frac_nnls'] = (cn[:-1] / max(cn[:-1].sum(), 1e-9)).tolist()
                out[f'{name}_chi2_nnls'] = chin
            # single-source chi2 for every physically possible host within 60"
            ch = {}
            for i in self.i_phys60:
                _, _, cc = lin([i])
                ch[int(self.S.gaia.source_id[i])] = cc
            out['single_host_chi2'] = ch
            _, _, chi_c3 = lin([self.ic3])
            out['dchi2_C3minusN'] = chi_c3 - chi2_N
        return out


def stack(Ds, sigs, scalar_weights=True, mask=None):
    Ds = np.asarray(Ds); sigs = np.asarray(sigs)
    if scalar_weights:
        if mask is None:
            mask = np.ones(Ds.shape[1:], bool)
        wk = 1.0 / np.array([np.median(s[mask] ** 2) for s in sigs])
        W = wk[:, None, None] * np.ones_like(Ds)
    else:
        W = 1.0 / sigs ** 2
    D = (W * Ds).sum(0) / W.sum(0)
    sig = np.sqrt((W ** 2 * sigs ** 2).sum(0)) / W.sum(0)
    return D, sig


def pointing_series(ctx, rmin=2.5, rmax=9.0):
    """Per-cadence pointing (dx,dy, in px) from stars other than target/neighbour, by linear
    regression of each frame onto [mean image, dM/dx, dM/dy, 1] over pixels rmin<r<rmax from the
    target (derivatives from the registered PSF model of the mean image)."""
    S = ctx.S
    g = S.gaia
    good = S.good
    mimg = S.flux[good].mean(0)
    # model gradient via PSF model of all T<16 stars in stamp with catalogue-free fluxes from mean image
    sel = np.where((g.tmag < 16) & (ctx.x > -3) & (ctx.x < S.nx + 2) & (ctx.y > -3) & (ctx.y < S.ny + 2))[0]
    h = 0.02
    A = np.column_stack([ctx.psf(ctx.x[i], ctx.y[i]).ravel() for i in sel] + [np.ones(S.nx * S.ny)])
    f, *_ = np.linalg.lstsq(A, mimg.ravel(), rcond=None)
    def model(dx, dy):
        mm = np.zeros((S.ny, S.nx))
        for k, i in enumerate(sel):
            mm += f[k] * ctx.psf(ctx.x[i] + dx, ctx.y[i] + dy)
        return mm
    gx = (model(h, 0) - model(-h, 0)) / (2 * h)
    gy = (model(0, h) - model(0, -h)) / (2 * h)
    yy, xx = np.mgrid[0:S.ny, 0:S.nx]
    r = np.hypot(xx - ctx.pT[0], yy - ctx.pT[1])
    pm = (r > rmin) & (r < rmax)
    X = np.column_stack([mimg[pm], gx[pm], gy[pm], np.ones(pm.sum())])
    Y = S.flux[:, pm].T  # npix x ncad
    Y = np.nan_to_num(Y)
    c, *_ = np.linalg.lstsq(X, Y, rcond=None)
    scale = c[0]
    with np.errstate(divide='ignore', invalid='ignore'):
        dxs = np.nan_to_num(c[1] / scale)
        dys = np.nan_to_num(c[2] / scale)
    return dxs, dys, scale
