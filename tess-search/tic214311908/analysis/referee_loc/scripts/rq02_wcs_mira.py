"""Referee WCS calibration and Mira apparent position (independent code).

Per sector, median direct image over good cadences:
 (1) Global fit: all Gaia stars with T<16.5 predicted in/near the cutout (Gaia+PM, header WCS, colour shift model)
     with a common pixel offset (dx,dy) free; free fluxes for T<14.5, catalogue fluxes for fainter stars; the Mira
     handled as a free-position/free-flux PRF but its core (r < RMASK px) masked; background plane free.
     sigma^2 = |data| + (FRAC*data)^2.  Variants: RMASK 2.0/3.0/4.0 px, FRAC 0.02/0.05, colour model curve/linear/none,
     and 'far' (pixels within 4 px of the Mira masked entirely and Mira wings also masked).
 (2) Mira position with the offset fixed (baseline variant): pixels within R=3 px of the Mira, Mira position+flux
     free, other stars fixed at their apparent positions with fluxes (M3V, pair, ref free; others catalogue);
     background constant.  Variants: PRF convolved with a Gaussian of free width (PSF broadening for a very red
     source); R = 2/4 px.  Also a 'centroid' estimate: flux-weighted centroid within 1.5 px of (data - other stars)
     compared with the same centroid of the PRF model placed at the Gaia position.
Outputs: out/wcs_mira.json
"""
import sys; sys.path.insert(0, 'tic214311908/analysis/referee_loc/scripts')
from rl_lib import *
from scipy.optimize import least_squares
from scipy.ndimage import gaussian_filter

g = gaia()
OUT = {}


def global_fit(S, O, psf, cmodel='curve', RMASK=2.5, FRAC=0.03, far=False):
    S.off = np.zeros(2)
    xy0 = S.xy(g, cmodel)
    ins = (xy0[:, 0] > -3) & (xy0[:, 0] < S.nx + 2) & (xy0[:, 1] > -3) & (xy0[:, 1] < S.ny + 2) & (g['T'].values < 16.5)
    idx = np.where(ins)[0]
    im = int(np.where(g.source_id.values[idx] == ID['Mira'])[0][0])
    freef = g['T'].values[idx] < 14.5
    yy, xx = np.mgrid[0:S.ny, 0:S.nx]
    mxy = xy0[idx[im]]
    rr = np.hypot(xx - mxy[0], yy - mxy[1])
    mask = rr > (4.0 if far else RMASK)
    def build(p):
        d, dm = p[:2], p[2:4]
        cols, fixed = [], np.zeros_like(O)
        for k, j in enumerate(idx):
            pos = xy0[j] + (d + dm if k == im else d)
            img = psf(*pos)
            if freef[k] or k == im: cols.append(img)
            else: fixed = fixed + g.F.values[j] * img
        cols += [np.ones_like(O), xx - 10.0, yy - 10.0]
        return np.array(cols), fixed
    def solve(p):
        A, fixed = build(p)
        w = 1 / np.sqrt(np.abs(O[mask]) + (FRAC * O[mask]) ** 2 + 4.0)
        c, *_ = np.linalg.lstsq(A[:, mask].T * w[:, None], (O - fixed)[mask] * w, rcond=None)
        return c, ((np.tensordot(c, A, 1) + fixed) - O)[mask] * w
    best = None
    for gx in (-0.15, 0, 0.15):
        for gy in (-0.15, 0, 0.15):
            r = least_squares(lambda p: solve(p)[1], [gx, gy, 0, 0], diff_step=1e-3, x_scale=0.05,
                              bounds=([-1, -1, -1.5, -1.5], [1, 1, 1.5, 1.5]))
            if best is None or r.cost < best.cost: best = r
    r = best
    dof = mask.sum() - len(solve(r.x)[0]) - 4
    cov = np.linalg.pinv(r.jac.T @ r.jac) * max(1.0, 2 * r.cost / dof)
    return dict(off=r.x[:2], err=np.sqrt(np.diag(cov))[:2], mira_rel=r.x[2:4], chi2r=2 * r.cost / dof, nstar=len(idx),
                nfree=int(freef.sum()))


def blurred(psf, x, y, s):
    img = psf(x, y)
    return img if s < 0.02 else gaussian_filter(img, s, mode='constant')


def mira_fit(S, O, psf, R=3.0, blur=False, cmodel='curve', nbscale=False):
    xy = S.xy(g, cmodel)
    xyn = S.xy(g, use_colour=False)
    jm = int(np.where(g.source_id == ID['Mira'])[0][0])
    mg = xyn[jm]
    yy, xx = np.mgrid[0:S.ny, 0:S.nx]
    mask = np.hypot(xx - mg[0], yy - mg[1]) <= R
    near = np.where((np.hypot(xy[:, 0] - mg[0], xy[:, 1] - mg[1]) < R + 5) & (g['T'].values < 18.5))[0]
    # reference star free; M3V + pair A + pair B either at catalogue flux (baseline) or scaled by ONE common free factor
    trio = [j for j in near if g.source_id.values[j] in (ID['M3V'], ID['pairA'], ID['pairB'])]
    freej = [j for j in near if g.source_id.values[j] == ID['ref']]
    fixed = sum(g.F.values[j] * psf(*xy[j]) for j in near if j != jm and j not in freej and (j not in trio or not nbscale))
    extra = [psf(*xy[j]) for j in freej]
    if nbscale:
        extra.append(sum(g.F.values[j] * psf(*xy[j]) for j in trio))
    def solve(p):
        s = p[2] if blur else 0.0
        A = np.array([blurred(psf, mg[0] + p[0], mg[1] + p[1], s)] + extra + [np.ones_like(O)])
        w = 1 / np.sqrt(np.abs(O[mask]) + (0.03 * O[mask]) ** 2)
        c, *_ = np.linalg.lstsq(A[:, mask].T * w[:, None], (O - fixed)[mask] * w, rcond=None)
        return c, (A[:, mask].T @ c - (O - fixed)[mask]) * w, A
    best = None
    for gx in np.linspace(-0.6, 0.6, 7):
        for gy in np.linspace(-0.6, 0.6, 7):
            p0 = [gx, gy] + ([0.3] if blur else [])
            c2 = (solve(p0)[1] ** 2).sum()
            if best is None or c2 < best[0]: best = (c2, p0)
    lb = [-1.5, -1.5] + ([0.0] if blur else []); ub = [1.5, 1.5] + ([1.5] if blur else [])
    r = least_squares(lambda p: solve(p)[1], best[1], diff_step=1e-3, x_scale=0.05, bounds=(lb, ub))
    c, res, A = solve(r.x)
    dof = mask.sum() - len(c) - len(r.x)
    cov = np.linalg.pinv(r.jac.T @ r.jac) * max(1.0, (res ** 2).sum() / dof)
    dEN = S.J @ r.x[:2]
    eEN = np.sqrt(np.diag(S.J @ cov[:2, :2] @ S.J.T))
    mod = np.tensordot(c, A, 1) + fixed
    rms = float(np.sqrt(np.mean(((O - mod)[mask] / O[mask]) ** 2)))
    return dict(dEN=dEN, err=eEN, rad=float(dEN @ S.u), tan=float(dEN @ np.array([S.u[1], -S.u[0]])), flux=float(c[0]),
                blur=float(r.x[2]) if blur else 0.0, chi2r=float((res ** 2).sum() / dof), resid_rms=rms,
                trio_scale=float(c[1 + len(freej)]) if nbscale else 1.0)


def centroid_shift(S, O, psf, cmodel='curve'):
    """Model-light check: flux-weighted centroid within 1.5 px of the Mira's Gaia position of (data - other stars - bkg)
    minus the same centroid of a pure PRF at the Mira's Gaia position.  Returns (dE,dN) arcsec."""
    xy = S.xy(g, cmodel); xyn = S.xy(g, use_colour=False)
    jm = int(np.where(g.source_id == ID['Mira'])[0][0]); mg = xyn[jm]
    yy, xx = np.mgrid[0:S.ny, 0:S.nx]
    others = sum(g.F.values[j] * psf(*xy[j]) for j in range(len(g))
                 if j != jm and g['T'].values[j] < 17 and np.hypot(*(xy[j] - mg)) < 8)
    ring = (np.hypot(xx - mg[0], yy - mg[1]) > 6) & (np.hypot(xx - mg[0], yy - mg[1]) < 9)
    bk = np.median((O - others)[ring])
    E = O - others - bk
    ap = np.hypot(xx - mg[0], yy - mg[1]) <= 1.5
    def cen(I):
        return np.array([(I * xx)[ap].sum() / I[ap].sum(), (I * yy)[ap].sum() / I[ap].sum()])
    model = psf(*mg)
    return S.J @ (cen(E) - cen(model))


SECL = [int(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else SECTORS
for sec in SECL:
    S = Sec(sec)
    psf = S.prf()
    O = np.median(S.flux[S.good], 0)
    res = {}
    for cm, RM in (('curve', 2.0), ('curve', 2.5), ('curve', 3.5), ('linear', 2.5), ('none', 2.5)):
        res[f'{cm}_R{RM}_F0.03'] = global_fit(S, O, psf, cm, RM, 0.03)
    res['curve_far'] = global_fit(S, O, psf, 'curve', 2.5, 0.03, far=True)
    res['curve_R2.5_F0.05'] = global_fit(S, O, psf, 'curve', 2.5, 0.05)
    base = res['curve_R2.5_F0.03']
    S.off = base['off']
    offs = np.array([v['off'] for v in res.values()])
    mf = {}
    for R in (2.0, 3.0, 4.0):
        for bl in (False, True):
            mf[f'R{R}_blur{int(bl)}'] = mira_fit(S, O, psf, R, bl)
    mf['R3.0_blur0_trioscale'] = mira_fit(S, O, psf, 3.0, False, nbscale=True)
    mf['R3.0_blur0_linear'] = mira_fit(S, O, psf, 3.0, False, cmodel='linear')
    cen = centroid_shift(S, O, psf)
    OUT[str(sec)] = dict(offset_px=base['off'], offset_err_px=base['err'], offset_variant_range_px=[offs.min(0), offs.max(0)],
                         global_variants=res, mira=mf, mira_centroid_dEN=cen, u=S.u, theta=S.theta,
                         colour_shift_pred_mira_curve=float(colour_shift(4.6837, 'curve')),
                         colour_shift_pred_m3v_curve=float(colour_shift(2.7561, 'curve')),
                         colour_shift_pred_m3v_linear=float(colour_shift(2.7561, 'linear')))
    b = mf['R3.0_blur0']; bb = mf['R3.0_blur1']
    print(f"S{sec}: WCS offset ({base['off'][0]:+.3f},{base['off'][1]:+.3f}) +- ({base['err'][0]:.3f},{base['err'][1]:.3f}) px, "
          f"variants x[{offs[:,0].min():+.3f},{offs[:,0].max():+.3f}] y[{offs[:,1].min():+.3f},{offs[:,1].max():+.3f}] chi2r {base['chi2r']:.1f} nstar {base['nstar']}")
    print(f"   Mira R3: dEN ({b['dEN'][0]:+.1f},{b['dEN'][1]:+.1f})+-({b['err'][0]:.1f},{b['err'][1]:.1f}) rad {b['rad']:+.1f} tan {b['tan']:+.1f} flux {b['flux']:.0f} rms {b['resid_rms']*100:.1f}%")
    print(f"   Mira R3 blur: dEN ({bb['dEN'][0]:+.1f},{bb['dEN'][1]:+.1f}) rad {bb['rad']:+.1f} tan {bb['tan']:+.1f} blur {bb['blur']:.2f}px rms {bb['resid_rms']*100:.1f}%; "
          f"centroid ({cen[0]:+.1f},{cen[1]:+.1f}) rad {cen @ S.u:+.1f}")
    for k, v in mf.items():
        print(f"     {k}: dEN ({v['dEN'][0]:+.1f},{v['dEN'][1]:+.1f}) rad {v['rad']:+.1f} tan {v['tan']:+.1f} flux {v['flux']:.0f} chi2r {v['chi2r']:.1f} blur {v['blur']:.2f} trio {v['trio_scale']:.2f}")
    for k, v in res.items():
        print(f"     WCS {k}: ({v['off'][0]:+.3f},{v['off'][1]:+.3f}) chi2r {v['chi2r']:.2f}")
jdump(OUT, f'{ROOT}/out/wcs_mira' + ('_' + sys.argv[1].replace(',', '_') if len(sys.argv) > 1 else '') + '.json')
