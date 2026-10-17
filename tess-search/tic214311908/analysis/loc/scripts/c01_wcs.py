"""WCS calibration on normal-colour field stars, and the Mira's apparent (colour-shifted) position.

Per sector, on the median direct image:
 (a) global fit: every Gaia star with T < 15.5 in/near the cutout as a TESS PRF at its WCS position + lateral-colour
     shift (l0lib.colour_shift, field-angle scaled) + a common offset d (free); the Mira's position free on its own;
     fluxes free (T < 15) or fixed to catalogue (15-15.5); background plane; pixels within 1.6 px of the Mira masked;
     per-pixel sigma^2 = O + (0.02*model)^2 (PRF-mismatch floor).
 (b) individual fits of isolated bright calibrators (T < 13.3, no star with dT < 2 within 2.5 px, > 3 px from the
     Mira, > 2 px from the edges): free position + flux + local constant, neighbours fixed from (a).
     Their offsets (fit - predicted, without d) give d and its empirical scatter; red calibrators test the colour curve.
Output: out/wcs_cal.json  {sector: offset_px, err_px, mira_offset_arcsec, ...}
"""
import sys
sys.path.insert(0, 'tic214311908/analysis/loc/scripts')
from l0lib import *
from scipy.optimize import least_squares

g = stars()
out = {}
if os.path.exists(f'{BASE}/out/wcs_cal.json'):
    os.remove(f'{BASE}/out/wcs_cal.json')
for sec in SECTORS:
    S = Sector(sec)
    psf = S.prf()
    O = np.median(S.flux[S.good], 0)
    xy = S.star_xy(g, colour=True)
    ins = (xy[:, 0] > -2) & (xy[:, 0] < S.nx + 1) & (xy[:, 1] > -2) & (xy[:, 1] < S.ny + 1) & (g['T'].values < 15.5)
    idx = np.where(ins)[0]
    im = int(np.where(g.source_id.values[idx] == IDS['Mira'])[0][0])
    free = g['T'].values[idx] < 15.0
    yy, xx = np.mgrid[0:S.ny, 0:S.nx]
    mxy = xy[idx[im]]
    mask = np.hypot(xx - mxy[0], yy - mxy[1]) > 1.6
    def model_parts(p):
        d = p[:2]; dm = p[2:4]
        cols, fixed = [], np.zeros_like(O)
        for k, j in enumerate(idx):
            pos = xy[j] + (dm if k == im else d)
            img = psf(pos[0], pos[1])
            if free[k]: cols.append(img)
            else: fixed += g.F.values[j] * img
        cols += [np.ones_like(O), xx - 10.0, yy - 10.0]
        return np.array(cols), fixed
    def solve(p, wts=None):
        A, fixed = model_parts(p)
        m = mask
        M = A[:, m].T; b = (O - fixed)[m]
        w = 1 / np.sqrt(np.abs(O[m]) + (0.02 * O[m]) ** 2) if wts is None else wts
        c, *_ = np.linalg.lstsq(M * w[:, None], b * w, rcond=None)
        mod = (A.T @ c).T + fixed if False else np.tensordot(c, A, axes=1) + fixed
        return c, mod, w
    def res(p):
        c, mod, w = solve(p)
        return (O - mod)[mask] * w
    r = least_squares(res, np.zeros(4), diff_step=1e-3, x_scale=0.05)
    c, mod, w = solve(r.x)
    chi2r = float((r.fun ** 2).sum() / (mask.sum() - len(c) - 4))
    Jc = r.jac; cov = np.linalg.pinv(Jc.T @ Jc) * max(chi2r, 1)
    d_glob = r.x[:2]; dm = r.x[2:4]
    # individual calibrators
    fl = dict(zip(idx[free], c[:free.sum()]))
    cal = []
    for k, j in enumerate(idx):
        if k == im or g['T'].values[j] > 13.3: continue
        x0, y0 = xy[j]
        if not (2 <= x0 <= S.nx - 3 and 2 <= y0 <= S.ny - 3): continue
        if np.hypot(x0 - mxy[0], y0 - mxy[1]) < 3.0: continue
        dd = np.hypot(xy[idx, 0] - x0, xy[idx, 1] - y0)
        if np.any((dd < 2.5) & (dd > 0.01) & (g['T'].values[idx] < g['T'].values[j] + 2)): continue
        # neighbours (all other stars) fixed from the global fit
        nb = np.zeros_like(O)
        for kk, jj in enumerate(idx):
            if jj == j: continue
            pos = xy[jj] + (dm if kk == im else d_glob)
            nb += (fl.get(jj, g.F.values[jj])) * psf(pos[0], pos[1])
        mm = np.hypot(xx - x0, yy - y0) <= 2.6
        def r1(p):
            img = psf(x0 + p[0], y0 + p[1])[mm]
            A = np.column_stack([img, np.ones(mm.sum())])
            b = (O - nb)[mm]; ww = 1 / np.sqrt(np.abs(O[mm]) + (0.02 * O[mm]) ** 2)
            cc, *_ = np.linalg.lstsq(A * ww[:, None], b * ww, rcond=None)
            return (A @ cc - b) * ww
        rr = least_squares(r1, d_glob.copy(), diff_step=1e-3, x_scale=0.05)
        cal.append(dict(source_id=int(g.source_id.values[j]), T=float(g['T'].values[j]), C=float(g.C.values[j]),
                        x=float(x0), y=float(y0), dx=float(rr.x[0]), dy=float(rr.x[1]), chi2=float((rr.fun ** 2).sum()), n=int(mm.sum())))
    cal = pd.DataFrame(cal)
    norm = cal[(cal.C > 0.5) & (cal.C < 1.8)]
    dmed = np.array([np.median(norm.dx), np.median(norm.dy)])
    dsd = np.array([1.4826 * np.median(np.abs(norm.dx - dmed[0])), 1.4826 * np.median(np.abs(norm.dy - dmed[1]))])
    # radial residual (arcsec) of each calibrator after removing the global offset: tests the colour curve
    resEN = (S.J @ (cal[['dx', 'dy']].values - d_glob[None, :]).T).T
    cal['rad_res'] = resEN @ S.u_rad
    cal['tan_res'] = resEN @ np.array([S.u_rad[1], -S.u_rad[0]])
    mEN = S.J @ (dm - d_glob)   # Mira apparent offset (arcsec) relative to its Gaia + colour-curve prediction
    # Mira offset relative to its plain Gaia position (no colour shift): add back the assumed colour shift
    csh = colour_shift(g.C.values[idx[im]], S.colscale) * S.u_rad
    mEN_raw = mEN + csh
    out[str(sec)] = dict(offset_px=d_glob.tolist(), offset_err_px=np.sqrt(np.diag(cov)[:2]).tolist(),
                         offset_px_calmedian=dmed.tolist(), cal_scatter_px=dsd.tolist(), n_cal=int(len(norm)),
                         offset_err_px_empirical=(dsd / np.sqrt(max(len(norm), 1))).tolist(),
                         chi2r_global=chi2r, mira_offset_vs_gaia_arcsec=mEN_raw.tolist(),
                         mira_offset_radial_arcsec=float(mEN_raw @ S.u_rad),
                         mira_offset_tangential_arcsec=float(mEN_raw @ np.array([S.u_rad[1], -S.u_rad[0]])),
                         mira_flux=float(c[list(idx[free]).index(idx[im])]),
                         calibrators=cal.round(4).to_dict('records'))
    print(f'S{sec}: global d=({d_glob[0]:+.3f},{d_glob[1]:+.3f}) px +- ({np.sqrt(cov[0,0]):.3f},{np.sqrt(cov[1,1]):.3f}); '
          f'calibrator median ({dmed[0]:+.3f},{dmed[1]:+.3f}) scatter ({dsd[0]:.3f},{dsd[1]:.3f}) n={len(norm)}; chi2r {chi2r:.1f}')
    print(f'       Mira apparent offset from Gaia (E,N) = ({mEN_raw[0]:+.1f},{mEN_raw[1]:+.1f})"  radial {mEN_raw @ S.u_rad:+.1f}" '
          f'(neg = away from boresight), tangential {mEN_raw @ np.array([S.u_rad[1], -S.u_rad[0]]):+.1f}"; Mira flux {out[str(sec)]["mira_flux"]:.0f} e/s')
    red = cal[cal.C > 1.8]
    for _, rw in red.iterrows():
        print(f'       red calibrator {rw.source_id} C={rw.C:.2f} T={rw["T"]:.1f}: radial residual {rw.rad_res:+.2f}" tangential {rw.tan_res:+.2f}" '
              f'(colour curve already applied: predicted shift {colour_shift(rw.C, S.colscale):+.2f}")')
    S.wcs_off = d_glob
to_json(out, f'{BASE}/out/wcs_cal.json')
