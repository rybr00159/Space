"""Check 2: is 'registration on field stars' independent of the target/neighbour pair, and does it matter?

(a) Fisher information for a global (x, y) registration shift, PRF scene of all Gaia stars (T < 18) with T-mag fluxes,
    shot noise from scene + background. Fraction of the information carried by the pair (target + neighbour) in
    the pixel sets the investigators used:
      loc_independent local fit: r <= 4 px around the target (pair included);
      loc_independent 'external': annulus 2.5-8 px;   psf_photometry: all pixels with r > 3 px.
(b) Registration error needed to turn a target host into the observed neighbour-like answer: 0.84*sep along u (the
    observed f_target ~ 0.16), vs the WCS-only PSF-free field-star offsets from r1_wcs_check.json projected on u.
"""
import json
import numpy as np
from rcommon import Cube, SECTORS, OUT, dump, TID, NID
from r3_pixel_ratio import prf_model

ZP = 0.86 * 15000.0
wcs = json.load(open(f'{OUT}/r1_wcs_check.json'))
res = {}
for s in SECTORS:
    c = Cube(s)
    prf = prf_model(c)
    g = c.g
    yy, xx = np.mgrid[0:c.ny, 0:c.nx]
    rr = np.hypot(xx - c.xt, yy - c.yt)
    ok = np.isfinite(g.x.values) & (g.x.values > -3) & (g.x.values < c.nx + 2) & (g.y.values > -3) & (g.y.values < c.ny + 2) \
        & (g.tmag.values < 18)
    scene = np.zeros((c.ny, c.nx)); gx_all = np.zeros_like(scene); gy_all = np.zeros_like(scene)
    gx_pair = np.zeros_like(scene); gy_pair = np.zeros_like(scene); pair_scene = np.zeros_like(scene)
    h = 0.02
    for j in np.where(ok)[0]:
        F = ZP * 10 ** (-0.4 * (g.tmag.values[j] - 10))
        x, y = g.x.values[j], g.y.values[j]
        M = prf(x, y)
        dMx = (prf(x + h, y) - prf(x - h, y)) / (2 * h); dMy = (prf(x, y + h) - prf(x, y - h)) / (2 * h)
        scene += F * M; gx_all += F * dMx; gy_all += F * dMy
        if g.source_id.values[j] in (TID, NID):
            gx_pair += F * dMx; gy_pair += F * dMy; pair_scene += F * M
    var = scene + 250.0          # shot-noise proxy (scene + sky floor), arbitrary units cancel in fractions
    out = {}
    for nm, m in (('loc_indep_local_r4', rr <= 4), ('loc_indep_external_2.5_8', (rr >= 2.5) & (rr <= 8)),
                  ('psfphot_field_r_gt_3', rr > 3)):
        I_all = np.array([[np.sum(a * b / var * m) for b in (gx_all, gy_all)] for a in (gx_all, gy_all)])
        I_pair = np.array([[np.sum(a * b / var * m) for b in (gx_pair, gy_pair)] for a in (gx_pair, gy_pair)])
        # information along u
        u = c.u
        out[nm] = dict(frac_info_pair_along_u=float((u @ I_pair @ u) / (u @ I_all @ u)),
                       frac_flux_pair=float(np.sum(pair_scene * m) / np.sum(scene * m)),
                       n_pix=int(m.sum()))
    off = np.array(wcs[str(s)]['field_star_resid_px_actual'])
    mean_off = off.mean(0); sem = off.std(0, ddof=1) / np.sqrt(len(off))
    needed = 0.84 * c.sep_px
    res[str(s)] = dict(info=out, sep_px=c.sep_px, registration_error_needed_px=needed,
                       wcs_only_fieldstar_offset_px=mean_off.tolist(), wcs_only_fieldstar_offset_sem_px=sem.tolist(),
                       wcs_only_offset_along_u_px=float(mean_off @ c.u))
    print(f"S{s}: pair info fraction along u: " + ', '.join(f"{k} {v['frac_info_pair_along_u']:.4f} (flux {v['frac_flux_pair']:.4f})" for k, v in out.items())
          + f" | needed registration error {needed:.3f} px; WCS-only field-star offset along u {mean_off @ c.u:+.3f} px", flush=True)
dump(res, f'{OUT}/r2_registration.json')
