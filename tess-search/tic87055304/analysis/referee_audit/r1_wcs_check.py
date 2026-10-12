"""Check 1: epoch-propagated Gaia positions versus the TESScut WCS, orientation/parity in every sector.

(a) WCS from ext 2 (APERTURE) vs the FLUX-column WCS in ext 1: positions of target and neighbour.
(b) Orientation: pixel-scale matrix, rotation, parity (sign of det), PA of the target->neighbour vector recovered
    by inverting the pixel vector through the WCS.
(c) Independent geometry from tess-point (spacecraft pointing model, no TESScut WCS): CCD col/row of both stars;
    the pixel vector target->neighbour must equal the WCS one (cutout x = CCD col - 1CRV4P, y = row - 2CRV4P).
(d) Field-star test with NO PSF model: flux-weighted centroids (3x3, local-background subtracted) of bright stars in
    the median out-of-transit image, compared with WCS-predicted positions under the actual WCS and under mirrored /
    transposed / rotated alternatives (mirror about the target).
"""
import numpy as np
from astropy.wcs import WCS
import tess_stars2px as tp
from rcommon import Cube, SECTORS, OUT, dump, TID, NID

res = {}
for s in SECTORS:
    c = Cube(s)
    r = {}
    # (a) ext-1 column WCS
    hdr = c.h1
    wc = WCS(naxis=2)
    wc.wcs.ctype = [hdr['1CTYP5'], hdr['2CTYP5']]
    wc.wcs.crpix = [hdr['1CRPX5'], hdr['2CRPX5']]
    wc.wcs.crval = [hdr['1CRVL5'], hdr['2CRVL5']]
    wc.wcs.cdelt = [hdr['1CDLT5'], hdr['2CDLT5']]
    wc.wcs.pc = [[hdr['11PC5'], hdr['12PC5']], [hdr['21PC5'], hdr['22PC5']]]
    g = c.g
    for nm, i in (('target', c.it), ('neighbour', c.ineb)):
        xa, ya = c.w.all_world2pix(g.ra_ep[i], g.dec_ep[i], 0)
        xb, yb = wc.all_world2pix(g.ra_ep[i], g.dec_ep[i], 0)
        r[f'{nm}_px_ext2'] = [float(xa), float(ya)]
        r[f'{nm}_px_ext1col'] = [float(xb), float(yb)]
    # (b) orientation
    M = c.w.pixel_scale_matrix * 3600.0   # d(E-ish intermediate, N)/d(x,y), arcsec/px
    r['pixel_scale_matrix_arcsec'] = M
    r['det_sign'] = float(np.sign(np.linalg.det(M)))
    # sky direction of +x and +y (E, N) and their PA (east of north)
    for k, nm in ((0, '+x'), (1, '+y')):
        e, n = M[0, k], M[1, k]
        r[f'PA_of_{nm}_deg'] = float(np.degrees(np.arctan2(e, n)) % 360)
    # recover PA of target->neighbour from pixel vector through WCS
    ra1, de1 = c.w.all_pix2world(c.xn, c.yn, 0)
    ra0, de0 = c.w.all_pix2world(c.xt, c.yt, 0)
    dE = (ra1 - ra0) * np.cos(np.radians(de0)) * 3600; dN = (de1 - de0) * 3600
    r['nb_offset_EN_arcsec_via_wcs'] = [float(dE), float(dN)]
    r['nb_PA_deg'] = float(np.degrees(np.arctan2(dE, dN)) % 360)
    r['nb_sep_arcsec'] = float(np.hypot(dE, dN))
    r['nb_vec_px'] = [c.xn - c.xt, c.yn - c.yt]
    r['nb_sep_px'] = c.sep_px
    r['epoch_jyear'] = c.epoch
    # Gaia-only offset at the epoch
    dEg = (g.ra_ep[c.ineb] - g.ra_ep[c.it]) * np.cos(np.radians(g.dec_ep[c.it])) * 3600
    dNg = (g.dec_ep[c.ineb] - g.dec_ep[c.it]) * 3600
    r['nb_offset_EN_arcsec_gaia_epoch'] = [float(dEg), float(dNg)]
    # (c) tess-point
    ids = [TID, NID]
    out = tp.tess_stars2px_function_entry(ids, [g.ra_ep[c.it], g.ra_ep[c.ineb]], [g.dec_ep[c.it], g.dec_ep[c.ineb]],
                                          trySector=s)
    sec, cam, ccd, col, row = out[3], out[4], out[5], out[6], out[7]
    r['tess_point'] = dict(cam=cam.tolist(), ccd=ccd.tolist(), col=col.tolist(), row=row.tolist())
    tp_vec = [float(col[1] - col[0]), float(row[1] - row[0])]
    r['tess_point_nb_vec_px'] = tp_vec
    r['tess_point_minus_wcs_vec_px'] = [tp_vec[0] - r['nb_vec_px'][0], tp_vec[1] - r['nb_vec_px'][1]]
    r['tess_point_target_minus_cutout_px'] = [float(col[0] - (c.col0 + c.xt)), float(row[0] - (c.row0 + c.yt))]
    r['header_cam_ccd'] = [c.cam, c.ccd]
    # (d) field stars, PSF-free centroids
    oot = c.good.copy()
    med = np.median(c.f[oot], axis=0)
    yy, xx = np.mgrid[0:c.ny, 0:c.nx]
    sel = (g.tmag < 15.0) & (g.x > 2) & (g.x < c.nx - 3) & (g.y > 2) & (g.y < c.ny - 3)
    sel &= ~g.source_id.isin([TID, NID])
    # isolation: no star brighter than (T+1.0) within 2.0 px
    stars = []
    for i in np.where(sel)[0]:
        d = np.hypot(g.x - g.x[i], g.y - g.y[i])
        if np.any((d < 2.0) & (d > 0) & (g.tmag < g.tmag[i] + 1.0)):
            continue
        # also skip if within 3 px of the target pair
        if np.hypot(g.x[i] - c.xt, g.y[i] - c.yt) < 3.0:
            continue
        stars.append(i)
    def centroid_at(xp, yp):
        ix, iy = int(round(xp)), int(round(yp))
        # brightest pixel within 1 px of the predicted position
        best = None
        for ddx in (-1, 0, 1):
            for ddy in (-1, 0, 1):
                jx, jy = ix + ddx, iy + ddy
                if 1 <= jx < c.nx - 1 and 1 <= jy < c.ny - 1:
                    if best is None or med[jy, jx] > best[0]:
                        best = (med[jy, jx], jx, jy)
        _, jx, jy = best
        ring = (np.abs(xx - jx) <= 3) & (np.abs(yy - jy) <= 3) & ~((np.abs(xx - jx) <= 1) & (np.abs(yy - jy) <= 1))
        b = np.median(med[ring])
        box = (np.abs(xx - jx) <= 1) & (np.abs(yy - jy) <= 1)
        w = np.clip(med - b, 0, None) * box
        return (w * xx).sum() / w.sum(), (w * yy).sum() / w.sum()
    meas = np.array([centroid_at(g.x[i], g.y[i]) for i in stars])
    pred = np.array([[g.x[i], g.y[i]] for i in stars])
    # alternatives: mirror about the target in x, in y; transpose about the target; 180-deg rotation about the target
    xt, yt = c.xt, c.yt
    alts = {
        'actual_wcs': pred,
        'mirror_x_about_target': np.column_stack([2 * xt - pred[:, 0], pred[:, 1]]),
        'mirror_y_about_target': np.column_stack([pred[:, 0], 2 * yt - pred[:, 1]]),
        'transpose_about_target': np.column_stack([xt + (pred[:, 1] - yt), yt + (pred[:, 0] - xt)]),
        'rot180_about_target': np.column_stack([2 * xt - pred[:, 0], 2 * yt - pred[:, 1]]),
        'array_transposed_[x,y]_swap': np.column_stack([pred[:, 1], pred[:, 0]]),
    }
    alt_res = {}
    for k, p in alts.items():
        # a PSF-free centroid in a 3x3 box is biased toward the box centre; compare with the same estimator applied
        # at the alternative position (i.e. re-measure there), and report the median distance
        mm = np.array([centroid_at(x, y) if (1 <= x < c.nx - 2 and 1 <= y < c.ny - 2) else (np.nan, np.nan)
                       for x, y in p])
        dd = np.hypot(mm[:, 0] - p[:, 0], mm[:, 1] - p[:, 1])
        # flux found at the alternative position (relative to the flux at the actual position)
        fl = np.array([med[int(round(y)), int(round(x))] if (0 <= round(x) < c.nx and 0 <= round(y) < c.ny) else np.nan
                       for x, y in p])
        alt_res[k] = dict(median_dist_px=float(np.nanmedian(dd)), rms_dist_px=float(np.sqrt(np.nanmean(dd**2))),
                          median_peak_flux=float(np.nanmedian(fl)), n=int(np.isfinite(dd).sum()))
    r['field_star_test'] = alt_res
    r['field_star_n'] = len(stars)
    r['field_star_resid_px_actual'] = (meas - pred).tolist()
    # background level for reference
    r['median_bkg_floor_e_s_px'] = float(np.percentile(med, 10))
    res[str(s)] = r
    print(s, 'cam/ccd', c.cam, c.ccd, 'det', r['det_sign'], 'PA(+x) %.1f PA(+y) %.1f' % (r['PA_of_+x_deg'], r['PA_of_+y_deg']),
          'nb vec px', np.round(r['nb_vec_px'], 3), 'PA %.1f sep %.2f"' % (r['nb_PA_deg'], r['nb_sep_arcsec']),
          '| tess-point vec', np.round(tp_vec, 3), 'diff', np.round(r['tess_point_minus_wcs_vec_px'], 4),
          '| ext1-ext2 target', np.round(np.array(r['target_px_ext1col']) - np.array(r['target_px_ext2']), 6))
    for k, v in alt_res.items():
        print('    %-30s median dist %.3f px  rms %.3f  peak flux %.0f  (n=%d)' % (k, v['median_dist_px'], v['rms_dist_px'],
                                                                              v['median_peak_flux'], v['n']))
dump(res, f'{OUT}/r1_wcs_check.json')
for s in SECTORS:
    r = res[str(s)]
    d = np.array(r['field_star_resid_px_actual'])
    print(s, 'n=%d' % len(d), 'mean resid (dx,dy) px', np.round(d.mean(0), 3), 'rms', np.round(d.std(0), 3),
          '| tess-point target minus cutout px', np.round(r['tess_point_target_minus_cutout_px'], 3))
