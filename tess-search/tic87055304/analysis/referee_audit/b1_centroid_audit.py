"""TASK B: audit of src/centroid.py (copied verbatim to centroid_copy.py, md5 1053541e...).

1. Replicate the vet run (vet/TIC87055304/report.json 'localisation'): same ephemeris and duration
   (P = 4.136190053729093, t0 = 1658.411589268281, dur = 0.12491872185913565 = the BLS duration, NOT the 0.1875 d T14),
   TIC J2000 position, 15x15 stamp (= central crop [3:18, 3:18] of the 21x21 TESScut cutouts; target pixel position
   agrees with report.json to < 0.001 px), QUALITY == 0, nboot = 300, 0.15 px floor.
2. Diagnose: per-window D floor (far-field median of each difference image), estimator variants (no clipping,
   floor-subtracted, both), and what the estimator returns for noiseless model difference images.
3. Injection: real cube at fake epochs (phases 0.25..0.75 step 0.05), transit (T14 0.1875 d, T12 0.045 d, flat bottom)
   of the real total amplitude injected as a PRF (and a Gaussian sigma 0.6 px) at the target or at the neighbour;
   centroid.py core logic run on the fake ephemeris exactly as in the vet.
   Reported: offset_arcsec, offset_sigma (as the code computes them), the projection on the target->neighbour axis
   (s_par = 0 target, 1 neighbour), and the 'lead's statistic' |c - x_N| / err_px (no floor) used for '>6 sigma'.
"""
import sys
import copy
import numpy as np
from scipy.special import erf
from astropy.wcs import WCS
sys.path.insert(0, '/home/user/work/tic87/referee_audit')
import centroid_copy as cc
from rcommon import Cube, SECTORS, OUT, dump, TIC_RA, TIC_DEC

P_V, T0_V, DUR_V = 4.136190053729093, 1658.411589268281, 0.12491872185913565
CROP = slice(3, 18)
PHASES = np.round(np.arange(0.25, 0.751, 0.05), 3)
T14, T12 = 0.1875, 0.045


class _H:
    def __init__(self, data=None, header=None):
        self.data, self.header = data, header


def make_hdul(c, flux=None):
    f = c.f if flux is None else flux
    hdr = copy.deepcopy(c.h2)
    hdr['CRPIX1'] -= 3; hdr['CRPIX2'] -= 3
    return [None, _H(dict(TIME=c.t, FLUX=f[:, CROP, CROP], QUALITY=c.q)), _H(header=hdr)]


def localise_h(h, ra, dec, P, t0, dur, nboot=300):
    """Body of centroid.localise() after the cutout download, verbatim logic."""
    r = cc.diff_image(h, P, t0, dur)
    if r is None:
        return None
    diffs, outs, w = r
    tx, ty = w.world_to_pixel_values(ra, dec)
    D = diffs.mean(0); O = outs.mean(0)
    cx, cy = cc.centroid(D, tx, ty); ox, oy = cc.centroid(O - np.median(O), tx, ty)
    rng = np.random.default_rng(0); bx, by = [], []
    for _ in range(nboot):
        k = rng.integers(0, len(diffs), len(diffs))
        x, y = cc.centroid(diffs[k].mean(0), tx, ty); bx.append(x); by.append(y)
    ex, ey = np.nanstd(bx), np.nanstd(by)
    dx, dy = cx - tx, cy - ty
    sx, sy = np.hypot(ex, 0.15), np.hypot(ey, 0.15)
    sig = np.hypot(dx / sx, dy / sy)
    return dict(n_transits=len(diffs), target_px=(float(tx), float(ty)), diff_centroid=(float(cx), float(cy)),
                offset_px=(float(dx), float(dy)), offset_arcsec=float(np.hypot(dx, dy) * 21.0),
                err_px=(float(ex), float(ey)), offset_sigma=float(sig), D=D, diffs=diffs, w=w)


def trap(dt):
    a = np.abs(dt); h = T14 / 2
    return np.clip((h - a) / T12, 0, 1)


def gauss_psf(shape, x, y, s):
    xe = np.arange(shape[1] + 1) - 0.5; ye = np.arange(shape[0] + 1) - 0.5
    cx = 0.5 * np.diff(erf((xe - x) / (np.sqrt(2) * s))); cy = 0.5 * np.diff(erf((ye - y) / (np.sqrt(2) * s)))
    return np.outer(cy, cx)


def summarize(L, c, xn15, yn15):
    tx, ty = L['target_px']; cx, cy = L['diff_centroid']; ex, ey = L['err_px']
    u = np.array([xn15 - tx, yn15 - ty]); sep = np.hypot(*u); u /= sep
    off = np.array([cx - tx, cy - ty])
    s_par = float(off @ u / sep)
    z_lead = float(np.hypot((cx - xn15) / ex, (cy - yn15) / ey))           # err_px only (no floor)
    z_floor = float(np.hypot((cx - xn15) / np.hypot(ex, 0.15), (cy - yn15) / np.hypot(ey, 0.15)))
    return dict(offset_px=L['offset_px'], offset_arcsec=L['offset_arcsec'], err_px=L['err_px'],
                offset_sigma=L['offset_sigma'], n_transits=L['n_transits'], s_par=s_par,
                z_vs_neighbour_errpx_only=z_lead, z_vs_neighbour_with_floor=z_floor,
                dist_to_neighbour_px=float(np.hypot(cx - xn15, cy - yn15)), sep_px=float(sep))


def main():
    import PRF
    res = {}
    for s in SECTORS:
        c = Cube(s)
        xn15, yn15 = c.xn - 3, c.yn - 3
        h = make_hdul(c)
        L = localise_h(h, TIC_RA, TIC_DEC, P_V, T0_V, DUR_V)
        r = dict(real=None)
        if L is not None:
            r['real'] = summarize(L, c, xn15, yn15)
            D = L['D']; diffs = L['diffs']
            tx, ty = L['target_px']
            yy, xx = np.mgrid[0:D.shape[0], 0:D.shape[1]]
            far = np.hypot(xx - tx, yy - ty) > 5
            floors = [float(np.median(d[far])) for d in diffs]
            r['real']['per_window_floor_e_s_px'] = floors
            r['real']['stack_floor_e_s_px'] = float(np.median(D[far]))
            ap = np.hypot(xx - tx, yy - ty) <= 3
            r['real']['n_ap_pix'] = int(ap.sum())
            r['real']['sum_D_aperture'] = float(D[ap].sum())
            r['real']['sum_clipped_D_aperture'] = float(np.clip(D, 0, None)[ap].sum())
            r['real']['n_negative_pix_in_aperture'] = int((D[ap] < 0).sum())
            # estimator variants on the same stacked D
            def cen(img, clip=True):
                m = ap
                w = (np.clip(img, 0, None) if clip else img) * m
                return (w * xx).sum() / w.sum(), (w * yy).sum() / w.sum()
            u = np.array([xn15 - tx, yn15 - ty]); sep = np.hypot(*u)
            var = {}
            for nm, img, clip in (('as_coded', D, True), ('no_clip', D, False),
                                  ('floor_subtracted', D - np.median(D[far]), True),
                                  ('floor_subtracted_no_clip', D - np.median(D[far]), False),
                                  ('per_window_floor_subtracted', np.mean([d - np.median(d[far]) for d in diffs], 0), True)):
                x_, y_ = cen(img, clip)
                var[nm] = dict(offset_px=[float(x_ - tx), float(y_ - ty)],
                               s_par=float(np.array([x_ - tx, y_ - ty]) @ u / sep ** 2))
            r['real']['estimator_variants'] = var
            # noiseless model difference images with the real amplitude
            A = max(float(D[ap].sum() - np.median(D[far]) * ap.sum()), 5.0)
            prf = PRF.TESS_PRF(c.cam, c.ccd, c.s, c.col0 + c.nx / 2, c.row0 + c.ny / 2, localdatadir=f'{OUT}/prf_files/')
            MT = prf.locate(c.xt, c.yt, (c.ny, c.nx))[CROP, CROP]; MN = prf.locate(c.xn, c.yn, (c.ny, c.nx))[CROP, CROP]
            nl = {}
            for nm, img in (('target_noiseless', A * MT), ('neighbour_noiseless', A * MN),
                            ('neighbour_plus_0.3e_s_pedestal', A * MN + 0.3), ('neighbour_plus_1e_s_pedestal', A * MN + 1.0)):
                x_, y_ = cen(img)
                nl[nm] = dict(offset_px=[float(x_ - tx), float(y_ - ty)],
                              s_par=float(np.array([x_ - tx, y_ - ty]) @ u / sep ** 2))
            r['real']['noiseless_models'] = nl
            r['real']['A_used_e_s'] = A
            print(f"S{s} REAL: n={L['n_transits']} offset {np.round(L['offset_px'],3)} px = {L['offset_arcsec']:.2f}\" "
                  f"err_px {np.round(L['err_px'],3)} sigma {L['offset_sigma']:.2f} s_par {r['real']['s_par']:.2f} "
                  f"| z_vs_nb(err_px only) {r['real']['z_vs_neighbour_errpx_only']:.1f} with floor {r['real']['z_vs_neighbour_with_floor']:.1f}"
                  f" | floors {np.round(floors,2)} | variants s_par: " +
                  ' '.join(f"{k}={v['s_par']:.2f}" for k, v in var.items()) +
                  ' | noiseless s_par: ' + ' '.join(f"{k}={v['s_par']:.2f}" for k, v in nl.items()), flush=True)
        else:
            print(f"S{s} REAL: None")
        # injections
        if L is not None:
            A = r['real']['A_used_e_s']
        else:
            A = 13.0
        prf = PRF.TESS_PRF(c.cam, c.ccd, c.s, c.col0 + c.nx / 2, c.row0 + c.ny / 2, localdatadir=f'{OUT}/prf_files/')
        models = {'prf': (prf.locate(c.xt, c.yt, (c.ny, c.nx)), prf.locate(c.xn, c.yn, (c.ny, c.nx))),
                  'gauss0.6': (gauss_psf((c.ny, c.nx), c.xt, c.yt, 0.6), gauss_psf((c.ny, c.nx), c.xn, c.yn, 0.6))}
        inj = {}
        for mname, (MT, MN) in models.items():
            for host, M in (('target', MT), ('neighbour', MN)):
                out = []
                for ph in PHASES:
                    t0f = T0_V + ph * P_V
                    ep = np.round((c.t - t0f) / P_V); dt = c.t - t0f - ep * P_V
                    flux = c.f - A * trap(dt)[:, None, None] * M[None]
                    Lf = localise_h(make_hdul(c, flux), TIC_RA, TIC_DEC, P_V, t0f, DUR_V)
                    if Lf is None or not np.isfinite(Lf['offset_sigma']):
                        out.append(None); continue
                    out.append(dict(phase=float(ph), **summarize(Lf, c, xn15, yn15)))
                inj[f'{mname}_{host}'] = out
                ok = [o for o in out if o]
                arr = lambda k: np.array([o[k] for o in ok])
                print(f"   S{s} inject {mname:8s} on {host:9s}: n={len(ok)}/{len(out)} offset\" med {np.median(arr('offset_arcsec')):.2f} "
                      f"[{np.min(arr('offset_arcsec')):.2f},{np.max(arr('offset_arcsec')):.2f}] sigma med {np.median(arr('offset_sigma')):.2f} "
                      f"max {np.max(arr('offset_sigma')):.2f}; s_par med {np.median(arr('s_par')):.2f} "
                      f"[{np.percentile(arr('s_par'),16):.2f},{np.percentile(arr('s_par'),84):.2f}]; "
                      f"z_vs_nb(err_px only) med {np.median(arr('z_vs_neighbour_errpx_only')):.1f} "
                      f"frac>6: {np.mean(arr('z_vs_neighbour_errpx_only')>6):.2f}", flush=True)
        r['injections'] = inj
        r['A_injected_e_s'] = A
        res[str(s)] = r
    dump(res, f'{OUT}/b1_centroid_audit.json')


if __name__ == '__main__':
    main()
