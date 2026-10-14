"""Empirical test for TESS lateral colour: colour-dependent PRF-fit positions of isolated field stars.

For each full calibrated FFI (two per sector):
 * Gaia DR3 stars (G 8.5-14.5, propagated to the FFI epoch) within +-BOX px of the target on the same CCD;
 * isolation: no Gaia source with dG < 4.5 within 3.5 px, and none with dG < 2 within 6 px;
 * 9x9-px stamp fit with the TESS PRF (PRF package, interpolated to a 7x7 grid of CCD positions across the box),
   free (x, y, flux) + planar background;
 * residual = fitted position - WCS(Gaia) position, in pixels.
Then residuals are modelled (robust, iteratively clipped) as
   r = poly2(x, y)  [registration / WCS errors]  +  m * (G - 12)  [magnitude term]
       + k_rad * (C - C0) * u_rad  +  k_tan * (C - C0) * u_tan
 with C = BP-RP, C0 = 0.786 (the target), u_rad = unit vector (pixel frame) from the star towards the camera
 boresight (tess-point), u_tan perpendicular.  k_rad < 0 means redder stars appear displaced away from the boresight.
Output: results/x10_ffi_<file>.npz (per-star table) and printed coefficient summary.
"""
import sys, os, glob, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from scipy.optimize import least_squares
import PRF
import tess_stars2px as t2p
from loclib import L, R, to_json, btjd_to_jyear

BOX = int(os.environ.get('X10_BOX', 350))
RA0, DE0 = 330.887948, -42.835511
C0 = 0.786
ISO_R = float(os.environ.get('X10_ISO_R', 3.5)); ISO_DG = float(os.environ.get('X10_ISO_DG', 4.5))
TAG = os.environ.get('X10_TAG', '')
gb = Table.read(f'{R}/gaia/cone3.5deg_G15.fits').to_pandas()
gall = Table.read(f'{R}/gaia/cone3deg_all.fits').to_pandas() if os.path.exists(f'{R}/gaia/cone3deg_all.fits') else None


def boresight(sec, cam):
    res = t2p.tess_stars2px_function_entry(1, RA0, DE0, trySector=sec)
    sc = res[-1]; i = list(sc.sectors).index(sec)
    return sc.camRa[cam - 1][i], sc.camDec[cam - 1][i]


def run(fn):
    t0 = time.time()
    with fits.open(fn) as h:
        hd0 = h[0].header; hd1 = h[1].header
        img = np.array(h[1].data, float)
        w = WCS(hd1)
    sec = int(os.path.basename(fn).split('-s0')[1][:3]); cam = int(hd1['CAMERA']); ccd = int(hd1['CCD'])
    tmid = 0.5 * (hd1['TSTART'] + hd1['TSTOP'])
    ep = btjd_to_jyear(tmid)
    dt = ep - 2016.0
    g = gb.copy()
    g['ra_e'] = g.ra + np.nan_to_num(g.pmra) * dt / 3.6e6 / np.cos(np.radians(g.dec))
    g['de_e'] = g.dec + np.nan_to_num(g.pmdec) * dt / 3.6e6
    x, y = w.all_world2pix(g.ra_e.values, g.de_e.values, 0)
    g['x'] = x; g['y'] = y
    xt, yt = w.all_world2pix([RA0], [DE0], 0); xt, yt = float(xt[0]), float(yt[0])
    # boresight pixel position (WCS extrapolation; only its direction is used)
    bra, bde = boresight(sec, cam)
    xb, yb = w.all_world2pix([bra], [bde], 0); xb, yb = float(xb[0]), float(yb[0])
    ny, nx = img.shape
    sel = (np.abs(g.x - xt) < BOX) & (np.abs(g.y - yt) < BOX) & (g.x > 44 + 8) & (g.x < 2091 - 8) & (g.y > 8) & (g.y < 2047 - 8)
    cand = g[sel & (g.phot_g_mean_mag > 8.5) & (g.phot_g_mean_mag < 14.5) & np.isfinite(g.bp_rp)].copy()
    # isolation (use the all-magnitude cone if available for faint neighbours)
    nb = gall if gall is not None else gb
    nbx, nby = w.all_world2pix(nb.ra.values, nb.dec.values, 0)
    nbg = nb.phot_g_mean_mag.values
    keep = []
    for i, r in cand.iterrows():
        d = np.hypot(nbx - r.x, nby - r.y)
        dg = nbg - r.phot_g_mean_mag
        close = (d > 0.02) & (d < 6.0)
        bad = np.any(close & (d < ISO_R) & (dg < ISO_DG)) or np.any(close & (dg < 2.0) & (d < 6.0 if ISO_R > 3 else d < 5.0))
        keep.append(not bad)
    cand = cand[np.array(keep)]
    # PRF grid
    xs = np.linspace(xt - BOX, xt + BOX, 7); ys = np.linspace(yt - BOX, yt + BOX, 7)
    sub = 's0001' if sec < 4 else 's0004'
    prfs = {}
    for i, xx in enumerate(xs):
        for j, yy in enumerate(ys):
            cc = float(np.clip(xx + 1, 50, 2085)); rr = float(np.clip(yy + 1, 6, 2040))
            prfs[(i, j)] = PRF.TESS_PRF(cam, ccd, sec, cc, rr, localdatadir=f'{L}/prf_files/{sub}/')
    rows = []
    hw = 4
    for _, r in cand.iterrows():
        i = int(np.argmin(np.abs(xs - r.x))); j = int(np.argmin(np.abs(ys - r.y)))
        P = prfs[(i, j)]
        ix, iy = int(round(r.x)), int(round(r.y))
        st = img[iy - hw:iy + hw + 1, ix - hw:ix + hw + 1]
        if st.shape != (9, 9) or not np.all(np.isfinite(st)):
            continue
        yy, xx = np.mgrid[0:9, 0:9]
        x0 = r.x - (ix - hw); y0 = r.y - (iy - hw)
        sig = np.sqrt(np.clip(st, 1, None) + 30.0 ** 2)

        def lin(p):
            M = P.locate(p[0], p[1], (9, 9))
            A = np.column_stack([M.ravel(), np.ones(81), (xx - 4).ravel(), (yy - 4).ravel()]) / sig.ravel()[:, None]
            c, *_ = np.linalg.lstsq(A, st.ravel() / sig.ravel(), rcond=None)
            return A, c

        def res(p):
            A, c = lin(p)
            return A @ c - st.ravel() / sig.ravel()
        try:
            o = least_squares(res, [x0, y0], bounds=([x0 - 1.5, y0 - 1.5], [x0 + 1.5, y0 + 1.5]), diff_step=1e-4,
                              x_scale=0.05, xtol=1e-10, ftol=1e-12)
        except Exception:
            continue
        A, c = lin(o.x)
        if c[0] <= 0:
            continue
        J = o.jac
        try:
            cov = np.linalg.inv(J.T @ J) * max((o.fun ** 2).sum() / (81 - 6), 1.0)
        except np.linalg.LinAlgError:
            continue
        peak = st.max()
        rows.append((int(r.source_id), r.x, r.y, o.x[0] + ix - hw, o.x[1] + iy - hw, np.sqrt(cov[0, 0]), np.sqrt(cov[1, 1]),
                     c[0], c[1], r.phot_g_mean_mag, r.bp_rp, r.ruwe, peak, (o.fun ** 2).sum() / (81 - 6)))
    a = np.array(rows, dtype=[('sid', 'i8'), ('xw', 'f8'), ('yw', 'f8'), ('xf', 'f8'), ('yf', 'f8'), ('ex', 'f8'), ('ey', 'f8'),
                              ('flux', 'f8'), ('bkg', 'f8'), ('G', 'f8'), ('C', 'f8'), ('ruwe', 'f8'), ('peak', 'f8'), ('chi2r', 'f8')])
    base = os.path.basename(fn).replace('.fits', '')
    np.savez_compressed(f'{R}/results/x10{TAG}_{base}.npz', stars=a, sec=sec, cam=cam, ccd=ccd, xt=xt, yt=yt, xb=xb, yb=yb, ep=ep)
    print(f'{base}: S{sec} cam{cam} ccd{ccd} target ({xt:.1f},{yt:.1f}) boresight px ({xb:.0f},{yb:.0f}) '
          f'{len(cand)} isolated candidates, {len(a)} fitted [{time.time()-t0:.0f}s]', flush=True)
    return a


if __name__ == '__main__':
    files = sys.argv[1:] or sorted(glob.glob(f'{R}/ffi/*_ffic.fits'))
    for f in files:
        run(f)
