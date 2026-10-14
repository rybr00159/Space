"""v3 independent test of TESS lateral colour with field stars in single calibrated FFIs.

For one FFI per sector (same CCD as the target):
  * Gaia DR3 (VizieR I/355, G < 18.5 within 2.7 deg) propagated to the FFI epoch, projected with the FFI's SIP WCS;
  * isolated stars G 9.5-14.5 within +-HW px of the target: no Gaia source with dG < 5 within 4 px and none with
    dG < 2.5 within 7 px;
  * 7x7-px stamp (the same box as the difference-image fits) fitted with the TESS PRF (PRF package; PRF objects on a
    100-px grid of CCD positions) with free (x, y, flux, constant);
  * residual = fitted - WCS(Gaia) position [px].
Output: f01_<ffi>.npz with the per-star table, the boresight sky direction, and the pixel->sky Jacobian at the target.
"""
import sys, os, glob, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from scipy.optimize import least_squares
import PRF
import tess_stars2px as t2p
from v3lib import V3, BASE, RA0, DE0, jyear

HW = int(os.environ.get('F01_HW', 430))
SB = 3   # stamp half-size -> 7x7
GMIN = float(os.environ.get('F01_GMIN', 9.5)); GMAX = float(os.environ.get('F01_GMAX', 15.5))
ISO_R1, ISO_DG1 = 3.0, 4.0
ISO_R2, ISO_DG2 = 6.0, 2.0
TAG = os.environ.get('F01_TAG', '')
gw = Table.read(f"{V3}/gaia/{os.environ.get('F01_CAT', 'wide_2p7deg_G18p5.fits')}").to_pandas()
RCMAX = float(os.environ.get('F01_RCMAX', 2.45))
gw = gw[np.isfinite(gw.ra) & np.isfinite(gw.dec)].reset_index(drop=True)


def boresight(sec):
    res = t2p.tess_stars2px_function_entry(1, RA0, DE0, trySector=sec)
    sc = res[-1]; cam = int(res[4][0]); i = list(sc.sectors).index(sec)
    return float(sc.camRa[cam - 1][i]), float(sc.camDec[cam - 1][i])


def pa_to(ra1, de1, ra2, de2):
    """Position angle (deg E of N) of point 2 seen from point 1 (great circle) and separation (deg)."""
    a1, d1, a2, d2 = map(np.radians, (ra1, de1, ra2, de2))
    y = np.sin(a2 - a1) * np.cos(d2)
    x = np.cos(d1) * np.sin(d2) - np.sin(d1) * np.cos(d2) * np.cos(a2 - a1)
    sep = np.degrees(np.arccos(np.clip(np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(a2 - a1), -1, 1)))
    return np.degrees(np.arctan2(y, x)) % 360, sep


def run(fn):
    t0 = time.time()
    with fits.open(fn) as h:
        img = np.array(h[1].data, float); hd = h[1].header
    w = WCS(hd)
    sec = int(os.path.basename(fn).split('-s0')[1][:3]); cam = int(hd['CAMERA']); ccd = int(hd['CCD'])
    ep = float(jyear(0.5 * (hd['TSTART'] + hd['TSTOP'])))
    dt = ep - 2016.0
    g = gw.copy()
    g['rae'] = g.ra + np.nan_to_num(g.pmra) * dt / 3.6e6 / np.cos(np.radians(g.dec))
    g['dee'] = g.dec + np.nan_to_num(g.pmdec) * dt / 3.6e6
    x, y = w.all_world2pix(g.rae.values, g.dee.values, 0)
    g['x'] = x; g['y'] = y
    xt, yt = [float(v[0]) for v in w.all_world2pix([RA0 + 38.5 * dt / 3.6e6 / np.cos(np.radians(DE0))], [DE0 - 3.13 * dt / 3.6e6], 0)]
    bra, bde = boresight(sec)
    near = g[(np.abs(g.x - xt) < HW + 10) & (np.abs(g.y - yt) < HW + 10)].reset_index(drop=True)
    rc = np.hypot((near.ra - RA0) * np.cos(np.radians(DE0)), near.dec - DE0)
    cand = near[(near.phot_g_mean_mag > GMIN) & (near.phot_g_mean_mag < GMAX) & np.isfinite(near.bp_rp) & (rc < RCMAX) &
                (np.abs(near.x - xt) < HW) & (np.abs(near.y - yt) < HW) &
                (near.x > 44 + 6) & (near.x < 2091 - 6) & (near.y > 6) & (near.y < 2047 - 6)]
    nx_, ny_, ng_ = near.x.values, near.y.values, near.phot_g_mean_mag.values
    keep = []
    for _, r in cand.iterrows():
        d = np.hypot(nx_ - r.x, ny_ - r.y); dg = ng_ - r.phot_g_mean_mag
        o = d > 0.01
        bad = np.any(o & (d < ISO_R1) & (dg < ISO_DG1)) or np.any(o & (d < ISO_R2) & (dg < ISO_DG2))
        keep.append(not bad)
    cand = cand[np.array(keep, bool)]
    # PRF grid
    sub = 's0001' if sec < 4 else 's0004'
    gx = np.arange(xt - HW, xt + HW + 1, 100.0); gy = np.arange(yt - HW, yt + HW + 1, 100.0)
    prfs = {}
    for i, cx in enumerate(gx):
        for j, cy in enumerate(gy):
            prfs[(i, j)] = PRF.TESS_PRF(cam, ccd, sec, float(np.clip(cx + 1, 46, 2090)), float(np.clip(cy + 1, 2, 2046)),
                                        localdatadir=f'{BASE}/prf_files/{sub}/')
    rows = []
    yy, xx = np.mgrid[0:2 * SB + 1, 0:2 * SB + 1]
    for _, r in cand.iterrows():
        Pm = prfs[(int(np.argmin(np.abs(gx - r.x))), int(np.argmin(np.abs(gy - r.y))))]
        ix, iy = int(round(r.x)), int(round(r.y))
        st = img[iy - SB:iy + SB + 1, ix - SB:ix + SB + 1]
        if st.shape != (2 * SB + 1, 2 * SB + 1) or not np.all(np.isfinite(st)):
            continue
        x0, y0 = r.x - (ix - SB), r.y - (iy - SB)
        sig = np.sqrt(np.clip(st, 1, None) * 1.0 + 20.0 ** 2)
        b = (st / sig).ravel()

        def lin(p):
            A = np.column_stack([Pm.locate(p[0], p[1], st.shape).ravel(), np.ones(st.size)]) / sig.ravel()[:, None]
            c, *_ = np.linalg.lstsq(A, b, rcond=None)
            return A, c

        def res(p):
            A, c = lin(p); return A @ c - b
        try:
            o = least_squares(res, [x0, y0], bounds=([x0 - 1.2, y0 - 1.2], [x0 + 1.2, y0 + 1.2]), diff_step=2e-4, x_scale=0.05,
                              xtol=1e-10, ftol=1e-12)
        except Exception:
            continue
        A, c = lin(o.x)
        if c[0] <= 0:
            continue
        chi2r = (o.fun ** 2).sum() / (st.size - 4)
        try:
            cov = np.linalg.inv(o.jac.T @ o.jac) * max(chi2r, 1.0)
        except np.linalg.LinAlgError:
            continue
        pa, sep = pa_to(r.rae, r.dee, bra, bde)
        rows.append((int(r.source_id), r.x, r.y, o.x[0] + ix - SB, o.x[1] + iy - SB, np.sqrt(cov[0, 0]), np.sqrt(cov[1, 1]),
                     c[0], c[1], r.phot_g_mean_mag, r.bp_rp, r.ruwe if np.isfinite(r.ruwe) else 99.0, st.max(), chi2r, pa, sep,
                     r.parallax if np.isfinite(r.parallax) else np.nan))
    a = np.array(rows, dtype=[('sid', 'i8'), ('xw', 'f8'), ('yw', 'f8'), ('xf', 'f8'), ('yf', 'f8'), ('ex', 'f8'), ('ey', 'f8'),
                              ('flux', 'f8'), ('bkg', 'f8'), ('G', 'f8'), ('C', 'f8'), ('ruwe', 'f8'), ('peak', 'f8'),
                              ('chi2r', 'f8'), ('pa_bs', 'f8'), ('theta', 'f8'), ('plx', 'f8')])
    # local pixel->sky Jacobian at the target
    def jac(xq, yq, hh=0.01):
        r0 = np.array(w.all_pix2world([[xq, yq]], 0))[0]
        cc = np.cos(np.radians(r0[1]))
        f = lambda rr: np.array([(rr[0] - r0[0]) * cc * 3600, (rr[1] - r0[1]) * 3600])
        return np.column_stack([f(np.array(w.all_pix2world([[xq + hh, yq]], 0))[0]) / hh,
                                f(np.array(w.all_pix2world([[xq, yq + hh]], 0))[0]) / hh])
    J = jac(xt, yt)
    pat, tht = pa_to(RA0, DE0, bra, bde)
    base = os.path.basename(fn).replace('.fits', '')
    np.savez_compressed(f'{V3}/results/f01{TAG}_{base}.npz', stars=a, sec=sec, cam=cam, ccd=ccd, xt=xt, yt=yt, ep=ep, J=J,
                        bra=bra, bde=bde, pa_t=pat, theta_t=tht)
    print(f'{base}: S{sec} cam{cam} ccd{ccd} target ({xt:.1f},{yt:.1f}) PA->boresight {pat:.1f} theta {tht:.2f} deg; '
          f'{len(cand)} isolated, {len(a)} fitted [{time.time() - t0:.0f}s]', flush=True)


if __name__ == '__main__':
    for f in (sys.argv[1:] or sorted(glob.glob(f'{V3}/ffi/*_ffic.fits'))):
        run(f)
