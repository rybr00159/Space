"""v4 field-star test of TESS lateral colour (colour-dependent PRF-fit positions) in single calibrated FFIs.

For each FFI (two per sector, same CCD as the target, different orbits):
  * Gaia DR3 G<17 within 3.6 deg (v4 query) propagated to the FFI epoch, projected with the FFI SIP WCS;
  * candidates 9.5 < G < 16.8, finite BP-RP, RUWE < 1.4, >= 6 px from CCD edges, within 2.6 deg of the target;
  * isolation from the deeper G<19 catalogue (v3 download, Gaia data only): contamination index
    sum_i 10^(-0.4 dG_i) exp(-d_i^2/2) (d in px, d<6) < 0.02 and no neighbour with dG < 0.5 within 8 px;
  * 7x7 stamp (same box as the difference-image fits) fitted with the TESS PRF (PRF package, objects on a 128-px grid)
    with free (x, y, flux, constant) by variable projection + Nelder-Mead;
  * residual = fitted - WCS(Gaia) position in px.
Output v4/results/f01_<ffi>.npz (per-star table, boresight, Jacobian at the target).
"""
import sys, os, glob, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from scipy.optimize import minimize
import PRF
import tess_stars2px as t2p
from v4lib import V4, BASE, RA0, DE0, jyear

SB = 3
GMIN, GMAX = 9.5, 16.8
RSEL = 2.6
CONTAM = float(os.environ.get('F01_CONTAM', 0.02))   # sum_i f_i/f * exp(-d_i^2/2) over neighbours (d in px)
wide = Table.read(f'{V4}/gaia/wide_3p6deg_G17.fits')
deep = Table.read(f'{BASE}/v3/gaia/wide_4p3deg_G19.fits')     # Gaia DR3 G<19 (used for isolation only)


def arr(t, c):
    x = t[c]
    return np.array(x.filled(np.nan) if hasattr(x, 'filled') else x, float)


def pa_sep(ra1, de1, ra2, de2):
    a1, d1, a2, d2 = map(np.radians, (ra1, de1, ra2, de2))
    y = np.sin(a2 - a1) * np.cos(d2)
    x = np.cos(d1) * np.sin(d2) - np.sin(d1) * np.cos(d2) * np.cos(a2 - a1)
    s = np.degrees(np.arccos(np.clip(np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(a2 - a1), -1, 1)))
    return np.degrees(np.arctan2(y, x)) % 360, s


def project(t, w, dty):
    ra, de = arr(t, 'ra'), arr(t, 'dec')
    de2 = de + np.nan_to_num(arr(t, 'pmdec')) * dty / 3.6e6
    ra2 = ra + np.nan_to_num(arr(t, 'pmra')) * dty / 3.6e6 / np.cos(np.radians(de))
    x, y = w.all_world2pix(ra2, de2, 0)
    return ra2, de2, x, y


def jac(w, xq, yq, h=0.05):
    r0 = np.array(w.all_pix2world([[xq, yq]], 0))[0]; c = np.cos(np.radians(r0[1]))
    cols = []
    for dx, dy in ((h, 0), (0, h)):
        a = np.array(w.all_pix2world([[xq + dx, yq + dy]], 0))[0]; b = np.array(w.all_pix2world([[xq - dx, yq - dy]], 0))[0]
        cols.append(np.array([((a[0] - b[0] + 180) % 360 - 180) * c * 3600, (a[1] - b[1]) * 3600]) / (2 * h))
    return np.column_stack(cols)


def run(fn):
    t0 = time.time()
    with fits.open(fn) as h:
        img = np.array(h[1].data, float); hd = h[1].header
    w = WCS(hd)
    base = os.path.basename(fn).replace('.fits', '')
    sec = int(base.split('-s0')[1][:3]); cam = int(hd['CAMERA']); ccd = int(hd['CCD'])
    ep = float(jyear(0.5 * (hd['TSTART'] + hd['TSTOP'])))
    texp = float(hd.get('EXPOSURE', 0.0)) * 86400.0 or 158.4
    dty = ep - 2016.0
    ra, de, x, y = project(wide, w, dty)
    _, _, xd, yd = project(deep, w, dty)
    gd = arr(deep, 'phot_g_mean_mag')
    xt, yt = [float(np.ravel(v)[0]) for v in w.all_world2pix([RA0 + 38.50 * dty / 3.6e6 / np.cos(np.radians(DE0))], [DE0 - 3.134 * dty / 3.6e6], 0)]
    res = t2p.tess_stars2px_function_entry(1, RA0, DE0, trySector=sec)
    sc = res[-1]; i = list(sc.sectors).index(sec)
    bra, bde = float(sc.camRa[cam - 1][i]), float(sc.camDec[cam - 1][i])
    G = arr(wide, 'phot_g_mean_mag'); C = arr(wide, 'bp_rp'); ruwe = arr(wide, 'ruwe'); plx = arr(wide, 'parallax')
    _, rt = pa_sep(RA0, DE0, arr(wide, 'ra'), arr(wide, 'dec'))
    cand = np.where((G > GMIN) & (G < GMAX) & np.isfinite(C) & (ruwe < 1.4) & (rt < RSEL) &
                    (x > 44 + 6) & (x < 2091 - 6) & (y > 6) & (y < 2047 - 6))[0]
    from scipy.spatial import cKDTree
    okd = np.isfinite(xd) & np.isfinite(yd)
    tree = cKDTree(np.column_stack([xd[okd], yd[okd]])); xdd, ydd, gdd = xd[okd], yd[okd], gd[okd]
    keep = []
    for k in cand:
        nn = np.array(tree.query_ball_point([x[k], y[k]], 8.0), int)
        dd = np.hypot(xdd[nn] - x[k], ydd[nn] - y[k]); dg = gdd[nn] - G[k]
        o = dd > 0.02
        contam = np.sum((10 ** (-0.4 * dg) * np.exp(-0.5 * dd ** 2))[o & (dd < 6)])
        bad = (contam > CONTAM) or np.any(o & (dd < 8) & (dg < 0.5))
        if not bad:
            keep.append(k)
    keep = np.array(keep, int)
    prfs = {}
    sub = 's0001' if sec < 4 else 's0004'
    def getprf(xx, yy):
        key = (int(xx // 128), int(yy // 128))
        if key not in prfs:
            cx = float(np.clip(key[0] * 128 + 64, 46, 2090)); cy = float(np.clip(key[1] * 128 + 64, 2, 2046))
            prfs[key] = PRF.TESS_PRF(cam, ccd, sec, cx + 1, cy + 1, localdatadir=f'{BASE}/prf_files/{sub}/')
        return prfs[key]
    rows = []
    for k in keep:
        Pm = getprf(x[k], y[k])
        ix, iy = int(round(x[k])), int(round(y[k]))
        st = img[iy - SB:iy + SB + 1, ix - SB:ix + SB + 1]
        if st.shape != (2 * SB + 1, 2 * SB + 1) or not np.all(np.isfinite(st)):
            continue
        x0, y0 = x[k] - (ix - SB), y[k] - (iy - SB)
        skyl = np.nanmedian(img[max(iy - 20, 0):iy + 21, max(ix - 20, 44):ix + 21])
        sig = np.sqrt(np.clip(st, 0, None) / texp + max(skyl, 0) / texp + 1.0)
        if not np.all(np.isfinite(sig)):
            continue
        b = (st / sig).ravel()
        def chi(p, ret=False):
            A = np.column_stack([Pm.locate(p[0], p[1], st.shape).ravel(), np.ones(st.size)]) / sig.ravel()[:, None]
            c, *_ = np.linalg.lstsq(A, b, rcond=None)
            r = A @ c - b
            return (float(r @ r), c) if ret else float(r @ r)
        # coarse 5x5 start grid (+-0.4 px) then Nelder-Mead
        try:
            best = min(((chi([x0 + gx, y0 + gy]), gx, gy) for gx in np.linspace(-0.4, 0.4, 5) for gy in np.linspace(-0.4, 0.4, 5)))
            o = minimize(chi, [x0 + best[1], y0 + best[2]], method='Nelder-Mead', options=dict(xatol=2e-4, fatol=1e-6, maxiter=300))
            c2, cf = chi(o.x, True)
        except (np.linalg.LinAlgError, ValueError):
            continue
        if cf[0] <= 0:
            continue
        # formal errors from finite-difference curvature
        hh = 0.01; H = np.zeros((2, 2))
        for a_ in range(2):
            for b_ in range(2):
                e1 = np.eye(2)[a_] * hh; e2 = np.eye(2)[b_] * hh
                H[a_, b_] = (chi(o.x + e1 + e2) - chi(o.x + e1 - e2) - chi(o.x - e1 + e2) + chi(o.x - e1 - e2)) / (4 * hh * hh)
        try:
            cov = 2 * np.linalg.inv(H) * max(c2 / (st.size - 4), 1.0)
        except np.linalg.LinAlgError:
            continue
        if not np.all(np.isfinite(cov)) or cov[0, 0] <= 0 or cov[1, 1] <= 0:
            continue
        pa, th = pa_sep(ra[k], de[k], bra, bde)
        rows.append((int(wide['source_id'][k]), x[k], y[k], o.x[0] + ix - SB, o.x[1] + iy - SB, np.sqrt(cov[0, 0]), np.sqrt(cov[1, 1]),
                     cf[0], cf[1], G[k], C[k], ruwe[k], plx[k], arr(wide, 'parallax_error')[k], st.max(), c2 / (st.size - 4), pa, th, rt[k]))
    a = np.array(rows, dtype=[('sid', 'i8'), ('xw', 'f8'), ('yw', 'f8'), ('xf', 'f8'), ('yf', 'f8'), ('ex', 'f8'), ('ey', 'f8'),
                              ('flux', 'f8'), ('bkg', 'f8'), ('G', 'f8'), ('C', 'f8'), ('ruwe', 'f8'), ('plx', 'f8'), ('eplx', 'f8'),
                              ('peak', 'f8'), ('chi2r', 'f8'), ('pa_bs', 'f8'), ('theta', 'f8'), ('rt', 'f8')])
    J = jac(w, xt, yt)
    pat, tht = pa_sep(RA0, DE0, bra, bde)
    np.savez_compressed(f'{V4}/results/f01_{base}.npz', stars=a, sec=sec, cam=cam, ccd=ccd, xt=xt, yt=yt, ep=ep, J=J, bra=bra, bde=bde,
                        pa_t=pat, theta_t=tht, texp=texp)
    print(f'{base}: S{sec} cam{cam} ccd{ccd} target ({xt:.1f},{yt:.1f}) PA->boresight {pat:.1f} theta {tht:.2f}; '
          f'{len(cand)} cand, {len(keep)} isolated, {len(a)} fitted [{time.time() - t0:.0f}s]', flush=True)


if __name__ == '__main__':
    for f in sys.argv[1:]:
        run(f)
