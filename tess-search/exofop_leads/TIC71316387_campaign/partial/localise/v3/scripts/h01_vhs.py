"""VISTA VHS (ESO Phase-3, SODA cutouts) companion search around the partly saturated target.

Per image (Ks: 2 pawprints + tile; J: 2 pawprints + tile; seeing 1.07-1.36"):
  * empirical PSF from unsaturated Gaia stars in the 6' cutout (raw peak < 20000 ADU, G < 17.3, no Gaia neighbour with
    dG < 5 within 6"), recentred with cubic-spline shifts, flux-normalised, median-combined;
  * target model:  F * ePSF(x0, y0) + b  [+ f * F * ePSF(x0 + dx, y0 + dy)],  fitted on a 31x31-px stamp with pixels
    whose raw value exceeds SATLIM (non-linear / saturated core) masked;  noise sigma^2 = raw / gain + rn^2;
  * scan: companion on a grid 0.6-4.0" (0.1" steps in x/y), linear f >= 0, dchi2 = chi2_single - chi2_pair;
  * null: the same scan for the unsaturated comparison stars (leave-one-out ePSF) -> dchi2 distribution of the best
    point per star, scaled to the target's S/N is NOT assumed; instead the detection threshold is the maximum dchi2 of
    the target scan restricted to an annulus where no companion is expected?  No: threshold = injection-calibrated:
  * injection: a companion (ePSF copy, flux ratio 10^(-0.4 dm)) added to the real target stamp at sep 0.8, 1.0, 1.25,
    1.5, 2, 3" x 8 PAs, dm = 1, 1.5, 2, 2.5, 3, 3.5, 4; recovered if the best-fit companion lies within 0.3" of the
    injected position and dchi2 exceeds the real (un-injected) target scan maximum.
Output: results/h01_vhs.json
"""
import sys, os, glob, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import shift as ndshift
from scipy.optimize import least_squares
from v3lib import V3, to_json
import warnings; warnings.filterwarnings('ignore')

SATLIM = float(os.environ.get('H01_SATLIM', 24000))
GAIN, RN = 4.2, 25.0
H = 15
EP = 2012.64
g = pd.read_csv(V3 + '/gaia/field_6arcmin.csv')


def load(f):
    h = fits.open(f); ext = [e for e in h if e.data is not None][0]
    d = ext.data.astype(float); w = WCS(ext.header)
    ra = g.ra + np.nan_to_num(g.pmra) * (EP - 2016) / 3.6e6 / np.cos(np.radians(g.dec))
    de = g.dec + np.nan_to_num(g.pmdec) * (EP - 2016) / 3.6e6
    x, y = w.all_world2pix(ra, de, 0)
    # pixel scale & orientation
    r0 = w.all_pix2world([[x[0], y[0]], [x[0] + 1, y[0]], [x[0], y[0] + 1]], 0)
    c = np.cos(np.radians(r0[0, 1]))
    J = np.array([[(r0[1, 0] - r0[0, 0]) * c * 3600, (r0[2, 0] - r0[0, 0]) * c * 3600],
                  [(r0[1, 1] - r0[0, 1]) * 3600, (r0[2, 1] - r0[0, 1]) * 3600]])
    return dict(d=d, x=x, y=y, J=J, hdr=ext.header, band=('Ks' if 'Ks' in str(ext.header.get('HIERARCH ESO INS FILT1 NAME', '')) or ext.header.get('FILTER', '') == 'Ks' else 'J'))


def stamp(d, x, y, n=H):
    ix, iy = int(round(x)), int(round(y))
    return d[iy - n:iy + n + 1, ix - n:ix + n + 1].copy(), x - ix, y - iy


def local_sky(d, x, y):
    ix, iy = int(round(x)), int(round(y))
    a = d[iy - 40:iy + 41, ix - 40:ix + 41]
    yy, xx = np.mgrid[-40:41, -40:41]
    r = np.hypot(xx, yy)
    return float(np.median(a[(r > 25) & (r < 40)]))


def build_epsf(I, idx, exclude=None):
    cube = []
    for i in idx:
        if i == exclude:
            continue
        st, fx, fy = stamp(I['d'], I['x'][i], I['y'][i], H + 3)
        st = st - local_sky(I['d'], I['x'][i], I['y'][i])
        st = ndshift(st, (-fy, -fx), order=3, mode='nearest')[3:-3, 3:-3]
        cube.append(st / st.sum())
    E = np.median(np.array(cube), axis=0)
    return E / E.sum()


def psf(E, dx, dy):
    return ndshift(E, (dy, dx), order=3, mode='constant', cval=0.0)


def fit_single(st, sig, mask, E):
    def res(p):
        m = psf(E, p[0], p[1])
        A = np.column_stack([m[mask], np.ones(mask.sum())]) / sig[mask][:, None]
        c, *_ = np.linalg.lstsq(A, st[mask] / sig[mask], rcond=None)
        return A @ c - st[mask] / sig[mask]
    r = least_squares(res, [0.0, 0.0], diff_step=1e-3)
    m = psf(E, *r.x)
    A = np.column_stack([m[mask], np.ones(mask.sum())]) / sig[mask][:, None]
    c, *_ = np.linalg.lstsq(A, st[mask] / sig[mask], rcond=None)
    return r.x, c, float((r.fun ** 2).sum())


class Scanner:
    """Pre-computed companion templates on a grid (arcsec, image axes) around a fixed primary position p0."""
    def __init__(self, E, p0, pix, rmin=0.5, rmax=4.0, step=0.1):
        grid = np.arange(-rmax, rmax + 1e-9, step)
        pts = [(a, b) for a in grid for b in grid if rmin <= np.hypot(a, b) <= rmax]
        self.pts = np.array(pts)
        self.T = np.array([psf(E, p0[0] + a / pix, p0[1] + b / pix) for a, b in pts])   # (K, ny, nx)
        self.m0 = psf(E, *p0)

    def run(self, st, sig, mask):
        w = 1 / sig[mask]; y = st[mask] * w
        m0 = self.m0[mask] * w; one = w
        T = self.T[:, mask] * w[None, :]
        # single model chi2 (primary fixed position, flux + background)
        A0 = np.column_stack([m0, one]); c0, *_ = np.linalg.lstsq(A0, y, rcond=None)
        chi0 = float(((A0 @ c0 - y) ** 2).sum())
        # pair model for every grid point via normal equations
        K = len(T)
        N = np.empty((K, 3, 3)); bv = np.empty((K, 3))
        N[:, 0, 0] = m0 @ m0; N[:, 0, 1] = N[:, 1, 0] = T @ m0; N[:, 0, 2] = N[:, 2, 0] = m0 @ one
        N[:, 1, 1] = (T * T).sum(1); N[:, 1, 2] = N[:, 2, 1] = T @ one; N[:, 2, 2] = one @ one
        bv[:, 0] = m0 @ y; bv[:, 1] = T @ y; bv[:, 2] = one @ y
        c = np.linalg.solve(N, bv[:, :, None])[:, :, 0]
        chi = (y @ y) - (c * bv).sum(1)
        d = chi0 - chi
        d[c[:, 1] <= 0] = 0
        k = int(np.argmax(d))
        return (float(d[k]), tuple(self.pts[k]), float(c[k, 1] / c[k, 0])), chi0, c0


def run(f):
    I = load(f)
    d = I['d']; x, y = I['x'], I['y']
    it = int(np.where(g.source_id.values == 6571343840265900544)[0][0])
    pix = float(np.sqrt(abs(np.linalg.det(I['J']))))
    # orientation: image +x / +y in sky (E, N) for converting offsets
    Jn = I['J'] / pix
    # PSF stars
    cand = []
    for i in range(len(g)):
        if i == it or g.phot_g_mean_mag[i] > 17.3:
            continue
        if not (H + 45 < x[i] < d.shape[1] - H - 45 and H + 45 < y[i] < d.shape[0] - H - 45):
            continue
        dd = np.hypot((g.ra - g.ra[i]) * np.cos(np.radians(g.dec[i])), g.dec - g.dec[i]) * 3600
        if np.any((dd > 0.01) & (dd < 6) & (g.phot_g_mean_mag - g.phot_g_mean_mag[i] < 5)):
            continue
        st, _, _ = stamp(d, x[i], y[i])
        if st.max() > 20000:
            continue
        cand.append(i)
    E = build_epsf(I, cand)
    sky = local_sky(d, x[it], y[it])
    raw, fx, fy = stamp(d, x[it], y[it])
    st = raw - sky
    sig = np.sqrt(np.clip(raw, 1, None) / GAIN + RN ** 2)
    yy, xx = np.mgrid[-H:H + 1, -H:H + 1]
    mask = (raw < SATLIM) & (np.hypot(xx, yy) <= H)
    nsat = int((~(raw < SATLIM) & (np.hypot(xx, yy) <= H)).sum())
    rsat = float(np.sqrt(nsat / np.pi) * pix)
    p0, c0, chi0 = fit_single(st, sig, mask, E)
    SC = Scanner(E, p0, pix, rmin=0.6, rmax=4.0, step=0.1)
    best, chi0f, c0f = SC.run(st, sig, mask)
    real_max = best[0]
    # null: comparison stars (leave-one-out ePSF), dchi2 maximum per star (normalised by chi2/dof to absorb S/N)
    nulls = []
    for i in cand:
        Ei = build_epsf(I, cand, exclude=i)
        rw, _, _ = stamp(d, x[i], y[i]); s_ = rw - local_sky(d, x[i], y[i]); sg = np.sqrt(np.clip(rw, 1, None) / GAIN + RN ** 2)
        mk = np.hypot(xx, yy) <= H
        pp, cc, ch = fit_single(s_, sg, mk, Ei)
        b, _, _ = Scanner(Ei, pp, pix, rmin=0.6, rmax=4.0, step=0.1).run(s_, sg, mk)
        nulls.append(dict(G=float(g.phot_g_mean_mag[i]), chi2r=ch / (mk.sum() - 4), dchi2_max=b[0], f=b[2], sep=float(np.hypot(*b[1])) if b[1] else None))
    # injections
    rng = np.random.default_rng(5)
    inj = []
    F = c0f[0]
    for sep in (0.8, 1.0, 1.25, 1.5, 2.0, 3.0):
        for dm in (1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5):
            rec = []
            for pa in np.arange(0, 360, 45.0) + rng.uniform(0, 45):
                dE, dN = sep * np.sin(np.radians(pa)), sep * np.cos(np.radians(pa))
                dpx = np.linalg.solve(Jn, np.array([dE, dN])) / 1.0   # arcsec in image axes
                fr = 10 ** (-0.4 * dm)
                comp = F * fr * psf(E, p0[0] + dpx[0] / pix * 1.0, p0[1] + dpx[1] / pix * 1.0)
                st2 = st + comp
                raw2 = raw + comp
                mk2 = (raw2 < SATLIM) & (np.hypot(xx, yy) <= H)
                sg2 = np.sqrt(np.clip(raw2, 1, None) / GAIN + RN ** 2)
                b, _, _ = SC.run(st2, sg2, mk2)
                ok = b[1] is not None and np.hypot(b[1][0] - dpx[0], b[1][1] - dpx[1]) < 0.35 and b[0] > real_max
                rec.append(bool(ok))
            inj.append(dict(sep=sep, dm=dm, frac=float(np.mean(rec))))
    # best real candidate in sky coordinates
    if best[1] is not None:
        sky_off = Jn @ np.array(best[1])
        bsep = float(np.hypot(*best[1])); bpa = float(np.degrees(np.arctan2(sky_off[0], sky_off[1])) % 360)
    else:
        bsep = bpa = None
    res = dict(file=os.path.basename(f), fwhm=float(I['hdr'].get('PSF_FWHM')), pix=pix, n_psf_stars=len(cand), satlim=SATLIM,
               masked_core_radius_arcsec=rsat, target_chi2r=chi0 / (mask.sum() - 4), real_dchi2_max=real_max, real_best_sep=bsep,
               real_best_pa=bpa, real_best_dm=(float(-2.5 * np.log10(best[2])) if best[2] else None), nulls=nulls, injections=inj)
    print(f"{res['file']}: FWHM {res['fwhm']:.2f}\" psf stars {len(cand)} masked core r~{rsat:.2f}\" chi2r {res['target_chi2r']:.2f}; "
          f"real best dchi2 {real_max:.1f} at {bsep}\" PA {bpa} dm {res['real_best_dm']}; null dchi2 max per star: "
          f"{[round(n['dchi2_max'], 1) for n in nulls]} (chi2r {[round(n['chi2r'], 2) for n in nulls]})", flush=True)
    for sep in (0.8, 1.0, 1.25, 1.5, 2.0, 3.0):
        print(f"   sep {sep}\": recovered fraction vs dm: " + ' '.join(f"{r['dm']}:{r['frac']:.2f}" for r in inj if r['sep'] == sep), flush=True)
    return res


if __name__ == '__main__':
    out = [run(f) for f in sorted(glob.glob(V3 + '/imaging/vhs_*.fits.fz'))]
    to_json(out, V3 + '/results/h01_vhs.json')
