"""SkyMapper DR2 short-survey images (unsaturated for the V = 11.9 target): search for an elongation / companion.

Per image:
  * comparison stars: Gaia DR3 G 10.0-14.5 in the 6' cutout, no Gaia source with dG < 6 within 10", peak < 30000 ADU,
    >= 20 px from the edge;
  * empirical PSF (ePSF): background-subtracted, flux-normalised 41x41 stamps of comparison stars, recentred with cubic
    spline shifts (centroids from a Moffat fit), median-combined (leave-one-out ePSF for each comparison star);
  * every star (target and comparisons) is fitted on a 25x25 stamp with
      single: F*ePSF(x0,y0) + b              and
      binary: F*[ePSF(x0,y0) + f*ePSF(x0+dx, y0+dy)] + b   (grid over sep 0.4-3.0", PA 0-360 deg, f >= 0 linear)
    statistic: dchi2 = chi2_single - chi2_binary, plus the best (sep, PA, f);
  * null distribution: comparison stars with their leave-one-out ePSF;
  * sensitivity: inject a companion (a shifted, scaled copy of the target's own stamp) into the target stamp,
    sep 0.5/0.75/1.0/1.5/2.0", 8 PAs, dm = 2.0/2.5/3.0/3.5, and record dchi2.
Output: results/x20_skymapper.json
"""
import sys, os, glob, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from scipy.ndimage import shift as ndshift
from scipy.optimize import least_squares
from loclib import R, to_json

IMG = f'{R}/imaging'
RA0, DE0 = 330.887948, -42.835511
PMRA, PMDE = 38.503, -3.134
gall = Table.read(f'{R}/gaia/cone3deg_all.fits').to_pandas()
H = 20           # half size of ePSF stamp
h = 12           # half size of fit stamp
yy, xx = np.mgrid[-h:h + 1, -h:h + 1]


def moffat(p, n):
    x0, y0, a, al, be, b = p
    yy_, xx_ = np.mgrid[-n:n + 1, -n:n + 1]
    r2 = (xx_ - x0) ** 2 + (yy_ - y0) ** 2
    return a * (1 + r2 / al ** 2) ** (-be) + b


def stamp(d, x, y, n):
    ix, iy = int(round(x)), int(round(y))
    return d[iy - n:iy + n + 1, ix - n:ix + n + 1].astype(float), x - ix, y - iy


def build_epsf(stamps):
    """stamps: list of (st, fx, fy, bkg) with st of size 2H+1 centred at integer pixel; fx,fy fractional offsets."""
    cube = []
    for st, fx, fy, bk in stamps:
        s = st - bk
        s = ndshift(s, (-fy, -fx), order=3, mode='nearest')
        s /= s[H - 10:H + 11, H - 10:H + 11].sum()
        cube.append(s)
    return np.median(np.array(cube), axis=0)


def psf_at(E, dx, dy):
    """ePSF shifted by (dx, dy) px, cropped to the fit stamp."""
    s = ndshift(E, (dy, dx), order=3, mode='constant', cval=0.0)
    return s[H - h:H + h + 1, H - h:H + h + 1]


def broaden(E, sx, sy):
    from scipy.ndimage import gaussian_filter
    if sx < 0.05 and sy < 0.05:
        return E
    return gaussian_filter(E, (max(sy, 1e-3), max(sx, 1e-3)), mode='constant')


def fit_star(st, E0, ivar, sc_px, J=None):
    """single (ePSF broadened by a free axis-aligned Gaussian) vs binary (same broadened PSF twice).
    J: pixel->sky (dE, dN) jacobian in arcsec/px for reporting the companion PA."""
    sw = np.sqrt(ivar.ravel())

    def lin(Ms):
        A = np.column_stack([m.ravel() for m in Ms] + [np.ones(st.size)]) * sw[:, None]
        c, *_ = np.linalg.lstsq(A, st.ravel() * sw, rcond=None)
        r = st.ravel() * sw - A @ c
        return c, (r ** 2).sum(), r

    def res1(p):
        E = broaden(E0, abs(p[2]), abs(p[3]))
        return lin([psf_at(E, p[0], p[1])])[2]
    o = least_squares(res1, [0.0, 0.0, 0.3, 0.3], diff_step=1e-3, x_scale=[0.05, 0.05, 0.1, 0.1])
    x0, y0, sx, sy = o.x[0], o.x[1], abs(o.x[2]), abs(o.x[3])
    c1 = (o.fun ** 2).sum()
    E = broaden(E0, sx, sy)
    P0 = psf_at(E, x0, y0)
    best = (c1, 0.0, 0.0, 0.0)
    for sarc in np.arange(0.5, 3.01, 0.1):
        sp = sarc / sc_px
        for pa in np.radians(np.arange(0, 360, 10)):
            dx = sp * np.cos(pa); dy = sp * np.sin(pa)
            c, chi, _ = lin([P0, psf_at(E, x0 + dx, y0 + dy)])
            if c[1] < 0:
                continue
            if chi < best[0]:
                best = (chi, dx, dy, c[1] / c[0])
    chi_b, bdx, bdy, f_b = best
    if f_b > 0:
        def res2(p):
            Eb = broaden(E0, abs(p[4]), abs(p[5]))
            return lin([psf_at(Eb, p[0], p[1]), psf_at(Eb, p[0] + p[2], p[1] + p[3])])[2]
        o2 = least_squares(res2, [x0 - f_b * bdx, y0 - f_b * bdy, bdx, bdy, sx, sy], diff_step=1e-3,
                           x_scale=[0.05, 0.05, 0.1, 0.1, 0.1, 0.1])
        Eb = broaden(E0, abs(o2.x[4]), abs(o2.x[5]))
        c, chi2b, _ = lin([psf_at(Eb, o2.x[0], o2.x[1]), psf_at(Eb, o2.x[0] + o2.x[2], o2.x[1] + o2.x[3])])
        if chi2b < chi_b and c[1] >= 0 and c[0] > 0:
            chi_b = chi2b; bdx, bdy = o2.x[2], o2.x[3]; f_b = c[1] / c[0]
    sky = (J @ np.array([bdx, bdy])) if J is not None else np.array([np.nan, np.nan])
    return dict(chi2_single=float(c1), dchi2=float(c1 - chi_b), sep=float(np.hypot(bdx, bdy) * sc_px), f=float(f_b),
                dE=float(sky[0]), dN=float(sky[1]), pa=float(np.degrees(np.arctan2(sky[0], sky[1])) % 360),
                sx=float(sx), sy=float(sy), npix=int(st.size))


def run(fn):
    hd = fits.open(fn)[0]
    d = hd.data.astype(float); hdr = hd.header; w = WCS(hdr)
    gain = float(hdr.get('GAIN', 1.0)); sat = float(hdr.get('SATURATE', 65000))
    filt = os.path.basename(fn)
    mjd = float(hdr['MJD-OBS']); ep = 2000 + (mjd - 51544.5) / 365.25
    sc_px = abs(hdr.get('CD2_2', 0.000138)) * 3600
    g = gall.copy()
    # pm only for the target (others: use 2016 positions; < 0.1" error for typical pm over 1-2 yr)
    x, y = w.all_world2pix(g.ra.values, g.dec.values, 0)
    g['x'] = x; g['y'] = y
    tx, ty = w.all_world2pix([[RA0 + PMRA * (ep - 2016) / 3.6e6 / np.cos(np.radians(DE0)), DE0 + PMDE * (ep - 2016) / 3.6e6]], 0)[0]
    ny, nx = d.shape
    bkg_all = np.median(d)
    r0 = np.array(w.all_pix2world([[tx, ty]], 0))[0]
    def _off(p):
        q = np.array(w.all_pix2world([p], 0))[0]
        return np.array([(q[0] - r0[0]) * np.cos(np.radians(r0[1])) * 3600, (q[1] - r0[1]) * 3600])
    J = np.column_stack([_off([tx + 1, ty]), _off([tx, ty + 1])])
    comps = []
    for i, r in g[(g.phot_g_mean_mag > 10.5) & (g.phot_g_mean_mag < 15.5)].iterrows():
        if not (H + 2 < r.x < nx - H - 3 and H + 2 < r.y < ny - H - 3):
            continue
        if np.hypot(r.x - tx, r.y - ty) < 30:
            continue
        dd = np.hypot(g.x - r.x, g.y - r.y) * sc_px
        dg = g.phot_g_mean_mag - r.phot_g_mean_mag
        if np.any((dd > 0.05) & (dd < 8) & (dg < 4)):
            continue
        st, fx, fy = stamp(d, r.x, r.y, H)
        if st.max() > 30000 or not np.all(np.isfinite(st)):
            continue
        comps.append((r, st, fx, fy))
    # refine centroids with a Moffat fit, background from stamp border
    stamps = []
    for r, st, fx, fy in comps:
        bk = np.median(np.r_[st[:3].ravel(), st[-3:].ravel(), st[:, :3].ravel(), st[:, -3:].ravel()])
        o = least_squares(lambda p: (moffat(p, H) - st).ravel(), [fx, fy, st.max() - bk, 4.0, 2.5, bk],
                          bounds=([-3, -3, 0, 0.5, 1.0, -1e4], [3, 3, 1e6, 30, 10, 1e5]))
        stamps.append((st, o.x[0], o.x[1], o.x[5], r, o.x))
    fwhm = [2 * s[5][3] * np.sqrt(2 ** (1 / s[5][4]) - 1) * sc_px for s in stamps]
    print(f'{filt}: {len(stamps)} comparison stars, median Moffat FWHM {np.median(fwhm):.2f}" (QAFWHM {hdr.get("QAFWHM")})', flush=True)
    if len(stamps) < 4:
        return None
    Eall = build_epsf([(s[0], s[1], s[2], s[3]) for s in stamps])
    noise_var = lambda st, bk: np.clip(st, 0, None) / gain + (np.std(np.r_[st[:2].ravel(), st[-2:].ravel()]) ** 2) + (0.01 * np.clip(st - bk, 0, None)) ** 2
    # comparison stars (null): leave-one-out ePSF
    null = []
    for k, s in enumerate(stamps):
        E = build_epsf([(q[0], q[1], q[2], q[3]) for j, q in enumerate(stamps) if j != k])
        st_full = s[0]
        st = st_full[H - h:H + h + 1, H - h:H + h + 1]
        # recentre ePSF frame on the star's integer pixel: fit offsets relative to stamp centre
        iv = 1.0 / noise_var(st, s[3])
        r = fit_star(st, E, iv, sc_px, J)
        r.update(G=float(s[4].phot_g_mean_mag), peak=float(st.max()))
        null.append(r)
    # target
    st_full, fx, fy = stamp(d, tx, ty, H)
    bk = np.median(np.r_[st_full[:3].ravel(), st_full[-3:].ravel(), st_full[:, :3].ravel(), st_full[:, -3:].ravel()])
    st = st_full[H - h:H + h + 1, H - h:H + h + 1]
    iv = 1.0 / noise_var(st, bk)
    tr = fit_star(st, Eall, iv, sc_px, J)
    tr.update(peak=float(st.max()))
    # sensitivity: inject a companion built from the target's own background-subtracted stamp (shifted copy)
    sens = []
    rng = np.random.default_rng(3)
    base = st_full - bk
    for sep in [0.5, 0.75, 1.0, 1.5, 2.0]:
        for dm in [2.0, 2.5, 3.0, 3.5]:
            ds = []
            for pa in np.radians(np.arange(0, 360, 45) + 10):
                dx, dy = np.linalg.solve(J, np.array([sep * np.sin(pa), sep * np.cos(pa)]))
                comp = ndshift(base, (dy, dx), order=3, mode='constant', cval=0.0) * 10 ** (-0.4 * dm)
                sti = (st_full + comp)[H - h:H + h + 1, H - h:H + h + 1]
                sti = sti + rng.normal(0, np.sqrt(np.clip(comp[H - h:H + h + 1, H - h:H + h + 1], 0, None) / gain))
                ri = fit_star(sti, Eall, 1.0 / noise_var(sti, bk), sc_px, J)
                inj_sky = J @ np.array([dx, dy])
                ok = np.hypot(ri['dE'] - inj_sky[0], ri['dN'] - inj_sky[1]) < 0.4
                ds.append((ri['dchi2'] - tr['dchi2'], bool(ok), ri['f']))
            sens.append(dict(sep=sep, dm=dm, ddchi2=[q[0] for q in ds], match=[q[1] for q in ds], f=[q[2] for q in ds]))
    nd = np.array([n['dchi2'] for n in null])
    thr = float(np.max(nd)) if len(nd) else np.nan
    print(f'   target: dchi2 = {tr["dchi2"]:.1f} (null comps: median {np.median(nd):.1f}, max {np.max(nd):.1f}, n={len(nd)}), '
          f'best sep {tr["sep"]:.2f}" PA {tr["pa"]:.0f} f {tr["f"]:.3f} broadening ({tr["sx"]:.2f},{tr["sy"]:.2f}) px; chi2_single/npix target {tr["chi2_single"]/tr["npix"]:.2f}, '
          f'comps median {np.median([n["chi2_single"]/n["npix"] for n in null]):.2f}', flush=True)
    for s_ in sens:
        print(f'      inj sep {s_["sep"]:.2f}" dm {s_["dm"]:.1f}: median extra dchi2 {np.median(s_["ddchi2"]):8.1f}  recovered at injected position: {np.mean(s_["match"]):.2f}  median f {np.median(s_["f"]):.3f} (true {10**(-0.4*s_["dm"]):.3f})')
    return dict(file=filt, mjd=mjd, fwhm_med=float(np.median(fwhm)), qafwhm=hdr.get('QAFWHM'), n_comp=len(stamps), target=tr,
                null=null, null_max=thr, sens=sens, scale=sc_px)


if __name__ == '__main__':
    files = sys.argv[1:] or [f for f in sorted(glob.glob(f'{IMG}/sm_*.fits')) if '20160926' not in f]
    out = []
    for f in files:
        r = run(f)
        if r:
            out.append(r)
    to_json(out, f'{R}/results/x20_skymapper.json')
