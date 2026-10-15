"""Own aperture photometry on the full 21x21 TESScut cutouts (all six sectors; S1/S28/S68 serve as a check of the
method against TESS-SPOC PDCSAP, S95/S105/S106 are the science light curves for those sectors).
 * time: TESScut TIME - TIMECORR(FFI) + target Roemer delay (a01_timefix.py, validated vs SPOC to ~1 ms)
 * sky: per-cadence robust plane fit to pixels free of stars (star mask from the median image, dilated, plus
   everything within 4 px of the target)
 * position: global (dx,dy) offset of the WCS fitted with the TESS PRF to the median sky-subtracted image, using all
   Gaia DR3 stars in the stamp with free fluxes for the 6 brightest and Gaia-predicted fluxes for the rest
 * apertures: pixel-centre circles r = 1.0..3.0 px; chosen by minimum 1-h out-of-eclipse scatter
 * contamination: fraction of aperture flux from Gaia neighbours, from the fitted PRF scene (target = Gaia source
   6571343840265900544; any unresolved companion is part of the 'target' flux by construction)
Output: data/tc_s{sec}.npz and results/phot.json"""
import numpy as np, json, glob, sys
from astropy.io import fits
from astropy.wcs import WCS
import pandas as pd
from scipy.optimize import least_squares
import PRF
V = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v5'
LOC = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/data/tesscut'
PRFDIR = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/prf_files'
GAIA = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/data/gaia_240as.csv'
P = 1.3036712; T0 = 3165.6697; T14 = 1.55 / 24
TID = 6571343840265900544

def ecl_mask(t, w=0.8):
    ph = ((t - T0) / P) % 1.0
    d1 = np.minimum(ph, 1 - ph) * P; d2 = np.abs(ph - 0.5) * P
    return (d1 < w * T14) | (d2 < w * T14)

def tmag(g, c):   # Stassun et al. 2019 TIC-8 relation G, Bp-Rp -> T
    return g - 0.00522555 * c ** 3 + 0.0891337 * c ** 2 - 0.633923 * c + 0.0324473

def scat1h(t, f, m):
    """robust scatter of 1-h bins of OOE flux after a 0.5-d running median"""
    tt, ff = t[~m], f[~m]
    rm = np.array([np.median(ff[np.abs(tt - x) < 0.25]) for x in tt])
    r = ff / rm - 1
    nb = np.floor((tt - tt[0]) * 24).astype(int)
    cnt = np.bincount(nb); s = np.bincount(nb, r)
    ok = cnt >= 0.6 * np.median(cnt[cnt > 0])
    b = s[ok] / cnt[ok]
    return 1.4826 * np.median(np.abs(b - np.median(b)))

def run(fn):
    h = fits.open(fn, memmap=True); sec = int(h[0].header['SECTOR'])
    cam, ccd = int(h[0].header['CAMERA']), int(h[0].header['CCD'])
    d = h[1].data
    t = np.array(d['TIME'], float); q = np.array(d['QUALITY']); F = np.array(d['FLUX'], float)
    tc = np.load(f'{V}/data/timecorr_tc_s{sec}.npy')
    ok = np.isfinite(t) & (t > 0) & (q == 0) & np.all(np.isfinite(F.reshape(len(t), -1)), axis=1)
    t, F = t[ok], F[ok]
    t = t + np.interp(t, tc[0], tc[1]) / 86400.0
    ny, nx = F.shape[1:]
    yy, xx = np.mgrid[:ny, :nx]
    w = WCS(h[2].header)
    g = pd.read_csv(GAIA); g['bp_rp'] = g['bp_rp'].fillna(1.2); g['T'] = tmag(g.phot_g_mean_mag, g.bp_rp)
    ep = 2000 + (np.median(t) + 2457000 - 2451545.0) / 365.25
    dty = ep - 2016.0
    ra = g.ra + g.pmra.fillna(0) / 3.6e6 / np.cos(np.deg2rad(g.dec)) * dty; de = g.dec + g.pmdec.fillna(0) / 3.6e6 * dty
    gx, gy = w.all_world2pix(ra.values, de.values, 0)
    itgt = int(np.where(g.source_id.values == TID)[0][0])
    # ---- star mask and sky
    med0 = np.median(F, axis=0)
    lo = np.percentile(med0, 30); rs = 1.4826 * np.median(np.abs(med0[med0 < np.percentile(med0, 60)] - lo))
    star = med0 > lo + 4 * rs
    from scipy.ndimage import binary_dilation
    star = binary_dilation(star, iterations=1)
    star |= np.hypot(xx - gx[itgt], yy - gy[itgt]) < 4.0
    skypix = ~star
    A = np.vstack([np.ones(skypix.sum()), xx[skypix] - nx / 2, yy[skypix] - ny / 2]).T
    sky = np.empty((len(t), 3))
    for i in range(len(t)):
        v = F[i][skypix]; m = np.ones(v.size, bool)
        for _ in range(3):
            c, *_ = np.linalg.lstsq(A[m], v[m], rcond=None)
            r = v - A @ c; s = 1.4826 * np.median(np.abs(r[m]))
            m = np.abs(r) < 3 * s
        sky[i] = c
    skyimg = sky[:, 0, None, None] + sky[:, 1, None, None] * (xx - nx / 2) + sky[:, 2, None, None] * (yy - ny / 2)
    img = F - skyimg
    med = np.median(img, axis=0)
    # ---- PRF scene fit to the median image
    colc = h[1].header['1CRV4P'] + gx[itgt] - (h[1].header['1CRPX4'] - 1)
    rowc = h[1].header['2CRV4P'] + gy[itgt] - (h[1].header['2CRPX4'] - 1)
    prf = PRF.TESS_PRF(cam, ccd, sec, colc, rowc, localdatadir=f'{PRFDIR}/' + ('s0001' if sec < 4 else 's0004'))
    inside = (gx > -0.5) & (gx < nx - 0.5) & (gy > -0.5) & (gy < ny - 0.5)
    idx = np.where(inside)[0]
    order = idx[np.argsort(g['T'].values[idx])]
    free = list(order[:6]); fixed = list(order[6:])
    if itgt not in free: free = [itgt] + free[:5]
    def scene(par, return_parts=False):
        dx, dy = par[0], par[1]; fl = par[2:2 + len(free)]; scale = par[2 + len(free)]
        mod = np.zeros((ny, nx)); parts = {}
        for j, i in enumerate(free):
            pm = prf.locate(gx[i] + dx, gy[i] + dy, (ny, nx)); mod += fl[j] * pm; parts[i] = fl[j] * pm
        for i in fixed:
            pm = prf.locate(gx[i] + dx, gy[i] + dy, (ny, nx)); fi = scale * 10 ** (-0.4 * (g['T'].values[i] - g['T'].values[itgt])); mod += fi * pm; parts[i] = fi * pm
        return (mod, parts) if return_parts else mod
    f0 = med.sum()
    p0 = [0.0, 0.0] + [f0 * 10 ** (-0.4 * (g['T'].values[i] - g['T'].values[itgt])) for i in free] + [f0]
    wts = 1 / np.sqrt(np.abs(med) + 50.0)
    res = least_squares(lambda p: ((scene(p) - med) * wts).ravel(), p0, x_scale='jac')
    mod, parts = scene(res.x, True)
    dx, dy = res.x[:2]
    xt, yt = gx[itgt] + dx, gy[itgt] + dy
    tgtimg = parts[itgt]; totimg = sum(parts.values())
    resid_frac = float(np.abs(med - mod)[np.hypot(xx - xt, yy - yt) < 3].sum() / med[np.hypot(xx - xt, yy - yt) < 3].sum())
    # ---- apertures
    m_ecl = ecl_mask(t)
    rr = np.hypot(xx - xt, yy - yt)
    best = None; trials = []
    for rad in (1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 3.0):
        ap = rr <= rad
        fl = img[:, ap].sum(1); f = fl / np.median(fl[~m_ecl])
        s = scat1h(t, f, m_ecl)
        contam = float(1 - tgtimg[ap].sum() / totimg[ap].sum())
        trials.append(dict(r=rad, npix=int(ap.sum()), scat1h_ppm=s * 1e6, contam=contam))
        if best is None or s < best[1]: best = (rad, s, ap, f, fl, contam)
    rad, s, ap, f, fl, contam = best
    fcorr = (f - contam) / (1 - contam)
    cx = (img * xx).reshape(len(t), -1)[:, ap.ravel()].sum(1) / fl
    cy = (img * yy).reshape(len(t), -1)[:, ap.ravel()].sum(1) / fl
    texp = float(np.median(np.diff(t)))
    np.savez(f'{V}/data/tc_s{sec}.npz', t=t, f=fcorr, f_raw=f, flux=fl, sky0=sky[:, 0], skyap=skyimg[:, ap].sum(1), cx=cx, cy=cy,
             ap=ap, contam=contam, texp=texp, xt=xt, yt=yt)
    out = dict(sector=sec, cam=cam, ccd=ccd, n=int(len(t)), texp_s=texp * 86400, wcs_offset_px=[float(dx), float(dy)],
               target_px=[float(xt), float(yt)], prf_resid_frac_r3=resid_frac, chosen_r=rad, contam=contam,
               scat1h_ppm=s * 1e6, sky_median_e_s_pix=float(np.median(sky[:, 0])), sky_in_ap_frac=float(np.median(skyimg[:, ap].sum(1) / fl)),
               trials=trials, n_free=len(free), target_flux_fit=float(res.x[2 + free.index(itgt)]))
    print(json.dumps({k: v for k, v in out.items() if k != 'trials'}), flush=True)
    return out

if __name__ == '__main__':
    allres = {}
    for fn in sorted(glob.glob(f'{LOC}/s*/*.fits')):
        o = run(fn); allres[o['sector']] = o
    json.dump(allres, open(f'{V}/results/phot.json', 'w'), indent=1)
