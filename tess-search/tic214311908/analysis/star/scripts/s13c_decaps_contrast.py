"""Contrast limits around TIC 214311908 from the DECaPS2 coadds (Saydjari+23; epoch ~2018.7; 0.262"/px; cutouts
from the Legacy Survey viewer, layers 'decaps2' (image) and 'decaps2-model' (crowdsource model of all catalogued
sources)). Sources that crowdsource detected are in the DECaPS2 catalogue (s05/s08). Here we test for sources it did
NOT model: in the residual image (image - model), PSF-matched-filter fluxes at separations 0.5-6" and 16 position
angles around the M3V (2018.7 position); the 5-sigma limit at each separation is 5 x the robust scatter over
position angles (positions within 1.2" of a Gaia source or a reliable DECaPS2 source (ndet_ok >= 2, finite i) excluded;
single-detection DECaPS2 entries around the M3V are treated as artefacts). Separations < 1" are not reliable (the
crowdsource model of the M3V can absorb a companion there). Empirical PSF from isolated
field stars (16 < i < 18.5) in the same cutout. Injection-recovery of point sources into the residual image checks
the calibration. Limits are expressed as delta-i and delta-z relative to the M3V and as Delta T (TESS ~ i/z).
Output: out/decaps_contrast.json, figs/decaps_resid.png"""
import json, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import shift as ndshift
B = 'tic214311908/analysis/star'
D = f'{B}/data/img/decaps'
a = fits.open(f'{D}/cut_decaps2.fits')[0]; m = fits.open(f'{D}/cut_decaps2-model.fits')[0]
w = WCS(a.header).celestial; scale = 0.262
ep = 2018.7; dt = ep - 2016.0
ra0 = 260.969674 + (-113.671152) * dt / 3.6e6 / np.cos(np.radians(-51.323821)); de0 = -51.323821 + 85.54355 * dt / 3.6e6
x0, y0 = [float(v) for v in w.all_world2pix(ra0, de0, 0)]
cat = pd.read_csv(f'{B}/data/decaps2_90.csv', usecols=['ra', 'dec', 'mean_mag_i', 'mean_mag_z', 'ndet_ok'])
cat = cat[(cat.ndet_ok >= 2) & np.isfinite(cat.mean_mag_i)].reset_index(drop=True)   # reliable catalogue sources only
cx, cy = w.all_world2pix(cat.ra.values, cat.dec.values, 0)
g = pd.read_csv(f'{B}/data/gaia_120.csv')
gx, gy = w.all_world2pix(g.ra.values + np.nan_to_num(g.pmra.values) * dt / 3.6e6 / np.cos(np.radians(g.dec.values)), g.dec.values + np.nan_to_num(g.pmdec.values) * dt / 3.6e6, 0)
out = {}
rng = np.random.default_rng(2)
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, axs = plt.subplots(1, 3, figsize=(13, 4.3))
for k, band in ((2, 'i'), (3, 'z')):
    img = a.data[k].astype(float); mod = m.data[k].astype(float); res = img - mod
    ny, nx = img.shape; yy, xx = np.mgrid[0:ny, 0:nx]
    # M3V total flux (nmgy) from catalogue magnitude
    F0 = 10 ** ((22.5 - (13.355 if band == 'i' else 12.761)) / 2.5)
    # empirical PSF from isolated stars
    mag = cat.mean_mag_i.values if band == 'i' else cat.mean_mag_z.values
    stamps = []
    for j in np.where((mag > 15.5) & (mag < 18.5))[0]:
        xj, yj = cx[j], cy[j]
        if not (12 < xj < nx - 13 and 12 < yj < ny - 13): continue
        dd = np.hypot(cx - xj, cy - yj); others = (dd > 0.5) & (dd < 4 / scale) & (mag < mag[j] + 3)
        if others.any() or np.hypot(xj - x0, yj - y0) < 8 / scale: continue
        st = img[int(round(yj)) - 10:int(round(yj)) + 11, int(round(xj)) - 10:int(round(xj)) + 11]
        st = ndshift(st, (round(yj) - yj, round(xj) - xj), order=3); st = st - np.median(img[int(round(yj)) - 15:int(round(yj)) + 16, int(round(xj)) - 15:int(round(xj)) + 16])
        if st.sum() > 0: stamps.append(st / st.sum())
    psf = np.median(np.array(stamps), axis=0); psf = np.clip(psf, 0, None); psf /= psf.sum()
    fwhm = 2 * np.sqrt(np.sum(psf > psf.max() / 2) / np.pi) * scale
    def mf_flux(R, xs, ys):
        xi, yi = int(round(xs)), int(round(ys)); P = ndshift(psf, (ys - yi, xs - xi), order=3)
        cut = R[yi - 10:yi + 11, xi - 10:xi + 11]
        return float(np.sum(cut * P) / np.sum(P * P))
    seps = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0]
    res_sep = []
    for s_ in seps:
        fl = []
        for th in np.arange(16) * np.pi / 8:
            xs, ys = x0 - s_ / scale * np.sin(th), y0 + s_ / scale * np.cos(th)   # PA from north through east (x increases to the west in this TAN WCS? checked below)
            dc = np.min(np.hypot(cx - xs, cy - ys)) * scale; dg = np.min(np.hypot(gx - xs, gy - ys)) * scale
            if dc < 1.2 and np.hypot(cx[np.argmin(np.hypot(cx - xs, cy - ys))] - x0, cy[np.argmin(np.hypot(cx - xs, cy - ys))] - y0) * scale > 0.3: continue
            if dg < 1.2 and np.hypot(gx[np.argmin(np.hypot(gx - xs, gy - ys))] - x0, gy[np.argmin(np.hypot(gx - xs, gy - ys))] - y0) * scale > 0.3: continue
            fl.append(mf_flux(res, xs, ys))
        fl = np.array(fl)
        if len(fl) < 6: res_sep.append(dict(sep=s_, n=int(len(fl)))); continue
        sd = 1.4826 * np.median(np.abs(fl - np.median(fl)))
        lim5 = 5 * sd + max(np.median(fl), 0)
        dm = -2.5 * np.log10(lim5 / F0)
        # injection check at this separation (PA 45 deg)
        th = np.pi / 4; xs, ys = x0 - s_ / scale * np.sin(th), y0 + s_ / scale * np.cos(th)
        Fi = lim5; xi, yi = int(round(xs)), int(round(ys)); P = ndshift(psf, (ys - yi, xs - xi), order=3) * Fi
        R2 = res.copy(); R2[yi - 10:yi + 11, xi - 10:xi + 11] += P
        rec = mf_flux(R2, xs, ys) - mf_flux(res, xs, ys)
        res_sep.append(dict(sep=s_, n=int(len(fl)), sd=float(sd), max_abs_snr=float(np.max(np.abs(fl - np.median(fl))) / sd), lim5_flux=float(lim5), delta_mag_5sig=float(dm), inj_recovery=float(rec / Fi)))
    out[band] = dict(psf_fwhm_arcsec=float(fwhm), n_psf_stars=len(stamps), M3V_flux_nmgy=F0, curve=res_sep)
    print(band, f'PSF FWHM {fwhm:.2f}" from {len(stamps)} stars;', ' '.join(f"{r['sep']}\":{r.get('delta_mag_5sig', float('nan')):.1f}(max|S/N| {r.get('max_abs_snr', float('nan')):.1f})" for r in res_sep))
    if band == 'i':
        s0 = int(round(8 / scale))
        for ax, im_, t_ in ((axs[0], img, 'DECaPS2 i image'), (axs[1], res, 'image - crowdsource model')):
            cut = im_[int(round(y0)) - s0:int(round(y0)) + s0 + 1, int(round(x0)) - s0:int(round(x0)) + s0 + 1]
            v = np.percentile(np.abs(cut), 99) if 'model' in t_ else np.percentile(cut, 99.5)
            ax.imshow(cut, origin='lower', cmap='gray_r' if 'image' == t_[-5:] else 'RdBu_r', vmin=-v if 'model' in t_ else np.percentile(cut, 5), vmax=v,
                      extent=(-s0 * scale, s0 * scale, -s0 * scale, s0 * scale))
            ax.set_title(t_ + ' (pixel axes, 16" box)', fontsize=9)
ax = axs[2]
for band, col in (('i', 'C0'), ('z', 'C3')):
    c = [r for r in out[band]['curve'] if 'delta_mag_5sig' in r]
    ax.plot([r['sep'] for r in c], [r['delta_mag_5sig'] for r in c], 'o-', color=col, label=f'DECaPS2 {band} (5 sigma, unmodelled sources)')
ax.axhline(2.5 * np.log10(1 / 0.0165), color='k', ls='--', label='max delta-T for a 100% eclipse')
ax.axhline(2.5 * np.log10(0.5 / 0.0165), color='k', ls=':', label='50% eclipse')
ax.invert_yaxis(); ax.set_xlabel('separation (arcsec)'); ax.set_ylabel('contrast (mag)'); ax.legend(fontsize=7); ax.set_ylim(9, 0)
fig.tight_layout(); fig.savefig(f'{B}/figs/decaps_resid.png', dpi=110)
json.dump(out, open(f'{B}/out/decaps_contrast.json', 'w'), indent=1)
