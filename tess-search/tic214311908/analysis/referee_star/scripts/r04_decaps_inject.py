"""Referee: honest injection-recovery for the DECaPS2 companion limits around the M3V.
The investigator (s13c) measured matched-filter fluxes in (image - crowdsource model) and 'checked' them by injecting
point sources INTO THE RESIDUAL and re-measuring with the same linear filter, which returns 1.0 by construction and
cannot test whether a companion would have been absorbed by the model of the M3V.
Here: empirical PSF from isolated field stars (Gaia positions propagated to 2018.7) in the same DECaPS2 cutout;
the M3V is fitted with a single PSF (flux, x, y, background plane free) in a 4.5" radius box -> residual.
Companions are injected INTO THE IMAGE (PSF x flux at separation s, 8 PAs) and the single-star fit is repeated,
so any absorption into the M3V fit is included.  Detection statistic: matched-filter flux at the injected position
in the new residual, in units of the robust scatter of the same statistic over 16 PAs at that separation in the
uninjected residual (same criterion as s13c).  5-sigma detection threshold; 90% of PAs detected -> limit.
Also a check of saturation: peak of the M3V vs PSF x catalogue flux.
Input data (downloaded by the investigator, used as data only): star/data/img/decaps/cut_decaps2.fits (i = plane 2, z = plane 3)
Output: out/r04_decaps_inject.json"""
import json, warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import shift as ndshift
from scipy.optimize import least_squares
R = 'tic214311908/analysis/referee_star'
IMG = 'tic214311908/analysis/star/data/img/decaps/cut_decaps2.fits'
h = fits.open(IMG)[0]; w = WCS(h.header).celestial; scale = 0.262
g = pd.read_csv(f'{R}/data/gaia_cone45.csv')
ep = 2018.7; dt = ep - 2016.0
gra = g.ra + np.nan_to_num(g.pmra) * dt / 3.6e6 / np.cos(np.radians(g.dec)); gde = g.dec + np.nan_to_num(g.pmdec) * dt / 3.6e6
gx, gy = w.all_world2pix(gra.values, gde.values, 0)
i0 = int(np.where(g.source_id == 5925373272850784768)[0][0]); x0, y0 = gx[i0], gy[i0]
G = g.phot_g_mean_mag.values
rng = np.random.default_rng(5)
out = {}
for k, band, mag0 in ((2, 'i', 13.355), (3, 'z', 12.761)):
    img = h.data[k].astype(float); ny, nx = img.shape
    F0 = 10 ** ((22.5 - mag0) / 2.5)
    # --- PSF from isolated stars (Gaia), stamps 31x31 px (8.1")
    hb = 15; stamps = []
    for j in range(len(g)):
        if not (15.0 < G[j] < 18.5): continue
        xj, yj = gx[j], gy[j]
        if not (hb + 2 < xj < nx - hb - 3 and hb + 2 < yj < ny - hb - 3): continue
        dd = np.hypot(gx - xj, gy - yj) * scale
        if np.any((dd > 0.05) & (dd < 3.5) & (G < G[j] + 3)) or np.hypot(xj - x0, yj - y0) * scale < 10: continue
        xi, yi = int(round(xj)), int(round(yj))
        st = img[yi - hb:yi + hb + 1, xi - hb:xi + hb + 1].copy()
        bg = np.median(np.r_[st[0], st[-1], st[:, 0], st[:, -1]]); st -= bg
        st = ndshift(st, (yi - yj, xi - xj), order=3)
        if st.sum() > 0: stamps.append(st / st[hb - 6:hb + 7, hb - 6:hb + 7].sum())
    psf = np.median(stamps, axis=0)
    psf /= psf[hb - 6:hb + 7, hb - 6:hb + 7].sum()   # normalised inside 13x13 px (3.4") box
    yy, xx = np.mgrid[0:2 * hb + 1, 0:2 * hb + 1]
    def P(dx, dy):  # PSF centred at (hb+dx, hb+dy)
        return ndshift(psf, (dy, dx), order=3, mode='constant')
    # --- fit region around the M3V
    rb = 17  # 4.5" half-box
    xi0, yi0 = int(round(x0)), int(round(y0))
    cut0 = img[yi0 - rb:yi0 + rb + 1, xi0 - rb:xi0 + rb + 1].copy()
    cy_, cx_ = np.mgrid[-rb:rb + 1, -rb:rb + 1]
    def model_star(p, shape_only=False):
        f, dx, dy, b0, bx, by = p
        st = np.zeros((2 * rb + 1, 2 * rb + 1))
        pp = P(dx + (x0 - xi0), dy + (y0 - yi0))      # 31x31 centred on hb
        st[rb - hb:rb + hb + 1, rb - hb:rb + hb + 1] = pp
        return f * st + b0 + bx * cx_ + by * cy_
    def fit(cut):
        p0 = [cut.sum(), 0, 0, 0, 0, 0]
        r = least_squares(lambda p: (cut - model_star(p)).ravel(), p0, x_scale=[cut.sum(), 0.2, 0.2, 1, 0.1, 0.1])
        return r.x, cut - model_star(r.x)
    p_star, res0 = fit(cut0)
    peak_obs = cut0.max(); peak_mod = model_star(p_star).max()
    sat_info = dict(peak_obs=float(peak_obs), peak_model=float(peak_mod), fitted_flux=float(p_star[0]), catalogue_flux=float(F0),
                    centre_resid_frac=float(res0[rb - 2:rb + 3, rb - 2:rb + 3].sum() / p_star[0]))
    def mf(resid, sx, sy):
        # matched filter at offset (sx, sy) px from the M3V position
        tpl = np.zeros_like(resid); pp = P(sx + (x0 - xi0), sy + (y0 - yi0))
        tpl[rb - hb:rb + hb + 1, rb - hb:rb + hb + 1] = pp
        return float(np.sum(resid * tpl) / np.sum(tpl * tpl))
    curve = []
    for s in (0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0):
        sp = s / scale
        th16 = np.arange(16) * np.pi / 8
        base = np.array([mf(res0, -sp * np.sin(t), sp * np.cos(t)) for t in th16])
        # exclude PAs within 1.2" of a Gaia source (other than the M3V)
        okpa = []
        for t in th16:
            px, py = x0 - sp * np.sin(t), y0 + sp * np.cos(t)
            dd = np.hypot(gx - px, gy - py) * scale; dd[i0] = 99
            okpa.append(dd.min() > 1.2)
        okpa = np.array(okpa); bb = base[okpa]
        sd = 1.4826 * np.median(np.abs(bb - np.median(bb))); med = np.median(bb)
        thr_naive = med + 5 * sd    # investigator-style 5 sigma threshold
        res_dm = {}
        for dm in np.arange(3.0, 7.01, 0.25):
            Fc = F0 * 10 ** (-0.4 * dm)
            det = []; recf = []
            for t in th16[okpa][::2]:
                sx, sy = -sp * np.sin(t), sp * np.cos(t)
                comp = np.zeros_like(cut0); pp = P(sx + (x0 - xi0), sy + (y0 - yi0)) * Fc
                comp[rb - hb:rb + hb + 1, rb - hb:rb + hb + 1] = pp
                pr, rr = fit(cut0 + comp)
                v = mf(rr, sx, sy); recf.append((v - mf(res0, sx, sy)) / Fc)
                det.append(v > thr_naive)
            res_dm[f'{dm:.2f}'] = dict(frac_det=float(np.mean(det)), recovered_fraction_med=float(np.median(recf)))
        dms = sorted(res_dm, key=float)
        lim = None
        for d_ in dms:
            if res_dm[d_]['frac_det'] >= 0.9: lim = float(d_)
            else: break
        # naive (no-absorption) limit like s13c
        naive = float(-2.5 * np.log10(max(thr_naive, 1e-6) / F0))
        curve.append(dict(sep=s, n_pa=int(okpa.sum()), sd=float(sd), naive_5sig_dm=naive, injected_90pct_dm=lim,
                          absorbed_fraction_at_dm4=1 - res_dm['4.00']['recovered_fraction_med'], detail=res_dm))
        print(f"{band} {s:4.2f}\": naive 5-sigma dm {naive:.2f}; image-injection 90% limit dm {lim}; recovered flux fraction at dm=4: {res_dm['4.00']['recovered_fraction_med']:.2f}", flush=True)
    out[band] = dict(n_psf_stars=len(stamps), saturation_check=sat_info, curve=curve)
    print(band, 'PSF stars', len(stamps), 'sat check', {k_: round(v_, 3) for k_, v_ in sat_info.items()})
# combine with the investigator's 5-sigma residual limits (s13c), corrected for the fraction of a companion's flux
# that survives the single-star fit (median recovered fraction over dm 3.5-5 and PAs)
inv = json.load(open('tic214311908/analysis/star/out/decaps_contrast.json'))
for band in ('i', 'z'):
    for c in out[band]['curve']:
        fr = np.median([c['detail'][d]['recovered_fraction_med'] for d in c['detail'] if 3.5 <= float(d) <= 5.0])
        c['recovered_fraction'] = float(fr)
        iv = [r for r in inv[band]['curve'] if abs(r['sep'] - c['sep']) < 1e-6 and 'delta_mag_5sig' in r]
        if iv:
            c['investigator_dm5'] = iv[0]['delta_mag_5sig']
            c['corrected_dm5'] = float(iv[0]['delta_mag_5sig'] - 2.5 * np.log10(1 / fr)) if fr > 0.02 else None
        print(band, c['sep'], 'recovered fraction', round(fr, 2), 'investigator dm5', c.get('investigator_dm5'), '-> corrected', c.get('corrected_dm5'))
json.dump(out, open(f'{R}/out/r04_decaps_inject.json', 'w'), indent=1)
