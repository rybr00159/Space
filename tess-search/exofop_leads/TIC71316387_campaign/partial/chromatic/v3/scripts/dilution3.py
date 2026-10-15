"""Dilution of the target's eclipse depth by resolved Gaia neighbours inside each survey aperture.

  * Gaia DR3 sources within 360" (data/gaia_360as.csv), positions propagated to the survey epoch.
  * Band fluxes relative to the target: G-band flux ratio x synthetic colour term (WASP-G or V-G) from
    the warped BT-Settl spectra at the Teff implied by each source's BP-RP (Mamajek sequence); sources
    without BP-RP get the median colour of the field.
  * PSF: circular Gaussian; aperture: circle of radius R centred on the target catalogue position.
  * SuperWASP: R and sigma calibrated on the WASP catalogue magnitudes of the two blended neighbours
    1SWASP J220337.10-425021.6 (46.7") and J220328.23-424954.8 (55.4"): their WASP_MAG are dominated
    by the target's PSF wing, so (WASP_MAG_nb - WASP_MAG_tgt) measures the PSF x aperture.
  Dilution D = sum_nb F_nb phi(d_nb) / (F_tgt phi(0));  true depth = observed depth x (1 + D).
"""
import sys, json
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/v3/scripts')
import lib3 as L
import sed3 as S

TID = 6571343840265900544
g = pd.read_csv(L.CH + '/data/gaia_360as.csv')
tab = S.mamajek()
keys = [k for k in tab if np.isfinite(tab[k].get('Bp-Rp', np.nan)) and 2900 < tab[k]['Teff'] < 7500]
BPRP = np.array([tab[k]['Bp-Rp'] for k in keys]); TE = np.array([tab[k]['Teff'] for k in keys])
o = np.argsort(BPRP); BPRP, TE = BPRP[o], TE[o]
_colcache = {}


def colour_term(bprp, band):
    """2.5log10(F_G/F_band) relative scale: returns m_band - m_G up to a constant (same for all stars)"""
    te = float(np.interp(np.clip(bprp, BPRP.min(), BPRP.max()), BPRP, TE))
    key = (round(te, -1), band)
    if key not in _colcache:
        m = S.star_model(teff=key[0], tab=tab)
        _colcache[key] = -2.5 * np.log10(S.band_flux(m['spec'], band) / S.band_flux(m['spec'], 'G'))
    return _colcache[key]


def phi(d, R, sig, n=400):
    """fraction of a circular Gaussian (sigma sig) centred at distance d that falls inside a circle of radius R"""
    r = (np.arange(n) + 0.5) / n * R
    dr = R / n
    from scipy.special import i0e
    # radial integral of the offset Gaussian: int_0^R r/s^2 exp(-(r^2+d^2)/2s^2) I0(r d/s^2) dr
    x = r * d / sig ** 2
    val = r / sig ** 2 * np.exp(-(r - d) ** 2 / (2 * sig ** 2)) * i0e(x)
    return float(np.sum(val) * dr)


def field(epoch):
    dt = epoch - 2016.0
    ra = g.ra + g.pmra.fillna(0) / 3.6e6 / np.cos(np.radians(g.dec)) * dt
    de = g.dec + g.pmdec.fillna(0) / 3.6e6 * dt
    it = int(np.where(g.source_id.values == TID)[0][0])
    sep = np.hypot((ra - ra[it]) * np.cos(np.radians(de[it])), de - de[it]).values * 3600
    return it, sep


def rel_flux(band):
    c = g.bp_rp.values.copy()
    c[~np.isfinite(c)] = np.nanmedian(c)
    m = g.phot_g_mean_mag.values + np.array([colour_term(x, band) for x in c])
    it = int(np.where(g.source_id.values == TID)[0][0])
    return 10 ** (-0.4 * (m - m[it]))


def dilution(band, R, sig, epoch, centre_sep=None):
    it, sep = field(epoch)
    F = rel_flux(band)
    ph = np.array([phi(s, R, sig) for s in sep])
    contrib = F * ph
    D = (contrib.sum() - contrib[it]) / contrib[it]
    top = np.argsort(-contrib)[:6]
    return float(D), [(int(g.source_id.values[i]), float(sep[i]), float(g.phot_g_mean_mag.values[i]), float(contrib[i] / contrib[it])) for i in top if i != it]


if __name__ == '__main__':
    out = {}
    # ---------------- SuperWASP calibration of R, sigma from blended neighbours
    PIX = 13.7
    wasp_mag = {'tgt': 12.06279, 'J220337.10': 12.66568, 'J220328.23': 13.14546}
    pos = {'J220337.10': (330.904608, -42.839339), 'J220328.23': (330.867631, -42.831897), 'tgt': (330.8878, -42.835533)}
    ep = 2007.3
    it, sep_t = field(ep)
    dt = ep - 2016.0
    ra = (g.ra + g.pmra.fillna(0) / 3.6e6 / np.cos(np.radians(g.dec)) * dt).values
    de = (g.dec + g.pmdec.fillna(0) / 3.6e6 * dt).values
    Fw = rel_flux('WASP')

    def apflux(ra0, de0, R, sig):
        s = np.hypot((ra - ra0) * np.cos(np.radians(de0)), de - de0) * 3600
        return np.sum(Fw * np.array([phi(x, R, sig) if x < R + 6 * sig else 0.0 for x in s]))
    print('SuperWASP aperture/PSF calibration from blended-neighbour catalogue magnitudes')
    cal = {}
    for Rpix in (2.5, 3.0, 3.5, 4.0, 4.5):
        R = Rpix * PIX

        def cost(sig):
            ft = apflux(*pos['tgt'], R, sig)
            c = 0
            for k in ('J220337.10', 'J220328.23'):
                dm_mod = -2.5 * np.log10(apflux(*pos[k], R, sig) / ft)
                c += (dm_mod - (wasp_mag[k] - wasp_mag['tgt'])) ** 2
            return c
        r = minimize_scalar(cost, bounds=(3, 60), method='bounded')
        sig = r.x
        ft = apflux(*pos['tgt'], R, sig)
        dms = {k: float(-2.5 * np.log10(apflux(*pos[k], R, sig) / ft)) for k in ('J220337.10', 'J220328.23')}
        D, top = dilution('WASP', R, sig, ep)
        cal[Rpix] = dict(R_arcsec=R, sigma_arcsec=float(sig), fwhm_pix=float(2.3548 * sig / PIX), rms_mag=float(np.sqrt(r.fun / 2)), dm_model=dms, D=D, top=top)
        print('  R=%.1f px: best sigma %.1f" (FWHM %.2f px), rms mismatch %.3f mag, model dm %s (obs %.3f, %.3f); dilution D = %.4f' % (
            Rpix, sig, 2.3548 * sig / PIX, np.sqrt(r.fun / 2), {k: round(v, 3) for k, v in dms.items()},
            wasp_mag['J220337.10'] - wasp_mag['tgt'], wasp_mag['J220328.23'] - wasp_mag['tgt'], D))
    out['wasp_cal'] = cal
    # adopted: R = 3.5 px (48", the SuperWASP DR1 photometry aperture of Pollacco et al. 2006 / Butters et al. 2010)
    Ds = [cal[k]['D'] for k in cal]
    out['wasp_D'] = dict(adopted=cal[3.5]['D'], min=float(min(Ds)), max=float(max(Ds)),
                         err=float(max(abs(cal[3.5]['D'] - min(Ds)), abs(cal[3.5]['D'] - max(Ds)))))
    print('WASP adopted D = %.4f (range over R 2.5-4.5 px with calibrated PSF: %.4f-%.4f)' % (cal[3.5]['D'], min(Ds), max(Ds)))
    print('  main contributors (R=3.5px):', [(round(s, 1), round(gm, 2), round(c, 4)) for _, s, gm, c in cal[3.5]['top']])
    # ---------------- ASAS-SN (8"/px, 2-px = 16" radius aperture, FWHM ~ 2 px), epoch 2016.5
    res = {}
    for R, fw in ((16.0, 15.0), (16.0, 20.0), (24.0, 16.0), (12.0, 16.0)):
        D, top = dilution('V', R, fw / 2.3548, 2016.5)
        res[f'R{R:.0f}_fwhm{fw:.0f}'] = dict(D=D, top=top)
        print('ASAS-SN V  R=%.0f" FWHM %.0f": D = %.4f; top %s' % (R, fw, D, [(round(s, 1), round(gm, 2), round(c, 4)) for _, s, gm, c in top[:3]]))
    out['asassn_D'] = dict(adopted=res['R16_fwhm15']['D'], err=float(max(abs(v['D'] - res['R16_fwhm15']['D']) for v in res.values())), all=res)
    # ---------------- ASAS-3 (15"/px), apertures of diameter 2-6 px; FWHM ~ 1.5-2.5 px; epoch 2005
    res = {}
    for k, diam in enumerate((2, 3, 4, 5, 6)):
        for fw in (1.5, 2.5):
            D, top = dilution('V', diam / 2 * 15.0, fw * 15.0 / 2.3548, 2005.0)
            res[f'ap{k}_fwhm{fw}'] = D
        print('ASAS-3 MAG_%d (diam %d px): D = %.4f (FWHM 1.5 px) .. %.4f (FWHM 2.5 px)' % (k, diam, res[f'ap{k}_fwhm1.5'], res[f'ap{k}_fwhm2.5']))
    out['asas3_D'] = res
    # ---------------- NEOWISE / KELT not needed (PSF fitting; resolved).  W1 contamination from the 18" star if unresolved:
    it, sep = field(2019.0)
    Fw1 = rel_flux('W1')
    j = np.argsort(sep)[1]
    print('W1 flux of the nearest neighbour (%.1f", G=%.2f) relative to target: %.4f (resolved by the 6" PSF; not a dilution term)' % (
        sep[j], g.phot_g_mean_mag.values[j], Fw1[j]))
    out['W1_nearest_rel'] = float(Fw1[j])
    json.dump(out, open(L.V3 + '/results/dilution3.json', 'w'), indent=1)
