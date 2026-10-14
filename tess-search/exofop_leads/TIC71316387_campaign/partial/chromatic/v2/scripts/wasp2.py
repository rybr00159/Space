"""SuperWASP DR1 (WASP-South cameras 227/228, May 2006 - May 2008) eclipse depths, v2.

File: 1SWASP J220333.07-425007.9 (NASA Exoplanet Archive DR1 FITS; RA_OBJ 330.8878, DEC_OBJ -42.835533,
0.2" from Gaia DR3 6571343840265900544 at the WASP epoch).  Times: TMID/86400 + JD_REF = HJD(UTC) ->
BJD_TDB.  Cleaning: cameras 227/228 (camera 221: 11 points, dropped); finite positive flux; segments
= night x camera, >= minpts points; iterative 4-sigma clipping on full-model residuals; per-segment
robust weights.  Depth model and errors: see lib2.run_survey (fake-phase null with real eclipses
modelled; injection-recovery at fake phases).
Variants: TAMFLUX2 (SysRem-detrended, FLAG==32 points; FLAG==0 points are un-detrended) vs raw FLUX2
(all flags), baseline order 0/1/2, each camera alone, each season alone, 3-sigma clipping.
Control stars: two WASP sources > 140" away (not blended with the target) measured at the same ephemeris.
"""
import sys, json, argparse
import numpy as np
from astropy.io import fits
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/v2/scripts')
import lib2 as L

D = L.CH + '/data/wasp/'
EXPT = 30.0 / 86400
TARGET = D + '1SWASP_J220333.07-425007.9.fits'
CONTROLS = [D + '1SWASP_J220321.80-425123.9.fits', D + '1SWASP_J220322.16-424809.1.fits']
BLENDED = [D + '1SWASP_J220337.10-425021.6.fits', D + '1SWASP_J220328.23-424954.8.fits']


def load(fn, col='TAMFLUX2', flag=32, cams=(227, 228)):
    with fits.open(fn) as h:
        d = h[1].data
        jdref = float(h[0].header['JD_REF'])
        cam = np.array([s[:3] for s in d['IMAGEID']]).astype(int)
        hjd = d['TMID'] / 86400.0 + jdref
        f = np.array(d[col], float)
        ok = np.isfinite(f) & (f > 0) & np.isin(cam, cams)
        if flag is not None:
            ok &= d['FLAG'] == flag
        t = L.C.hjd_utc_to_bjd_tdb(hjd[ok])
        f = f[ok]; cam = cam[ok]
    seg = L.segments(t, 0.3, cam)
    med = np.array([np.median(f[seg == k]) for k in range(seg.max() + 1)])
    return t, f / med[seg], seg, cam, hjd[ok]


def summary(r):
    return ('%-22s n=%5d nseg=%3d in-deep %4d (%2d ev) in-shallow %4d (%2d ev) | deep %6.0f +- %4.0f (formal %4.0f beta %.2f null mean %5.0f) p_up=%.4f | '
            'shallow %6.0f +- %4.0f (null mean %5.0f) | gain(6000) %.3f' % (
                r['name'], r['n'], r['nseg'], r['n_in_deep'], r['n_deep_events'], r['n_in_shallow'], r['n_shallow_events'],
                r['deep'] * 1e6, r['deep_err'] * 1e6, r['deep_err_formal'] * 1e6, r['deep_beta'], r['deep_null_mean'] * 1e6,
                r['deep_p_null_upper'], r['shallow'] * 1e6, r['shallow_err'] * 1e6, r['shallow_null_mean'] * 1e6,
                r['injection']['0.006']['deep_gain']))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--nfake', type=int, default=500)
    ap.add_argument('--ninj', type=int, default=200)
    ap.add_argument('--quick', action='store_true')
    a = ap.parse_args()
    shp = L.load_shapes()
    out = {}
    t, f, seg, cam, hjd = load(TARGET)
    print('BJD_TDB - HJD_UTC: median %.1f s (range %.1f..%.1f)' % (np.median(t - hjd) * 86400, np.min(t - hjd) * 86400, np.max(t - hjd) * 86400))
    print('time span BJD %.2f - %.2f; cycles from T0: %.0f .. %.0f' % (t.min(), t.max(), (t.min() - L.T0) / L.P, (t.max() - L.T0) / L.P))
    r, (tt, ff, ww, sg) = L.run_survey('TAM_o1', t, f, seg, shp, EXPT, 1, a.nfake, a.ninj)
    out['default'] = r; print(summary(r), flush=True)
    np.savez(L.V2 + '/results/wasp2_default_cleaned.npz', t=tt, f=ff, w=ww, seg=sg)
    # phase scan of the deep-eclipse amplitude (ephemeris check in WASP data itself)
    grid = np.linspace(-0.05, 0.05, 201)
    sc = []
    for g in grid:
        aa, ee, rr = L.measure(tt, ff, ww, sg, {'d': dict(shp['d'], phc=shp['d']['phc'] + g), 's': dict(shp['s'], phc=shp['s']['phc'] + g)}, EXPT, 1)
        sc.append((g, aa[0], aa[1], ee[0], np.sum(ww * rr ** 2)))
    sc = np.array(sc)
    out['phase_scan'] = sc.tolist()
    ib = np.argmin(sc[:, 4])
    print('phase scan: chi2 min at dphase %.4f (%.1f min); amp there %.0f ppm; amp at 0: %.0f' % (sc[ib, 0], sc[ib, 0] * L.P * 1440, sc[ib, 1] * 1e6, sc[100, 1] * 1e6))
    if not a.quick:
        nf, ni = 300, 100
        var = {}
        for nm, kw in (('TAM_o0', dict(order=0)), ('TAM_o2', dict(order=2)), ('TAM_clip3', dict(nclip=3.0))):
            rr, _ = L.run_survey(nm, t, f, seg, shp, EXPT, nfake=nf, ninj=ni, **kw); var[nm] = rr; print(summary(rr), flush=True)
        t2, f2, seg2, cam2, _ = load(TARGET, col='FLUX2', flag=None)
        for nm, o in (('RAW_o1', 1), ('RAW_o2', 2)):
            rr, _ = L.run_survey(nm, t2, f2, seg2, shp, EXPT, o, nf, ni); var[nm] = rr; print(summary(rr), flush=True)
        for c in (227, 228):
            m = cam == c
            _, sgc = np.unique(seg[m], return_inverse=True)
            rr, _ = L.run_survey(f'TAM_cam{c}', t[m], f[m], sgc, shp, EXPT, 1, nf, ni); var[f'cam{c}'] = rr; print(summary(rr), flush=True)
        for nm, (lo, hi) in (('2006', (2453800, 2454100)), ('2007', (2454100, 2454450)), ('2008', (2454450, 2454700))):
            m = (t > lo) & (t < hi)
            if m.sum() < 500:
                continue
            _, sgc = np.unique(seg[m], return_inverse=True)
            rr, _ = L.run_survey(f'TAM_season{nm}', t[m], f[m], sgc, shp, EXPT, 1, nf, ni); var['season' + nm] = rr; print(summary(rr), flush=True)
        out['variants'] = var
        ctl = {}
        for fn in CONTROLS + BLENDED:
            tc, fc, sgc, _, _ = load(fn)
            nm = fn.split('/')[-1][:-5]
            rr, _ = L.run_survey(nm, tc, fc, sgc, shp, EXPT, 1, nf, 60); ctl[nm] = rr; print(summary(rr), flush=True)
        out['controls'] = ctl
    for k in list(out.get('variants', {})) :
        out['variants'][k].pop('null_samples', None)
    for k in list(out.get('controls', {})):
        out['controls'][k].pop('null_samples', None)
    json.dump(out, open(L.V2 + '/results/wasp2.json', 'w'), indent=1)
