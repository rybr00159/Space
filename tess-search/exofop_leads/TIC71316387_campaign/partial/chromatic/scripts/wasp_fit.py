"""SuperWASP DR1 (WASP-South, 2006-2008) eclipse depths for 1SWASP J220333.07-425007.9.

Cleaning: cameras 227/228 only (221 has 11 points); FLAG==32 (SysRem/'TAMUZ'-corrected points;
FLAG==0 points have TAMFLUX2==FLUX2, i.e. uncorrected, 3% of data) ; positive finite fluxes;
segments = night x camera with >= minpts points; iterative 4-sigma clipping on residuals of the
full model; per-segment robust weights.
Variants: TAMFLUX2 (default) vs raw FLUX2; baseline order 0/1/2; camera 227 vs 228 alone;
season 2006 vs 2007.
"""
import sys, json
import numpy as np
from astropy.io import fits
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/scripts')
import chromlib as C
import survey_fit as SF

OUT = SF.OUT
EXPT = 30.0 / 86400   # WASP exposure 30 s


def load(fn, col='TAMFLUX2', flag=32, cams=(227, 228)):
    h = fits.open(fn)
    d = h[1].data
    cam = np.array([s[:3] for s in d['IMAGEID']]).astype(int)
    hjd = d['TMID'] / 86400 + float(h[0].header['JD_REF'])
    f = np.array(d[col], float)
    ok = np.isfinite(f) & (f > 0) & np.isin(cam, cams)
    if flag is not None:
        ok &= d['FLAG'] == flag
    t = C.hjd_utc_to_bjd_tdb(hjd[ok])
    f = f[ok]
    # normalise each night x camera segment by its median -> fractional flux (baselines still free)
    seg = C.segments(t, 0.3, cam[ok])
    med = np.array([np.median(f[seg == k]) for k in range(seg.max() + 1)])
    f = f / med[seg]
    return t, f, cam[ok], hjd[ok], np.array(d['CCDX'][ok]), np.array(d['CCDY'][ok])


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--file', default=OUT + '/data/wasp/1SWASP_J220333.07-425007.9.fits')
    ap.add_argument('--tag', default='wasp_target')
    ap.add_argument('--nfake', type=int, default=400)
    ap.add_argument('--ninj', type=int, default=150)
    ap.add_argument('--variants', action='store_true')
    a = ap.parse_args()

    t, f, cam, hjd, cx, cy = load(a.file)
    print('BJD-HJD offset (s): median %.1f' % (np.median(t - hjd) * 86400))
    seg = C.segments(t, 0.3, cam)
    results = {}
    res, (tt, ff, ww, sg) = SF.run(a.tag + '_TAM_o1', t, f, seg, EXPT, order=1, nfake=a.nfake, ninj=a.ninj)
    results['default'] = res
    print('%s: n=%d nseg=%d in-deep=%d (%d events) in-shallow=%d (%d events); deep %.0f +- %.0f (null %.0f, formal %.0f) ppm, p_null=%.4f; shallow %.0f +- %.0f ppm p_null=%.4f' % (
        res['name'], res['n'], res['nseg'], res['n_in_deep'], res['n_deep_events'], res['n_in_shallow'], res['n_shallow_events'],
        res['deep'] * 1e6, res['deep_err_null'] * 1e6, res['null_deep_mean'] * 1e6, res['deep_err_formal'] * 1e6, res['deep_null_p'],
        res['shallow'] * 1e6, res['shallow_err_null'] * 1e6, res['shallow_null_p']), flush=True)
    for k, v in res['injection'].items():
        print('   inj %s: deep rec %.0f +- %.0f (ratio %.3f); shallow rec %.0f +- %.0f (ratio %.3f)' % (
            k, v['deep_rec_mean'] * 1e6, v['deep_rec_sd'] * 1e6, v['deep_ratio'], v['shallow_rec_mean'] * 1e6, v['shallow_rec_sd'] * 1e6, v['shallow_ratio']), flush=True)
    # ephemeris check: phase scan
    scan = SF.phase_scan(tt, ff, ww, sg, EXPT, 1, np.linspace(-0.08, 0.08, 161))
    results['phase_scan'] = scan.tolist()
    ib = np.argmax(scan[:, 1])
    print('   phase scan: max deep amplitude %.0f ppm at dphase %.4f (%.1f min)' % (scan[ib, 1] * 1e6, scan[ib, 0], scan[ib, 0] * C.P * 1440))
    # chi2 minimum on fine grid -> timing offset with error from delta-chi2 = 1 (formal, weights rescaled)
    chi = scan[:, 3]
    results['phase_scan_chi2_min_dphase'] = float(scan[np.argmin(chi), 0])
    np.savez(OUT + f'/results/{a.tag}_cleaned.npz', t=tt, f=ff, w=ww, seg=sg)

    if a.variants:
        var = {}
        for name, kw in [('TAM_o0', dict(order=0)), ('TAM_o2', dict(order=2))]:
            r, _ = SF.run(a.tag + '_' + name, t, f, seg, EXPT, nfake=200, ninj=60, **kw)
            var[name] = r
        # raw flux
        t2, f2, cam2, *_ = load(a.file, col='FLUX2')
        seg2 = C.segments(t2, 0.3, cam2)
        for name, kw in [('RAW_o1', dict(order=1)), ('RAW_o2', dict(order=2))]:
            r, _ = SF.run(a.tag + '_' + name, t2, f2, seg2, EXPT, nfake=200, ninj=60, **kw)
            var[name] = r
        # cameras separately
        for c in (227, 228):
            m = cam == c
            r, _ = SF.run(a.tag + f'_cam{c}', t[m], f[m], seg[m], EXPT, nfake=200, ninj=60)
            var[f'cam{c}'] = r
        # seasons
        for yr, (lo, hi) in {'2006': (2453700, 2454100), '2007_08': (2454100, 2454700)}.items():
            m = (t > lo) & (t < hi)
            r, _ = SF.run(a.tag + f'_season{yr}', t[m], f[m], seg[m], EXPT, nfake=200, ninj=60)
            var['season' + yr] = r
        # clipping threshold
        r, _ = SF.run(a.tag + '_clip3', t, f, seg, EXPT, nfake=200, ninj=60, nclip=3.0)
        var['clip3'] = r
        results['variants'] = var
        for k, r in var.items():
            print('  variant %-10s n=%5d deep %6.0f +- %4.0f (null mean %5.0f)  shallow %6.0f +- %4.0f  inj6000 ratio %.3f' % (
                k, r['n'], r['deep'] * 1e6, r['deep_err_null'] * 1e6, r['null_deep_mean'] * 1e6, r['shallow'] * 1e6,
                r['shallow_err_null'] * 1e6, r['injection']['0.006']['deep_ratio']), flush=True)
    # strip big null arrays from variants for size
    json.dump(results, open(OUT + f'/results/{a.tag}.json', 'w'), indent=1, default=float)
