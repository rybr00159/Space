"""v3 SuperWASP DR1 analysis of 1SWASP J220333.07-425007.9 (= TIC 71316387).

File from the NASA Exoplanet Archive SuperWASP DR1 holdings (WASP-South, cameras 221/227/228,
2006-05-07 .. 2008-05-29).  RA_OBJ/DEC_OBJ 330.8878/-42.835533 = 0.17" from the Gaia DR3 position.
Times: TMID/86400 + JD_REF = HJD(UTC) (archive docs) -> BJD_TDB.
Default: TAMFLUX2 (SysRem-corrected), FLAG==32 (FLAG==0 points have TAMFLUX2 == FLUX2, i.e. not
corrected), cameras 227+228 (221 has 11 points); segments = night x camera; linear baseline per
segment fitted simultaneously with the two eclipse amplitudes (TESS template, frozen shape).
"""
import sys, json, argparse
import numpy as np
from astropy.io import fits
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/v3/scripts')
import lib3 as L

D = L.CH + '/data/wasp/'
EXPT = 30.0 / 86400
TARGET = '1SWASP_J220333.07-425007.9'
UNBLENDED = ['1SWASP_J220321.80-425123.9', '1SWASP_J220322.16-424809.1']      # 145", 135" away
BLENDED = ['1SWASP_J220337.10-425021.6', '1SWASP_J220328.23-424954.8']        # 46.7", 55.4" away


def load(name, col='TAMFLUX2', flag=32, cams=(227, 228)):
    with fits.open(D + name + '.fits') as h:
        d = h[1].data
        jdref = float(h[0].header['JD_REF'])
        cam = np.array([s[:3] for s in d['IMAGEID']]).astype(int)
        hjd = d['TMID'] / 86400.0 + jdref
        f = np.array(d[col], float)
        ok = np.isfinite(f) & (f > 0) & np.isin(cam, cams)
        if flag is not None:
            ok &= d['FLAG'] == flag
        hjd, f, cam = hjd[ok], f[ok], cam[ok]
    t = L.hjd_utc_to_bjd_tdb(hjd)
    seg = L.segments(t, 0.25, cam)
    med = np.array([np.median(f[seg == k]) for k in range(seg.max() + 1)])
    return t, f / med[seg], seg, cam, hjd


def summ(nm, r):
    i6 = r['injection'].get('0.006', {})
    return ('%-26s n=%5d nseg=%3d in-d %4d (%2d ev) in-s %4d (%2d ev) | deep %6.0f +- %4.0f (formal %4.0f beta %.2f nullmean %5.0f) p>=%.4f | '
            'shallow %6.0f +- %4.0f (nullmean %5.0f) | gain6000 d %.3f s %.3f' % (
                nm, r['n'], r['nseg'], r['n_in_deep'], r['ev_deep'], r['n_in_shallow'], r['ev_shallow'],
                r['deep'] * 1e6, r['deep_err'] * 1e6, r['deep_formal'] * 1e6, r['deep_beta'], r['deep_null_mean'] * 1e6,
                r['deep_p_ge'], r['shallow'] * 1e6, r['shallow_err'] * 1e6, r['shallow_null_mean'] * 1e6,
                i6.get('deep_gain', np.nan), i6.get('shallow_gain', np.nan)))


def strip(r):
    r = dict(r); r.pop('null_samples_deep', None); r.pop('null_samples_shallow', None); return r


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--nfake', type=int, default=1000)
    ap.add_argument('--ninj', type=int, default=200)
    a = ap.parse_args()
    T = json.load(open(L.V3 + '/results/tess3.json'))
    E3 = json.load(open(L.V3 + '/results/ephem3.json'))
    shp_nom, sod = T['shapes'], T['s_over_d']
    shp = L.survey_shapes(shp_nom, E3)          # refined TESS ephemeris (P' = %.8f)
    out = dict(shapes=shp, shapes_nominal=shp_nom)
    t, f, seg, cam, hjd = load(TARGET)
    print('BJD_TDB-HJD_UTC median %.1f s; span BJD %.2f..%.2f; cycles %.0f..%.0f; ephemeris 1-sigma at mid-span %.2f min' % (
        np.median(t - hjd) * 86400, t.min(), t.max(), (t.min() - L.T0) / L.P, (t.max() - L.T0) / L.P,
        abs((np.median(t) - L.T0) / L.P) * L.P_ERR * 1440))
    r, cl = L.run_survey(t, f, seg, shp, EXPT, 1, a.nfake, a.ninj, s_over_d=sod)
    out['default'] = r
    print(summ('TAM_o1 (default)', r), flush=True)
    np.savez(L.V3 + '/data/wasp3_default_clean.npz', t=cl['t'], f=cl['f'], w=cl['w'], seg=cl['seg'])

    # time-offset scan of both amplitudes in the WASP data alone (ephemeris check)
    tt, ff, ww, pr = cl['t'], cl['f'], cl['w'], cl['proj']
    grid = np.linspace(-0.04, 0.04, 161)
    sc = []
    for g in grid:
        sh = {k: dict(v, phc=v['phc'] + g) for k, v in shp.items()}
        aa, ee, rr, _ = L.fit_amps(tt, ff, ww, None, [L.tmpl(tt, sh, 'd', 0, EXPT), L.tmpl(tt, sh, 's', 0, EXPT)], proj=pr)
        sc.append((g, aa[0], aa[1], ee[0], ee[1], float(np.sum(ww * rr ** 2))))
    sc = np.array(sc)
    out['phase_scan'] = sc.tolist()
    ib = np.argmax(sc[:, 1] / sc[:, 3])
    ic = np.argmin(sc[:, 5])
    print('offset scan: max deep S/N at dphase %+.4f (%+.1f min), amp %.0f; chi2 min at %+.4f (%+.1f min); amp at 0 %.0f' % (
        sc[ib, 0], sc[ib, 0] * L.P * 1440, sc[ib, 1] * 1e6, sc[ic, 0], sc[ic, 0] * L.P * 1440, sc[80, 1] * 1e6))

    nf, ni = 500, 100
    var = {}
    runs = [('TAM_o0', dict(order=0)), ('TAM_o2', dict(order=2)), ('TAM_clip3', dict(nclip=3.0)), ('TAM_clip6', dict(nclip=6.0)),
            ('TAM_minpts40', dict(minpts=40))]
    for nm, kw in runs:
        kw2 = dict(order=1); kw2.update(kw)
        rr, _ = L.run_survey(t, f, seg, shp, EXPT, nfake=nf, ninj=ni, s_over_d=sod, **kw2)
        var[nm] = strip(rr); print(summ(nm, rr), flush=True)
    t2, f2, seg2, cam2, _ = load(TARGET, col='FLUX2', flag=None)
    for nm, o in (('RAW_o1', 1), ('RAW_o2', 2)):
        rr, _ = L.run_survey(t2, f2, seg2, shp, EXPT, o, nf, ni, s_over_d=sod)
        var[nm] = strip(rr); print(summ(nm, rr), flush=True)
    for c in (227, 228):
        m = cam == c
        rr, _ = L.run_survey(t[m], f[m], seg[m], shp, EXPT, 1, nf, ni, s_over_d=sod)
        var[f'TAM_cam{c}'] = strip(rr); print(summ(f'TAM_cam{c}', rr), flush=True)
    for nm, (lo, hi) in (('2006', (2453800, 2454100)), ('2007', (2454100, 2454450)), ('2008', (2454450, 2454700))):
        m = (t > lo) & (t < hi)
        if m.sum() < 1000:
            print('season', nm, 'only', m.sum(), 'points; skipped'); continue
        rr, _ = L.run_survey(t[m], f[m], seg[m], shp, EXPT, 1, nf, ni, s_over_d=sod)
        var['TAM_season' + nm] = strip(rr); print(summ('TAM_season' + nm, rr), flush=True)
    # airmass/hour-angle proxy: drop the first and last 0.04 d of each night (twilight/high airmass)
    keep = np.ones(len(t), bool)
    for k in np.unique(seg):
        m = seg == k
        keep[m] = (t[m] > t[m].min() + 0.04) & (t[m] < t[m].max() - 0.04)
    rr, _ = L.run_survey(t[keep], f[keep], seg[keep], shp, EXPT, 1, nf, ni, s_over_d=sod)
    var['TAM_trim_night_ends'] = strip(rr); print(summ('TAM_trim_night_ends', rr), flush=True)
    # nominal (task) ephemeris instead of the refined one
    rr, _ = L.run_survey(t, f, seg, shp_nom, EXPT, 1, nf, ni, s_over_d=sod)
    var['TAM_nominal_ephemeris'] = strip(rr); print(summ('TAM_nominal_ephemeris', rr), flush=True)
    out['variants'] = var
    ctl = {}
    for nm in UNBLENDED + BLENDED:
        tc, fc, sgc, _, _ = load(nm)
        rr, _ = L.run_survey(tc, fc, sgc, shp, EXPT, 1, nf, 60, s_over_d=sod)
        ctl[nm] = strip(rr); print(summ(nm, rr), flush=True)
    out['controls'] = ctl
    json.dump(out, open(L.V3 + '/results/wasp3.json', 'w'), indent=1)
    print('saved')
