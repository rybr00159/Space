"""v3 TESS reference: eclipse template shape, per-sector depths (fake-phase nulls), per-sector
eclipse times and an independent check of the ephemeris precision.

Data
  SPOC/TESS-SPOC PDCSAP (HLSP) S1 (1800 s), S28 (600 s), S68 (200 s): QUALITY==0; crowding-corrected
     by SPOC (CROWDSAP 0.970-0.979) -> depth relative to the target + anything unresolved in the TIC.
  TESScut aperture photometry from the lcmodel/v5 sub-campaign S95, S105, S106 (200 s): f_raw with
     the PRF-model contamination fraction removed, times corrected by its timecorr arrays.
Segmentation: quarter-orbit chunks with boundaries at phases (dphi + 0.125 + k/4), i.e. each real or
fake eclipse sits at the centre of its own chunk; linear baseline per chunk, fitted simultaneously.
"""
import sys, json, glob
import numpy as np
from astropy.io import fits
from scipy.optimize import minimize
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/v3/scripts')
import lib3 as L

SPOC = L.BASE + '/exofop/agents/ctoi_vet/recent/data/mastDownload/HLSP/'
V5 = L.V3 + '/data/'   # frozen copies of lcmodel/v5/data/tc_s*.npz (times already corrected there)


def load_all():
    data = {}
    for sec, ex in ((1, 1800.), (28, 600.), (68, 200.)):
        fn = glob.glob(SPOC + f'*0071316387-s{sec:04d}*/*lc.fits')[0]
        with fits.open(fn) as h:
            assert int(h[0].header['TICID']) == 71316387
            d = h[1].data
            t = d['TIME'] + 2457000.0
            f = np.array(d['PDCSAP_FLUX'], float)
            ok = np.isfinite(t) & np.isfinite(f) & (d['QUALITY'] == 0)
            data[sec] = dict(t=t[ok], f=f[ok] / np.nanmedian(f[ok]), ex=ex / 86400, src='SPOC PDCSAP',
                             crowdsap=float(h[1].header['CROWDSAP']))
    for sec in (95, 105, 106):
        z = np.load(V5 + f'tc_s{sec}.npz')
        t = z['t'] + 2457000.0          # a02_phot.py already applied its timecorr
        f = np.array(z['f_raw'], float); c = float(z['contam'])
        ok = np.isfinite(t) & np.isfinite(f)
        f = f[ok] / np.median(f[ok]); f = (f - c) / (1 - c)
        data[sec] = dict(t=t[ok], f=f, ex=float(z['texp']), src='TESScut aperture (lcmodel v5), PRF contamination removed', contam=c)
    return data


def chunks(t, dphi=0.0):
    """quarter-orbit chunk id with boundaries at phases dphi+0.125+k/4"""
    x = (t - L.T0) / L.P - dphi - 0.125
    return np.floor(x * 4).astype(int)


def shp_from(par):
    return {'d': dict(T=par[0], p=par[1], q=par[2], phc=par[3]), 's': dict(T=par[4], p=par[5], q=par[6], phc=0.5 + par[7])}


def fit_sector(d, shp, dphi=0.0, w=None, extra_real=False):
    t, f, ex = d['t'], d['f'], d['ex']
    seg = chunks(t, dphi)
    cols = [L.tmpl(t, shp, 'd', dphi, ex), L.tmpl(t, shp, 's', dphi, ex)]
    if extra_real:
        cols += [L.tmpl(t, shp, 'd', 0, ex), L.tmpl(t, shp, 's', 0, ex)]
    if w is None:
        w = np.ones(len(t))
    a, e, r, pr = L.fit_amps(t, f, w, seg, cols, 1)
    return a, e, r


if __name__ == '__main__':
    data = load_all()
    # --- clipping and white-noise weights per sector (iterate with a provisional shape)
    shp0 = shp_from([0.068, 2.0, 1.0, 0.0, 0.068, 2.0, 1.0, 0.0])
    for sec, d in data.items():
        keep = np.ones(len(d['t']), bool)
        for _ in range(5):
            dd = dict(t=d['t'][keep], f=d['f'][keep], ex=d['ex'])
            a, e, r = fit_sector(dd, shp0)
            s = L.mad_std(r)
            bad = np.abs(r) > 5 * s
            if not bad.any():
                break
            idx = np.where(keep)[0]; keep[idx[bad]] = False
        d['t'], d['f'] = d['t'][keep], d['f'][keep]
        d['sig'] = s
        d['w'] = np.ones(len(d['t'])) / s ** 2
        print('S%d n=%d point rms %.0f ppm' % (sec, len(d['t']), s * 1e6), flush=True)

    # --- template shape: joint fit to the 200 s + 600 s sectors
    SHAPE_SECS = (28, 68, 95, 105, 106)

    def chi2(par):
        T1, p1, q1, d1, T2, p2, q2, d2 = par
        if not (0.04 < T1 < 0.12 and 0.04 < T2 < 0.12 and 0.5 < p1 < 10 and 0.5 < p2 < 10 and 0.3 < q1 < 10
                and 0.3 < q2 < 10 and abs(d1) < 0.01 and abs(d2) < 0.01):
            return 1e30
        shp = shp_from(par)
        return sum(np.sum(data[s]['w'] * fit_sector(data[s], shp, w=data[s]['w'])[2] ** 2) for s in SHAPE_SECS)

    best = None
    for st in ([0.068, 2.0, 1.0, 0, 0.068, 2.0, 1.0, 0], [0.064, 1.3, 1.6, 0, 0.066, 1.3, 1.5, 0], [0.075, 3.5, 3.0, 0, 0.075, 3.5, 3.0, 0]):
        r = minimize(chi2, st, method='Nelder-Mead', options=dict(xatol=1e-6, fatol=1e-4, maxiter=6000, maxfev=6000))
        r = minimize(chi2, r.x, method='Nelder-Mead', options=dict(xatol=1e-6, fatol=1e-4, maxiter=6000, maxfev=6000))
        print('shape start', st[:3], '-> chi2 %.2f' % r.fun, np.round(r.x, 5), flush=True)
        if best is None or r.fun < best.fun:
            best = r
    shp = shp_from(best.x)
    ntot = sum(len(data[s]['t']) for s in SHAPE_SECS)
    print('shape chi2 %.1f for %d points' % (best.fun, ntot))
    for k in 'ds':
        s = shp[k]
        # FWHM-like and "effective" duration where template > 0.01
        x = np.linspace(-1, 1, 20001); y = L.vshape(x, s['p'], s['q'])
        fw = (x[y >= 0.5].max() - x[y >= 0.5].min()) * s['T'] / 2 * 24
        print(k, 'T=%.4f d (%.3f h) p=%.3f q=%.3f FWHM %.3f h; centre offset %.2f min' % (
            s['T'], s['T'] * 24, s['p'], s['q'], fw, (s['phc'] - (0 if k == 'd' else 0.5)) * L.P * 1440))
    out = dict(shapes=shp, shape_chi2=float(best.fun), shape_npts=int(ntot), sectors={})

    # --- per-sector depths with fake-phase null (chunks re-centred on every fake pair; real eclipses modelled)
    rng = np.random.default_rng(5)
    for sec in sorted(data):
        d = data[sec]
        a, e, r = fit_sector(d, shp, w=d['w'])
        beta_r = L.mad_std(r * np.sqrt(d['w']))
        fk, fe = [], []
        for dphi in L.fake_offsets(300, rng):
            aa, ee, rr = fit_sector(d, shp, dphi, w=d['w'], extra_real=True)
            ns = [np.sum(L.tmpl(d['t'], shp, k, dphi, d['ex']) > 0.5) for k in 'ds']
            if min(ns) < 3:
                continue
            fk.append(aa[:2]); fe.append(ee[:2])
        fk, fe = np.array(fk), np.array(fe)
        beta = np.std(fk / fe, axis=0, ddof=1)
        row = dict(src=d['src'], exptime_s=d['ex'] * 86400, n=int(len(d['t'])), rms=float(d['sig']),
                   deep=float(a[0]), shallow=float(a[1]), deep_formal=float(e[0]), shallow_formal=float(e[1]),
                   null_mean=fk.mean(0).tolist(), null_sd=fk.std(0, ddof=1).tolist(), beta=beta.tolist(), n_fake=len(fk),
                   crowdsap=d.get('crowdsap'), contam=d.get('contam'))
        row['deep_err'] = float(max(row['null_sd'][0], beta[0] * e[0]))
        row['shallow_err'] = float(max(row['null_sd'][1], beta[1] * e[1]))
        out['sectors'][str(sec)] = row
        print('S%-3d deep %5.0f +- %3.0f (formal %3.0f null sd %3.0f mean %4.0f beta %.2f) shallow %5.0f +- %3.0f (null mean %4.0f) s/d %.3f' % (
            sec, a[0] * 1e6, row['deep_err'] * 1e6, e[0] * 1e6, row['null_sd'][0] * 1e6, row['null_mean'][0] * 1e6, beta[0],
            a[1] * 1e6, row['shallow_err'] * 1e6, row['null_mean'][1] * 1e6, a[1] / a[0]), flush=True)

    def combine(secs, key):
        v = np.array([out['sectors'][str(s)][key] for s in secs])
        e = np.array([out['sectors'][str(s)][key + '_err'] for s in secs])
        w = 1 / e ** 2
        m = np.sum(w * v) / np.sum(w)
        chi = float(np.sum(w * (v - m) ** 2))
        es = np.sqrt(1 / np.sum(w))
        return dict(mean=float(m), err=float(es * np.sqrt(max(1, chi / (len(v) - 1)))), err_stat=float(es), chi2=chi,
                    dof=len(v) - 1, sd=float(np.std(v, ddof=1)), values=v.tolist(), errs=e.tolist())

    for nm, secs in (('spoc', (1, 28, 68)), ('tesscut', (95, 105, 106)), ('all', (1, 28, 68, 95, 105, 106))):
        for key in ('deep', 'shallow'):
            c = combine(secs, key); out[nm + '_' + key] = c
            print('%-8s %-8s %5.0f +- %3.0f ppm (chi2 %.1f/%d, sd %.0f)' % (nm, key, c['mean'] * 1e6, c['err'] * 1e6, c['chi2'], c['dof'], c['sd'] * 1e6))
    # adopted TESS reference: all six sectors; add half the SPOC-vs-TESScut difference as a systematic
    for key in ('deep', 'shallow'):
        al, sp, tc = out['all_' + key], out['spoc_' + key], out['tesscut_' + key]
        sysf = 0.5 * abs(sp['mean'] - tc['mean'])
        out['adopted_' + key] = dict(mean=al['mean'], err=float(np.hypot(al['err'], sysf)), err_stat=al['err'], sys=float(sysf))
        print('ADOPTED %s %.0f +- %.0f ppm (stat %.0f, sys %.0f)' % (key, al['mean'] * 1e6, np.hypot(al['err'], sysf) * 1e6, al['err'] * 1e6, sysf * 1e6))
    sd = out['adopted_shallow']['mean'] / out['adopted_deep']['mean']
    out['s_over_d'] = float(sd)
    print('s/d = %.3f' % sd)

    # --- per-sector eclipse times (shape fixed, time offset scanned), errors by event bootstrap
    def sector_time(d, which, grid):
        chi = []
        for g in grid:
            sh = {k: dict(v) for k, v in shp.items()}
            sh['d']['phc'] += g / L.P; sh['s']['phc'] += g / L.P
            # only the eclipse under study is shifted; the other stays at the fitted value
            if which == 'd':
                sh['s']['phc'] = shp['s']['phc']
            else:
                sh['d']['phc'] = shp['d']['phc']
            a, e, r = fit_sector(d, sh, w=d['w'])
            chi.append(np.sum(d['w'] * r ** 2))
        chi = np.array(chi)
        i = np.argmin(chi)
        j = slice(max(i - 6, 0), min(i + 7, len(grid)))
        c = np.polyfit(grid[j], chi[j], 2)
        tbest = -c[1] / (2 * c[0])
        return tbest, 1 / np.sqrt(c[0]), chi
    grid = np.linspace(-6, 6, 121) / 1440.0
    times = {}
    for sec in sorted(data):
        d = data[sec]
        E = np.round((np.median(d['t']) - L.T0) / L.P)
        rowt = {}
        for which in 'ds':
            tb, te, chi = sector_time(d, which, grid)
            # bootstrap over events (resample eclipse cycles with replacement) for an empirical error
            cyc = np.round((d['t'] - L.T0) / L.P - (0 if which == 'd' else 0.5)).astype(int)
            uc = np.unique(cyc)
            bs = []
            for b in range(60):
                pick = rng.choice(uc, len(uc), replace=True)
                idx = np.concatenate([np.where(cyc == c)[0] for c in pick])
                db = dict(t=d['t'][idx], f=d['f'][idx], ex=d['ex'], w=d['w'][idx])
                # duplicated points -> sort to keep chunks coherent
                o = np.argsort(db['t'], kind='stable'); db = {k: (v[o] if hasattr(v, '__len__') else v) for k, v in db.items()}
                tb2, _, _ = sector_time(db, which, grid[::2])
                bs.append(tb2)
            rowt[which] = dict(dt_min=float(tb * 1440), err_formal_min=float(te * 1440), err_boot_min=float(np.std(bs, ddof=1) * 1440), E=float(E))
        times[str(sec)] = rowt
        print('S%-3d E=%6.0f  deep dt %+6.2f +- %.2f (boot %.2f) min   shallow dt %+6.2f +- %.2f (boot %.2f) min' % (
            sec, E, rowt['d']['dt_min'], rowt['d']['err_formal_min'], rowt['d']['err_boot_min'],
            rowt['s']['dt_min'], rowt['s']['err_formal_min'], rowt['s']['err_boot_min']), flush=True)
    out['times'] = times
    json.dump(out, open(L.V3 + '/results/tess3.json', 'w'), indent=1)
    np.savez(L.V3 + '/data/tess3_clean.npz', **{f's{s}_{k}': data[s][k] for s in data for k in ('t', 'f', 'w')},
             **{f's{s}_ex': data[s]['ex'] for s in data})
