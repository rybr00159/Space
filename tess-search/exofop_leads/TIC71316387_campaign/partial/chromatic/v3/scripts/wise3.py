"""NEOWISE-R (2014-2024) + AllWISE multi-epoch (2010) single-exposure W1/W2 profile-fit photometry.

Baseline: constant per visit (a visit = ~1 d of exposures every ~95 min).  Fake-phase null (5000 offsets,
real eclipses modelled), end-to-end injection (re-cleaned), event jackknife, quality-cut variants.
"""
import sys, json
import numpy as np
import pandas as pd
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/v3/scripts')
import lib3 as L

T = json.load(open(L.V3 + '/results/tess3.json'))
E3 = json.load(open(L.V3 + '/results/ephem3.json'))
shp = L.survey_shapes(T['shapes'], E3)
sod = T['s_over_d']
EX = 7.7 / 86400
DATA = L.CH + '/data/wise/'
out = {}


def load(cut='strict', allwise=True):
    n = pd.read_csv(DATA + 'neowise.csv')
    n['src'] = 'NEOWISE'
    if allwise:
        a = pd.read_csv(DATA + 'allwise_mep.csv').rename(columns={'w1mpro_ep': 'w1mpro', 'w1sigmpro_ep': 'w1sigmpro', 'w2mpro_ep': 'w2mpro',
                                                                 'w2sigmpro_ep': 'w2sigmpro'})
        a['qual_frame'] = 10 * (a.qi_fact.values > 0); a['src'] = 'AllWISE'
        n = pd.concat([n, a], ignore_index=True)
    m = np.isfinite(n.w1mpro.values) & (n.nb.values == 1) & (n.cc_flags.astype(str).isin(['0', '0000']))
    if cut == 'strict':
        m &= (n.qual_frame.values > 0) & (n.saa_sep.values > 5) & (n.moon_masked.astype(str).isin(['0', '00', '0000']))
    elif cut == 'loose':
        m &= (n.qual_frame.values > 0)
    n = n[m]
    t = L.jd_utc_to_bjd_tdb(n.mjd.values + 2400000.5)
    return n, t


def run(n, t, band, nfake=5000, seed=1, order=0):
    fl = 10 ** (-0.4 * (n[band + 'mpro'].values - 10.3))
    ok = np.isfinite(fl)
    t, fl = t[ok], fl[ok]
    seg = L.segments(t, 5.0)
    med = np.array([np.median(fl[seg == j]) for j in range(seg.max() + 1)])
    f = fl / med[seg]
    r, cl = L.run_survey(t, f, seg, shp, EX, order, nfake=nfake, ninj=200, inj_depths=(0.006, 0.015), s_over_d=sod, minpts=5, floor=3e-3, seed=seed)
    return r, cl, (t, f, seg)


for cut in ('strict', 'loose', 'none'):
    for aw in (True, False):
        n, t = load(cut, aw)
        for band in ('w1', 'w2'):
            r, cl, raw = run(n, t, band, nfake=3000 if (cut, aw) != ('strict', True) else 5000)
            key = f'{band.upper()}_{cut}_{"allwise" if aw else "neowiseonly"}'
            z = (r['deep'] - r['deep_null_mean']) / r['deep_err']
            print('%-26s n=%3d nvis=%2d in-d %2d (%d ev) in-s %2d | deep %6.0f +- %5.0f (z=%.2f, p>=%.4f, null kurt %.2f) | shallow %6.0f +- %5.0f' % (
                key, r['n'], r['nseg'], r['n_in_deep'], r['ev_deep'], r['n_in_shallow'], r['deep'] * 1e6, r['deep_err'] * 1e6, z, r['deep_p_ge'],
                r['deep_null_kurt'], r['shallow'] * 1e6, r['shallow_err'] * 1e6), flush=True)
            if (cut, aw) == ('strict', True):
                out[key + '_null_deep'] = r['null_samples_deep']
                # leave-one-event-out jackknife on the deep amplitude
                tc, fc, wc, sgc = cl['t'], cl['f'], cl['w'], cl['seg']
                s = L.tmpl(tc, shp, 'd', 0, EX)
                cyc = np.round((tc - shp['d']['t0']) / shp['d']['per'] - shp['d']['phc']).astype(int)
                evs = np.unique(cyc[s > 0.2])
                jk = []
                for c in evs:
                    m = ~((cyc == c) & (s > 0.01))
                    pr = L.Projector(tc[m], wc[m], sgc[m], 0)
                    a, e = L.measure_at(tc[m], fc[m], wc[m], pr, shp, EX, 0.0)
                    jk.append(a[0])
                jk = np.array(jk)
                r['jackknife_deep'] = jk.tolist()
                r['jackknife_sd'] = float(np.sqrt((len(jk) - 1) / len(jk) * np.sum((jk - jk.mean()) ** 2)))
                print('     leave-one-event-out deep amplitudes: %s ; jackknife sd %.0f ppm' % (np.round(jk * 1e6).astype(int).tolist(), r['jackknife_sd'] * 1e6))
                # end-to-end injection (re-cleaned, paired)
                rng = np.random.default_rng(9)
                t0_, f0_, seg0_ = raw
                for D in (0.006, 0.015):
                    rec = L.e2e_injection(t0_, f0_, seg0_, shp, EX, 0, L.fake_offsets(120, rng), D, sod, minpts=5, floor=3e-3)
                    ok = np.isfinite(rec[:, 2])
                    r[f'e2e_{D}'] = dict(rec_mean=float(np.nanmean(rec[:, 0])), rec_sd=float(np.nanstd(rec[:, 0])),
                                         gain_deep=float(np.nanmean(rec[:, 2])), gain_deep_sem=float(np.nanstd(rec[:, 2]) / np.sqrt(ok.sum())),
                                         gain_shallow=float(np.nanmean(rec[:, 3])))
                    print('     e2e injection %.3f: recovered %.0f +- %.0f ppm, gain deep %.3f +- %.3f shallow %.3f' % (
                        D, r[f'e2e_{D}']['rec_mean'] * 1e6, r[f'e2e_{D}']['rec_sd'] * 1e6, r[f'e2e_{D}']['gain_deep'], r[f'e2e_{D}']['gain_deep_sem'],
                        r[f'e2e_{D}']['gain_shallow']), flush=True)
                out[key] = {k: v for k, v in r.items() if not k.startswith('null_samples')}
            else:
                out[key] = {k: v for k, v in r.items() if not k.startswith('null_samples')}
json.dump(out, open(L.V3 + '/results/wise3.json', 'w'), indent=1)
print('saved')
