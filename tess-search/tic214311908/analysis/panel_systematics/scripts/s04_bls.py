"""Blind box search of the raw-pixel aperture light curves (s03 'blind'), own code.
Bins to 30 min, per-sector robust normalisation, astropy BoxLeastSquares over P = 0.8-20 d, durations 1.0-3.0 h.
Reports the top peaks (power = delta-chi2 'likelihood'), the rank of 5.6915 d, and the peak SDE.
Output: out/s04_bls_{ap}.json
"""
import sys, numpy as np
sys.path.insert(0, 'tic214311908/analysis/panel_systematics/scripts')
from ps_lib import ROOT, jdump, P as P0
from astropy.timeseries import BoxLeastSquares

ap = sys.argv[1] if len(sys.argv) > 1 else 'box'
kn = sys.argv[2] if len(sys.argv) > 2 else '0.75'
d = np.load(f'{ROOT}/out/s03_lc_blind_k{kn}.npz')
t, f, s = d['t'], d[ap], d['sec']
ok = np.isfinite(f)
t, f, s = t[ok], f[ok], s[ok]
# 30-min bins per sector, errors from the per-sector robust scatter of binned points
TB, FB, EB = [], [], []
for sec in np.unique(s):
    m = s == sec
    tb = np.floor(t[m] / (0.5 / 24))
    u, inv = np.unique(tb, return_inverse=True)
    n = np.bincount(inv)
    fb = np.bincount(inv, f[m]) / n
    tt = np.bincount(inv, t[m]) / n
    keep = n >= max(1, 0.5 * np.median(n))
    fb, tt = fb[keep], tt[keep]
    r = 1.4826 * np.median(np.abs(fb - np.median(fb)))
    clip = np.abs(fb - np.median(fb)) < 6 * r
    TB.append(tt[clip]); FB.append(fb[clip]); EB.append(np.full(clip.sum(), r))
t, f, e = map(np.concatenate, (TB, FB, EB))
bls = BoxLeastSquares(t, f / 1e3 + 1, dy=e / 1e3)
periods = np.exp(np.arange(np.log(0.8), np.log(20.0), 0.06 / 2500 / 3))  # df ~ q-limited over a 2550-d baseline
durs = np.array([1.0, 1.4, 1.8, 2.4, 3.0]) / 24
print('n periods', len(periods), 'npts', len(t), flush=True)
r = bls.power(periods, durs, objective='likelihood')
pw = np.array(r.power)
# positive (dimming) only: depth sign
dep = np.array(r.depth)
pw = np.where(dep > 0, pw, 0)
o = np.argsort(pw)[::-1]
top = []
for i in o:
    if all(abs(np.log(periods[i] / p['P'])) > 0.002 for p in top):
        top.append(dict(P=float(periods[i]), power=float(pw[i]), depth_ppt=float(dep[i] * 1e3), dur_h=float(r.duration[i] * 24),
                        t0=float(r.transit_time[i])))
    if len(top) >= 15:
        break
med = np.median(pw); sd = 1.4826 * np.median(np.abs(pw - med))
sde = (top[0]['power'] - np.mean(pw)) / np.std(pw)
near = np.abs(periods / P0 - 1) < 2e-4
out = dict(ap=ap, knots=kn, npts=len(t), nper=len(periods), top=top, sde_top=float(sde),
           power_at_P=float(pw[near].max()), rank_of_P=int(np.sum(pw > pw[near].max())),
           second_ratio=float(top[1]['power'] / top[0]['power']))
for p in top[:10]:
    print('%.5f %.1f depth %.1f ppt dur %.1f h' % (p['P'], p['power'], p['depth_ppt'], p['dur_h']))
print('SDE', round(sde, 1), 'power at P', round(out['power_at_P'], 1), 'rank', out['rank_of_P'])
np.save(f'{ROOT}/out/s04_pgram_{ap}_k{kn}.npy', np.c_[periods[::20], pw.reshape(-1)[::20]].astype(np.float32))
jdump(out, f'{ROOT}/out/s04_bls_{ap}_k{kn}.json')
