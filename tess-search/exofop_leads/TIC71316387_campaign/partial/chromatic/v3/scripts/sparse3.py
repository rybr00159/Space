"""Sparse archival surveys: ASAS-SN V (photometry database, 2014-2018), ASAS-3 V (2000-2009) and
NEOWISE-R W1/W2 single exposures (2014-2024).

These have ~1 point per night (ASAS-SN, ASAS-3) or ~12-15 points per ~1-d visit (NEOWISE), so the
baseline is per camera x season (ASAS-SN, ASAS-3; linear in time) or per visit (NEOWISE; constant).
Depth estimate, fake-phase null (real eclipses modelled), paired injection-recovery as in lib3.run_survey.
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
DATA = L.CH + '/data/'
out = {}


def season_seg(t, cam=None, gap=60.0):
    return L.segments(t, gap, cam)


def report(nm, r):
    print('%-18s n=%4d nseg=%3d in-d %3d (%2d ev) in-s %3d | deep %6.0f +- %5.0f (null mean %5.0f, beta %.2f) p>=%.3f | shallow %6.0f +- %5.0f | rms %.4f' % (
        nm, r['n'], r['nseg'], r['n_in_deep'], r['ev_deep'], r['n_in_shallow'], r['deep'] * 1e6, r['deep_err'] * 1e6,
        r['deep_null_mean'] * 1e6, r['deep_beta'], r['deep_p_ge'], r['shallow'] * 1e6, r['shallow_err'] * 1e6, r['rms_med']), flush=True)


def strip(r):
    r = dict(r); r.pop('null_samples_deep'); r.pop('null_samples_shallow'); return r


# ------------------------------------------------------------------ ASAS-SN V (photometry database)
d = np.loadtxt(DATA + 'asassn_web/asassn_photdb_V.txt')
cam = np.loadtxt(DATA + 'asassn_web/asassn_photdb_V_cam.txt', dtype=str)
hjd, mag, merr, flux, ferr = d.T
ok = np.isfinite(flux) & (flux > 0) & (merr < 0.05)
t = L.hjd_utc_to_bjd_tdb(hjd[ok]); f = flux[ok]; c = cam[ok]
seg = season_seg(t, c)
med = np.array([np.median(f[seg == k]) for k in range(seg.max() + 1)])
f = f / med[seg]
EX = 3 * 90.0 / 86400 + 120 / 86400      # 3 dithered 90-s exposures spread over ~6.5 min
for nm, cams in (('ASASSN_V_all', ('bf', 'bh')), ('ASASSN_V_bf', ('bf',)), ('ASASSN_V_bh', ('bh',))):
    m = np.isin(c, cams)
    r, cl = L.run_survey(t[m], f[m], seg[m], shp, EX, 1, nfake=2000, ninj=200, inj_depths=(0.003, 0.006), s_over_d=sod, minpts=5, floor=1e-3)
    if r is None:
        print('%-18s fewer than 3 points inside either eclipse: no measurement' % nm); out[nm] = None; continue
    out[nm] = strip(r); report(nm, r)
r, cl = L.run_survey(t, f, seg, shp, EX, 0, nfake=2000, ninj=100, s_over_d=sod, minpts=5, floor=1e-3)
out['ASASSN_V_all_o0'] = strip(r); report('ASASSN_V_all_o0', r)

# ------------------------------------------------------------------ ASAS-3 V
# parse the file with all datasets (asas3_220333-4250.0.txt); select datasets of the target designation
rows = []
cur = None
for line in open(DATA + 'asas3/asas3_220333-4250.0.txt'):
    if line.startswith('#dataset'):
        cur = line.split(';')[1].strip()
        continue
    if line.startswith('#desig'):
        desig = line.split('=')[1].strip()
        continue
    if line.startswith('#') or not line.strip():
        continue
    p = line.split()
    rows.append([float(p[0])] + [float(x) for x in p[1:11]] + [p[11], cur, desig])
A = pd.DataFrame(rows, columns=['hjd'] + [f'm{k}' for k in range(5)] + [f'e{k}' for k in range(5)] + ['grade', 'field', 'desig'])
A = A[(A.desig == '220333-4250.1') & (A.grade.isin(['A', 'B']))]
print('ASAS-3: %d A/B points of 220333-4250.1 in fields %s' % (len(A), A.field.value_counts().to_dict()))
best = None
for k in range(5):
    mm = A[f'm{k}'].values
    g = mm < 29
    s = L.mad_std(mm[g])
    print('  aperture MAG_%d: robust scatter %.4f mag (n=%d)' % (k, s, g.sum()))
    if best is None or s < best[1]:
        best = (k, s)
k = best[0]
mm = A[f'm{k}'].values; g = mm < 29
t = L.hjd_utc_to_bjd_tdb(A.hjd.values[g] + 2450000.0)
f = 10 ** (-0.4 * (mm[g] - np.median(mm[g])))
seg = season_seg(t, A.field.values[g])
med = np.array([np.median(f[seg == j]) for j in range(seg.max() + 1)])
f = f / med[seg]
r, cl = L.run_survey(t, f, seg, shp, 180 / 86400, 1, nfake=2000, ninj=100, inj_depths=(0.006, 0.02), s_over_d=sod, minpts=5, floor=3e-3)
r['aperture'] = int(k)
out['ASAS3_V'] = strip(r); report('ASAS3_V (MAG_%d)' % k, r)

# ------------------------------------------------------------------ NEOWISE W1 / W2
n = pd.read_csv(DATA + 'wise/neowise.csv')
n = n[(n.qual_frame > 0) & (n.saa_sep > 5) & (n.moon_masked == 0) & (n.cc_flags.astype(str) == '0') & (n.nb == 1)]
t = L.jd_utc_to_bjd_tdb(n.mjd.values + 2400000.5)
seg = L.segments(t, 5.0)
print('NEOWISE: %d points in %d visits after quality cuts' % (len(t), seg.max() + 1))
for band in ('w1', 'w2'):
    fl = 10 ** (-0.4 * (n[band + 'mpro'].values - 10.3))
    ok = np.isfinite(fl)
    med = np.array([np.median(fl[ok][seg[ok] == j]) for j in range(seg.max() + 1)])
    ff = fl[ok] / med[seg[ok]]
    r, cl = L.run_survey(t[ok], ff, seg[ok], shp, 7.7 / 86400, 0, nfake=2000, ninj=100, inj_depths=(0.006, 0.012), s_over_d=sod, minpts=5, floor=3e-3)
    out['NEOWISE_' + band.upper()] = strip(r); report('NEOWISE_' + band.upper(), r)
json.dump(out, open(L.V3 + '/results/sparse3.json', 'w'), indent=1)
print('saved')
