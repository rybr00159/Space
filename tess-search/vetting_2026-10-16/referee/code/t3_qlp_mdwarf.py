"""TIC 214311906 vs TIC 214311908 (M dwarf 18.9" W): QLP fractional dip depth
at the 5.69149-d ephemeris in each sector, plus a global shape test (flat vs V)
on my own TESScut photometry centred on the M dwarf."""
import sys, json, glob
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/five_referee/code')
import numpy as np
from astropy.io import fits
from scipy.optimize import least_squares
from reflib import *

P, T0, T14, T12 = 5.69149, 3840.913, 1.55 / 24, 0.25 / 24
OUT = {}


def readqlp(path):
    d = fits.open(path)[1].data
    q = d['QUALITY'] == 0
    t = np.array(d['TIME'][q], float); f = np.array(d['SAP_FLUX'][q], float)
    m = np.isfinite(t) & np.isfinite(f)
    return t[m], f[m] / np.median(f[m])


def joint(t, f, T14=T14, T12=T12):
    texp = np.median(np.diff(t))
    n0, n1 = int(np.floor((t.min() - T0) / P)), int(np.ceil((t.max() - T0) / P))
    blocks = []
    for n in range(n0, n1 + 1):
        c = T0 + n * P
        s = np.abs(t - c) < 0.3
        if s.sum() > 8 and np.sum(np.abs(t[s] - c) < T14 / 2) >= 2 and np.sum(t[s] < c - T14 / 2) >= 3 and np.sum(t[s] > c + T14 / 2) >= 3:
            blocks.append((c, s))
    nW = len(blocks)
    X = np.zeros((0, 1 + 3 * nW)); Y = np.zeros(0)
    for k, (c, s) in enumerate(blocks):
        dt = t[s] - c
        sh = 1 - supersample(trapezoid, t[s], texp, t0=c, depth=1.0, T14=T14, T12=T12)
        Xk = np.zeros((s.sum(), 1 + 3 * nW)); Xk[:, 0] = -sh
        Xk[:, 1 + 3 * k] = 1; Xk[:, 2 + 3 * k] = dt; Xk[:, 3 + 3 * k] = dt ** 2
        X = np.vstack([X, Xk]); Y = np.concatenate([Y, f[s]])
    b, *_ = np.linalg.lstsq(X, Y, rcond=None)
    r = Y - X @ b
    cov = np.linalg.pinv(X.T @ X) * np.var(r, ddof=X.shape[1])
    return b[0], np.sqrt(cov[0, 0]), nW, np.sum(r ** 2)


for sec in (12, 39, 66, 93, 103, 104):
    row = {}
    for tic in (214311906, 214311908, 214311910):
        if tic == 214311906:
            p = lc_path(tic, sec, 'QLP')
        else:
            gg = glob.glob(f'{BASE}/data/dl_qlp_{tic}/*s{sec:04d}-{tic:016d}*llc.fits')
            if not gg:
                continue
            p = gg[0]
        t, f = readqlp(p)
        if sec in (103, 104):
            t = t + np.nanmedian(fits.open(cut_path(214311906, sec))[1].data['TIMECORR'])
        d, e, n, _ = joint(t, f)
        row[tic] = (d * 1e6, e * 1e6, n)
    OUT[f'S{sec}'] = row
    print(f"S{sec:3d}: QLP SAP depth (ppm) " + '   '.join(f'TIC{k}: {v[0]:6.0f}+-{v[1]:4.0f}' for k, v in row.items()))

# shape test on TESScut photometry centred on the M dwarf (r=1 px), all sectors combined, relative chi2 for V vs flat
ts, fs = [], []
for sec in (12, 39, 66, 93, 103, 104):
    cut = load_cut(214311906, sec)
    g = gaia_table(214311906, btjd_to_jyear(np.median(cut['time'])))
    a = g[g.source_id == 5925373272850784768].iloc[0]
    t, f, e, ex = aperture_lc(214311906, sec, r=1.0, center=(a['ra_ep'], a['dec_ep']))
    ts.append(t); fs.append(f)
sig = []
for t, f in zip(ts, fs):
    d, e, n, ss = joint(t, f)
    sig.append(np.sqrt(ss / max(1, len(t))))
res = {}
for name, (t14, t12) in {'flat_T12=0.15h': (1.55 / 24, 0.15 / 24), 'flat_T12=0.3h': (1.55 / 24, 0.3 / 24), 'V_1.55h': (1.55 / 24, 0.775 / 24),
                         'V_2.2h': (2.2 / 24, 1.1 / 24), 'V_1.8h': (1.8 / 24, 0.9 / 24), 'box_1.2h': (1.2 / 24, 0.05 / 24), 'flat_1.8h_T12=0.3': (1.8 / 24, 0.3 / 24), 'flat_1.3h_T12=0.2': (1.3 / 24, 0.2 / 24)}.items():
    chi = 0; dd = []
    for (t, f), sg in zip(zip(ts, fs), sig):
        d, e, n, c2 = joint(t, f, T14=t14, T12=t12)
        chi += c2 / sg ** 2; dd.append((d, e))
    dd = np.array(dd)
    w = 1 / dd[:, 1] ** 2
    res[name] = dict(chi2=chi, depth_ppm=np.sum(dd[:, 0] * w) / np.sum(w) * 1e6, err_ppm=1 / np.sqrt(np.sum(w)) * 1e6)
best = min(v['chi2'] for v in res.values())
for k, v in res.items():
    print(f'   shape {k:16s}: delta-chi2 = {v["chi2"]-best:6.1f}   depth {v["depth_ppm"]:.0f}+-{v["err_ppm"]:.0f} ppm (r=1 px on M dwarf)')
OUT['shape'] = res
save_json(OUT, f'{BASE}/results/t3_qlp_mdwarf.json')
