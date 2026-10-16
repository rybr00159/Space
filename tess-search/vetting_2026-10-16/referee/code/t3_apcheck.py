"""TIC 214311906: model-light check of the localisation - absolute dip depth
(e-/s) in small apertures centred on the Mira, on the two M dwarfs, on the
faint stars ~17" W, and on mirror positions. Joint fit over all events of a sector."""
import sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/five_referee/code')
import numpy as np
from reflib import *

TIC = 214311906
P, T0, T14, T12 = 5.69149, 3840.913, 1.55 / 24, 0.25 / 24
OUT = {}
pos = {'Mira': None, 'MdwarfA_18.9W': 5925373272850784768, 'MdwarfB_7W': 5925373272850782208, 'G17.5_16.7W': 5925373272827260544}


def abs_depth(t, f):
    texp = np.median(np.diff(t))
    n0, n1 = int(np.floor((t.min() - T0) / P)), int(np.ceil((t.max() - T0) / P))
    blocks = []
    for n in range(n0, n1 + 1):
        c = T0 + n * P
        s = np.abs(t - c) < 0.25
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
    return b[0], np.sqrt(cov[0, 0])


tot = {}
for sec in (12, 39, 66, 93, 103, 104):
    cut = load_cut(TIC, sec)
    g = gaia_table(TIC, btjd_to_jyear(np.median(cut['time'])))
    tgt = g.iloc[0]
    rec = {}
    centres = {}
    for name, sid in pos.items():
        r = tgt if sid is None else g[g.source_id == sid].iloc[0]
        centres[name] = (r['ra_ep'], r['dec_ep'])
    # mirror of M dwarf A through the Mira
    rA = g[g.source_id == pos['MdwarfA_18.9W']].iloc[0]
    centres['mirror_of_A'] = (2 * tgt['ra_ep'] - rA['ra_ep'], 2 * tgt['dec_ep'] - rA['dec_ep'])
    centres['N_of_Mira_19'] = (tgt['ra_ep'], tgt['dec_ep'] + 19 / 3600.)
    centres['S_of_Mira_19'] = (tgt['ra_ep'], tgt['dec_ep'] - 19 / 3600.)
    line = f'S{sec:3d}:'
    for name, c in centres.items():
        t, f, e, ex = aperture_lc(TIC, sec, r=0.75, center=c)
        fa = f * ex['median']  # e-/s
        d, de = abs_depth(t, fa)
        rec[name] = (d, de, ex['npix'])
        tot.setdefault(name, []).append((d, de))
        line += f'  {name}: {d:5.1f}+-{de:4.1f}'
    OUT[f'S{sec}'] = rec
    print(line + '  (e-/s, r=0.75px apertures)')
print('weighted means over sectors (e-/s):')
OUT['mean'] = {}
for name, v in tot.items():
    v = np.array(v)
    w = 1 / v[:, 1] ** 2
    m = np.sum(v[:, 0] * w) / np.sum(w)
    s = 1 / np.sqrt(np.sum(w))
    OUT['mean'][name] = (m, s)
    print(f'   {name:14s} {m:6.2f} +- {s:4.2f}')
save_json(OUT, f'{BASE}/results/t3_apcheck.json')
