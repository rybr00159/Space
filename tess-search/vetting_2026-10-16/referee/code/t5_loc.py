"""TIC 61141412: difference-image localisation of the 6.53-d dips (event times
from the quadratic ephemeris fitted in t5_etv), with injection-on-target check."""
import sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/five_referee/code')
import numpy as np
from reflib import *
from diffimg import *

TIC = 61141412
ETV = json.load(open(f'{BASE}/results/t5_etv.json'))
c2 = ETV['tesscut_r2_win0.25_poly1_detNone']['deg2']['coef']
T14 = 1.0 / 24
rng = np.random.default_rng(11)
OUT = {}
rows = []
for sec in (13, 27, 67, 94, 104):
    cut = load_cut(TIC, sec)
    t = cut['time']
    n0, n1 = int(np.floor((t.min() - c2[0]) / c2[1])) - 1, int(np.ceil((t.max() - c2[0]) / c2[1])) + 1
    tcs = np.array([c2[0] + c2[1] * n + c2[2] * n * n for n in range(n0, n1 + 1)])
    tcs = tcs[(tcs > t.min() + 0.1) & (tcs < t.max() - 0.1)]
    evs, direct = event_diff_images(cut, tcs, T14, fin=0.35, gap=0.8, wout=0.15)
    if len(evs) < 2:
        continue
    D, E = combine(evs)
    g = gaia_table(TIC, btjd_to_jyear(np.median(t)))
    tgt = g.iloc[0]
    xt, yt = [float(v) for v in cut['wcs'].world_to_pixel_values(tgt['ra_ep'], tgt['dec_ep'])]
    prf = PRFModel(cut)
    derr = np.sqrt((0.002 * np.abs(direct)) ** 2 + 5 ** 2)
    ft = fit_prf(direct, derr, prf, xt, yt, box=2)
    fd = fit_prf(D, E, prf, xt, yt, box=2)
    o = pix_to_EN(cut, fd['x'], fd['y'], ft['x'], ft['y'])
    bs = []
    for b in range(150):
        idx = rng.integers(0, len(evs), len(evs))
        Db, Eb = combine(evs, idx)
        fb = fit_prf(Db, Eb, prf, fd['x'], fd['y'], box=2)
        bs.append(pix_to_EN(cut, fb['x'], fb['y'], ft['x'], ft['y']))
    bs = np.array(bs)
    rsd = 1.4826 * np.median(np.abs(bs - np.median(bs, axis=0)), axis=0)
    # injection of the same deficit on the target (direct image as template) at phase +0.5
    tmpl = direct - np.median(direct)
    tci = tcs + c2[1] / 2
    tci = tci[(tci > t.min() + 0.1) & (tci < t.max() - 0.1)]
    evi, di = event_diff_images(cut, tci, T14, fin=0.35, gap=0.8, wout=0.15, inject=(tmpl, fd['A'] / ft['A'] / 0.8))
    oi = (np.nan, np.nan)
    if len(evi) >= 2:
        Di, Ei = combine(evi)
        fdi = fit_prf(di, derr, prf, xt, yt, box=2)
        fi = fit_prf(Di, Ei, prf, xt, yt, box=2)
        oi = pix_to_EN(cut, fi['x'], fi['y'], fdi['x'], fdi['y'])
    rows.append((sec, o[0], o[1], rsd[0], rsd[1], oi[0], oi[1], fd['A'] / ft['A']))
    print(f'S{sec:3d}: n={len(evs)} offset ({o[0]:+5.1f},{o[1]:+5.1f})" robust-boot sd ({rsd[0]:.1f},{rsd[1]:.1f}); injection-on-target recovers ({oi[0]:+.1f},{oi[1]:+.1f})"; depth(PRF) {fd["A"]/ft["A"]*1e6:.0f} ppm')
r = np.array(rows)
w = 1 / np.maximum(r[:, 3:5], 0.3) ** 2
d = r[:, 1:3] - np.nan_to_num(r[:, 5:7])
mE = np.sum(d[:, 0] * w[:, 0]) / np.sum(w[:, 0]); mN = np.sum(d[:, 1] * w[:, 1]) / np.sum(w[:, 1])
sE = 1 / np.sqrt(np.sum(w[:, 0])); sN = 1 / np.sqrt(np.sum(w[:, 1]))
print(f'combined (offset - injection bias): dE={mE:+.1f}+-{sE:.1f}", dN={mN:+.1f}+-{sN:.1f}"')
OUT = dict(per_sector=r.tolist(), combined=[mE, sE, mN, sN])
save_json(OUT, f'{BASE}/results/t5_loc.json')
