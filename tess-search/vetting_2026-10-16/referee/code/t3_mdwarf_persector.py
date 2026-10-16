"""TIC 214311906 field: per-sector depth S/N at the 5.69149-d ephemeris and an
independent per-sector BLS (1-13 d) on TESScut photometry centred on the M dwarf
TIC 214311908 (r=1 px), after biweight detrending."""
import sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/five_referee/code')
import numpy as np
from astropy.timeseries import BoxLeastSquares
from wotan import flatten
from reflib import *

P, T0, T14, T12 = 5.69149, 3840.913, 1.6 / 24, 0.3 / 24
OUT = {}
for sec in (12, 39, 66, 93, 103, 104):
    cut = load_cut(214311906, sec)
    g = gaia_table(214311906, btjd_to_jyear(np.median(cut['time'])))
    a = g[g.source_id == 5925373272850784768].iloc[0]
    t, f, e, ex = aperture_lc(214311906, sec, r=1.0, center=(a['ra_ep'], a['dec_ep']))
    fl = flatten(t, f, method='biweight', window_length=0.5, break_tolerance=0.5)
    m = np.isfinite(fl)
    t, fl = t[m], fl[m]
    s = robust_std(fl)
    # depth at ephemeris (all events, box of T14)
    ph = ((t - T0) / P + 0.5) % 1 - 0.5
    inn = np.abs(ph * P) < 0.4 * T14
    d = 1 - np.mean(fl[inn]); de = s / np.sqrt(inn.sum())
    bls = BoxLeastSquares(t, fl, dy=np.full(len(t), s))
    pg = bls.autopower([0.05, 0.065, 0.08], minimum_period=1.0, maximum_period=min(13.0, np.ptp(t) / 2), frequency_factor=1.0)
    k = np.argmax(pg.depth_snr)
    top = np.argsort(pg.depth_snr)[::-1]
    # rank of the true period among distinct peaks
    rank, seen = None, []
    for i in top:
        p = pg.period[i]
        if any(abs(p - q) / q < 0.01 for q in seen):
            continue
        seen.append(p)
        if abs(p - P) / P < 0.01 or abs(p - 2 * P) / (2 * P) < 0.01 or abs(p - P / 2) / (P / 2) < 0.01:
            rank = len(seen)
            break
        if len(seen) > 20:
            break
    OUT[sec] = dict(depth_ppm=d * 1e6, err_ppm=de * 1e6, n_in=int(inn.sum()), bls_best_P=float(pg.period[k]), bls_best_snr=float(pg.depth_snr[k]), rank_of_P=rank)
    print(f'S{sec:3d}: depth at ephemeris {d*1e6:6.0f}+-{de*1e6:4.0f} ppm (S/N {d/de:4.1f}); BLS best P={pg.period[k]:.4f} d (SNR {pg.depth_snr[k]:.1f}); rank of 5.69 d (or 2x, 1/2x) among distinct peaks: {rank}')
save_json(OUT, f'{BASE}/results/t3_mdwarf_persector.json')
