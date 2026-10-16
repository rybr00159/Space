"""TIC 61332759: difference images phased on the 0.7403-d 'transit' ephemeris
(in-dip minus local out-of-dip), PRF-fit position vs Gaia target/neighbour."""
import sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/five_referee/code')
import numpy as np
from reflib import *
from diffimg import *

TIC = 61332759
NBR = 6711181657530577280
P2, T02 = 0.7402990, 3148.78573
T14 = 2.8 / 24
OUT = {}
for sec in (13, 27, 67, 94, 104):
    cut = load_cut(TIC, sec)
    t = cut['time']
    if sec == 104:
        pass  # TESScut times are barycentric; ephemeris from QLP (S104 offset ~6 min) irrelevant for a 2.8-h dip
    n0, n1 = int(np.floor((t.min() - T02) / P2)), int(np.ceil((t.max() - T02) / P2))
    tcs = T02 + np.arange(n0, n1 + 1) * P2
    tcs = tcs[(tcs > t.min() + 0.2) & (tcs < t.max() - 0.2)]
    evs, direct = event_diff_images(cut, tcs, T14, fin=0.3, gap=0.6, wout=0.12)
    D, E = combine(evs)
    g = gaia_table(TIC, btjd_to_jyear(np.median(t)))
    tgt = g.iloc[0]; nb = g[g.source_id == NBR].iloc[0]
    xt, yt = [float(v) for v in cut['wcs'].world_to_pixel_values(tgt['ra_ep'], tgt['dec_ep'])]
    xn, yn = [float(v) for v in cut['wcs'].world_to_pixel_values(nb['ra_ep'], nb['dec_ep'])]
    prf = PRFModel(cut)
    derr = np.sqrt((0.002 * np.abs(direct)) ** 2 + 5 ** 2)
    ft = fit_prf(direct, derr, prf, xt, yt, box=2)
    ny, nx = D.shape
    yy, xx = np.mgrid[:ny, :nx]
    snr = D / E
    cen = np.hypot(xx - xt, yy - yt) < 5
    j = np.unravel_index(np.argmax(np.where(cen, snr, -np.inf)), D.shape)
    fd = fit_prf(D, E, prf, j[1], j[0], box=2)
    o = pix_to_EN(cut, fd['x'], fd['y'], ft['x'], ft['y'])
    # S/N of the deficit at the target pixel and at the neighbour pixel
    it = (int(round(yt)), int(round(xt))); inb = (int(round(yn)), int(round(xn)))
    rec = dict(n_events=len(evs), fit_offset_from_target=list(o), gaia_nbr_offset=[nb['dE'], nb['dN']],
               snr_target_pix=float(snr[it]), snr_nbr_pix=float(snr[inb]), peak_pix=[int(j[1]), int(j[0])],
               frac_deficit_target_pix_ppm=float(D[it] / direct[it] * 1e6), frac_deficit_nbr_pix_ppm=float(D[inb] / direct[inb] * 1e6))
    OUT[f'S{sec}'] = rec
    print(f"S{sec}: n={len(evs)} dip-phased diff image PRF fit at ({o[0]:+.1f},{o[1]:+.1f})\" from target; Gaia neighbour at ({nb['dE']:+.1f},{nb['dN']:+.1f})\"; "
          f"S/N at target pix {snr[it]:.1f}, at neighbour pix {snr[inb]:.1f}; fractional deficit target pix {D[it]/direct[it]*1e6:.0f} ppm, nbr pix {D[inb]/direct[inb]*1e6:.0f} ppm")
save_json(OUT, f'{BASE}/results/t2_dipdiff.json')
