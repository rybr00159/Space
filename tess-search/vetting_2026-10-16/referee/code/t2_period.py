"""TIC 61332759 neighbour (Gaia DR3 6711181657530577280): periodogram of its own
TESScut aperture light curve per sector; 1.4806 d vs 0.7403 d (fundamental vs
harmonic power), and Gaia DR3 variability classification."""
import sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/five_referee/code')
import numpy as np
from astropy.timeseries import LombScargle
from astroquery.gaia import Gaia
from reflib import *

TIC = 61332759
NBR = 6711181657530577280
OUT = {}
for sec in (13, 27, 67, 94, 104):
    cut = load_cut(TIC, sec)
    g = gaia_table(TIC, btjd_to_jyear(np.median(cut['time'])))
    nb = g[g.source_id == NBR].iloc[0]
    t, f, e, ex = aperture_lc(TIC, sec, r=1.0, center=(nb['ra_ep'], nb['dec_ep']))
    # remove slow trends per segment with a 3rd-order polynomial
    seg = np.concatenate([[0], np.cumsum(np.diff(t) > 0.8)])
    r = f.copy()
    for s in np.unique(seg):
        m = seg == s
        r[m] = f[m] - np.polyval(np.polyfit(t[m] - t[m].mean(), f[m], 3), t[m] - t[m].mean())
    fr = np.linspace(0.1, 4.0, 20000)
    pw = LombScargle(t, r).power(fr)
    k = np.argmax(pw)
    # 2-harmonic fit at 1.4806 d: compare odd/even minima depth
    P = 1.480598
    ph = ((t - 3148.7857) / P) % 1
    X = np.vstack([np.ones_like(t), np.sin(2 * np.pi * ph), np.cos(2 * np.pi * ph), np.sin(4 * np.pi * ph), np.cos(4 * np.pi * ph)]).T
    b, *_ = np.linalg.lstsq(X, r, rcond=None)
    A1, A2 = np.hypot(b[1], b[2]), np.hypot(b[3], b[4])
    resid = r - X @ b
    eA = np.std(resid) * np.sqrt(2 / len(t))
    OUT[sec] = dict(LS_peak_P=1 / fr[k], A1_ppm=A1 * 1e6, A2_ppm=A2 * 1e6, A_err_ppm=eA * 1e6)
    print(f'S{sec:3d}: neighbour LS peak P={1/fr[k]:.4f} d; at P=1.4806: fundamental semi-amp {A1*1e6:.0f} ppm, 2nd harmonic {A2*1e6:.0f} ppm (+-{eA*1e6:.0f})')
q = f"SELECT source_id, best_class_name, best_class_score FROM gaiadr3.vari_classifier_result WHERE source_id={NBR}"
try:
    t = Gaia.launch_job(q).get_results()
    print('Gaia DR3 vari_classifier_result:', [(int(x['source_id']), str(x['best_class_name']), float(x['best_class_score'])) for x in t])
    OUT['gaia_class'] = [(str(x['best_class_name']), float(x['best_class_score'])) for x in t]
except Exception as e:
    print('gaia class query failed', e)
save_json(OUT, f'{BASE}/results/t2_period.json')
