"""TIC 214311906 and TIC 55829034: per-sector reality/periodicity of the dips,
fractional vs absolute (e-/s) depths, per-sector BLS, odd/even, phase 0.5,
and the physical (stellar density / orbit size) argument.
Usage: python -I t34_events.py <TIC>"""
import sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/five_referee/code')
import numpy as np
from astropy.timeseries import BoxLeastSquares
from wotan import flatten
from reflib import *

TIC = int(sys.argv[1])
CFG = {214311906: dict(P=5.69149, T0=3840.913, T14=1.55 / 24, T12=0.25 / 24, secs=(12, 39, 66, 93, 103, 104)),
       55829034: dict(P=0.616146, T0=3837.8667, T14=2.0 / 24, T12=1.0 / 24, secs=(13, 66, 93, 104))}[TIC]
P, T0, T14, T12 = CFG['P'], CFG['T0'], CFG['T14'], CFG['T12']
OUT = dict(cfg=CFG, sectors={})
QLP_SHIFT = {103: None, 104: None}  # filled from TESScut TIMECORR below


def joint_depth(t, f, centres, win, texp, absolute=False):
    blocks = []
    for c in centres:
        s = np.abs(t - c) < win
        if s.sum() > 8 and np.sum(np.abs(t[s] - c) < T14 / 2) >= 2 and np.sum(t[s] < c - T14 / 2) >= 3 and np.sum(t[s] > c + T14 / 2) >= 3:
            blocks.append((c, s))
    nW = len(blocks)
    if nW == 0:
        return np.nan, np.nan, 0, []
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
    # per-event depths (separately)
    per = []
    for k, (c, s) in enumerate(blocks):
        dt = t[s] - c
        sh = 1 - supersample(trapezoid, t[s], texp, t0=c, depth=1.0, T14=T14, T12=T12)
        Xk = np.vstack([-sh, np.ones_like(dt), dt, dt ** 2]).T
        bk, *_ = np.linalg.lstsq(Xk, f[s], rcond=None)
        rk = f[s] - Xk @ bk
        ck = np.linalg.pinv(Xk.T @ Xk) * np.var(rk, ddof=4)
        per.append((c, bk[0], np.sqrt(ck[0, 0])))
    return b[0], np.sqrt(cov[0, 0]), nW, per


for sec in CFG['secs']:
    rec = {}
    cut = load_cut(TIC, sec)
    tc_corr = None
    for name, getter in [('CUT_r1.5', lambda: aperture_lc(TIC, sec, r=1.5)), ('CUT_r2.5', lambda: aperture_lc(TIC, sec, r=2.5)),
                         ('QLP_SAP', lambda: load_lc(TIC, sec, 'QLP', 'SAP_FLUX'))]:
        r = getter()
        if r is None:
            continue
        t, f, e, ex = r
        if name == 'QLP_SAP' and sec in (103, 104):
            from astropy.io import fits as _f
            tcorr = np.nanmedian(_f.open(cut_path(TIC, sec))[1].data['TIMECORR'])
            t = t + tcorr
        texp = np.median(np.diff(t))
        n0, n1 = int(np.floor((t.min() - T0) / P)), int(np.ceil((t.max() - T0) / P))
        cen = T0 + np.arange(n0, n1 + 1) * P
        win = max(0.25, 2.5 * T14)
        d, de, nW, per = joint_depth(t, f, cen, win, texp)
        # absolute depth in flux units (e-/s for CUT; QLP SAP normalised so use median)
        absd = d * ex['median']
        absde = de * ex['median']
        # half-phase and odd/even
        d5, de5, n5, _ = joint_depth(t, f, cen + P / 2, win, texp)
        # per-sector BLS on detrended data (independent period search)
        fl = flatten(t, f, method='biweight', window_length=max(0.5, 4 * T14), break_tolerance=0.5)
        m = np.isfinite(fl)
        bls = BoxLeastSquares(t[m], fl[m], dy=np.full(m.sum(), robust_std(fl[m])))
        pmax = min(15.0, np.ptp(t) / 2)
        pg = bls.autopower([T14 * 0.8, T14, T14 * 1.3], minimum_period=0.3 if TIC == 55829034 else 1.0, maximum_period=pmax, frequency_factor=1.0)
        k = np.argmax(pg.depth_snr)
        # where does the input period rank?
        kp = np.argmin(np.abs(pg.period - P))
        # check harmonics/subharmonics of P near the top
        top = np.argsort(pg.depth_snr)[::-1][:5]
        rec[name] = dict(depth_ppm=d * 1e6, depth_err_ppm=de * 1e6, n_events=nW, abs_depth=absd, abs_depth_err=absde, median_flux=ex['median'],
                         half_phase_ppm=d5 * 1e6, half_phase_err_ppm=de5 * 1e6,
                         per_event=[dict(tc=c, depth_ppm=dd * 1e6, err_ppm=ee * 1e6) for c, dd, ee in per],
                         bls_best_P=float(pg.period[k]), bls_best_snr=float(pg.depth_snr[k]), bls_snr_at_P=float(pg.depth_snr[kp]),
                         bls_top5=[(float(pg.period[i]), float(pg.depth_snr[i])) for i in top])
        print(f"S{sec:3d} {name:8s} median flux {ex['median']:9.1f}  depth {d*1e6:6.0f}+-{de*1e6:4.0f} ppm (n={nW})  abs {absd:7.2f}+-{absde:5.2f}  "
              f"half-phase {d5*1e6:6.0f}+-{de5*1e6:4.0f}  BLS best P={pg.period[k]:.4f} (SNR {pg.depth_snr[k]:.1f}), SNR at P {pg.depth_snr[kp]:.1f}")
        if name == 'CUT_r2.5':
            print('        per-event depths (ppm): ' + ', '.join(f'{c:.2f}:{dd*1e6:.0f}+-{ee*1e6:.0f}' for c, dd, ee in per))
    OUT['sectors'][sec] = rec
save_json(OUT, f'{BASE}/results/t34_events_{TIC}.json')
