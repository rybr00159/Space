"""Pointing-jitter audit and per-transit sign test.

(a) Pointing: per-cadence (dx,dy) of the stamp from stars at 2.5-9 px from the target (linear
    regression of each frame on the mean image and its PSF-model gradients; pipeline.pointing_series).
    For each real transit window the same baseline operator as for the images gives the
    in-transit pointing shift delta = baseline(in) - mean(in).  A pointing shift delta produces a
    target-dominated dipole that moves the fitted difference-image position by about
    -delta * F_target / A (pixels); this predicted spurious offset is projected on the T->N axis.
(b) Momentum dumps (QUALITY bit 32) and other flagged cadences inside each window.
(c) Single-transit calibration: for every fake window, inject on the target / on the neighbour
    and fit the single-window difference image -> P(s_par > 0.5 | host) for one transit.
    The probability that ALL real transits have s_par > 0.5 under each hypothesis follows.
Outputs: jitter_pertransit.json, fig_jitter_pertransit.png
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pipeline import *

real = json.load(open(f'{HERE}/diff_prf_main.json'))
reg = json.load(open(f'{HERE}/registration.json'))
out = {}
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
logp_all_T = 0.0; logp_all_N = 0.0
cols = dict(zip(SECTORS, ['tab:blue', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown']))
for sec in SECTORS:
    ctx = SectorContext(sec, psfkind='prf', reg_variant='local_r4', rfit=3.0)
    S = ctx.S
    dxs, dys, sc = pointing_series(ctx)
    FT = reg[str(sec)]['variants']['prf:local_r4']['F_target']
    A = real[str(sec)]['stacked']['free_amp']
    wins = ctx.real_windows()
    vT = ctx.jac @ ctx.vecTN
    rows = []
    for w, p in zip(wins, real[str(sec)]['per_transit']):
        cube = np.stack([dxs, dys], axis=1)[:, :, None]  # (ncad, 2, 1)
        dlt, _, rms = difference_image(S.t, cube, w)
        delta = -dlt[:, 0]          # mean(in) - baseline(in): pointing shift of in-transit frames, px
        spur_pix = -delta * FT / A  # predicted spurious shift of the fitted position (px)
        spur_sky = ctx.jac @ spur_pix
        s_spur = float(spur_sky @ vT / (vT @ vT))
        idx = np.arange(len(S.t))
        inwin = (S.t > w['tc'] - OOT_OUT) & (S.t < w['tc'] + OOT_OUT)
        dumps = S.t[inwin & ((S.q & 32) > 0)]
        flagged = int(np.sum(inwin & (S.q != 0)))
        rows.append(dict(epoch=w['epoch'], delta_px=delta.tolist(), delta_rms_px=rms[:, 0].tolist(),
                         spurious_s_par=s_spur, s_par=p['s_par'], s_par_corrected=p['s_par'] - s_spur,
                         dumps_in_window_dt=(dumps - w['tc']).tolist(), n_flagged_in_window=flagged))
    # (c) single-transit calibration
    fw = ctx.fake_windows()
    shape_mean = np.mean([ctx.in_window_shape_mean(w) for w in wins])
    amp = A / shape_mean
    sT, sN = [], []
    for w in fw:
        D, sg, _ = ctx.diff(w)
        aT = ctx.analyse(D + ctx.inject_diff(w, ctx.pT, amp), sg, free=True, extra=False)
        aN = ctx.analyse(D + ctx.inject_diff(w, ctx.pN, amp), sg, free=True, extra=False)
        sT.append(aT['s_par']); sN.append(aN['s_par'])
    sT = np.array(sT); sN = np.array(sN)
    pT = float(np.mean(sT > 0.5)); pN = float(np.mean(sN > 0.5))
    s_real = np.array([r['s_par'] for r in rows])
    k = int(np.sum(s_real > 0.5)); n = len(s_real)
    from scipy.stats import binom
    # probability of >=k of n above 0.5 under each hypothesis (binomial with per-transit rate)
    pT_eff = max(pT, 0.5 / len(sT))  # floor: half a count if none of the fakes exceeded 0.5
    P_T = float(binom.sf(k - 1, n, pT_eff)); P_N = float(binom.sf(k - 1, n, pN))
    logp_all_T += np.log10(P_T); logp_all_N += np.log10(P_N)
    out[str(sec)] = dict(F_target=FT, A=A, leverage_FT_over_A=FT / A, transits=rows,
                         pointing_series_rms_px=[float(np.std(dxs[S.good])), float(np.std(dys[S.good]))],
                         single_transit_cal=dict(n_fake=len(fw), P_spar_gt_half_T=pT, P_spar_gt_half_N=pN,
                                                 sT_mean=float(np.mean(sT)), sT_std=float(np.std(sT)),
                                                 sN_mean=float(np.mean(sN)), sN_std=float(np.std(sN))),
                         n_real=n, n_real_gt_half=k, P_ge_k_T=P_T, P_ge_k_N=P_N)
    print(f"S{sec}: F_T/A={FT/A:.0f}; per-transit delta (mpx): " +
          ' '.join(f"({r['delta_px'][0]*1e3:+.2f},{r['delta_px'][1]*1e3:+.2f})" for r in rows) +
          f"\n   spurious s_par: " + ' '.join(f"{r['spurious_s_par']:+.2f}" for r in rows) +
          f"\n   real s_par:     " + ' '.join(f"{r['s_par']:+.2f}" for r in rows) +
          f"\n   dumps in windows: {[len(r['dumps_in_window_dt']) for r in rows]}"
          f"\n   single-transit cal: H_T s_par={np.mean(sT):+.2f}+-{np.std(sT):.2f} P(>0.5)={pT:.3f};  H_N {np.mean(sN):+.2f}+-{np.std(sN):.2f} P(>0.5)={pN:.3f}"
          f"\n   real: {k}/{n} transits with s_par>0.5 -> P(H_T)={P_T:.2e}  P(H_N)={P_N:.2f}", flush=True)
    ax[0].hist(sT, np.linspace(-2, 3, 41), histtype='step', color=cols[sec], ls='--')
    ax[0].hist(sN, np.linspace(-2, 3, 41), histtype='step', color=cols[sec])
    ax[0].plot(s_real, np.full(n, -0.5 - SECTORS.index(sec) * 0.6), 'v', color=cols[sec], ms=8, label=f'S{sec} real transits')
    ax[1].plot([r['spurious_s_par'] for r in rows], s_real, 'o', color=cols[sec], label=f'S{sec}')
out['combined'] = dict(log10_P_all_T=logp_all_T, log10_P_all_N=logp_all_N)
print('combined sign test: log10 P(H_T) =', logp_all_T, ' log10 P(H_N) =', logp_all_N)
ax[0].set_xlabel('single-transit s_par'); ax[0].set_ylabel('fake windows')
ax[0].set_title('single-transit calibration: dashed = injected on target, solid = on neighbour')
ax[0].legend(fontsize=7)
ax[1].axhline(0, color='r'); ax[1].axhline(1, color='m'); ax[1].axvline(0, color='k', lw=0.5)
ax[1].set_xlabel('predicted spurious s_par from measured in-transit pointing shift'); ax[1].set_ylabel('measured per-transit s_par')
ax[1].set_title('pointing-jitter audit (PRF)'); ax[1].legend(fontsize=7)
plt.tight_layout(); plt.savefig(f'{HERE}/fig_jitter_pertransit.png', dpi=110); plt.close()
to_json(out, f'{HERE}/jitter_pertransit.json')
