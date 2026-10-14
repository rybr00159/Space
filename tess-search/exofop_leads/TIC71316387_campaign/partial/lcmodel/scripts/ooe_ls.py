"""Lomb-Scargle of the out-of-eclipse light curve (eclipses masked, per-segment Legendre trend with ~4-d knots removed),
per sector and for the S105+S106 pair, to check whether the strongest periodicity in 0.4-6 d is at P_EB."""
import numpy as np, sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/scripts')
import ooe_var as O, lclib as L
from astropy.timeseries import LombScargle
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
P = O.P
sets = [(1, O.load_spoc_col(1, 'PDCSAP_FLUX')), (28, O.load_spoc_col(28, 'PDCSAP_FLUX')), (68, O.load_spoc_col(68, 'PDCSAP_FLUX')),
        (95, L.load_tc(95)), (105, L.load_tc(105)), (106, L.load_tc(106))]
out = {}
fig, axs = plt.subplots(6, 1, figsize=(8, 12), sharex=True)
freq = np.linspace(1 / 6.0, 1 / 0.4, 20000)
for ax, (sec, lc) in zip(axs, sets):
    t, f, seg = O.prep(lc)
    A = O.design_trend(t, seg, 0, dper=4.0)
    c, *_ = np.linalg.lstsq(A, f - 1, rcond=None); r = f - 1 - A @ c
    pw = LombScargle(t, r).power(freq)
    j = np.argmax(pw); pk = 1 / freq[j]
    # power at P_EB and P_EB/2
    jP = np.argmin(np.abs(freq - 1 / P)); jP2 = np.argmin(np.abs(freq - 2 / P))
    fap = LombScargle(t, r).false_alarm_probability(pw[j])
    out[sec] = dict(peak_period=float(pk), peak_power=float(pw[j]), power_at_PEB=float(pw[jP]), power_at_PEB2=float(pw[jP2]), fap_peak=float(fap),
                    rank_PEB=int(np.sum(pw > pw[jP])))
    print(f'S{sec}: highest peak at {pk:.4f} d (power {pw[j]:.4f}, FAP {fap:.1e}); power at P_EB {pw[jP]:.4f}; at P_EB/2 {pw[jP2]:.4f}', flush=True)
    ax.plot(1 / freq, pw, 'k-', lw=0.6); ax.axvline(P, color='r', ls=':'); ax.axvline(P / 2, color='r', ls=':'); ax.set_ylabel(f'S{sec}')
axs[-1].set_xscale('log'); axs[-1].set_xlabel('period (d)')
plt.tight_layout(); plt.savefig(L.HERE + '/plots/ooe_ls.png', dpi=100)
json.dump(out, open(L.HERE + '/results/ooe_ls.json', 'w'), indent=1)
