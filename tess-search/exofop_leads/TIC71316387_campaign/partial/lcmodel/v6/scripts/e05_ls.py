"""Lomb-Scargle periodograms of the eclipse-masked out-of-eclipse light curves (after removing a 3-d-knot spline
trend per orbit segment, as in e01_ooe.py), 0.2-2.5 c/d.
Is the out-of-eclipse modulation at P_EB?  Per sector (and S105+S106 combined):
  * period of the highest peak in 1.0-1.7 d and of the global peak, with the frequency resolution 1/T;
  * empirical significance of the power at f_EB from the SAME periodogram: fraction of local maxima at frequencies
    away from f_EB and its harmonics/subharmonics (|f - j f_EB/2| > 3/T) with power >= the power at f_EB, and the
    number of independent frequencies, giving p_any = 1 - (1 - frac)^1 (chance that noise peaks reach that power at
    one pre-specified frequency) -- the frequency is fixed in advance (P_EB), so no look-elsewhere factor applies;
  * chance that the global peak falls within +-0.5/T of f_EB if peaks were randomly placed: (1/T)/(band width).
(The v5 'circular-shift' null was wrong: rolling a periodic residual series preserves its periodicity.)"""
import numpy as np, sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6/lib')
import data as D
from astropy.timeseries import LombScargle
import importlib.util
spec = importlib.util.spec_from_file_location('ooe', '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6/scripts/e01_ooe.py')
sys.argv = [sys.argv[0], '3.0']
O = importlib.util.module_from_spec(spec); spec.loader.exec_module(O)
P = D.P0; fEB = 1 / P

def detrend(t, f):
    B = O.spline_basis(t, 3.0)
    c, *_ = np.linalg.lstsq(B, f, rcond=None)
    r = f / (B @ c) - 1
    s = 1.4826 * np.median(np.abs(r)); m = np.abs(r) < 5 * s
    return t[m], r[m]

out = {}
sets = {'SPOC-PDC1': [D.load_spoc(1)], 'SPOC-PDC28': [D.load_spoc(28)], 'SPOC-PDC68': [D.load_spoc(68)], 'TC95': [D.load_tc(95)],
        'TC105': [D.load_tc(105)], 'TC106': [D.load_tc(106)], 'TC105+106': [D.load_tc(105), D.load_tc(106)]}
fmin, fmax = 0.2, 2.5
freq = np.linspace(fmin, fmax, 20000)
for name, lcs in sets.items():
    T, R = [], []
    for lc in lcs:
        t, f = O.prep(D.clip(lc)); t, r = detrend(t, f); T.append(t); R.append(r)
    t = np.concatenate(T); r = np.concatenate(R)
    Tspan = t.max() - t.min(); res_f = 1 / Tspan
    ls = LombScargle(t, r); pw = ls.power(freq)
    pEB = float(ls.power(np.array([fEB]))[0])
    near = (1 / freq > 1.0) & (1 / freq < 1.7)
    j = np.argmax(np.where(near, pw, 0)); Ppk = 1 / freq[j]
    jall = np.argmax(pw)
    # local maxima away from f_EB multiples (j f_EB / 2, j = 1..6)
    lm = np.where((pw[1:-1] > pw[:-2]) & (pw[1:-1] > pw[2:]))[0] + 1
    away = np.ones(lm.size, bool)
    for jj in range(1, 7):
        away &= np.abs(freq[lm] - jj * fEB / 2) > 3 * res_f
    frac = float(np.mean(pw[lm][away] >= pEB))
    n_indep = int((fmax - fmin) / res_f)
    def amp(Pf):
        ph = 2 * np.pi * t / Pf; X = np.vstack([np.ones_like(t), np.cos(ph), np.sin(ph)]).T
        c, *_ = np.linalg.lstsq(X, r, rcond=None); return float(np.hypot(c[1], c[2]))
    out[name] = dict(T_span_d=float(Tspan), freq_resolution=float(res_f), P_peak_1to1p7=float(Ppk), P_global_peak=float(1 / freq[jall]),
                     P_global_peak_minus_PEB_in_res=float((freq[jall] - fEB) / res_f), power_at_PEB=pEB, power_global_peak=float(pw[jall]),
                     n_local_max_away=int(away.sum()), frac_local_max_away_ge_PEB=frac, max_power_away=float(pw[lm][away].max()),
                     n_indep_freq=n_indep, p_global_peak_within_half_res_by_chance=float(res_f / (fmax - fmin)),
                     amp_at_PEB_ppm=amp(P) * 1e6, amp_at_peak_ppm=amp(Ppk) * 1e6)
    o = out[name]
    print(f"{name:10s} global peak P={o['P_global_peak']:.4f} d ((f-fEB)/res {o['P_global_peak_minus_PEB_in_res']:+.2f}); power at P_EB {pEB:.4f}, "
          f"max power of {o['n_local_max_away']} peaks away from f_EB harmonics {o['max_power_away']:.4f}, fraction >= P_EB power {frac:.3f}; "
          f"amp at P_EB {o['amp_at_PEB_ppm']:.0f} ppm", flush=True)
# combined: number of independent sectors whose global peak lies within 0.5 resolution of f_EB, and binomial chance
secs = ['SPOC-PDC1', 'SPOC-PDC28', 'SPOC-PDC68', 'TC95', 'TC105', 'TC106']
k = sum(abs(out[s]['P_global_peak_minus_PEB_in_res']) < 0.5 for s in secs)
pch = np.mean([out[s]['p_global_peak_within_half_res_by_chance'] for s in secs])
from scipy.stats import binom
out['combined'] = dict(n_sectors=len(secs), n_global_peak_at_PEB=int(k), p_single=float(pch), p_binomial_ge_k=float(binom.sf(k - 1, len(secs), pch)))
print('sectors with the global peak within 0.5 resolution of f_EB:', k, 'of', len(secs), '; chance per sector', round(pch, 4), '; binomial p(>=k)', out['combined']['p_binomial_ge_k'])
json.dump(out, open(f'{D.V}/results/ooe_ls.json', 'w'), indent=1)
