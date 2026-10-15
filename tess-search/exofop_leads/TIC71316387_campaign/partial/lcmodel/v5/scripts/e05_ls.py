"""Lomb-Scargle periodograms of the eclipse-masked out-of-eclipse light curves (after removing a 3-d-knot spline
trend per orbit segment, as in e01_ooe.py): peak period near P_EB per sector and for the contiguous S105+S106
pair, with an empirical false-alarm level from 300 phase-scrambled-by-segment (circularly shifted) realisations.
Also the first-harmonic semi-amplitude at the peak and at exactly P_EB."""
import numpy as np, sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v5/lib')
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v5/scripts')
import data as D
from astropy.timeseries import LombScargle
import importlib.util
spec = importlib.util.spec_from_file_location('ooe', '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v5/scripts/e01_ooe.py')
sys.argv = [sys.argv[0], '3.0']
O = importlib.util.module_from_spec(spec); spec.loader.exec_module(O)
P = D.P0
rng = np.random.default_rng(8)
def detrend(t, f):
    B = O.spline_basis(t, 3.0)
    c, *_ = np.linalg.lstsq(B, f, rcond=None)
    r = f / (B @ c) - 1
    s = 1.4826 * np.median(np.abs(r)); m = np.abs(r) < 5 * s
    return t[m], r[m]
out = {}
sets = {'SPOC-PDC1': [D.load_spoc(1)], 'SPOC-PDC28': [D.load_spoc(28)], 'SPOC-PDC68': [D.load_spoc(68)], 'TC95': [D.load_tc(95)],
        'TC105': [D.load_tc(105)], 'TC106': [D.load_tc(106)], 'TC105+106': [D.load_tc(105), D.load_tc(106)]}
freq = np.linspace(1 / 5.0, 1 / 0.4, 6000)
for name, lcs in sets.items():
    T, R = [], []
    for lc in lcs:
        t, f = O.prep(D.clip(lc)); t, r = detrend(t, f); T.append(t); R.append(r)
    t = np.concatenate(T); r = np.concatenate(R)
    ls = LombScargle(t, r); pw = ls.power(freq)
    near = (1 / freq > 1.0) & (1 / freq < 1.7)
    j = np.argmax(np.where(near, pw, 0)); Ppk = 1 / freq[j]
    jall = np.argmax(pw)
    # empirical FAP: circular shifts of the residuals within the time series destroy the phase coherence with P
    mx = []
    for k in range(300):
        rr = np.roll(r, rng.integers(r.size // 10, r.size - r.size // 10))
        mx.append(LombScargle(t, rr).power(freq[near]).max())
    A = np.sqrt(4 * pw[j] * np.var(r) / 2) if False else None
    # semi-amplitude of a sinusoid at P_EB and at the peak
    def amp(Pf):
        ph = 2 * np.pi * t / Pf; X = np.vstack([np.ones_like(t), np.cos(ph), np.sin(ph)]).T
        c, *_ = np.linalg.lstsq(X, r, rcond=None); return float(np.hypot(c[1], c[2]))
    hw = 1 / (t.max() - t.min())   # frequency resolution
    out[name] = dict(P_peak=float(Ppk), P_peak_unc_d=float(Ppk ** 2 * hw / 2), power_peak=float(pw[j]), power_null_99=float(np.percentile(mx, 99)),
                     power_null_max=float(np.max(mx)), P_global_peak=float(1 / freq[jall]), amp_at_peak_ppm=amp(Ppk) * 1e6, amp_at_PEB_ppm=amp(P) * 1e6,
                     PEB_minus_peak_in_resolution=float((1 / P - 1 / Ppk) / hw))
    o = out[name]
    print(f"{name:10s} peak P={Ppk:.4f} d (+-{o['P_peak_unc_d']:.3f}, global peak {o['P_global_peak']:.3f}) power {pw[j]:.4f} vs scrambled 99% {o['power_null_99']:.4f} max {o['power_null_max']:.4f}; "
          f"amp at peak {o['amp_at_peak_ppm']:.0f} ppm, at P_EB {o['amp_at_PEB_ppm']:.0f} ppm; (f_EB - f_peak)/resolution = {o['PEB_minus_peak_in_resolution']:+.2f}", flush=True)
json.dump(out, open(f'{D.V}/results/ooe_ls.json', 'w'), indent=1)
