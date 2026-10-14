"""Likelihood of the target's Gaia IPD signature as a function of companion separation and contrast, from the
SOAR speckle binaries (empirical), plus the harmonic-phase / PA calibration.

Target: G = 11.78, RUWE = 1.026, ipd_frac_multi_peak = 41, BP blended 47/47, RP blended 48/48,
        ipd_gof_harmonic_amplitude = 0.029, phase = 18.24 deg, ipd_frac_odd_win = 0.

'Match' definition (primary): IPDfmp in [28, 56], RUWE < 1.4, BP blended fraction >= 0.8.
Variants: without the blend condition; RUWE < 1.2.
Likelihood per (sep, dm) bin = fraction of speckle binaries in that bin (G 9.5-13.5, Gaia-unresolved) matching,
with a Bayesian (Jeffreys) interval. Single controls give the false-alarm rate.
"""
import numpy as np, pandas as pd, json
from scipy.stats import beta
O = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v2/speckle'
R = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v2'
S = pd.read_csv(f'{O}/calib_sample_v06.csv')
B = S[(S.kind == 'binary') & (S.Gmag > 9.5) & (S.Gmag < 13.5) & ~S.gaia_resolved.astype(bool)].copy()
C = S[(S.kind == 'single') & (S.Gmag > 9.5) & (S.Gmag < 13.5) & (S.dm1 >= 3.5)].copy()
out = {}


def jeff(k, n):
    return [float(beta.ppf(0.16, k + 0.5, n - k + 0.5)), float(beta.ppf(0.84, k + 0.5, n - k + 0.5))]


defs = {
    'fmp+ruwe+blend': lambda d: (d.IPDfmp >= 28) & (d.IPDfmp <= 56) & (d.RUWE < 1.4) & (d.bfrac >= 0.8),
    'fmp+ruwe': lambda d: (d.IPDfmp >= 28) & (d.IPDfmp <= 56) & (d.RUWE < 1.4),
    'fmp+ruwe1.2+blend': lambda d: (d.IPDfmp >= 28) & (d.IPDfmp <= 56) & (d.RUWE < 1.2) & (d.bfrac >= 0.8),
}
sepb = [0.03, 0.15, 0.3, 0.45, 0.6, 0.8, 1.0, 1.3, 2.0]
dmb = [0, 1.5, 2.5, 3.5, 4.5, 7]
for name, f in defs.items():
    mB = f(B); mC = f(C)
    print(f'\n== {name}: singles {mC.sum()}/{len(C)} = {mC.mean():.4f} (68% {np.round(jeff(mC.sum(), len(C)),4)}); binaries {mB.sum()}/{len(B)}')
    grid = []
    print('   sep\\dm   ' + ''.join(f'{dmb[j]:.1f}-{dmb[j+1]:.1f}     ' for j in range(len(dmb) - 1)))
    for i in range(len(sepb) - 1):
        line = f'  {sepb[i]:.2f}-{sepb[i+1]:.2f} '
        for j in range(len(dmb) - 1):
            s = (B.sep >= sepb[i]) & (B.sep < sepb[i + 1]) & (B.dm >= dmb[j]) & (B.dm < dmb[j + 1])
            n = int(s.sum()); k = int((s & mB).sum())
            grid.append(dict(sep=[sepb[i], sepb[i + 1]], dm=[dmb[j], dmb[j + 1]], n=n, k=k, frac=(k / n if n else None),
                             ci68=(jeff(k, n) if n else None)))
            line += f' {k:2d}/{n:3d}={k/n:4.2f}' if n else '    -        '
        print(line)
    sel = B[mB]
    out[name] = dict(singles_k=int(mC.sum()), singles_n=int(len(C)), singles_rate=float(mC.mean()), singles_ci68=jeff(mC.sum(), len(C)),
                     n_match=int(len(sel)), grid=grid,
                     sep_pct=[float(q) for q in np.percentile(sel.sep, [5, 16, 50, 84, 95])] if len(sel) else None,
                     dm_pct=[float(q) for q in np.percentile(sel.dm, [5, 16, 50, 84, 95])] if len(sel) else None,
                     matches=sel[['sep', 'dm', 'Gmag', 'RUWE', 'IPDfmp', 'bfrac', 'IPDgofha', 'IPDgofhp', 'pa']].round(3).to_dict('records'))
    if len(sel):
        print(f'   matched binaries: n={len(sel)} sep 5/16/50/84/95% {np.round(out[name]["sep_pct"],2)}  dm {np.round(out[name]["dm_pct"],2)}')

# harmonic phase calibration: d' = (phase - PA - 90) wrapped to [-90, 90)
H = B[(B.sep > 0.2) & (B.sep < 2.0) & np.isfinite(B.pa)].copy()
H['dp'] = ((H.IPDgofhp - (H.pa % 180) - 90) + 90) % 180 - 90
print('\nharmonic phase: d\' = phase - PA - 90 (mod 180)')
ph = []
for lo, hi in [(0, 0.015), (0.015, 0.025), (0.025, 0.035), (0.035, 0.05), (0.05, 0.1), (0.1, 1)]:
    s = H[(H.IPDgofha >= lo) & (H.IPDgofha < hi)]
    if len(s) < 5:
        continue
    f30 = float((np.abs(s.dp) < 30).mean()); f45 = float((np.abs(s.dp) < 45).mean())
    print(f'   gofha {lo:.3f}-{hi:.3f}: n={len(s):4d}  frac |d\'|<30 {f30:.2f} (random 0.33)  <45 {f45:.2f} (random 0.50)')
    ph.append(dict(amp=[lo, hi], n=int(len(s)), f30=f30, f45=f45))
s = H[(H.IPDgofha >= 0.02) & (H.IPDgofha < 0.04) & (H.IPDfmp >= 20)]
print(f'   gofha 0.02-0.04 & IPDfmp>=20: n={len(s)} frac |d\'|<30 {(np.abs(s.dp)<30).mean():.2f}  <45 {(np.abs(s.dp)<45).mean():.2f}')
out['phase'] = dict(bins=ph, target_like=dict(n=int(len(s)), f30=float((np.abs(s.dp) < 30).mean()), f45=float((np.abs(s.dp) < 45).mean())),
                    target_phase=18.24, implied_pa_mod180=(18.24 + 90) % 180)
json.dump(out, open(f'{R}/results/v07_gaia_likelihood.json', 'w'), indent=1)
