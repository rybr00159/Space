"""Validation of the two-sphere EB model (lib/ebmodel.py):
 1. occ_frac vs batman (p<1), on z grids that avoid batman's singular points z=|1-p| and z=p (where batman's elliptic
    integrals can fail), with 16 Gauss-Legendre nodes (as used in the fits) and 64 nodes.
 2. occ_frac vs brute-force 2-D pixel integration for p>1 (occulter larger than occulted star) and for total occultations.
 3. total-occultation depth equals the light fraction L2 = k^2 J/(1+k^2 J) times (1-L3).
 4. eccentric orbit: secondary-eclipse phase offset ~ (2/pi) e cos w (1+1/sin^2 i) and duration ratio ~ (1+e sin w)/(1-e sin w)
    (Sterne / Winn 2010 approximations).
 5. supersampling error for 1800/600/200-s cadences with NSUB used in fits (9/3/2)."""
import numpy as np, sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4/lib')
from ebmodel import occ_frac, eb_flux
from batman import _quadratic_ld
OUT = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4/results/validate_model.json'
res = {}
for NG in (16, 64):
    xg, wg = np.polynomial.legendre.leggauss(NG)
    worst = 0.0
    for u1, u2 in [(0.45, 0.2), (0.3, 0.3), (0.6, 0.05)]:
        for p in [0.05, 0.1, 0.3, 0.5, 0.6, 0.75, 0.9]:
            zz = np.linspace(0, 1 + p, 503)
            zz = zz[(np.abs(zz - abs(1 - p)) > 1e-3) & (np.abs(zz - p) > 1e-3)]
            fb = _quadratic_ld._quadratic_ld(zz, p, u1, u2, 1)
            fm = np.array([1 - occ_frac(z, p, u1, u2, xg, wg) for z in zz])
            worst = max(worst, float(np.max(np.abs(fb - fm))))
    res[f'batman_maxdiff_NG{NG}'] = worst
    print(f'NG={NG}: max |batman - ours| over p in [0.05,0.9], 3 LD laws = {worst:.2e}')
def brute(z, p, u1, u2, N=4000):
    x = np.linspace(-1, 1, N); X, Y = np.meshgrid(x, x); R2 = X ** 2 + Y ** 2; inn = R2 < 1
    mu = np.sqrt(np.clip(1 - R2, 0, 1)); I = (1 - u1 * (1 - mu) - u2 * (1 - mu) ** 2) * inn
    cov = ((X - z) ** 2 + Y ** 2) < p * p
    return (I * cov).sum() / I.sum()
xg, wg = np.polynomial.legendre.leggauss(16)
bw = 0.0
for p, z in [(1.5, 1.2), (1.5, 2.0), (2.0, 1.5), (1.2, 0.5), (1.2, 2.1), (2.04, 1.3), (2.04, 2.9), (0.5, 0.2), (0.5, 0.9)]:
    a = occ_frac(z, p, 0.3, 0.3, xg, wg); b = brute(z, p, 0.3, 0.3)
    bw = max(bw, abs(a - b)); print(f'p={p} z={z}: ours {a:.6f} brute {b:.6f}')
res['brute_maxdiff'] = bw
# 3. total occultation depth
k, J, L3 = 0.49, 0.5, 0.97
t = np.linspace(-0.01, 0.01, 101) + 1.3036712 / 2
f = eb_flux(t, 0.0, 1.3036712, 0.165, k, 0.0, J, L3, 0.0, 0.0, 0.4, 0.2, 0.3, 0.3, 1e-6, 1, xg, wg)
L2 = k * k * J / (1 + k * k * J)
res['total_occ_depth'] = [float(1 - f.min()), float((1 - L3) * L2)]
print('total-occultation depth', 1 - f.min(), 'expected (1-L3) L2 =', (1 - L3) * L2)
# 4. eccentric timing and durations
P = 1.3036712
for ecw, esw in [(0.01, 0.0), (0.0, 0.05), (0.003, -0.03)]:
    tt = np.linspace(-0.2 * P, 0.8 * P, 400001)
    ff = eb_flux(tt, 0.0, P, 0.2, 0.5, 0.02, 0.5, 0.0, ecw, esw, 0.4, 0.2, 0.3, 0.3, 1e-6, 1, xg, wg)
    prim = (tt > -0.1) & (tt < 0.1) & (ff < 1 - 1e-9); sec = (tt > P / 2 - 0.15) & (tt < P / 2 + 0.15) & (ff < 1 - 1e-9)
    tp = 0.5 * (tt[prim].min() + tt[prim].max()); ts = 0.5 * (tt[sec].min() + tt[sec].max())
    dp = tt[prim].max() - tt[prim].min(); ds = tt[sec].max() - tt[sec].min()
    sini2 = 1 - 0.02 ** 2
    pred_off = 2 / np.pi * ecw * (1 + 1 / sini2) / 2 * P   # Sterne: dt = P/pi * ecosw * (1 + csc^2 i)
    print(f'ecosw={ecw} esinw={esw}: secondary - primary - P/2 = {(ts-tp-P/2)*1440:.3f} min (approx {pred_off*1440:.3f}); '
          f'Tsec/Tprim = {ds/dp:.4f} (approx (1+esinw)/(1-esinw) = {(1+esw)/(1-esw):.4f})')
    res[f'ecc_{ecw}_{esw}'] = dict(dt_min=float((ts - tp - P / 2) * 1440), approx=float(pred_off * 1440), dur_ratio=float(ds / dp), approx_ratio=float((1 + esw) / (1 - esw)))
# 5. supersampling
tt = np.linspace(-0.06, 0.06, 2001)
for texp, ns in [(1800 / 86400, 9), (600 / 86400, 3), (200 / 86400, 2)]:
    worst = 0
    for par in [(0.165, 0.49, 0.01, 0.5, 0.975), (0.255, 0.79, 0.195, 0.42, 0.91), (0.25, 0.1, 0.2, 0.3, 0.0)]:
        rsum, k, cosi, J, L3 = par
        ref = eb_flux(tt, 0., P, rsum, k, cosi, J, L3, 0., 0., 0.3, 0.3, 0.3, 0.3, texp, 61, xg, wg)
        m = eb_flux(tt, 0., P, rsum, k, cosi, J, L3, 0., 0., 0.3, 0.3, 0.3, 0.3, texp, ns, xg, wg)
        worst = max(worst, np.abs(m - ref).max())
    res[f'nsub_err_ppm_{int(texp*86400)}s'] = float(worst * 1e6)
    print(f'cadence {texp*86400:.0f} s with NSUB={ns}: max supersampling error {worst*1e6:.1f} ppm')
json.dump(res, open(OUT, 'w'), indent=1)
