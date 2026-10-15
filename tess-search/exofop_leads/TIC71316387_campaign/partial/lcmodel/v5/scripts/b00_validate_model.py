"""Validation of lib/ebm.py:
 (1) occ() vs batman's Mandel-Agol quadratic routine for p < 1;
 (2) occ() vs brute-force 2-D integration for p > 1 and p < 1;
 (3) full eccentric light curve vs a brute-force 2-D two-disc sky integration (positions from an independent
     Kepler solver), including the secondary-eclipse timing;
 (4) exposure integration (cumulative-trapezoid) vs direct supersampling."""
import numpy as np, sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v5/lib')
import ebm
V = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v5'
res = {}
# (1) batman
import batman
from batman import _quadratic_ld
zs = np.linspace(0.0, 1.6, 400)
worst = 0
for p in (0.05, 0.3, 0.6, 0.9, 0.99):
    for u1, u2 in ((0.3, 0.25), (0.5, 0.1), (0.1, 0.5)):
        fb = _quadratic_ld._quadratic_ld(np.ascontiguousarray(zs), p, u1, u2, 1)
        fm = 1 - np.array([ebm.occ(z, p, u1, u2, ebm.XG, ebm.WG) for z in zs])
        worst = max(worst, np.max(np.abs(fb - fm)))
res['occ_vs_batman_maxabs'] = float(worst)
print('occ vs batman (p<1): max |diff| =', worst)
# (2) brute force
def brute(z, p, u1, u2, n=3000):
    x = (np.arange(n) + 0.5) / n * 2 - 1
    X, Y = np.meshgrid(x, x)
    r2 = X ** 2 + Y ** 2; ins = r2 < 1
    mu = np.sqrt(np.clip(1 - r2, 0, 1))
    I = np.where(ins, 1 - u1 * (1 - mu) - u2 * (1 - mu) ** 2, 0)
    cov = (X - z) ** 2 + Y ** 2 < p * p
    return (I * cov).sum() / I.sum()
worst2 = 0; rows = []
for p in (0.5, 1.3, 1.9, 2.6):
    for z in np.linspace(max(0, p - 1.05), p + 0.98, 7):
        b = brute(z, p, 0.35, 0.25); m = ebm.occ(z, p, 0.35, 0.25, ebm.XG, ebm.WG)
        worst2 = max(worst2, abs(b - m)); rows.append((p, z, b, m))
res['occ_vs_brute_maxabs'] = float(worst2)
print('occ vs brute-force 2-D (incl. p>1): max |diff| =', worst2, '(brute-force grid error ~1e-4)')
# (3) full light curve vs brute-force sky integration
P = 1.3036712
par = dict(rsum=0.22, k=0.62, b0=0.70, J=0.45, ecosw=0.004, esinw=-0.03, u1a=0.35, u2a=0.25, u1b=0.40, u2b=0.25)
def sky_positions(t, P, a1, e, w, inc):
    # independent: relative orbit of star 2 around star 1, conjunction (star 2 in front) at t=0
    nu0 = np.pi / 2 - w
    E0 = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(nu0 / 2)); M0 = E0 - e * np.sin(E0)
    M = 2 * np.pi * t / P + M0
    E = M.copy()
    for _ in range(100): E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    X = np.cos(E) - e; Y = np.sqrt(1 - e * e) * np.sin(E)
    # orbital plane -> sky: rotate by w, then incline; observer along +Z
    xo = X * np.cos(w) - Y * np.sin(w); yo = X * np.sin(w) + Y * np.cos(w)
    sx = xo; sy = yo * np.cos(inc); sz = yo * np.sin(inc)   # sz>0 : star 2 towards observer
    return sx, sy, sz
def brute_lc(t, p, n=900):
    e = np.hypot(p['ecosw'], p['esinw']); w = np.arctan2(p['esinw'], p['ecosw'])
    r1 = p['rsum'] / (1 + p['k']); r2 = p['rsum'] * p['k'] / (1 + p['k'])
    inc = np.arccos(p['b0'] * p['rsum'])
    sx, sy, sz = sky_positions(t, P, 1.0, e, w, inc)
    out = []
    g = (np.arange(n) + 0.5) / n * 2 - 1
    X, Y = np.meshgrid(g, g)
    def disc(R, u1, u2):
        r2_ = X ** 2 + Y ** 2; mu = np.sqrt(np.clip(1 - r2_, 0, 1))
        return np.where(r2_ < 1, 1 - u1 * (1 - mu) - u2 * (1 - mu) ** 2, 0)
    I1 = disc(1, p['u1a'], p['u2a']); I2 = disc(1, p['u1b'], p['u2b'])
    F1 = I1.sum() * r1 ** 2; F2 = I2.sum() * r2 ** 2
    s2 = p['J'] * F1 / F2 * p['k'] ** 2   # scale star-2 brightness so that L2/L1 = k^2 J
    Ltot = F1 + s2 * F2
    for i in range(t.size):
        if sz[i] > 0:   # star 2 in front of star 1: hidden part of star 1 (grid on star 1)
            cov = (X * r1 - sx[i]) ** 2 + (Y * r1 - sy[i]) ** 2 < r2 ** 2
            out.append(1 - (I1 * cov).sum() * r1 ** 2 / Ltot)
        else:
            cov = (X * r2 + sx[i]) ** 2 + (Y * r2 + sy[i]) ** 2 < r1 ** 2
            out.append(1 - s2 * (I2 * cov).sum() * r2 ** 2 / Ltot)
    return np.array(out)
dsec = ebm.secondary_offset(P, par['ecosw'], par['esinw'])
t = np.concatenate([np.linspace(-0.06, 0.06, 41), dsec + np.linspace(-0.06, 0.06, 41)])
fb = brute_lc(t, par)
fm = ebm.flux_grid(t, P, par['rsum'], par['k'], par['b0'], par['J'], par['ecosw'], par['esinw'], par['u1a'], par['u2a'], par['u1b'], par['u2b'], ebm.XG, ebm.WG)
res['lc_vs_brute_maxabs'] = float(np.max(np.abs(fb - fm))); res['lc_depths'] = [float(1 - fm[:41].min()), float(1 - fm[41:].min())]
# timing of brute-force secondary minimum (parabola through lowest points) vs Sterne offset
j = 41 + np.argmin(fb[41:]); c = np.polyfit(t[j - 3:j + 4], fb[j - 3:j + 4], 2); tmin_b = -c[1] / (2 * c[0])
res['secondary_offset_model_minus_halfP_min'] = float((dsec - P / 2) * 1440)
res['secondary_min_brute_minus_model_s'] = float((tmin_b - dsec) * 86400)
res['secondary_offset_approx_2Pecosw_over_pi_min'] = float(2 * P * par['ecosw'] / np.pi * (1 + 1 / np.sin(np.arccos(par['b0'] * par['rsum'])) ** 2) / 2 * 1440)
print('LC vs brute-force sky integration: max |diff| =', res['lc_vs_brute_maxabs'], 'depths', res['lc_depths'])
print('secondary offset (model) - P/2 =', res['secondary_offset_model_minus_halfP_min'], 'min; brute-force minimum - model =', res['secondary_min_brute_minus_model_s'], 's')
# (4) exposure integration
gm = ebm.GridModel(P)
tt = np.linspace(-0.2, 0.2, 300)
for texp in (1800 / 86400, 600 / 86400, 200 / 86400):
    a = gm.eval(tt, np.full(tt.size, texp), P, par['rsum'], par['k'], par['b0'], par['J'], par['ecosw'], par['esinw'], par['u1a'], par['u2a'], par['u1b'], par['u2b'])
    ns = 201; sub = (np.arange(ns) + 0.5) / ns - 0.5
    tsub = (tt[:, None] + texp * sub[None, :]).ravel()
    b = ebm.flux_grid(tsub, P, par['rsum'], par['k'], par['b0'], par['J'], par['ecosw'], par['esinw'], par['u1a'], par['u2a'], par['u1b'], par['u2b'], ebm.XG, ebm.WG).reshape(tt.size, ns).mean(1)
    res[f'smear_maxabs_texp{int(texp*86400)}s'] = float(np.max(np.abs(a - b)))
    print(f'exposure integration texp={texp*86400:.0f}s: max |grid - supersampled| = {np.max(np.abs(a-b)):.2e}')
json.dump(res, open(f'{V}/results/validate_model.json', 'w'), indent=1)
