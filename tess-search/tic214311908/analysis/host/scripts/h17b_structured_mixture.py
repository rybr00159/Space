"""h17b: structured (identifiable) population mixture for the twin-vs-metal-rich question.
The free 2D Gaussian mixture of h17 absorbs everything into 'singles' (not identifiable). Physical structure used here:
for single stars metallicity moves a star along a line in the plane (x = CMD offset, y = colour residual at fixed
BP-RP): y = c + r x + e, e ~ N(0, s_y) (s_y = intrinsic + photometric scatter); the x-distribution of singles is a free
3-component 1D Gaussian mixture. A near-equal companion moves a star by (dx, dy) = (-0.753, ~0): twins therefore sit
r*0.753 below the singles' line at the same x. Unequal binaries (q 0.6-0.93; offsets from the Mamajek MS at fixed
combined colour, weight ~ q^0.5) form a third group; a broad outlier component absorbs the rest.
Free: 1D GMM (8), c, r, s_y, w_twin, w_inter, w_out. ML with L-BFGS-B + Nelder-Mead; bootstrap (12 resamples, L-BFGS-B from the ML solution) for errors.
The M3V's P(twin | x, y) uses its errors (Ks error enters x and y with opposite signs for the 2MASS plane).
Planes: (dMKs, r(G-Ks)) and the 2MASS-free (dMG, r(G-W1)); samples as in h17 (RUWE < 1.1, clean IPD) plus a
v_tan 30-70 km/s subsample. Also the empirical singles' slope from the singles-dominated range |x| < 0.4.
Output: out/h17b_structured_mixture.json, figs/h17b_structured_mixture.png"""
import json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from scipy.optimize import minimize
from scipy.stats import norm
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
B = 'tic214311908/analysis/host'
h13 = json.load(open(f'{B}/out/h13_colour_feh.json')); R = {k: np.array(v) for k, v in h13['ridges'].items()}
h17 = json.load(open(f'{B}/out/h17_population_mixture.json'))
F = pd.read_csv(f'{B}/data/field_mdwarfs.csv')
x = F.bp_rp; fx = 1.162004 + 0.011464 * x + 0.049255 * x ** 2 - 0.005879 * x ** 3
F['Cs'] = (F.phot_bp_rp_excess_factor - fx) / (0.0059898 + 8.817481e-12 * F.phot_g_mean_mag ** 7.618399)
F['dMK'] = F.ks_m + 5 * np.log10(F.parallax / 100) - np.polyval(R['MK'], F.bp_rp)
F['dMG'] = F.phot_g_mean_mag + 5 * np.log10(F.parallax / 100) - np.polyval(R['MG'], F.bp_rp)
F['rGK'] = F.phot_g_mean_mag - F.ks_m - np.polyval(R['GK'], F.bp_rp)
F['rGW1'] = F.phot_g_mean_mag - F.w1mpro - np.polyval(np.array(h17['GW1_ridge']), F.bp_rp)
F['vtan'] = 4.74047 * np.hypot(F.pmra, F.pmdec) / F.parallax
base = (np.abs(F.bp_rp - 2.756) < 0.08) & (F.ruwe < 1.1) & (F.ipd_frac_multi_peak == 0) & (F.ipd_gof_harmonic_amplitude < 0.03) & \
       (F.parallax / F.parallax_error > 20) & (F.ph_qual == 'AAA') & (F.ks_msigcom < 0.05) & (F.Cs.abs() < 3)
wok = F.wise_cc.astype(str).str[:2].eq('00') & F.wise_ph_qual.astype(str).str[:2].eq('AA') & (F.wise_ext == 0)
T = h17['target']
# binary offsets (q) - recomputed as in h17 (Mamajek MS); twin group q 0.95-1 ~ (-0.753, 0)
S_ = 'tic214311908/analysis/star'
tab = []
for line in open(f'{S_}/data/ext/EEM_dwarf.txt'):
    p = line.split()
    if len(p) > 30 and p[0][0] == 'M' and p[0].endswith('V'):
        try: tab.append((float(p[30]), float(p[11]), float(p[13]), float(p[22]), float(p[23])))
        except Exception: pass
tab = np.array(sorted(tab)); ms = np.linspace(0.09, 0.6, 3000)
cI = np.interp(ms, tab[:, 0], tab[:, 1]); gI = np.interp(ms, tab[:, 0], tab[:, 2]); kI = np.interp(ms, tab[:, 0], tab[:, 3]); w1I = kI - np.interp(ms, tab[:, 0], tab[:, 4])
add = lambda a, b: -2.5 * np.log10(10 ** (-0.4 * a) + 10 ** (-0.4 * b))
rK = lambda col: np.interp(col, cI[::-1], kI[::-1]); rG = lambda col: np.interp(col, cI[::-1], gI[::-1]); rW = lambda col: np.interp(col, cI[::-1], w1I[::-1])
def offs(q):
    best = None
    for m1 in np.linspace(0.2, 0.45, 801):
        i1 = np.argmin(np.abs(ms - m1)); i2 = np.argmin(np.abs(ms - q * m1))
        BP1, RP1 = gI[i1] + 0.562 * cI[i1], gI[i1] - 0.438 * cI[i1]; BP2, RP2 = gI[i2] + 0.562 * cI[i2], gI[i2] - 0.438 * cI[i2]
        cc = add(BP1, BP2) - add(RP1, RP2)
        if best is None or abs(cc - 2.756) < best[0]:
            Gc, Kc, Wc = add(gI[i1], gI[i2]), add(kI[i1], kI[i2]), add(w1I[i1], w1I[i2])
            best = (abs(cc - 2.756), Kc - rK(cc), Gc - rG(cc), (Gc - Kc) - (rG(cc) - rK(cc)), (Gc - Wc) - (rG(cc) - rW(cc)))
    return best[1:]
qT = np.array([0.95, 0.975, 1.0]); qI = np.linspace(0.6, 0.925, 8); wI = qI ** 0.5; wI /= wI.sum()
OT = np.array([offs(q) for q in qT]); OI = np.array([offs(q) for q in qI])
print('offsets (dK, dG, dGK, dGW1): q=1', OT[-1].round(3), 'q=0.8', np.array(offs(0.8)).round(3), 'q=0.6', OI[0].round(3))
def fit(X, Y, ex, ey, cxy, plane, label, nboot=12, r_fixed=None):
    jx, jy = (0, 2) if plane == 'K' else (1, 3)
    shT = OT[:, [jx, jy]]; shI = OI[:, [jx, jy]]
    Xm, Xs = X.mean(), X.std(); Ym, Ys = Y.mean(), Y.std()
    def unpack(p):
        w = np.exp(p[0:3]); w = w / w.sum(); mu = p[3:6]; sg = np.exp(p[6:9])
        c, r, sy = p[9], (p[10] if r_fixed is None else r_fixed), np.exp(p[11])
        e = np.exp(p[12:15]); tot = 1 + e.sum(); wT, wIn, wO = e / tot
        return w, mu, sg, c, r, sy, wT, wIn, wO
    def comps(Xp, Yp, p, ex2=0.0, ey2=0.0, cv=0.0):
        w, mu, sg, c, r, sy, wT, wIn, wO = unpack(p)
        def single(xx, yy):
            # p(x) * p(y | x) with measurement errors: joint Gaussian per GMM component
            out = 0
            for k in range(3):
                vx = sg[k] ** 2 + ex2; mx = mu[k]; my = c + r * mu[k]
                vy = r ** 2 * sg[k] ** 2 + sy ** 2 + ey2; cxy_ = r * sg[k] ** 2 + cv
                det = vx * vy - cxy_ ** 2; dx = xx - mx; dy = yy - my
                out = out + w[k] * np.exp(-0.5 * (vy * dx * dx - 2 * cxy_ * dx * dy + vx * dy * dy) / det) / (2 * np.pi * np.sqrt(det))
            return out
        s = single(Xp, Yp)
        t = np.mean([single(Xp - d[0], Yp - d[1]) for d in shT], axis=0)
        i = np.sum([ww * single(Xp - d[0], Yp - d[1]) for ww, d in zip(wI, shI)], axis=0)
        o = norm.pdf(Xp, Xm, 3 * Xs) * norm.pdf(Yp, Ym, 3 * Ys)
        return (1 - wT - wIn - wO) * s, wT * t, wIn * i, wO * o
    def nll(p, XX, YY):
        a, b, c, d = comps(XX, YY, p); return -np.sum(np.log(a + b + c + d + 1e-300))
    # singles' slope from the singles-dominated range
    k0 = np.abs(X) < 0.4; r0 = np.polyfit(X[k0], Y[k0], 1)[0]
    p0 = np.r_[0, 0, 0, np.percentile(X, [25, 50, 75]), np.log([0.2, 0.2, 0.3]), np.median(Y), r0, np.log(0.03), np.log([0.02, 0.03, 0.01])]
    def solve(XX, YY, pstart):
        r1 = minimize(nll, pstart, args=(XX, YY), method='L-BFGS-B')
        r1 = minimize(nll, r1.x, args=(XX, YY), method='Nelder-Mead', options=dict(maxiter=8000, maxfev=8000, xatol=1e-5, fatol=1e-6))
        return r1
    best = None
    for wt0 in (-4.0, -2.5):
        pp = p0.copy(); pp[12] = wt0
        r_ = solve(X, Y, pp)
        if best is None or r_.fun < best.fun: best = r_
    p = best.x; w, mu, sg, c, r, sy, wT, wIn, wO = unpack(p)
    def target_post(pp, tx, ty):
        a, b, c_, d = comps(np.array([tx]), np.array([ty]), pp, ex ** 2, ey ** 2, cxy)
        tot = a + b + c_ + d
        return float(a[0] / tot[0]), float(b[0] / tot[0]), float(c_[0] / tot[0]), float(b[0] / a[0])
    Ps, Pt, Pi, O = target_post(p, Tx, Ty)
    bs = []
    rng = np.random.default_rng(4)
    for _ in range(nboot):
        ii = rng.integers(0, len(X), len(X))
        rb = minimize(nll, p, args=(X[ii], Y[ii]), method='L-BFGS-B')
        u = unpack(rb.x); bs.append((u[6], u[4], target_post(rb.x, Tx, Ty)[1], target_post(rb.x, Tx, Ty)[3]))
    bs = np.array(bs)
    near = np.abs(X - Tx) < 0.08
    res = dict(n=int(len(X)), slope_singles_r=float(r), slope_from_absx_lt_0p4=float(r0), s_y=float(sy), w_twin=float(wT), w_twin_16_84=[float(np.percentile(bs[:, 0], 16)), float(np.percentile(bs[:, 0], 84))],
               w_inter=float(wIn), w_out=float(wO), target_x=float(Tx), target_y=float(Ty),
               P_single=Ps, P_twin=Pt, P_inter=Pi, odds_twin_single=O, P_twin_16_84=[float(np.percentile(bs[:, 2], 16)), float(np.percentile(bs[:, 2], 84))],
               odds_16_84=[float(np.percentile(bs[:, 3], 16)), float(np.percentile(bs[:, 3], 84))],
               singles_line_y_at_target_x=float(c + r * Tx), twins_line_y_at_target_x=float(c + r * (Tx + 0.753)),
               near_target=dict(n=int(near.sum()), median_y=float(np.median(Y[near])), frac_y_ge_target=float(np.mean(Y[near] >= Ty))))
    print(f'[{label}] n={len(X)} r={r:+.3f} (|x|<0.4: {r0:+.3f}) s_y={sy:.3f} w_twin={wT:.4f} [{res["w_twin_16_84"][0]:.4f},{res["w_twin_16_84"][1]:.4f}] w_inter={wIn:.4f} | '
          f'target y={Ty:.3f}: singles line {res["singles_line_y_at_target_x"]:.3f}, twins line {res["twins_line_y_at_target_x"]:.3f} | P(twin)={Pt:.3f} '
          f'[{res["P_twin_16_84"][0]:.3f},{res["P_twin_16_84"][1]:.3f}] P(inter)={Pi:.3f} odds T:S={O:.3f}', flush=True)
    return res, (X, Y, p, comps)
out = {}
eK = 0.026
runs = [('K', base, 'dMK', 'rGK', T['dMK'], T['rGK'], np.hypot(eK, 0.036), np.hypot(eK, 0.01), -eK ** 2, 'dMK-rGK'),
        ('K', base & F.vtan.between(30, 70), 'dMK', 'rGK', T['dMK'], T['rGK'], np.hypot(eK, 0.036), np.hypot(eK, 0.01), -eK ** 2, 'dMK-rGK vtan30-70'),
        ('K', base, 'dMK', 'rGK', T['dMK'] + 0.04, T['rGK'] - 0.04, np.hypot(eK, 0.036), np.hypot(eK, 0.01), -eK ** 2, 'dMK-rGK Ks+0.04'),
        ('W', base & wok, 'dMG', 'rGW1', T['dMG'], T['rGW1'], 0.036, np.hypot(0.03, 0.01), 0.0, 'dMG-rGW1'),
        ('W', base & wok & F.vtan.between(30, 70), 'dMG', 'rGW1', T['dMG'], T['rGW1'], 0.036, np.hypot(0.03, 0.01), 0.0, 'dMG-rGW1 vtan30-70'),
        ('W', base & wok, 'dMG', 'rGW1', T['dMG'], T['rGW1'] - 0.03, 0.036, np.hypot(0.03, 0.01), 0.0, 'dMG-rGW1 W1+0.03')]
keep = {}
for plane, sel, xc, yc, Tx, Ty, ex, ey, cxy, lab in runs:
    d = F[sel & np.isfinite(F[xc]) & np.isfinite(F[yc])]
    res, kk = fit(d[xc].values, d[yc].values, ex, ey, cxy, plane, lab)
    out[lab] = res; keep[lab] = kk
json.dump(out, open(f'{B}/out/h17b_structured_mixture.json', 'w'), indent=1)
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for a_, lab, (tx, ty) in ((ax[0], 'dMK-rGK', (T['dMK'], T['rGK'])), (ax[1], 'dMG-rGW1', (T['dMG'], T['rGW1']))):
    X, Y, p, comps = keep[lab]; s, t, i, o = comps(X, Y, p); pt = t / (s + t + i + o)
    sc = a_.scatter(X, Y, c=pt, s=4, cmap='coolwarm', vmin=0, vmax=1)
    r_ = out[lab]['slope_singles_r']; c_ = out[lab]['singles_line_y_at_target_x'] - r_ * tx
    xx = np.linspace(-1.5, 1.0, 10); a_.plot(xx, c_ + r_ * xx, 'k-', lw=1, label='singles line'); a_.plot(xx, c_ + r_ * (xx + 0.753), 'k--', lw=1, label='twins line')
    a_.scatter([tx], [ty], marker='*', s=250, c='k', label='M3V'); a_.set_xlim(1.0, -1.5); a_.set_ylim(-0.2, 0.3); a_.legend(fontsize=8)
    a_.set_xlabel(lab.split('-')[0]); a_.set_ylabel(lab.split('-')[1]); plt.colorbar(sc, ax=a_, label='P(twin)')
plt.tight_layout(); plt.savefig(f'{B}/figs/h17b_structured_mixture.png', dpi=110)
