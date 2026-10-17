"""h17: empirical twin-vs-metal-rich test from the field population (no literature priors, no [Fe/H] labels).
Field M dwarfs (h03) at the M3V's colour (|BP-RP - 2.756| < 0.08), selected like the M3V: RUWE < 1.1,
ipd_frac_multi_peak = 0, ipd_gof_harmonic_amplitude < 0.03, plx/err > 20, 2MASS AAA, ks_msigcom < 0.05, |C*| < 3.
Plane: x = dMKs (offset from the CMD ridge, h13) and y = r(G-Ks) (G-Ks residual at fixed BP-RP, h13) - metallicity
moves stars along a correlated direction (metal-rich: brighter AND redder in G-Ks), an equal-mass companion moves them
by (-0.753, 0). Alternative 2MASS-free plane: x = dMG, y = r(G-W1) (AllWISE W1, cc_flags 0000, ph_qual AA).
Model: singles = 2-component 2D Gaussian mixture; binaries = the singles' distribution shifted by the binary offsets
(dx(q), dy(q)) from the Mamajek MS at fixed combined BP-RP, for q in a 'twin' group (0.95-1) and an 'intermediate'
group (0.6-0.95, weight ~ q^0.5); free weights w_twin, w_int; plus a broad outlier component. Maximum likelihood (EM-like
direct optimisation). Per-star twin probabilities; the M3V's P(twin | x, y) with its measurement errors (Ks error enters
x and y with opposite signs). Also: posterior predictive check (fraction of stars brighter than the M3V), a v_tan
30-70 km/s kinematic subsample, and Ks shifted by +0.04 mag (possible Mira diffraction-spike contamination).
Output: out/h17_population_mixture.json, figs/h17_population_mixture.png"""
import json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from scipy.optimize import minimize
from scipy.stats import multivariate_normal as mvn
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
B = 'tic214311908/analysis/host'
S_ = 'tic214311908/analysis/star'
h13 = json.load(open(f'{B}/out/h13_colour_feh.json')); R = {k: np.array(v) for k, v in h13['ridges'].items()}
F = pd.read_csv(f'{B}/data/field_mdwarfs.csv')
x = F.bp_rp; fx = 1.162004 + 0.011464 * x + 0.049255 * x ** 2 - 0.005879 * x ** 3
F['Cs'] = (F.phot_bp_rp_excess_factor - fx) / (0.0059898 + 8.817481e-12 * F.phot_g_mean_mag ** 7.618399)
F['MK'] = F.ks_m + 5 * np.log10(F.parallax / 100); F['MG'] = F.phot_g_mean_mag + 5 * np.log10(F.parallax / 100)
F['dMK'] = F.MK - np.polyval(R['MK'], F.bp_rp); F['dMG'] = F.MG - np.polyval(R['MG'], F.bp_rp)
F['rGK'] = F.phot_g_mean_mag - F.ks_m - np.polyval(R['GK'], F.bp_rp)
GW1r = None
F['GW1'] = F.phot_g_mean_mag - F.w1mpro
F['vtan'] = 4.74047 * np.hypot(F.pmra, F.pmdec) / F.parallax
base = (np.abs(F.bp_rp - 2.756) < 0.08) & (F.ruwe < 1.1) & (F.ipd_frac_multi_peak == 0) & (F.ipd_gof_harmonic_amplitude < 0.03) & \
       (F.parallax / F.parallax_error > 20) & (F.ph_qual == 'AAA') & (F.ks_msigcom < 0.05) & (F.Cs.abs() < 3)
wok = F.wise_cc.astype(str).str[:2].eq('00') & F.wise_ph_qual.astype(str).str[:2].eq('AA') & (F.wise_ext == 0)
# G-W1 ridge (median, from the wider sample)
wide = (F.ruwe < 1.2) & wok & F.bp_rp.between(2.4, 3.15) & (F.parallax / F.parallax_error > 20)
pts = [(c0, np.median(F.GW1[wide & (np.abs(F.bp_rp - c0) < 0.025)])) for c0 in np.arange(2.42, 3.13, 0.02)]
pts = np.array(pts); cGW1 = np.polyfit(pts[:, 0], pts[:, 1], 3)
F['rGW1'] = F.GW1 - np.polyval(cGW1, F.bp_rp)
# target
plx = 15.549
T = dict(bp_rp=2.756149, G=14.281043, K=10.482, W1=10.332 - 0.014)
T['MK'] = T['K'] + 5 * np.log10(plx / 100); T['MG'] = T['G'] + 5 * np.log10(plx / 100)
T['dMK'] = T['MK'] - np.polyval(R['MK'], T['bp_rp']); T['dMG'] = T['MG'] - np.polyval(R['MG'], T['bp_rp'])
T['rGK'] = T['G'] - T['K'] - np.polyval(R['GK'], T['bp_rp']); T['rGW1'] = T['G'] - T['W1'] - np.polyval(cGW1, T['bp_rp'])
print('target', {k: round(v, 3) for k, v in T.items()})
# binary offsets from the Mamajek MS (mass, BP-RP, M_G, M_Ks, Ks-W1)
tab = []
for line in open(f'{S_}/data/ext/EEM_dwarf.txt'):
    p = line.split()
    if len(p) > 30 and p[0][0] == 'M' and p[0].endswith('V'):
        try: tab.append((float(p[30]), float(p[11]), float(p[13]), float(p[22]), float(p[23])))
        except Exception: pass
tab = np.array(sorted(tab)); ms = np.linspace(0.09, 0.6, 3000)
cI = np.interp(ms, tab[:, 0], tab[:, 1]); gI = np.interp(ms, tab[:, 0], tab[:, 2]); kI = np.interp(ms, tab[:, 0], tab[:, 3]); w1I = kI - np.interp(ms, tab[:, 0], tab[:, 4])
def offs(q):
    best = None
    for m1 in np.linspace(0.2, 0.45, 801):
        i1 = np.argmin(np.abs(ms - m1)); i2 = np.argmin(np.abs(ms - q * m1))
        G1, G2 = gI[i1], gI[i2]; K1, K2 = kI[i1], kI[i2]; c1, c2 = cI[i1], cI[i2]; W1_, W2_ = w1I[i1], w1I[i2]
        BP1, RP1, BP2, RP2 = G1 + 0.562 * c1, G1 - 0.438 * c1, G2 + 0.562 * c2, G2 - 0.438 * c2
        add = lambda a, b: -2.5 * np.log10(10 ** (-0.4 * a) + 10 ** (-0.4 * b))
        cc = add(BP1, BP2) - add(RP1, RP2)
        if best is None or abs(cc - c1 - 0) < 9 and abs(cc - 2.756) < best[0]:
            Gc, Kc, Wc = add(G1, G2), add(K1, K2), add(W1_, W2_)
            rK = lambda col: np.interp(col, cI[::-1], kI[::-1]); rG = lambda col: np.interp(col, cI[::-1], gI[::-1]); rW = lambda col: np.interp(col, cI[::-1], w1I[::-1])
            best = (abs(cc - 2.756), Kc - rK(cc), Gc - rG(cc), (Gc - Kc) - (rG(cc) - rK(cc)), (Gc - Wc) - (rG(cc) - rW(cc)))
    return dict(dK=best[1], dG=best[2], dGK=best[3], dGW1=best[4])
qs_t = np.linspace(0.95, 1.0, 3); qs_i = np.linspace(0.6, 0.925, 6)
OT = [offs(q) for q in qs_t]; OI = [offs(q) for q in qs_i]
wi = qs_i ** 0.5; wi /= wi.sum()
print('twin offsets q=1:', OT[-1], ' q=0.95:', OT[0], ' q=0.8:', offs(0.8), ' q=0.6:', OI[0])
def g2(X, m, S):
    a, b, c = S[0, 0], S[0, 1], S[1, 1]; det = a * c - b * b
    dx = X[:, 0] - m[0]; dy = X[:, 1] - m[1]
    return np.exp(-0.5 * (c * dx * dx - 2 * b * dx * dy + a * dy * dy) / det) / (2 * np.pi * np.sqrt(det))
mvn_pdf = lambda X, m, S: g2(np.atleast_2d(X), np.asarray(m), np.asarray(S))
def run(xcol, ycol, sel, Tx, Ty, Terr, label):
    d = F[sel & np.isfinite(F[xcol]) & np.isfinite(F[ycol])]
    X = np.column_stack([d[xcol].values, d[ycol].values])
    kx, ky = ('dK', 'dGK') if xcol == 'dMK' else ('dG', 'dGW1')
    shT = np.array([[o[kx], o[ky]] for o in OT]); shI = np.array([[o[kx], o[ky]] for o in OI])
    def unpack(p):
        m1 = p[0:2]; m2 = p[2:4]
        L1 = np.array([[np.exp(p[4]), 0], [p[5], np.exp(p[6])]]); L2 = np.array([[np.exp(p[7]), 0], [p[8], np.exp(p[9])]])
        a = 1 / (1 + np.exp(-p[10])); wT = np.exp(p[11]); wI = np.exp(p[12]); wO = np.exp(p[13]); tot = 1 + wT + wI + wO
        return m1, m2, L1 @ L1.T, L2 @ L2.T, a, wT / tot, wI / tot, wO / tot
    Xm = X.mean(0); Xc = np.cov(X.T) * 9
    def comp_pdf(Xp, p, extra=np.zeros((2, 2))):
        m1, m2, S1, S2, a, wT, wI, wO = unpack(p)
        single = a * mvn_pdf(Xp, m1, S1 + extra) + (1 - a) * mvn_pdf(Xp, m2, S2 + extra)
        twin = np.mean([a * mvn_pdf(Xp, m1 + s, S1 + extra) + (1 - a) * mvn_pdf(Xp, m2 + s, S2 + extra) for s in shT], axis=0)
        inter = np.sum([w * (a * mvn_pdf(Xp, m1 + s, S1 + extra) + (1 - a) * mvn_pdf(Xp, m2 + s, S2 + extra)) for w, s in zip(wi, shI)], axis=0)
        outl = mvn_pdf(Xp, Xm, Xc + extra)
        ws = 1 - wT - wI - wO
        return ws * single, wT * twin, wI * inter, wO * outl
    def nll(p):
        a, b, c, e = comp_pdf(X, p); return -np.sum(np.log(a + b + c + e + 1e-300))
    sd = X.std(0)
    p0 = np.r_[np.median(X, 0), np.median(X, 0) + [0.3, -0.01], np.log(sd[0] * 0.6), 0.0, np.log(sd[1] * 0.6), np.log(sd[0] * 0.8), 0.0, np.log(sd[1] * 0.8), 0.0, np.log(0.03), np.log(0.05), np.log(0.01)]
    best = None
    for trial in range(2):
        pp = p0 + (np.random.default_rng(trial).normal(0, 0.05, len(p0)) if trial else 0)
        r = minimize(nll, pp, method='L-BFGS-B')
        r = minimize(nll, r.x, method='Nelder-Mead', options=dict(maxiter=6000, maxfev=6000, xatol=1e-5, fatol=1e-5))
        if best is None or r.fun < best.fun: best = r
    p = best.x; m1, m2, S1, S2, a, wT, wI, wO = unpack(p)
    # target posterior
    E = np.array([[Terr[0] ** 2, Terr[2]], [Terr[2], Terr[1] ** 2]])
    s_, t_, i_, o_ = comp_pdf(np.array([[Tx, Ty]]), p, E)
    tot = s_ + t_ + i_ + o_
    res = dict(n=int(len(X)), w_twin=float(wT), w_inter=float(wI), w_out=float(wO), single_means=[m1.tolist(), m2.tolist()], single_frac1=float(a),
               target=[float(Tx), float(Ty)], P_single=float(s_[0] / tot[0]), P_twin=float(t_[0] / tot[0]), P_inter=float(i_[0] / tot[0]), P_outlier=float(o_[0] / tot[0]),
               odds_twin_vs_single=float(t_[0] / s_[0]), nll=float(best.fun))
    # population check: among stars with x within 0.08 of target x, fraction assigned twin
    sp, tp, ip, op = comp_pdf(X, p); tt = sp + tp + ip + op
    nb = np.abs(X[:, 0] - Tx) < 0.08
    res['near_target_x'] = dict(n=int(nb.sum()), mean_P_twin=float(np.mean(tp[nb] / tt[nb])), mean_P_single=float(np.mean(sp[nb] / tt[nb])),
                                median_y=float(np.median(X[nb, 1])), frac_y_ge_target=float(np.mean(X[nb, 1] >= Ty)))
    res['frac_brighter_than_target'] = float(np.mean(X[:, 0] <= Tx))
    # model-predicted fraction brighter than target (posterior predictive): sample
    rng = np.random.default_rng(3); N = 200000
    comp = rng.choice(4, N, p=[1 - wT - wI - wO, wT, wI, wO])
    ss = np.where(rng.random(N) < a, 0, 1)
    Z = np.where(ss[:, None] == 0, rng.multivariate_normal(m1, S1, N), rng.multivariate_normal(m2, S2, N))
    Z[comp == 1] += shT[rng.integers(0, len(shT), (comp == 1).sum())]
    Z[comp == 2] += shI[rng.choice(len(shI), (comp == 2).sum(), p=wi)]
    Z[comp == 3] = rng.multivariate_normal(Xm, Xc, (comp == 3).sum())
    res['model_frac_brighter_than_target'] = float(np.mean(Z[:, 0] <= Tx))
    res['model_frac_brighter_by_component'] = {n_: float(np.mean((Z[:, 0] <= Tx) & (comp == i))) for i, n_ in enumerate(['single', 'twin', 'inter', 'outl'])}
    print(f'[{label}] n={len(X)} w_twin={wT:.4f} w_inter={wI:.4f} w_out={wO:.4f} | target P(single)={res["P_single"]:.3f} P(twin)={res["P_twin"]:.3f} '
          f'P(inter)={res["P_inter"]:.3f} odds T:S={res["odds_twin_vs_single"]:.3f} | frac brighter obs {res["frac_brighter_than_target"]:.4f} model {res["model_frac_brighter_than_target"]:.4f} {res["model_frac_brighter_by_component"]}')
    print('    near target x:', res['near_target_x'])
    return res, (X, sp / tt, tp / tt, p)
out = dict(target=T)
eK = 0.026; eG = 0.003
runs = {}
runs['dMK_rGK_ruwe1.1'], keep = run('dMK', 'rGK', base, T['dMK'], T['rGK'], (np.hypot(eK, 0.036), np.hypot(eK, 0.01), -eK ** 2), 'dMK-rGK')
runs['dMK_rGK_ruwe1.1_Ks+0.04'], _ = run('dMK', 'rGK', base, T['dMK'] + 0.04, T['rGK'] - 0.04, (np.hypot(eK, 0.036), np.hypot(eK, 0.01), -eK ** 2), 'dMK-rGK, Ks 0.04 fainter')
runs['dMK_rGK_vtan30-70'], _ = run('dMK', 'rGK', base & F.vtan.between(30, 70), T['dMK'], T['rGK'], (np.hypot(eK, 0.036), np.hypot(eK, 0.01), -eK ** 2), 'dMK-rGK, vtan 30-70')
eW = 0.03
runs['dMG_rGW1_ruwe1.1'], keepW = run('dMG', 'rGW1', base & wok, T['dMG'], T['rGW1'], (0.036, np.hypot(eW, 0.01), 0.0), 'dMG-rGW1 (no 2MASS)')
runs['dMG_rGW1_vtan30-70'], _ = run('dMG', 'rGW1', base & wok & F.vtan.between(30, 70), T['dMG'], T['rGW1'], (0.036, np.hypot(eW, 0.01), 0.0), 'dMG-rGW1, vtan 30-70')
out['runs'] = runs; out['GW1_ridge'] = cGW1.tolist(); out['binary_offsets_q1'] = OT[-1]
json.dump(out, open(f'{B}/out/h17_population_mixture.json', 'w'), indent=1, default=float)
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for a_, (X, ps, pt, p), lab, t in ((ax[0], keep, 'dMKs', ('dMK', 'rGK')), (ax[1], keepW, 'dMG', ('dMG', 'rGW1'))):
    sc = a_.scatter(X[:, 0], X[:, 1], c=pt, s=4, cmap='coolwarm', vmin=0, vmax=1)
    a_.scatter([T[t[0]]], [T[t[1]]], marker='*', s=250, c='k'); a_.set_xlabel(t[0]); a_.set_ylabel(t[1]); a_.set_xlim(-1.6, 1.2); a_.set_ylim(-0.25, 0.3)
    a_.invert_xaxis(); plt.colorbar(sc, ax=a_, label='P(twin)')
plt.tight_layout(); plt.savefig(f'{B}/figs/h17_population_mixture.png', dpi=110)
