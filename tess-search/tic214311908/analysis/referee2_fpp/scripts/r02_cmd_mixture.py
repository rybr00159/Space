"""r02: empirical calibration of the over-luminosity likelihood (the main driver of P(unresolved companion)).
The analyst (fpp/f03) assumed single-star offsets eps ~ N(0, sigma_s = 0.30 mag) and argued the bright side of the
comparison sample sets sigma_s. But the bright side of any field sample contains unresolved binaries, so its width is
not the single-star width. Here I fit a mixture model to the M_Ks offsets of Gaia-clean comparison M dwarfs at the
target colour:
    p(x) = (1 - f) S(x) + f * [S convolved with the binary brightening delta(q)]
    S(x) = (1 - t) SplitNormal(mu, sig_b (bright side), sig_f (faint side)) + t N(mu_t, sig_t)   (metal-poor/bad tail)
    delta(q) = dM_bin(q) - slope * dColour(q)  (companion light brightens M_Ks and reddens BP-RP; offsets are measured
               from the ridge at the observed colour; primary = 0.27 Msun, the ridge star at BP-RP 2.756)
    q ~ (1 - F_twin) U(0.075/m1, 1) + F_twin U(0.95, 1)
Fitted by maximum likelihood (truncated to -1.8 < x < 1.0), with bootstrap errors. Outputs the single-star bright-side
width and the population posterior P(binary | x = target offset) for Gaia-clean stars.
Input: referee_star/data/mdw_cmd.csv (comparison sample built in phase 1). Output: out/r02_cmd_mixture.json"""
import numpy as np, pandas as pd, sys
from scipy import optimize, stats
sys.path.insert(0, 'tic214311908/analysis/referee2_fpp/scripts')
from rlib import *
rng = np.random.default_rng(2)
d = pd.read_csv(f'{B}/referee_star/data/mdw_cmd.csv')
d = d[np.isfinite(d.ks_m) & (d.ph_qual.astype(str).str[2] == 'A') & (d.ks_msigcom < 0.05) & np.isfinite(d.ruwe)].copy()
d['MKs'] = d.ks_m + 5*np.log10(d.parallax/100.)
# ---- ridge: KDE mode of MKs in colour bins (RUWE < 1.1), straight-line fit
cb = np.arange(2.55, 2.951, 0.05); modes = []; cens = []
for lo in cb[:-1]:
    s = d[(d.bp_rp >= lo) & (d.bp_rp < lo + 0.05) & (d.ruwe < 1.1)].MKs.values
    xs = np.linspace(np.percentile(s, 1), np.percentile(s, 99), 3000)
    modes.append(xs[np.argmax(stats.gaussian_kde(s, 0.12)(xs))]); cens.append(np.median(d[(d.bp_rp >= lo) & (d.bp_rp < lo + 0.05)].bp_rp))
slope, icpt = np.polyfit(cens, modes, 1)
ridge = lambda c: icpt + slope*c
d['off'] = d.MKs - ridge(d.bp_rp)
x_t = MKS_OBS - ridge(2.756149)
print(f'ridge slope {slope:.3f} mag/mag; ridge at 2.756: {ridge(2.756149):.3f}; target offset {x_t:.3f}')
# ---- binary brightening table for a 0.27 Msun ridge primary (q grid)
m1 = float(mass_of_mks(ridge(2.756149)))
qg = np.linspace(0.075/m1, 1.0, 120)
m2 = qg*m1
mk1 = mks_of_mass(m1); mk2 = mks_of_mass(m2)
dM = -2.5*np.log10(1 + 10**(-0.4*(mk2 - mk1)))
bp = lambda m, mk: mk + colour_offsets(m, 'M_BP'); rp = lambda m, mk: mk + colour_offsets(m, 'M_RP')
BP = -2.5*np.log10(10**(-0.4*bp(m1, mk1)) + 10**(-0.4*bp(m2, mk2)))
RP = -2.5*np.log10(10**(-0.4*rp(m1, mk1)) + 10**(-0.4*rp(m2, mk2)))
dcol = (BP - RP) - (bp(m1, mk1) - rp(m1, mk1))
delta = dM - slope*dcol
print('m1', round(m1, 3), 'delta at q=1/0.8/0.6/0.4:', [round(float(np.interp(q, qg, delta)), 3) for q in (1, 0.8, 0.6, 0.4)],
      'colour shift at q=0.6:', round(float(np.interp(0.6, qg, dcol)), 3))

def pdf_S(x, mu, sb, sf, t, mt, st):
    a = 2/(np.sqrt(2*np.pi)*(sb + sf))
    core = np.where(x < mu, a*np.exp(-0.5*((x - mu)/sb)**2), a*np.exp(-0.5*((x - mu)/sf)**2))
    return (1 - t)*core + t*stats.norm.pdf(x, mt, st)
def pdf_B(x, par, Ft):
    mu, sb, sf, t, mt, st = par
    w = (1 - Ft)/(1 - qg[0])*np.ones_like(qg) + Ft*(qg >= 0.95)/0.05
    w = w/np.trapezoid(w, qg)
    return np.trapezoid(pdf_S(x[:, None] - delta[None, :], mu, sb, sf, t, mt, st)*w[None, :], qg, axis=1)
LO, HI = -1.8, 1.0
xn = np.linspace(LO, HI, 300)
def nll(theta, x, Ft, fix_f=None):
    mu, lsb, lsf, lt, mt, lst, lf = theta
    f = 1/(1 + np.exp(-lf)) if fix_f is None else fix_f
    par = (mu, np.exp(lsb), np.exp(lsf), 1/(1 + np.exp(-lt)), mt, np.exp(lst))
    p = (1 - f)*pdf_S(x, *par) + f*pdf_B(x, par, Ft)
    pn = (1 - f)*pdf_S(xn, *par) + f*pdf_B(xn, par, Ft)
    norm = np.trapezoid(pn, xn)
    return -np.sum(np.log(np.clip(p/norm, 1e-300, None)))
def fit(x, Ft, fix_f=None):
    best = None
    for lf0 in (-2.0, -0.7):
        for sb0 in (0.2,):
            th0 = [0.0, np.log(sb0), np.log(0.25), -2.0, 0.8, np.log(0.6), lf0]
            r = optimize.minimize(nll, th0, args=(x, Ft, fix_f), method='Nelder-Mead',
                                  options=dict(maxiter=3000, xatol=1e-3, fatol=1e-3))
            if best is None or r.fun < best.fun: best = r
    mu, lsb, lsf, lt, mt, lst, lf = best.x
    f = 1/(1 + np.exp(-lf)) if fix_f is None else fix_f
    par = (mu, np.exp(lsb), np.exp(lsf), 1/(1 + np.exp(-lt)), mt, np.exp(lst))
    pS = (1 - f)*pdf_S(np.array([x_t]), *par)[0]; pB = f*pdf_B(np.array([x_t]), par, Ft)[0]
    # breakdown of the binary posterior at x_t by q (twin q>0.9 vs others)
    w = (1 - Ft)/(1 - qg[0])*np.ones_like(qg) + Ft*(qg >= 0.95)/0.05; w = w/np.trapezoid(w, qg)
    integ = pdf_S(x_t - delta, *par)*w
    frac_twin = float(np.trapezoid(integ*(qg >= 0.9), qg)/np.trapezoid(integ, qg))
    frac_q_0p6_0p9 = float(np.trapezoid(integ*(qg >= 0.6)*(qg < 0.9), qg)/np.trapezoid(integ, qg))
    return dict(nll=float(best.fun), mu=float(mu), sig_bright=float(par[1]), sig_faint=float(par[2]), tail_frac=float(par[3]),
                tail_mu=float(mt), tail_sig=float(par[5]), f_bin=float(f), P_binary_at_target=float(pB/(pS + pB)),
                binary_post_frac_q_ge_0p9=frac_twin, binary_post_frac_q_0p6_0p9=frac_q_0p6_0p9,
                single_density_ratio_target_vs_mode=float(pdf_S(np.array([x_t]), *par)[0]/pdf_S(np.array([mu]), *par)[0]))
samples = {
    'clean_ruwe1.15': (np.abs(d.bp_rp - 2.756) < 0.10) & (d.ruwe < 1.15) & (d.ipd_frac_multi_peak == 0) & (d.ipd_gof_harmonic_amplitude < 0.04) & (d.non_single_star == 0),
    'ruwe1.1_only': (np.abs(d.bp_rp - 2.756) < 0.10) & (d.ruwe < 1.1),
}
out = dict(method=__doc__.split('\n')[0], ridge=dict(slope=float(slope), intercept=float(icpt), MKs_ridge_at_target=float(ridge(2.756149))),
           target_offset=float(x_t), m1_ridge=m1, results={})
for name, m in samples.items():
    x = d.off[m].values
    x = x[(x > LO) & (x < HI)]
    emp = dict(n=int(m.sum()), n_in_range=int(len(x)), frac_le_target=float(np.mean(d.off[m].values <= x_t)),
               pct=np.percentile(d.off[m].values, [1, 2.5, 5, 16, 50, 84]).tolist())
    res = {'empirical': emp}
    for Ft in (0.0, 0.1, 0.25):
        res[f'F_twin={Ft}'] = fit(x, Ft)
    # model without binaries (single Gaussian bright side fitted to everything): what the analyst effectively assumed
    res['no_binaries'] = fit(x, 0.1, fix_f=1e-9)
    # bootstrap of the baseline (F_twin 0.1)
    bs = []
    for i in range(10):
        xb = rng.choice(x, len(x), replace=True)
        r = fit(xb, 0.1); bs.append([r['sig_bright'], r['f_bin'], r['P_binary_at_target']])
    bs = np.array(bs)
    res['bootstrap_F_twin0.1'] = dict(sig_bright_16_50_84=np.percentile(bs[:, 0], [16, 50, 84]).tolist(),
                                      f_bin_16_50_84=np.percentile(bs[:, 1], [16, 50, 84]).tolist(),
                                      P_binary_16_50_84=np.percentile(bs[:, 2], [16, 50, 84]).tolist())
    out['results'][name] = res
    print(name, json.dumps(res, default=float)[:1500], flush=True)
jdump(out, 'out/r02_cmd_mixture.json')
