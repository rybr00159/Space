"""r02: empirical calibration of the over-luminosity likelihood, using BOTH the M_Ks offset x (from the ridge at fixed
BP-RP) and the G-Ks colour residual r at fixed BP-RP (a test the fpp analysis did not use).
Why: a near-equal twin brightens M_Ks by 0.75 mag but leaves the colours of a normal star (r ~ 0), whereas a metal-rich
single star is both over-luminous and red in G-Ks at fixed BP-RP (the comparison sample shows r rising steadily from
-0.04 at x = +0.6 to +0.07 at x = -0.8). The target has x = -0.84 and r = +0.097.
Model for Gaia-clean comparison M dwarfs at the target colour (|BP-RP - 2.756| < 0.1, RUWE < 1.15, ipd_frac_multi_peak
= 0, gof_harmonic < 0.04, non_single_star = 0), fitted by maximum likelihood on -1.6 < x < 1.0:
  single:  x ~ (1 - t) SplitNormal(mu, sig_b, sig_f) + t N(mt, st);  r | x ~ N(a + c x + c2 x^2, s_r)
  binary:  (x, r) = single (x_s, r_s) shifted by (dx(q), dr(q)) of an unresolved companion (Mann+19 Ks budget,
           Mamajek colours; primary 0.253 Msun = the ridge star); q ~ (1 - F_twin) U(qmin, 1) + F_twin U(0.95, 1)
  p = (1 - f_bin) single + f_bin binary
sig_b (single-star bright-side width) is profiled on a grid because it is degenerate with f_bin in x alone.
Outputs for every sig_b: fitted parameters, profile log-likelihood, and the population posterior at the target:
P(binary | x), P(binary | x, r), and the q breakdown. Output: out/r02_cmd_colour.json"""
import numpy as np, pandas as pd, sys, json
from scipy import optimize, stats
sys.path.insert(0, 'tic214311908/analysis/referee2_fpp/scripts')
from rlib import *
d = pd.read_csv(f'{B}/referee_star/data/mdw_cmd.csv')
d = d[np.isfinite(d.ks_m) & (d.ph_qual.astype(str).str[2] == 'A') & (d.ks_msigcom < 0.05) & np.isfinite(d.ruwe)].copy()
d['MKs'] = d.ks_m + 5*np.log10(d.parallax/100.)
d['GK'] = d.phot_g_mean_mag - d.ks_m
# ridge (KDE modes in colour bins, RUWE < 1.1)
cb = np.arange(2.55, 2.951, 0.05); modes, cens = [], []
for lo in cb[:-1]:
    s = d[(d.bp_rp >= lo) & (d.bp_rp < lo + 0.05) & (d.ruwe < 1.1)].MKs.values
    xs = np.linspace(np.percentile(s, 1), np.percentile(s, 99), 3000)
    modes.append(xs[np.argmax(stats.gaussian_kde(s, 0.12)(xs))]); cens.append(np.median(d[(d.bp_rp >= lo) & (d.bp_rp < lo + 0.05)].bp_rp))
SL, IC = np.polyfit(cens, modes, 1)
d['x'] = d.MKs - (IC + SL*d.bp_rp)
core = d[(d.ruwe < 1.1) & (np.abs(d.x) < 0.25)]
CG = np.polyfit(core.bp_rp, core.GK, 1)
d['r'] = d.GK - np.polyval(CG, d.bp_rp)
TBPRP, TG, TK = 2.756149, 14.281043, 10.482
x_t = MKS_OBS - (IC + SL*TBPRP); r_t = (TG - TK) - np.polyval(CG, TBPRP)
print(f'ridge slope {SL:.3f}; target x = {x_t:.3f}, r = {r_t:.3f}; G-Ks slope {CG[0]:.3f}', flush=True)
# binary shifts for the ridge primary
m1 = float(mass_of_mks(IC + SL*TBPRP))
qg = np.linspace(0.075/m1, 1.0, 50)
k1 = mks_of_mass(m1); k2 = mks_of_mass(qg*m1)
fl = lambda mm: 10**(-0.4*mm)
def comp(col):
    a1 = k1 + colour_offsets(m1, col); a2 = k2 + colour_offsets(qg*m1, col)
    return -2.5*np.log10(fl(a1) + fl(a2)) - a1
dK = -2.5*np.log10(fl(k1) + fl(k2)) - k1
dG, dBP, dRP = comp('M_G'), comp('M_BP'), comp('M_RP')
dcol = dBP - dRP
DX = dK - SL*dcol
DR = (dG - dK) - CG[0]*dcol
sel = (np.abs(d.bp_rp - 2.756) < 0.10) & (d.ruwe < 1.15) & (d.ipd_frac_multi_peak == 0) & (d.ipd_gof_harmonic_amplitude < 0.04) & (d.non_single_star == 0)
LO, HI = -1.6, 1.0
S = d[sel & (d.x > LO) & (d.x < HI)]
X, Rr = S.x.values, S.r.values
xn = np.linspace(LO, HI, 260)
def px_single(x, mu, sb, sf, t, mt, st):
    a = 2/(np.sqrt(2*np.pi)*(sb + sf))
    c = np.where(x < mu, a*np.exp(-0.5*((x - mu)/sb)**2), a*np.exp(-0.5*((x - mu)/sf)**2))
    return (1 - t)*c + t*np.exp(-0.5*((x - mt)/st)**2)/(np.sqrt(2*np.pi)*st)
def pr_given_x(r, x, a, c, c2, sr):
    return np.exp(-0.5*((r - a - c*x - c2*x*x)/sr)**2)/(np.sqrt(2*np.pi)*sr)
def qweights(Ft):
    w = (1 - Ft)/(1 - qg[0])*np.ones_like(qg) + Ft*(qg >= 0.95)/0.05
    return w/np.trapezoid(w, qg)
def unpack(th, sb):
    mu, lsf, lt, mt, lst, a, c, c2, lsr, lf = th
    return dict(mu=mu, sb=sb, sf=np.exp(lsf), t=1/(1 + np.exp(-lt)), mt=mt, st=np.exp(lst), a=a, c=c, c2=c2, sr=np.exp(lsr), f=1/(1 + np.exp(-lf)))
def dens(x, r, p, Ft, use_r=True):
    xs = (p['mu'], p['sb'], p['sf'], p['t'], p['mt'], p['st'])
    ps = px_single(x, *xs)*(pr_given_x(r, x, p['a'], p['c'], p['c2'], p['sr']) if use_r else 1)
    w = qweights(Ft)
    xx = x[:, None] - DX[None, :]
    pb = px_single(xx, *xs)*(pr_given_x(r[:, None] - DR[None, :], xx, p['a'], p['c'], p['c2'], p['sr']) if use_r else 1)
    pb = np.trapezoid(pb*w[None, :], qg, axis=1)
    return (1 - p['f'])*ps, p['f']*pb
def norm_x(p, Ft):
    s, b = dens(xn, np.zeros_like(xn), p, Ft, use_r=False)
    return np.trapezoid(s + b, xn)
def nll(th, sb, Ft):
    p = unpack(th, sb)
    s, b = dens(X, Rr, p, Ft)
    return -np.sum(np.log(np.clip(s + b, 1e-300, None))) + len(X)*np.log(norm_x(p, Ft))
BOUNDS = [(-0.3, 0.3), (np.log(0.15), np.log(1.0)), (-6, -0.5), (0.2, 2.0), (np.log(0.15), np.log(1.5)), (-0.1, 0.1), (-0.2, 0.1),
          (-0.1, 0.1), (np.log(0.015), np.log(0.1)), (-7, 0.5)]
def fit(sb, Ft, starts=((0.0, np.log(0.35), -2.5, 0.9, np.log(0.5), 0.0, -0.06, 0.0, np.log(0.035), -2.5),
                        (0.0, np.log(0.35), -2.5, 0.9, np.log(0.5), 0.0, -0.06, 0.0, np.log(0.035), -1.0))):
    best = None
    for th0 in starts:
        r = optimize.minimize(nll, th0, args=(sb, Ft), method='L-BFGS-B', bounds=BOUNDS, options=dict(maxiter=500))
        if best is None or r.fun < best.fun: best = r
    p = unpack(best.x, sb)
    s1, b1 = dens(np.array([x_t]), np.array([r_t]), p, Ft)
    s0, b0 = dens(np.array([x_t]), np.array([r_t]), p, Ft, use_r=False)
    w = qweights(Ft); xs = (p['mu'], p['sb'], p['sf'], p['t'], p['mt'], p['st'])
    integ = px_single(x_t - DX, *xs)*pr_given_x(r_t - DR, x_t - DX, p['a'], p['c'], p['c2'], p['sr'])*w
    return dict(nll=float(best.fun), params={k: float(v) for k, v in p.items()},
                P_binary_given_x=float(b0[0]/(s0[0] + b0[0])), P_binary_given_x_r=float(b1[0]/(s1[0] + b1[0])),
                frac_q_ge_0p9_given_x_r=float(np.trapezoid(integ*(qg >= 0.9), qg)/np.trapezoid(integ, qg)),
                frac_q_0p6_0p9_given_x_r=float(np.trapezoid(integ*(qg >= 0.6)*(qg < 0.9), qg)/np.trapezoid(integ, qg)),
                r_expected_single_at_xt=float(p['a'] + p['c']*x_t + p['c2']*x_t**2))
out = dict(method=__doc__.split('\n')[0], n=int(len(X)), ridge=dict(slope=float(SL), intercept=float(IC)), GK_fit=CG.tolist(),
           target=dict(x=float(x_t), r=float(r_t)), m1_ridge=m1, shifts=dict(q=qg.tolist(), dx=DX.tolist(), dr=DR.tolist()), profiles={})
SBG = [0.12, 0.16, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45]
for Ft in (0.0, 0.1, 0.25):
    rows = []
    for sb in SBG:
        rr = fit(sb, Ft); rr['sig_b'] = sb; rows.append(rr)
        print(f"F_twin {Ft} sig_b {sb:.2f}: nll {rr['nll']:.2f} f_bin {rr['params']['f']:.3f} P(bin|x) {rr['P_binary_given_x']:.3f} "
              f"P(bin|x,r) {rr['P_binary_given_x_r']:.3f} q>=0.9 {rr['frac_q_ge_0p9_given_x_r']:.2f} r_single(xt) {rr['r_expected_single_at_xt']:.3f} "
              f"c {rr['params']['c']:.3f} sr {rr['params']['sr']:.3f}", flush=True)
    nl = np.array([r_['nll'] for r_ in rows]); wts = np.exp(-(nl - nl.min())); wts /= wts.sum()
    out['profiles'][f'F_twin={Ft}'] = dict(rows=rows, profile_weights=wts.tolist(),
                                           marg_P_binary_given_x=float(np.sum(wts*[r_['P_binary_given_x'] for r_ in rows])),
                                           marg_P_binary_given_x_r=float(np.sum(wts*[r_['P_binary_given_x_r'] for r_ in rows])),
                                           best_sig_b=SBG[int(np.argmin(nl))])
    print('F_twin', Ft, 'marginal P(bin|x) %.3f, P(bin|x,r) %.3f, best sig_b %.2f' % (out['profiles'][f'F_twin={Ft}']['marg_P_binary_given_x'],
          out['profiles'][f'F_twin={Ft}']['marg_P_binary_given_x_r'], out['profiles'][f'F_twin={Ft}']['best_sig_b']), flush=True)
jdump(out, 'out/r02_cmd_colour.json')
