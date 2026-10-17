"""r07: diagnostics and figures. (1) PTP-weighted properties of the residual companion-hosted planet (radius, V/I depth
ratio, separation) in the referee-2 baseline; (2) analyst-CMD vs referee-2-CMD configuration posteriors side by side;
(3) figures: figs/r07_cmd_colour.png (comparison sample in the (x, r) plane with the target, the single-star metallicity
trend and the binary tracks), figs/r07_fpp_draws.png (FPP / FPP_strict / NFPP distributions of the 1500 joint draws).
Output: out/r07_diagnostics.json"""
import sys, json, numpy as np, pandas as pd
sys.path.insert(0, 'tic214311908/analysis/referee2_fpp/scripts')
import r05_fpp as F
from rlib import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
out = {}
for name, extra in (('referee2_baseline', {}), ('analyst_CMD_form', dict(cmd='gauss', sig_s=0.30, xobs=-0.753)),
                    ('referee2_INV', dict(ds='inv_beta'))):
    p = dict(F.BASE); p.update(extra)
    r = F.evaluate(p, detail=True)
    # PTP weights
    CP = F.CP
    Lc = CP['L_ipd']*CP['L_ast']*CP['L_dec']*CP['L_wide']
    Lcmd = F.cmd_like(p, F.DX_C, F.DR_C)
    q = CP['q']; qmin = 0.075/CP['m1']
    wq = ((1 - p['F_twin'])/1.0*np.ones_like(q) + p['F_twin']*(q >= 0.95)/0.05)*(1 - qmin)
    tab = F.build(p['ds'], p['ecc'], p['alpha'])
    Z1, Z2 = F.ZCACHE[(p['ds'], p['ecc'], p['alpha'])]
    Sa = np.where(CP['a'] <= 10, p['S_close'], np.where(CP['a'] >= 200, 1.0, p['S_close'] + (1 - p['S_close'])*np.log(CP['a']/10)/np.log(20)))
    occ = lambda cl, M: p['f_' + cl]/F.LNP*(np.asarray(M)/0.37)**p['beta_' + cl]
    p2 = sum(occ(cl, CP['m2'])*Z2[cl] for cl in ('sub', 'N', 'S', 'J'))*Sa*CP['stab2']
    w = wq*Lc*Lcmd*p2*(CP['dT'] > p['dT_twin'])
    ratioVI = CP['DcV']/CP['DcI']
    Rp = np.sqrt(0.0157/CP['Dc'])*1.03*CP['R2']*RSUN/REARTH
    pct = lambda v: [float(x) for x in np.interp([0.16, 0.5, 0.84], np.cumsum(w[np.argsort(v)])/w.sum(), np.sort(v))]
    out[name] = dict(FPP=r['FPP'], FPP_strict=r['FPP_strict'], FPP_system=r['FPP_system'], NFPP=r['NFPP'], config=r['config'],
                     post=r['post'], PTP_VoverI_16_50_84=pct(ratioVI), PTP_Rp_16_50_84=pct(Rp), PTP_sep_16_50_84=pct(CP['sep']),
                     PTP_a_16_50_84=pct(CP['a']), PTP_m2_16_50_84=pct(CP['m2']), PTP_dT_16_50_84=pct(CP['dT']))
    print(name, json.dumps(out[name], default=float)[:900])
jdump(out, 'out/r07_diagnostics.json')
# ---- figure 1: (x, r) plane
J = json.load(open(f'{R2}/out/r02_cmd_colour.json'))
d = pd.read_csv(f'{B}/referee_star/data/mdw_cmd.csv')
d = d[np.isfinite(d.ks_m) & (d.ph_qual.astype(str).str[2] == 'A') & (d.ks_msigcom < 0.05) & np.isfinite(d.ruwe)].copy()
d['x'] = d.ks_m + 5*np.log10(d.parallax/100.) - (J['ridge']['intercept'] + J['ridge']['slope']*d.bp_rp)
d['r'] = d.phot_g_mean_mag - d.ks_m - np.polyval(J['GK_fit'], d.bp_rp)
s = d[(np.abs(d.bp_rp - 2.756) < 0.1) & (d.ruwe < 1.15) & (d.ipd_frac_multi_peak == 0) & (d.ipd_gof_harmonic_amplitude < 0.04) & (d.non_single_star == 0)]
best = [r_ for r_ in J['profiles']['F_twin=0.1']['rows'] if abs(r_['sig_b'] - 0.40) < 1e-6][0]['params']
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
ax[0].scatter(s.x, s.r, s=3, c='0.6', label='Gaia-clean comparison M dwarfs (BP-RP 2.66-2.86)')
xx = np.linspace(-1.4, 1.0, 50)
ax[0].plot(xx, best['a'] + best['c']*xx + best['c2']*xx**2, 'C0', lw=2, label='single-star trend (fit)')
sh = J['shifts']; dx = np.array(sh['dx']); dr = np.array(sh['dr'])
for x0 in (-0.1, 0.0, 0.1):
    ax[0].plot(x0 + dx, best['a'] + best['c']*x0 + dr, 'C1', lw=1, alpha=0.7)
ax[0].plot([], [], 'C1', label='binary tracks (q = 0.3 -> 1), normal primary')
ax[0].plot(J['target']['x'], J['target']['r'], 'r*', ms=16, label='TIC 214311908')
ax[0].set_xlim(-1.5, 1.0); ax[0].set_ylim(-0.25, 0.3); ax[0].invert_xaxis()
ax[0].set_xlabel('M_Ks offset from ridge at fixed BP-RP (mag; left = over-luminous)'); ax[0].set_ylabel('G-Ks residual at fixed BP-RP (mag)')
ax[0].legend(fontsize=7, loc='upper left')
for ft, c in (('F_twin=0.0', 'C0'), ('F_twin=0.1', 'C1'), ('F_twin=0.25', 'C2')):
    rows = J['profiles'][ft]['rows']; sb = [r_['sig_b'] for r_ in rows]; nl = np.array([r_['nll'] for r_ in rows])
    ax[1].plot(sb, 2*(nl - nl.min()), 'o-', color=c, label=f'{ft}: -2 dlnL')
    ax[1].plot(sb, [100*r_['P_binary_given_x_r'] for r_ in rows], 's--', color=c, alpha=0.6, label=f'{ft}: 100 x P(binary | x, r)')
ax[1].axvline(0.30, color='k', ls=':', label='analyst sigma_s = 0.30')
ax[1].set_xlabel('single-star bright-side width sigma_b (mag)'); ax[1].set_ylim(0, 110); ax[1].legend(fontsize=7)
plt.tight_layout(); plt.savefig(f'{R2}/figs/r07_cmd_colour.png', dpi=120)
# ---- figure 2: draws
R = pd.read_csv(f'{R2}/out/r05_random.csv')
fig, ax = plt.subplots(1, 3, figsize=(12, 3.6))
for a_, col, thr in zip(ax, ('FPP', 'FPP_strict', 'NFPP'), (0.015, 0.015, 1e-3)):
    v = np.log10(np.clip(R[col].values, 1e-12, None))
    a_.hist(v, 50, color='0.5'); a_.axvline(np.log10(thr), color='r', ls='--'); a_.set_xlabel(f'log10 {col}')
    a_.set_title(f'{col}: median {np.median(R[col]):.2g}, P(<thr) {np.mean(R[col] < thr):.2f}', fontsize=9)
plt.tight_layout(); plt.savefig(f'{R2}/figs/r07_fpp_draws.png', dpi=120)
print('figures written')
