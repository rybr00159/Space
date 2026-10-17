"""h07: What does the Gaia DR3 RV (16.28 +- 2.22 km/s, 21 transits, rv_method_used = 2 'faint', grvs_mag 12.55,
rv_expected_sig_to_noise 14.7, 12 deblended transits) say about an SB2 twin?
 (1) Empirical: rv_error percentile of the M3V among comparison M dwarfs with rv_method_used = 2 and similar
     grvs_mag (+-0.25), rv_expected_sig_to_noise (+-25%), rv_nb_transits (15-30), RUWE < 1.2, and near-normal CMD
     position (|dMKs| < 0.25). Same for over-luminous (dMKs < -0.5) RUWE < 1.2 stars (twins + metal-rich) and for
     RUWE > 2 stars (astrometric binaries) - does binarity inflate the faint-star RV error empirically?
     (rv_chisq_pvalue / rv_renormalised_gof / rv_amplitude_robust are not published for faint stars.)
 (2) Simulation: for faint stars DR3 measures the RV from the combined (summed) cross-correlation function of all
     transits and derives the error from the CCF peak height and curvature (Zucker 2003). Single star: Gaussian CCF of
     width s_w (11/15/20 km/s). Equal-mass twin (2 x 0.264 Msun; flux ratio 1 in G_RVS): per-transit CCF =
     0.5 G(v-v1) + 0.5 G(v-v2); combined over 21 transits in 8 visibility periods spanning 974 d. Orbits: a log-uniform
     0.01-10 au, isotropic inclination, e = 0 for P < 10 d else uniform 0-0.8 (capped by e_max(P)), random phase.
     Likelihood ratio twin/single = E_orbits[p(log e_obs - log R)] / p(log e_obs) with p the KDE of log(rv_error) of normal comparison stars. Error inflation R = sigma_bin/sigma_single from Zucker's formula at the highest peak of the combined CCF
     (sigma ~ 1/sqrt(C * |C''|)); a 0.5 km/s floor in quadrature (option). RV offset of the chosen peak from systemic.
     P(observed error <= 2.217 | twin, a) = mean over orbits of F_single(2.217 / R), with F_single the empirical CDF
     of rv_error for normal comparison stars (step 1). 'Flagged' = simulated error above the 95th percentile of the
     normal comparison stars' errors (3.6 km/s).
Output: out/h07_gaia_rv.json, figs/h07_gaia_rv.png"""
import json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
B = 'tic214311908/analysis/host'
rng = np.random.default_rng(7)
T = dict(rv_err=2.217367, grvs=12.548391, sn=14.695826, ntr=21, teff_t=3500.)
m = pd.read_csv(f'{B}/data/field_mdwarfs.csv')
m['MK'] = m.ks_m + 5 * np.log10(m.parallax / 100)
# CMD ridge (KDE-mode in 0.04-mag colour bins, quadratic) as in h10; here a simple robust quadratic of the mode
from scipy.stats import gaussian_kde
s0 = m[(m.ruwe < 1.4) & (m.ph_qual == 'AAA') & (m.bp_rp > 2.45) & (m.bp_rp < 3.05)]
pts = []
for x0 in np.arange(2.50, 3.01, 0.02):
    k = np.abs(s0.bp_rp - x0) < 0.02
    yy = s0.MK[k].values; kde = gaussian_kde(yy, bw_method=0.25); gr = np.linspace(np.percentile(yy, 1), np.percentile(yy, 99), 400)
    pts.append((x0, gr[np.argmax(kde(gr))]))
pts = np.array(pts); cr = np.polyfit(pts[:, 0], pts[:, 1], 2)
m['dMK'] = m.MK - np.polyval(cr, m.bp_rp)
c2 = (m.rv_method_used == 2) & np.isfinite(m.radial_velocity_error) & (np.abs(m.grvs_mag - T['grvs']) < 0.25) & \
     (np.abs(m.rv_expected_sig_to_noise / T['sn'] - 1) < 0.25) & m.rv_nb_transits.between(15, 30) & (m.rv_template_teff <= 3700) & (m.ph_qual == 'AAA')
groups = {'normal': c2 & (m.ruwe < 1.2) & (m.dMK.abs() < 0.25), 'overlum_lowruwe': c2 & (m.ruwe < 1.2) & (m.dMK < -0.5),
          'overlum_ruwe_lt1.1': c2 & (m.ruwe < 1.1) & (m.dMK < -0.5) & (m.dMK > -1.0), 'ruwe_gt2': c2 & (m.ruwe > 2)}
out = {'target': T, 'ridge_coeffs': cr.tolist()}
emp = {}
for g, k in groups.items():
    e = m.radial_velocity_error[k].values
    emp[g] = dict(n=int(k.sum()), median=float(np.median(e)) if k.sum() else None, p16=float(np.percentile(e, 16)) if k.sum() else None,
                  p84=float(np.percentile(e, 84)) if k.sum() else None, p95=float(np.percentile(e, 95)) if k.sum() else None,
                  target_pct=float(np.mean(e <= T['rv_err'])) if k.sum() else None,
                  frac_gt_1p5x_median=float(np.mean(e > 1.5 * np.median(e))) if k.sum() else None)
    print(g, emp[g])
out['empirical'] = emp
# deblended-transit dependence (target has 12/21 deblended)
k = groups['normal']; fr = m.rv_nb_deblended_transits[k] / m.rv_nb_transits[k]
emp['normal_deblended_frac_gt_0.3'] = dict(n=int((fr > 0.3).sum()), median=float(np.median(m.radial_velocity_error[k][fr > 0.3])) if (fr > 0.3).sum() else None,
                                           target_pct=float(np.mean(m.radial_velocity_error[k][fr > 0.3] <= T['rv_err'])) if (fr > 0.3).sum() else None)
print('normal with >30% deblended transits', emp['normal_deblended_frac_gt_0.3'])
e_norm = np.sort(m.radial_velocity_error[groups['normal']].values)
def Fsingle(x):
    return np.searchsorted(e_norm, x, side='right') / len(e_norm)
kde_log = gaussian_kde(np.log(e_norm), bw_method=0.2)
def Lratio(R, floor=0.0):
    # likelihood ratio twin/single for the observed error: density of log(err) shifted by log R (with optional floor)
    y = np.log(T['rv_err'])
    if floor > 0:
        s = np.sqrt(np.clip((T['rv_err'] ** 2 - floor ** 2) / R ** 2 + floor ** 2, 1e-6, None)); ys = np.log(s)
        # Jacobian d ln s / d ln e_obs
        jac = (T['rv_err'] ** 2 / R ** 2) / s ** 2
        return float(np.mean(kde_log(ys) * jac) / kde_log(y)[0])
    return float(np.mean(kde_log(y - np.log(R))) / kde_log(y)[0])
# ---------- simulation
Mtot = 2 * 0.264
def gaia_epochs():
    vp = np.sort(rng.uniform(0, 974.27, 8))
    t = []
    nper = rng.multinomial(21 - 16, np.ones(8) / 8) + 2       # 21 transits: >=2 per visibility period
    for v, n in zip(vp, nper):
        for j in range(n):
            t.append(v + (j // 2) * 0.25 + (j % 2) * 106.5 / 1440)
    return np.array(t)
def kepler(M, e):
    E = M.copy()
    for _ in range(40): E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    return E
vgrid = np.arange(-150, 150.01, 0.25)
def sim(a, sw, nsim=400, floor=0.0):
    P = np.sqrt(a ** 3 / Mtot) * 365.25
    R = np.zeros(nsim); off = np.zeros(nsim)
    G0 = np.exp(-0.5 * vgrid ** 2 / sw ** 2)
    C0 = 1.0; curv0 = C0 / sw ** 2
    for i in range(nsim):
        t = gaia_epochs()
        e = 0.0 if P < 10 else min(rng.uniform(0, 0.8), 1 - (P / 2.0) ** (-2. / 3.))
        inc = np.arccos(rng.uniform(-1, 1)); om = rng.uniform(0, 2 * np.pi); tp = rng.uniform(0, P)
        K = (2 * np.pi * a * 1.496e8 / (P * 86400)) / np.sqrt(1 - e ** 2) * np.sin(inc) / 2   # each component, equal mass
        E = kepler(2 * np.pi * (t - tp) / P % (2 * np.pi), e)
        nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
        v1 = K * (np.cos(nu + om) + e * np.cos(om)); v2 = -v1
        C = np.zeros_like(vgrid)
        for a1, a2 in zip(v1, v2):
            C += 0.5 * np.exp(-0.5 * (vgrid - a1) ** 2 / sw ** 2) + 0.5 * np.exp(-0.5 * (vgrid - a2) ** 2 / sw ** 2)
        C /= len(t)
        j = np.argmax(C); j = min(max(j, 2), len(vgrid) - 3)
        h = vgrid[1] - vgrid[0]
        cpp = (C[j + 1] - 2 * C[j] + C[j - 1]) / h ** 2
        sig_ratio = np.sqrt(curv0 * C0 / max(-cpp * C[j], 1e-12))
        R[i] = sig_ratio; off[i] = vgrid[j]
    return R, off
agrid = np.logspace(-2, 1, 16)
res = {}
for sw in (11.0, 15.0, 20.0):
    rows = []
    for a in agrid:
        R, off = sim(a, sw)
        # error with optional floor: sigma_obs = sqrt((R*s)^2 + f^2); invert for s: P(sigma_obs <= 2.217)
        Ppass = float(np.mean(Fsingle(T['rv_err'] / R)))
        f = 0.5; s_lim = np.sqrt(np.clip((T['rv_err'] ** 2 - f ** 2) / R ** 2 + f ** 2, 0, None))
        Ppass_floor = float(np.mean(Fsingle(s_lim)))
        es = rng.choice(e_norm, len(R)); p95 = np.percentile(e_norm, 95)
        flag95 = float(np.mean(R * es > p95)); flag_obs = float(np.mean(R * es > T['rv_err']))
        rows.append(dict(a_au=float(a), P_days=float(np.sqrt(a ** 3 / Mtot) * 365.25), R_median=float(np.median(R)), R_p84=float(np.percentile(R, 84)),
                         frac_flagged_gt_p95_of_singles=flag95, frac_error_gt_observed=flag_obs,
                         frac_R_gt_1p5=float(np.mean(R > 1.5)), P_pass=Ppass, P_pass_floor0p5=Ppass_floor, L_ratio=Lratio(R), L_ratio_floor0p5=Lratio(R, 0.5),
                         rv_offset_abs_median=float(np.median(np.abs(off))), frac_offset_gt_5=float(np.mean(np.abs(off) > 5))))
        print(f'sw {sw}: a {a:.3f} au P {rows[-1]["P_days"]:.1f} d  R_med {rows[-1]["R_median"]:.2f} R84 {rows[-1]["R_p84"]:.2f}  '
              f'P(err<=obs) {Ppass:.3f} (floor {Ppass_floor:.3f}) Lratio {rows[-1]["L_ratio"]:.3f} ({rows[-1]["L_ratio_floor0p5"]:.3f}) |offset|>5 km/s {rows[-1]["frac_offset_gt_5"]:.2f}', flush=True)
    res[f'sw{int(sw)}'] = rows
out['simulation'] = res
out['P_pass_single'] = float(Fsingle(T['rv_err']))
json.dump(out, open(f'{B}/out/h07_gaia_rv.json', 'w'), indent=1)
fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
for g, k in groups.items():
    if k.sum() > 5:
        ax[0].hist(m.radial_velocity_error[k], bins=np.arange(0, 12, 0.25), histtype='step', density=True, label=f'{g} (n={k.sum()})')
ax[0].axvline(T['rv_err'], c='k', ls='--', label='M3V'); ax[0].set_xlabel('DR3 radial_velocity_error (km/s)'); ax[0].legend(fontsize=7)
for key, rows in res.items():
    ax[1].plot([r['a_au'] for r in rows], [r['L_ratio'] for r in rows], label=key)
ax[1].set_xscale('log'); ax[1].set_xlabel('twin semi-major axis (au)'); ax[1].set_ylabel('L(RV error | twin) / L(single)'); ax[1].legend()
plt.tight_layout(); plt.savefig(f'{B}/figs/h07_gaia_rv.png', dpi=110)
