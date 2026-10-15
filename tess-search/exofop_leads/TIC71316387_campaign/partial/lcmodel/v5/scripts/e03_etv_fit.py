"""Fit the eclipse timings: linear ephemeris (T0, P, secondary offset), quadratic (+ Q n^2), cubic, and per-sector
binned O-C.  Per-event errors = formal x z-scale (z-scale from fake-phase injection-recovery for that product and
eclipse type).  Events > 5 sigma from the linear fit are rejected (iteratively).  Empirical null for the
quadratic term: 20000 synthetic timing sets with the real epochs, residuals drawn from the injection-recovery error
distribution of each product/eclipse type (bootstrap of normalised errors), fitted identically.
Usage: e03_etv_fit.py <etv json> [label]"""
import numpy as np, sys, json
V = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v5'
fn = sys.argv[1]; label = sys.argv[2] if len(sys.argv) > 2 else 'default'
E = json.load(open(fn))
stats = E.get('stats') or json.load(open(f'{V}/results/etv_default.json'))['stats']
inj = E.get('injections') or json.load(open(f'{V}/results/etv_default.json'))['injections']
P = E['P']; T0 = E['T0']
ev = [e for e in E['events'] if np.isfinite(e['delta_min']) and np.isfinite(e['sig_formal_min'])]
n = np.array([e['n'] for e in ev], float); ph = np.array([e['ph'] for e in ev]); key = np.array([e['key'] for e in ev])
tobs = np.array([e['t_obs'] for e in ev])
zs = np.array([stats[f"{e['key']}_ph{e['ph']}"]['zscale'] for e in ev])
sig = np.array([e['sig_formal_min'] for e in ev]) * np.maximum(zs, 1.0) / 1440.0
yr = 2018.0 + (tobs - 1325.3) / 365.25 + (2018.57 - 2018.0)   # approximate decimal year (S1 starts 2018-07-25)

def design(kind, n, ph):
    cols = [np.ones_like(n), n, (ph == 0.5).astype(float)]
    if kind in ('quad', 'cubic'): cols.append(n ** 2 / 1e6)
    if kind == 'cubic': cols.append(n ** 3 / 1e9)
    if kind == 'quad_sep':   # separate quadratic for primaries and secondaries
        cols += [np.where(ph == 0, n ** 2, 0) / 1e6, np.where(ph == 0.5, n ** 2, 0) / 1e6]
    if kind == 'lin_sep_P':  # separate period for secondaries (apsidal-motion-like divergence)
        cols.append(np.where(ph == 0.5, n, 0.0))
    return np.vstack(cols).T

def wls(A, y, s):
    Aw = A / s[:, None]; yw = y / s
    c, *_ = np.linalg.lstsq(Aw, yw, rcond=None)
    cov = np.linalg.inv(Aw.T @ Aw)
    r = y - A @ c
    return c, cov, float(((r / s) ** 2).sum()), r

keep = np.ones(n.size, bool)
for it in range(5):
    A = design('lin', n[keep], ph[keep]); c, cov, chi, r = wls(A, tobs[keep], sig[keep])
    rall = tobs - design('lin', n, ph) @ c
    new = np.abs(rall / sig) < 5
    if np.array_equal(new, keep): break
    keep = new
nrej = int((~keep).sum())
nn, pp, tt, ss, kk = n[keep], ph[keep], tobs[keep], sig[keep], key[keep]
res = dict(label=label, n_events=int(n.size), n_rejected_5sigma=nrej, rejected=[dict(key=str(key[i]), n=int(n[i]), ph=float(ph[i]), oc_min=float(rall[i] * 1440), sig_min=float(sig[i] * 1440)) for i in np.where(~keep)[0]])
fits = {}
for kind in ('lin', 'quad', 'cubic', 'quad_sep', 'lin_sep_P'):
    A = design(kind, nn, pp); c, cov, chi, r = wls(A, tt, ss)
    k = A.shape[1]; N = tt.size
    fits[kind] = dict(chi2=chi, dof=N - k, bic=chi + k * np.log(N), coef=c.tolist(), err=np.sqrt(np.diag(cov)).tolist())
lin = fits['lin']; quad = fits['quad']
res['fits'] = fits
res['ephemeris'] = dict(T0_btjd=lin['coef'][0], T0_err_min=lin['err'][0] * 1440, P_d=lin['coef'][1], P_err_s=lin['err'][1] * 86400,
                        sec_offset_min=lin['coef'][2] * 1440, sec_offset_err_min=lin['err'][2] * 1440,
                        sec_minus_halfP_min=lin['coef'][2] * 1440 + E['dsec_d'] * 1440 - 0.5 * P * 1440 if False else None)
# secondary offset relative to exact half period: tobs includes template dsec; recompute from raw
A = design('lin', nn, pp); c, cov, chi, r = wls(A, tt, ss)
# T(n, sec) = c0 + c1 n + c2 where the secondary events have n counted from the primary of the same cycle and
# t_obs = T0 + n P + dsec + delta  -> c2 = (dsec + <delta_sec>) - <delta_prim> = secondary - primary - 0 * P;  minus P/2:
res['ephemeris']['sec_minus_halfP_min'] = float((c[2] - 0.5 * c[1]) * 1440)
res['ephemeris']['sec_minus_halfP_err_min'] = float(np.sqrt(cov[2, 2] + 0.25 * cov[1, 1] - cov[1, 2]) * 1440)
# quadratic term: Q in d/cycle^2 ->  dP/dt, and quadratic O-C amplitude over the baseline
Q = quad['coef'][3] / 1e6; Qe = quad['err'][3] / 1e6
span_cyc = nn.max() - nn.min()
res['quadratic'] = dict(Q_d_per_cycle2=Q, Q_err=Qe, dPdt_s_per_yr=2 * Q / P * 365.25 * 86400, dPdt_err_s_per_yr=2 * Qe / P * 365.25 * 86400,
                        curvature_min_over_baseline=Q * (span_cyc / 2) ** 2 * 1440, curvature_err_min=Qe * (span_cyc / 2) ** 2 * 1440,
                        c_min_per_yr2=Q / P ** 2 * 365.25 ** 2 * 1440, c_err_min_per_yr2=Qe / P ** 2 * 365.25 ** 2 * 1440,
                        dchi2_lin_minus_quad=lin['chi2'] - quad['chi2'], dBIC_lin_minus_quad=lin['bic'] - quad['bic'])
# empirical null for the quadratic term via bootstrap of injection-recovery normalised errors
pool = {}
for x in inj:
    if not np.isfinite(x['rec_min']): continue
    k_ = f"{x['key']}_ph{x['ph']}"
    pool.setdefault(k_, []).append((x['rec_min'] - x['inj_min'] - stats[k_]['bias_min']) / (x['sig_formal_min'] * max(stats[k_]['zscale'], 1.0)))
pool = {k_: np.array(v) for k_, v in pool.items()}
rng = np.random.default_rng(3)
kk_ph = np.array([f'{a}_ph{b}' for a, b in zip(kk, pp)])
Aq = design('quad', nn, pp); Al = design('lin', nn, pp)
nul_dchi = []; nul_Q = []
for it in range(20000):
    e = np.empty(tt.size)
    for k_ in np.unique(kk_ph):
        m = kk_ph == k_; e[m] = rng.choice(pool[k_], m.sum()) * ss[m]
    y = Al @ c + e
    cl, _, chl, _ = wls(Al, y, ss); cq, _, chq, _ = wls(Aq, y, ss)
    nul_dchi.append(chl - chq); nul_Q.append(cq[3] / 1e6)
nul_dchi = np.array(nul_dchi); nul_Q = np.array(nul_Q)
res['quadratic']['null_p_dchi2'] = float(np.mean(nul_dchi >= res['quadratic']['dchi2_lin_minus_quad']))
res['quadratic']['null_Q_sd'] = float(np.std(nul_Q)); res['quadratic']['null_Q_sd_over_formal'] = float(np.std(nul_Q) / Qe)
res['quadratic']['null_Q_95_abs'] = float(np.percentile(np.abs(nul_Q), 95))
# linear-fit chi2 calibration against the null
nul_chil = []
for it in range(4000):
    e = np.empty(tt.size)
    for k_ in np.unique(kk_ph):
        m = kk_ph == k_; e[m] = rng.choice(pool[k_], m.sum()) * ss[m]
    nul_chil.append(wls(Al, Al @ c + e, ss)[2])
res['lin_chi2_null_median'] = float(np.median(nul_chil)); res['lin_chi2_null_p'] = float(np.mean(np.array(nul_chil) >= lin['chi2']))
# per-sector binned O-C (relative to the linear fit), weighted mean per sector and eclipse type
rlin = tt - Al @ c
sect = {}
for k_ in np.unique(kk):
    for p_ in (0.0, 0.5):
        m = (kk == k_) & (pp == p_)
        if m.sum() == 0: continue
        w = 1 / ss[m] ** 2; mu = (rlin[m] * w).sum() / w.sum(); er = 1 / np.sqrt(w.sum())
        chi_in = float((((rlin[m] - mu) / ss[m]) ** 2).sum())
        sect[f'{k_}_ph{p_}'] = dict(n=int(m.sum()), oc_min=float(mu * 1440), err_min=float(er * 1440), chi2_within=chi_in, mean_epoch=float(nn[m].mean()),
                                   err_scatter_min=float(np.std(rlin[m]) / np.sqrt(m.sum()) * 1440))
    m = kk == k_
    w = 1 / ss[m] ** 2; mu = (rlin[m] * w).sum() / w.sum()
    sect[f'{k_}_both'] = dict(n=int(m.sum()), oc_min=float(mu * 1440), err_min=float(1 / np.sqrt(w.sum()) * 1440), mean_epoch=float(nn[m].mean()))
res['per_sector'] = sect
res['events_used'] = [dict(key=str(a), n=int(b), ph=float(c_), oc_min=float(d_ * 1440), sig_min=float(e_ * 1440)) for a, b, c_, d_, e_ in zip(kk, nn, pp, rlin, ss)]
json.dump(res, open(f'{V}/results/etvfit_{label}.json', 'w'), indent=1)
print(f"events {n.size}, rejected {nrej}")
for k_, f_ in fits.items():
    print(f"{k_:10s} chi2 {f_['chi2']:.1f} dof {f_['dof']} BIC {f_['bic']:.1f}")
e_ = res['ephemeris']; q_ = res['quadratic']
print(f"P = {e_['P_d']:.8f} +- {e_['P_err_s']:.3f} s; T0 = {e_['T0_btjd']:.6f} +- {e_['T0_err_min']:.2f} min; sec - P/2 = {e_['sec_minus_halfP_min']:+.2f} +- {e_['sec_minus_halfP_err_min']:.2f} min")
print(f"quadratic: c = {q_['c_min_per_yr2']:+.4f} +- {q_['c_err_min_per_yr2']:.4f} min/yr^2; dP/dt = {q_['dPdt_s_per_yr']:+.4f} +- {q_['dPdt_err_s_per_yr']:.4f} s/yr; "
      f"dchi2 {q_['dchi2_lin_minus_quad']:.2f} null p {q_['null_p_dchi2']:.3f}; null Q sd / formal {q_['null_Q_sd_over_formal']:.2f}")
print(f"linear chi2 {lin['chi2']:.1f} (null median {res['lin_chi2_null_median']:.1f}, p {res['lin_chi2_null_p']:.3f})")
for k_, v in sect.items():
    print(f"  {k_:14s} N={v['n']:3d} O-C {v['oc_min']:+6.2f} +- {v['err_min']:.2f} min" + (f" (scatter-based {v['err_scatter_min']:.2f}; chi2_within {v['chi2_within']:.1f})" if 'chi2_within' in v else ''))
