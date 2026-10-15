"""Ellipsoidal-variation limit and what it implies for f_EB and the mass ratio (v6).
(1) Per-sector ellipsoidal amplitude ell = -a2 (cos 4 pi phi, minima at the eclipses) from e01_ooe.py for the default
    product set (SPOC PDCSAP S1/S28/S68, own TESScut S95/S105/S106), each with its empirical (fake-period) noise sd.
    The sector-to-sector scatter exceeds the noise (spots evolve; their second harmonic leaks into cos 4 pi phi), so
    the sectors are combined with a DerSimonian-Laird random-effects model (between-sector variance tau^2 estimated
    from the data).  One-sided 95% upper limit UL = max(mu, 0) + t_0.95(k-1) sigma_RE (Student t, 5 dof: 2.015).  (Fixed-effect and unweighted
    combinations are reported for comparison.)
(2) Predicted DILUTED ellipsoidal semi-amplitude for a light-curve solution:
      A = f_EB * alpha * sin^2 i * [ q r1^3 L1 + r2^3 L2 / q ],   r_j = R_j / a,  q = M2/M1,
    alpha = 0.15 (15 + u)(1 + tau) / (3 - u) (Morris & Naftilan 1993) with u = 0.5, tau = 0.4 (cool dwarfs, TESS band;
    alpha = 1.30; tau in 0.2-0.6 changes alpha by +-14%, so A is quoted with alpha x 1.14 as the conservative upper).
    For each solution along the geometric L3 profile: the largest q allowed by the UL (with alpha x 0.86, i.e. the
    most permissive), compared with the mass ratio a main-sequence secondary needs to have the fitted surface
    brightness ratio J (ms.m2_from_J with M1 from the eclipse density at that q, solved self-consistently).
(3) For MS-constrained chains: predicted A distribution and the fraction above the UL.
Usage: e06_ellipsoidal.py [ms chain tag ...]"""
import numpy as np, sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6/lib')
import ms, derive
V = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6'
O = json.load(open(f'{V}/results/ooe_ks3.json'))
labs = ['SPOC-PDC1', 'SPOC-PDC28', 'SPOC-PDC68', 'TC95', 'TC105', 'TC106']
sel = [next(o for o in O['per_product'] if o['label'] == l) for l in labs]
ell = np.array([-o['a2'] / o['inj_recovery_ell'] for o in sel]); err = np.array([o['null_sd_a2'] for o in sel])
w = 1 / err ** 2; mu_fe = (w * ell).sum() / w.sum(); e_fe = 1 / np.sqrt(w.sum()); Q = ((ell - mu_fe) ** 2 * w).sum(); k = len(ell)
tau2 = max(0.0, (Q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum()))
wr = 1 / (err ** 2 + tau2); mu_re = (wr * ell).sum() / wr.sum(); e_re = 1 / np.sqrt(wr.sum())
from scipy.stats import t as student_t
UL_norm = max(mu_re, 0) + 1.645 * e_re
UL = max(mu_re, 0) + student_t.ppf(0.95, k - 1) * e_re   # adopted: Student-t (k-1 dof) because tau is estimated from only k sectors
res = dict(per_sector_ell_ppm=dict(zip(labs, (ell * 1e6).round(1).tolist())), per_sector_err_ppm=dict(zip(labs, (err * 1e6).round(1).tolist())),
           fixed_effect_ppm=[mu_fe * 1e6, e_fe * 1e6], Q=Q, dof=k - 1, tau_ppm=np.sqrt(tau2) * 1e6,
           random_effects_ppm=[mu_re * 1e6, e_re * 1e6], UL95_onesided_ppm=UL * 1e6, UL95_onesided_normal_ppm=UL_norm * 1e6,
           unweighted_ppm=[float(ell.mean() * 1e6), float(ell.std(ddof=1) / np.sqrt(k) * 1e6)])
print(f"ellipsoidal (diluted, ppm): per sector {np.round(ell*1e6,0)} +- {np.round(err*1e6,0)}")
print(f"  fixed effect {mu_fe*1e6:+.1f} +- {e_fe*1e6:.1f} (Q = {Q:.1f} for {k-1} dof); random effects tau = {np.sqrt(tau2)*1e6:.0f} ppm -> {mu_re*1e6:+.1f} +- {e_re*1e6:.1f}; "
      f"one-sided 95% UL {UL*1e6:.0f} ppm (Student-t; normal: {UL_norm*1e6:.0f}); unweighted {ell.mean()*1e6:+.1f} +- {ell.std(ddof=1)/np.sqrt(k)*1e6:.1f}")
res['A1_ppm'] = {o['label']: round(o['A1'] * 1e6, 1) for o in sel}; res['A1_null95_ppm'] = {o['label']: round(o['null_A1_95'] * 1e6, 1) for o in sel}
ALPHA = 0.15 * (15 + 0.5) * (1 + 0.4) / (3 - 0.5)
def a_ell(fEB, rsum, k_, b0, J, q, alpha=ALPHA):
    r1 = rsum / (1 + k_); r2 = rsum * k_ / (1 + k_); sin2i = 1 - (b0 * rsum) ** 2
    L2 = k_ * k_ * J / (1 + k_ * k_ * J); L1 = 1 - L2
    return fEB * alpha * (q * r1 ** 3 * L1 + r2 ** 3 * L2 / q) * sin2i
def q_ms_selfconsistent(rsum, k_, J, P=1.3036712):
    """MS mass ratio implied by J with M1 from the eclipse density rho1(r1, q) (fixed-point iteration)"""
    r1 = rsum / (1 + k_); qq = 0.5; m1 = np.nan
    for _ in range(30):
        m1 = ms.mass_from_density(derive.rho1(r1, qq))
        if not np.isfinite(m1): return np.nan, np.nan
        m2 = ms.m2_from_J(J, m1)
        if not np.isfinite(m2): return np.nan, m1
        qq = 0.5 * qq + 0.5 * min(1.0, m2 / m1)
    return qq, m1
prof = json.load(open(f'{V}/results/profile_L3_default.json'))['profile']
chi0 = min(r['chi2'] for r in prof)
rows = []
qs = np.linspace(0.01, 1.0, 991)
for r in prof:
    p = r['p']; fe = 1 - p['L3']
    A_q = a_ell(fe, p['rsum'], p['k'], p['b0'], p['J'], qs, ALPHA * 0.86)
    okq = qs[A_q <= UL]
    qmax = float(okq.max()) if okq.size else 0.0
    qms, m1 = q_ms_selfconsistent(p['rsum'], p['k'], p['J'])
    pred_ms = float(a_ell(fe, p['rsum'], p['k'], p['b0'], p['J'], qms) * 1e6) if np.isfinite(qms) else np.nan
    rows.append(dict(L3=p['L3'], k=p['k'], b0=p['b0'], rsum=p['rsum'], J=p['J'], dchi2=r['chi2'] - chi0, q_max_allowed=qmax, q_ms_from_J=qms,
                     M1_ms_from_density=m1, A_pred_at_q_ms_ppm=pred_ms, excluded_if_MS=bool(np.isfinite(qms) and qms > qmax)))
    print(f"L3={p['L3']:.3f} k={p['k']:.3f} b0={p['b0']:.2f} rsum={p['rsum']:.3f} J={p['J']:.2f} dchi2={r['chi2']-chi0:6.1f} | ellipsoidal UL allows q < {qmax:.3f}; "
          f"MS secondary with this J needs q = {qms:.2f} (M1 = {m1:.2f}) -> A_pred {pred_ms:6.1f} ppm {'EXCLUDED' if rows[-1]['excluded_if_MS'] else 'ok'}")
res['profile_predictions'] = rows
for tag in sys.argv[1:]:
    import emcee, msfit
    meta = json.load(open(f'{V}/results/chain_{tag}_meta.json'))
    rd = emcee.backends.HDFBackend(f'{V}/results/chain_{tag}.h5', read_only=True)
    fl = rd.get_chain(discard=rd.iteration // 2, thin=20, flat=True)
    S = {n: fl[:, i] for i, n in enumerate(meta['free'])}
    for kk, v in meta['fixed'].items():
        if kk not in S: S[kk] = np.full(fl.shape[0], v)
    A = []
    for j in range(fl.shape[0]):
        rs, kk_, JJ, a = msfit.geom(S['M1'][j], S['q'][j], S['dR1'][j], S['dR2'][j], S['dJ'][j], S['P'][j])
        if 'L3' in meta['free']:
            fe = 1 - S['L3'][j]
        else:
            fe = derive.fEB_if_bound(S['M1'][j], S['q'][j] * S['M1'][j]) * 10 ** (-0.4 * S['dmag'][j])
        A.append(a_ell(fe, rs, kk_, S['b0'][j], JJ, S['q'][j]))
    A = np.array(A) * 1e6
    res[f'chain_{tag}_pred_ppm'] = dict(quantiles_2p5_50_97p5=np.percentile(A, [2.5, 50, 97.5]).tolist(), frac_above_UL=float(np.mean(A > UL * 1e6)))
    print(f'chain {tag}: predicted diluted A_ell (ppm) 2.5/50/97.5%: {np.round(np.percentile(A, [2.5, 50, 97.5]), 1)}, fraction above UL {np.mean(A > UL*1e6):.3f}')
json.dump(res, open(f'{V}/results/ellipsoidal.json', 'w'), indent=1)
