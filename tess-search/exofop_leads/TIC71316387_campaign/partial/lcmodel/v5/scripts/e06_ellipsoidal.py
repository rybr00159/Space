"""Ellipsoidal-variation limit and what it implies for f_EB and the mass ratio.
(1) Combined ellipsoidal amplitude (-a2, cos 4 pi phi with minima at the eclipses) from the per-sector fits of the
    default product set (e01_ooe.py), weighted by the fake-period null sd and inflated by sqrt(chi2/dof) of the
    sector-to-sector scatter (spot-induced second harmonics); 95% upper limit.
(2) Predicted diluted amplitude for each light-curve solution along the L3 profile and for the MS posterior:
    A = f_EB * sum_j alpha_j q_j (R_j/a)^3 sin^2 i (L_j/L_EB),  alpha = 0.15 (15+u)(1+tau)/(3-u)
    (Morris & Naftilan 1993), u = 0.5, tau = 0.4 (cool dwarfs, TESS band; alpha = 1.30; tau = 0.2-0.6 changes alpha
    by -14/+14%), q_1 = M2/M1, q_2 = M1/M2.
Usage: e06_ellipsoidal.py [ms chain tag]"""
import numpy as np, sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v5/lib')
V = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v5'
O = json.load(open(f'{V}/results/ooe_ks3.json'))
labs = ['SPOC-PDC1', 'SPOC-PDC28', 'SPOC-PDC68', 'TC95', 'TC105', 'TC106']
sel = [o for o in O['per_product'] if o['label'] in labs]
ell = np.array([-o['a2'] / o['inj_recovery_ell'] for o in sel]); err = np.array([o['null_sd_a2'] for o in sel])
w = 1 / err ** 2; mu = (w * ell).sum() / w.sum(); e = 1 / np.sqrt(w.sum()); chi = ((ell - mu) ** 2 * w).sum(); dof = len(ell) - 1
e_infl = e * max(1, np.sqrt(chi / dof))
ul95 = max(mu, 0) + 1.96 * e_infl
res = dict(per_sector_ell_ppm=dict(zip(labs, (ell * 1e6).round(1).tolist())), per_sector_err_ppm=(err * 1e6).round(1).tolist(),
           ell_weighted_ppm=mu * 1e6, ell_err_formal_ppm=e * 1e6, chi2=chi, dof=dof, ell_err_inflated_ppm=e_infl * 1e6, ell_UL95_ppm=ul95 * 1e6,
           ell_unweighted_mean_ppm=float(ell.mean() * 1e6), ell_unweighted_sem_ppm=float(ell.std(ddof=1) / np.sqrt(len(ell)) * 1e6))
print(f"ellipsoidal (diluted): weighted {mu*1e6:+.1f} +- {e*1e6:.1f} ppm (chi2 {chi:.1f}/{dof}) -> inflated +- {e_infl*1e6:.1f}; 95% UL {ul95*1e6:.0f} ppm; "
      f"unweighted {ell.mean()*1e6:+.1f} +- {ell.std(ddof=1)/np.sqrt(len(ell))*1e6:.1f}")
# first harmonic (spot) amplitude per sector, for the record
res['A1_ppm'] = {o['label']: round(o['A1'] * 1e6, 1) for o in sel}; res['A1_null95_ppm'] = {o['label']: round(o['null_A1_95'] * 1e6, 1) for o in sel}
ALPHA = 0.15 * (15 + 0.5) * (1 + 0.4) / (3 - 0.5)
def a_ell(fEB, rsum, k, b0, J, q):
    r1 = rsum / (1 + k); r2 = rsum * k / (1 + k); sin2i = 1 - (b0 * rsum) ** 2
    L2 = k * k * J / (1 + k * k * J); L1 = 1 - L2
    return fEB * ALPHA * (q * r1 ** 3 * L1 + r2 ** 3 * L2 / q) * sin2i
prof = json.load(open(f'{V}/results/profile_L3_default.json'))['profile']
chi0 = min(r['chi2'] for r in prof)
rows = []
for r in prof:
    p = r['p']; fe = 1 - p['L3']
    pred = {f'q={q_}': a_ell(fe, p['rsum'], p['k'], p['b0'], p['J'], q_) * 1e6 for q_ in (0.05, 0.1, 0.3, 0.5, 1.0)}
    qmax = None
    for q_ in np.linspace(0.02, 1.0, 491):
        if a_ell(fe, p['rsum'], p['k'], p['b0'], p['J'], q_) > ul95: qmax = q_; break
    rows.append(dict(L3=p['L3'], k=p['k'], dchi2=r['chi2'] - chi0, pred_ppm=pred, q_above_which_excluded=qmax))
    print(f"L3={p['L3']:.3f} k={p['k']:.3f} dchi2={r['chi2']-chi0:5.1f}  predicted A_ell (ppm) " + ' '.join(f"{k_} {v:6.1f}" for k_, v in pred.items()) +
          f"  -> excluded for q > {qmax if qmax is None else round(qmax, 3)}")
res['profile_predictions'] = rows
if len(sys.argv) > 1:
    tag = sys.argv[1]
    import emcee, msfit, ms
    meta = json.load(open(f'{V}/results/chain_{tag}_meta.json'))
    rd = emcee.backends.HDFBackend(f'{V}/results/chain_{tag}.h5', read_only=True)
    fl = rd.get_chain(discard=rd.iteration // 2, thin=20, flat=True)
    S = {n: fl[:, i] for i, n in enumerate(meta['free'])}
    A = []
    for j in range(fl.shape[0]):
        rs, kk, JJ, a = msfit.geom(S['M1'][j], S['q'][j], S['dR1'][j], S['dR2'][j], S['dJ'][j], S['P'][j])
        L3 = S['L3'][j] if 'L3' in S else np.nan
        A.append(a_ell(1 - L3, rs, kk, S['b0'][j], JJ, S['q'][j]))
    A = np.array(A) * 1e6
    res['ms_chain_pred_ppm'] = dict(tag=tag, quantiles_2p5_50_97p5=np.percentile(A, [2.5, 50, 97.5]).tolist(), frac_above_UL=float(np.mean(A > ul95 * 1e6)))
    print('MS chain predicted A_ell (ppm) 2.5/50/97.5%:', np.round(np.percentile(A, [2.5, 50, 97.5]), 1), 'fraction above UL', np.mean(A > ul95 * 1e6))
json.dump(res, open(f'{V}/results/ellipsoidal.json', 'w'), indent=1)
