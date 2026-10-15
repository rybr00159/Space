"""Posterior summaries and derived quantities from an emcee chain (geometric or MS-constrained).
Geometric chain: r1, r2, i, f_EB = 1 - L3, undiluted eclipse depths, T2/T1 (blackbody, TESS band) and, with a
mass-ratio prior q ~ U(0.2, 1) [or q from the MS surface-brightness relation], rho1 and the MS mass of star 1.
MS chain: M1, M2, R1, R2, Teff1/2, f_EB, the photometric distance of the pair implied by f_EB, and the ratio of
the bound-hypothesis light fraction to the fitted one.
Usage: d01_derive.py <chain tag> [burn fraction] [thin]"""
import numpy as np, sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6/lib')
import emcee, ebm, ms, derive, msfit
V = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v6'
tag = sys.argv[1]; burn = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5; thin = int(sys.argv[3]) if len(sys.argv) > 3 else 10
meta = json.load(open(f'{V}/results/chain_{tag}_meta.json'))
rd = emcee.backends.HDFBackend(f'{V}/results/chain_{tag}.h5', read_only=True)
nst = rd.iteration
ch = rd.get_chain(discard=int(burn * nst), thin=thin)        # (steps, walkers, dim)
lp = rd.get_log_prob(discard=int(burn * nst), thin=thin)
try:
    tau = rd.get_autocorr_time(discard=int(burn * nst), tol=0)
except Exception:
    tau = np.full(ch.shape[2], np.nan)
flat = ch.reshape(-1, ch.shape[2]); lpf = lp.reshape(-1)
free = meta['free']; fixed = meta['fixed']
S = {n: flat[:, i] for i, n in enumerate(free)}
for k, v in fixed.items():
    if k not in S: S[k] = np.full(flat.shape[0], v)
rng = np.random.default_rng(4)
N = flat.shape[0]
def q(x, p=(2.5, 16, 50, 84, 97.5)):
    return [float(v) for v in np.percentile(x, p)]
out = dict(tag=tag, nsteps=int(nst), nwalkers=int(ch.shape[1]), n_samples=int(N), burn=burn, thin=thin,
           tau=dict(zip(free, [float(t) for t in tau])), max_tau=float(np.nanmax(tau)), steps_over_tau=float(nst / np.nanmax(tau)),
           acceptance=float(np.mean(rd.accepted / nst)), best_lp=float(lpf.max()))
if 'M1' in S:   # MS-constrained chain
    G = {k: np.empty(N) for k in ('rsum', 'k', 'J', 'a_rsun', 'L3')}
    for j in range(N):
        rs, kk, JJ, a = msfit.geom(S['M1'][j], S['q'][j], S['dR1'][j], S['dR2'][j], S['dJ'][j], S['P'][j])
        G['rsum'][j], G['k'][j], G['J'][j], G['a_rsun'][j] = rs, kk, JJ, a
        if meta['mode'] == 'bound':
            G['L3'][j] = 1 - derive.fEB_if_bound(S['M1'][j], S['q'][j] * S['M1'][j]) * 10 ** (-0.4 * S['dmag'][j])
        else:
            G['L3'][j] = S['L3'][j]
    S.update(G)
    M1 = S['M1']; M2 = S['q'] * M1
    R1 = ms.radius(M1) * (1 + S['dR1']); R2 = ms.radius(M2) * (1 + S['dR2'])
    S['M2'] = M2; S['R1'] = R1; S['R2'] = R2; S['Teff1'] = ms.teff(M1); S['Teff2'] = ms.teff(M2)
    S['T2/T1_ms'] = S['Teff2'] / S['Teff1']
    S['T2/T1_bb'] = derive.T2_from_J_bb(S['J'], S['Teff1']) / S['Teff1']
    fEB = 1 - S['L3']; S['f_EB'] = fEB
    S['fEB_bound_pred'] = derive.fEB_if_bound(M1, M2)
    S['ratio_fEB_fit_over_bound'] = fEB / S['fEB_bound_pred']
    S['d_EB_pc'] = derive.distance_from_fEB(fEB, M1, M2)
    S['MT_pair'] = derive.pair_MT(M1, M2)
    S['dmag_implied'] = -2.5 * np.log10(fEB / S['fEB_bound_pred'])   # pair brighter (<0) / fainter (>0) than a MS pair at 297 pc
    S['MKs_pair'] = -2.5 * np.log10(10 ** (-0.4 * ms.MKs(M1)) + 10 ** (-0.4 * ms.MKs(M2)))
    S['dKs_pair_minus_Gstar_if_bound'] = S['MKs_pair'] + derive.MU_G - 10.341   # Ks contrast of the pair vs the blend's Ks = 10.341
    S['dT_pair_minus_blend'] = -2.5 * np.log10(fEB)
    S['M1_plus_M2'] = M1 + M2
    S['a_au'] = S['a_rsun'] * 0.00465047
    S['rho1'] = M1 / R1 ** 3
    out['spt1_median'] = str(ms.spt(float(np.median(M1)))); out['spt2_median'] = str(ms.spt(float(np.median(M2))))
else:
    S['f_EB'] = 1 - S['L3']
qq = rng.uniform(0.2, 1.0, N)
r1 = S['rsum'] / (1 + S['k']); r2 = S['rsum'] * S['k'] / (1 + S['k'])
S['r1'] = r1; S['r2'] = r2; S['inc_deg'] = np.degrees(np.arccos(S['b0'] * S['rsum']))
S['b_prim_R1'] = S['b0'] * S['rsum'] / r1 * (1 - (S['ecosw'] ** 2 + S['esinw'] ** 2)) / (1 + S['esinw'])
S['e'] = np.hypot(S['ecosw'], S['esinw'])
if 'M1' not in S:
    S['rho1_qU0.2-1'] = derive.rho1(r1, qq)
    S['M1_ms_from_rho1'] = ms.mass_from_density(S['rho1_qU0.2-1'])
    # self-consistent q: MS secondary mass from J and the MS primary mass from density (2 iterations)
    m1 = ms.mass_from_density(derive.rho1(r1, 0.6))
    for it in range(3):
        m2 = ms.m2_from_J(S['J'], np.nan_to_num(m1, nan=0.5)); qs = np.clip(m2 / m1, 0.05, 1.0)
        m1 = ms.mass_from_density(derive.rho1(r1, qs))
    S['M1_ms_selfconsistent'] = m1; S['q_ms_selfconsistent'] = qs
    S['k_ms_selfconsistent'] = ms.radius(qs * m1) / ms.radius(m1)
    S['T2/T1_bb_T1=3900'] = derive.T2_from_J_bb(S['J'], 3900.0) / 3900.0
# undiluted depths (EB alone), by evaluating the model at conjunctions on a thinned subsample
gm = ebm.GridModel(1.3036722)
sub = rng.choice(N, min(N, 1500), replace=False)
dp, ds = [], []
for j in sub:
    dsec = ebm.secondary_offset(S['P'][j], S['ecosw'][j], S['esinw'][j])
    tt = np.array([0.0, dsec]) + np.linspace(-0.004, 0.004, 9)[:, None]
    f = ebm.flux_grid(tt.ravel(), S['P'][j], S['rsum'][j], S['k'][j], S['b0'][j], S['J'][j], S['ecosw'][j], S['esinw'][j],
                      S['u1a'][j], S['usa'][j] - S['u1a'][j], S['u1b'][j], S['usb'][j] - S['u1b'][j], ebm.XG, ebm.WG).reshape(9, 2)
    dp.append(1 - f[:, 0].min()); ds.append(1 - f[:, 1].min())
dp = np.array(dp); ds = np.array(ds)
out['undiluted_depth_primary'] = q(dp); out['undiluted_depth_secondary'] = q(ds)
out['diluted_depth_primary_ppm'] = q(dp * (1 - S['L3'][sub]) * 1e6); out['diluted_depth_secondary_ppm'] = q(ds * (1 - S['L3'][sub]) * 1e6)
out['depth_ratio_sec_over_prim'] = q(ds / dp)
out['quantiles_2.5_16_50_84_97.5'] = {k: q(v[np.isfinite(v)]) for k, v in S.items() if np.ndim(v) == 1 and np.isfinite(v).sum() > 10}
# correlations among the key geometric parameters
keys = [k for k in ('L3', 'k', 'b0', 'rsum', 'J', 'esinw', 'u1a', 'u1b') if k in S]
C = np.corrcoef(np.vstack([S[k] for k in keys]))
out['corr'] = dict(keys=keys, matrix=np.round(C, 3).tolist())
if 'M1' in S:
    for thr in (0.5, 1.0, 2.0):
        out[f'P(|log10 ratio fEB fit/bound| < {thr} dex)'] = float(np.mean(np.abs(np.log10(S['ratio_fEB_fit_over_bound'])) < thr))
    out['P(d_EB within 297 pc x/÷1.5)'] = float(np.mean(np.abs(np.log(S['d_EB_pc'] / 297.4)) < np.log(1.5)))
json.dump(out, open(f'{V}/results/derived_{tag}.json', 'w'), indent=1)
qd = out['quantiles_2.5_16_50_84_97.5']
print(f"{tag}: steps {nst}, walkers {ch.shape[1]}, max tau {out['max_tau']:.0f} (steps/tau {out['steps_over_tau']:.1f}), acc {out['acceptance']:.2f}")
for k in ('L3', 'f_EB', 'k', 'b0', 'rsum', 'r1', 'r2', 'inc_deg', 'J', 'esinw', 'ecosw', 'u1a', 'u1b', 'M1', 'q', 'M2', 'R1', 'R2', 'Teff1', 'Teff2', 'T2/T1_ms',
          'T2/T1_bb', 'dR1', 'dR2', 'dJ', 'rho1', 'rho1_qU0.2-1', 'M1_ms_from_rho1', 'M1_ms_selfconsistent', 'q_ms_selfconsistent', 'k_ms_selfconsistent', 'T2/T1_bb_T1=3900',
          'fEB_bound_pred', 'ratio_fEB_fit_over_bound', 'd_EB_pc', 'MT_pair', 'dmag', 'dmag_implied', 'MKs_pair', 'dKs_pair_minus_Gstar_if_bound',
          'dT_pair_minus_blend', 'M1_plus_M2', 'a_rsun', 'e', 'b_prim_R1'):
    if k in qd: print(f"  {k:24s} " + '  '.join(f'{v:.4g}' for v in qd[k]))
print('  undiluted depth prim', [round(v, 4) for v in out['undiluted_depth_primary']], 'sec', [round(v, 4) for v in out['undiluted_depth_secondary']])
print('  diluted depth prim ppm', [round(v) for v in out['diluted_depth_primary_ppm']], 'sec', [round(v) for v in out['diluted_depth_secondary_ppm']])
print('  corr', keys); print(np.round(C, 2))
for k, v in out.items():
    if k.startswith('P('): print(' ', k, round(v, 3))
