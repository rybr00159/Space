"""Posterior summaries and derived quantities from an emcee chain (geo or ms).
Distance-independent: r1 = R1/a, r2, b1 = cos i / r1, i, f_EB = 1 - L3, undiluted eclipse depths (model with L3 = 0, SPOC
scale), J, T2/T1 (ms: from the MS relations; geo: from J by inverting the smooth MS surface-brightness relation at the
primary mass implied by rho1, and by a blackbody in a 600-1000 nm top-hat band as a cross-check),
rho1 = 3 pi/(G P^2) r1^-3 / (1+q): ms -> q from the model; geo -> two mass-ratio priors (q implied by k on the MS; q ~ U(0.1,1)).
Bound test (ms only; uses the Gaia distance only at this final step): m_T(EB) = Tmag - 2.5 log10 f_EB versus M_T(pair) from
the smooth MS relations (+0.10 mag scatter, A_T = 0.02 +- 0.02):  D_EB / D_Gaia, and f_EB predicted if bound.
Usage: s12_derived.py <chain tag> [burn fraction] [thin]"""
import numpy as np, sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4/lib')
import fit3, lclib as L, msfit as M, ms_smooth as MS
from ebmodel import eb_flux
FD = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4'
tag = sys.argv[1]; burn = float(sys.argv[2]) if len(sys.argv) > 2 else 0.4; thin = int(sys.argv[3]) if len(sys.argv) > 3 else 10
z = np.load(f'{FD}/results/chain_{tag}.npz', allow_pickle=True)
ch = z['chain']; lnp = z['lnp']; names = list(z['names'])
b0 = int(burn * ch.shape[0])
X = ch[b0::thin].reshape(-1, ch.shape[2]); LP = lnp[b0::thin].reshape(-1)
ok = np.isfinite(LP); X, LP = X[ok], LP[ok]
S = {n: X[:, i] for i, n in enumerate(names)}
kind = 'geo' if 'rsum' in names else 'ms'
G = 6.674e-8; Psec = L.P0 * 86400; RHO_SUN = 1.408; KR = 3 * np.pi / (G * Psec ** 2)
rng = np.random.default_rng(1)
N = X.shape[0]
if kind == 'ms':
    M1 = S['M1']; q = S['q']; M2 = M1 * q
    R1 = MS.radius(M1) * S['fR']; R2 = MS.radius(M2) * S['fR']
    a = M.a_over_R(M1 + M2, L.P0); rsum = (R1 + R2) / a; k = R2 / R1
    J = MS.SB_T(M2) / MS.SB_T(M1) * np.exp(S['dlnJ'])
else:
    rsum, k, J = S['rsum'], S['k'], S['J']
r1 = rsum / (1 + k); r2 = rsum - r1; cosi = S['bn'] * rsum; inc = np.degrees(np.arccos(cosi)); b1 = cosi / r1
fEB = 1 - S['L3']
D = dict(L3=S['L3'], fEB=fEB, r1=r1, r2=r2, k=k, J=J, rsum=rsum, b1=b1, inc=inc, ecosw=S['ecosw'], esinw=S['esinw'],
         e=np.hypot(S['ecosw'], S['esinw']), s95=S['s95'], s105=S['s105'], s106=S['s106'],
         total_frac=(S['bn'] <= np.abs(1 - k) / (1 + k)).astype(float))
def mass_of_rho_vec(rho, fR=1.04):
    mg = np.exp(np.linspace(np.log(0.09), np.log(1.25), 3000)); rg = mg / (fR * MS.radius(mg)) ** 3
    return np.interp(rho, rg[::-1], mg[::-1])
if kind == 'ms':
    rho1 = KR / r1 ** 3 / (1 + q) / RHO_SUN
    T1 = MS.teff(M1); T2 = MS.teff(M2)
    D.update(M1=M1, M2=M2, q=q, R1=R1, R2=R2, a_Rsun=a, fR=S['fR'], dlnJ=S['dlnJ'], rho1_sun=rho1, rho1_check_sun=M1 / R1 ** 3,
             T1=T1, T2=T2, T2T1=T2 / T1)
else:
    # q implied by k on the MS (iterate with rho), and a flat prior
    qk = np.clip(k, 0.1, 1.0)
    for _ in range(30):
        m1 = mass_of_rho_vec(KR / r1 ** 3 / (1 + qk) / RHO_SUN)
        qn = np.array([np.interp(kk, MS.radius(np.linspace(0.08, m, 300)) / MS.radius(m), np.linspace(0.08, m, 300) / m) if kk < 1 else 1.0
                       for kk, m in zip(k, m1)])
        qk = 0.5 * (qk + qn)
    qu = rng.uniform(0.1, 1.0, N)
    rho_k = KR / r1 ** 3 / (1 + qk) / RHO_SUN; rho_u = KR / r1 ** 3 / (1 + qu) / RHO_SUN
    m1k = mass_of_rho_vec(rho_k); m1u = mass_of_rho_vec(rho_u)
    # T2/T1 from J: SB_T(M2) = J SB_T(M1) at M1 = m1k
    mg = np.exp(np.linspace(np.log(0.08), np.log(1.25), 3000)); sb = MS.SB_T(mg)
    m2J = np.interp(J * MS.SB_T(m1k), sb, mg)
    T2T1_ms = MS.teff(m2J) / MS.teff(m1k)
    # blackbody, 600-1000 nm top hat
    lam = np.linspace(600e-9, 1000e-9, 200); hc_k = 1.4388e-2
    def bb(T): return np.trapezoid(1 / lam ** 5 / (np.exp(hc_k / (lam * T[:, None])) - 1), lam, axis=1)
    Tg = np.linspace(2000, 7000, 501); bg = bb(Tg)
    T1g = MS.teff(m1k); T2bb = np.interp(J * np.interp(T1g, Tg, bg), bg, Tg)
    D.update(q_MSk=qk, rho1_sun_qMSk=rho_k, rho1_sun_qflat=rho_u, M1_from_rho_qMSk=m1k, M1_from_rho_qflat=m1u,
             T2T1_from_J_MS=T2T1_ms, T2T1_from_J_bb=T2bb / T1g, M2_from_J=m2J)
# undiluted and diluted depths, T14, for 300 random samples
idx = rng.choice(N, min(300, N), replace=False)
dp, ds, dpu, dsu, t14, t23 = [], [], [], [], [], []
for i in idx:
    u1a, u2a = fit3.q2u(S['q1a'][i], S['q2a'][i]); u1b, u2b = fit3.q2u(S['q1b'][i], S['q2b'][i])
    tt = np.linspace(-0.07, 0.07, 2801); sh = L.P0 / 2 + 2 * L.P0 / np.pi * S['ecosw'][i]
    for L3v, store in ((S['L3'][i], (dp, ds)), (0.0, (dpu, dsu))):
        args = (rsum[i], k[i], cosi[i], J[i], L3v, S['ecosw'][i], S['esinw'][i], u1a, u2a, u1b, u2b, 1e-6, 1, fit3.XG, fit3.WG)
        fp = eb_flux(tt, 0.0, L.P0, *args); fs = eb_flux(tt + sh, 0.0, L.P0, *args)
        store[0].append(1 - fp.min()); store[1].append(1 - fs.min())
    w = tt[fp < 1 - 1e-9]; t14.append((w.max() - w.min()) * 24)
    fl = fs.min(); w2 = tt[fs < fl + 1e-7 * (1 - fl)]; t23.append((w2.max() - w2.min()) * 24 if w2.size > 1 else 0.0)
D2 = dict(depth_prim_ppm=np.array(dp) * 1e6, depth_sec_ppm=np.array(ds) * 1e6, depth_prim_undiluted=np.array(dpu), depth_sec_undiluted=np.array(dsu),
          depth_ratio=np.array(ds) / np.array(dp), T14_h=np.array(t14), T23_sec_flat_h=np.array(t23))
if kind == 'geo':
    # distance-independent masses for the bound test: M1 from rho1 (q from k on the MS), M2 = q M1
    M1 = np.clip(m1k, 0.09, 1.25); M2 = np.clip(qk * m1k, 0.08, None)
    D.update(M1=M1, M2=M2)
if True:
    MTEB = -2.5 * np.log10(10 ** (-0.4 * MS.MTess(M1)) + 10 ** (-0.4 * MS.MTess(M2))) + rng.normal(0, 0.10, N)
    Tmag = 11.358 + rng.normal(0, 0.006, N); mTEB = Tmag - 2.5 * np.log10(fEB)
    AT = np.clip(rng.normal(0.02, 0.02, N), 0, None)
    plx = rng.normal(3.3622 + 0.025, np.hypot(0.017, 0.015), N); DG = 1000 / plx
    DEB = 10 ** ((mTEB - MTEB - AT + 5) / 5)
    MTtot = Tmag - 5 * np.log10(DG / 10) - AT
    fpred = 10 ** (-0.4 * (MTEB - MTtot))
    MK = -2.5 * np.log10(10 ** (-0.4 * MS.MKs(M1)) + 10 ** (-0.4 * MS.MKs(M2))); KEB = MK + 5 * np.log10(DG / 10)
    KG = -2.5 * np.log10(np.clip(10 ** (-0.4 * 10.341) - 10 ** (-0.4 * KEB), 1e-30, None))
    MGp = -2.5 * np.log10(10 ** (-0.4 * MS.MG(M1)) + 10 ** (-0.4 * MS.MG(M2))); GEB = MGp + 5 * np.log10(DG / 10) + 1.4 * AT
    GG = -2.5 * np.log10(np.clip(10 ** (-0.4 * 11.775) - 10 ** (-0.4 * GEB), 1e-30, None))
    TG = -2.5 * np.log10(np.clip(10 ** (-0.4 * Tmag) - 10 ** (-0.4 * mTEB), 1e-30, None))
    D.update(MT_EB=MTEB, mT_EB=mTEB, D_EB_pc=DEB, D_ratio=DEB / DG, fEB_pred_bound=fpred, dT_EB_minus_G=mTEB - TG,
             Ks_EB_bound=KEB, dKs_bound=KEB - KG, G_EB_bound=GEB, dG_bound=GEB - GG)
# predicted ellipsoidal semi-amplitude (total light), Morris & Naftilan leading term, alpha = 1.30 (u=0.5, tau=0.4)
ALPHA = 0.15 * (15 + 0.5) * (1 + 0.4) / (3 - 0.5)
qq = D['q'] if 'q' in D else D['q_MSk']
k2J = k ** 2 * J; L2f = k2J / (1 + k2J); L1f = 1 - L2f
D['Aell_pred_ppm'] = 1e6 * fEB * ALPHA * (1 - cosi ** 2) * (L1f * qq * r1 ** 3 + L2f * r2 ** 3 / qq)
ooe = json.load(open(f'{FD}/results/ooe.json'))['default']
# importance weights from the OOE ellipsoidal measurement (conservative sector-mean version, and joint version)
W_sec = np.exp(-0.5 * ((D['Aell_pred_ppm'] - ooe['Aell_sectors']) / ooe['Aell_sectors_err']) ** 2)
W_joint = np.exp(-0.5 * ((D['Aell_pred_ppm'] - ooe['Aell_joint']) / ooe['Aell_joint_err']) ** 2)
def wpct(x, w, qs=(2.5, 16, 50, 84, 97.5)):
    o = np.argsort(x); c = np.cumsum(w[o]); c /= c[-1]; return [float(np.interp(q / 100, c, x[o])) for q in qs]
def summ(x):
    return [float(v) for v in np.percentile(x, [2.5, 16, 50, 84, 97.5])]
out = {k_: summ(v) for k_, v in D.items()}; out.update({k_: summ(v) for k_, v in D2.items()})
out['with_ellipsoidal_weight_sectormean'] = {k_: wpct(np.asarray(v, float), W_sec) for k_, v in D.items() if np.ndim(v) == 1 and len(v) == N}
out['with_ellipsoidal_weight_joint'] = {k_: wpct(np.asarray(v, float), W_joint) for k_, v in D.items() if np.ndim(v) == 1 and len(v) == N}
out['ellipsoidal_weight_ESS_frac'] = [float(W_sec.sum() ** 2 / (W_sec ** 2).sum() / N), float(W_joint.sum() ** 2 / (W_joint ** 2).sum() / N)]
print('ellipsoidal reweighting (sector-mean): ESS fraction %.3f; L3 %s; M1 %s' % (out['ellipsoidal_weight_ESS_frac'][0],
      np.round(out['with_ellipsoidal_weight_sectormean']['L3'], 4), np.round(out['with_ellipsoidal_weight_sectormean'].get('M1', [np.nan]), 3)))
out['n_samples'] = int(N); out['lnp_max'] = float(LP.max()); out['kind'] = kind; out['burn'] = burn; out['thin'] = thin
for k_, v in out.items():
    if isinstance(v, list): print(f'{k_:22s} median {v[2]:.5g}   68%: [{v[1]:.5g}, {v[3]:.5g}]   95%: [{v[0]:.5g}, {v[4]:.5g}]')
if True:
    out['P_DEB_gt_1.5DG'] = float(np.mean(D['D_ratio'] > 1.5)); out['P_DEB_gt_2DG'] = float(np.mean(D['D_ratio'] > 2.0))
    out['P_DEB_lt_0.67DG'] = float(np.mean(D['D_ratio'] < 0.667))
    print('P(D_EB > 1.5 D_G) = %.4f; P(D_EB > 2 D_G) = %.4f; P(D_EB < D_G/1.5) = %.4f' % (out['P_DEB_gt_1.5DG'], out['P_DEB_gt_2DG'], out['P_DEB_lt_0.67DG']))
C = {}
for a_, b_ in [('L3', 'k'), ('L3', 'b1'), ('L3', 'r1'), ('k', 'b1'), ('r1', 'b1'), ('L3', 'J')] + ([('M1', 'L3'), ('M1', 'fR'), ('M1', 'dlnJ')] if kind == 'ms' else []):
    C[f'{a_}-{b_}'] = float(np.corrcoef(D[a_], D[b_])[0, 1])
out['correlations'] = C; print('correlations:', {k_: round(v, 2) for k_, v in C.items()})
json.dump(out, open(f'{FD}/results/post_{tag}.json', 'w'), indent=1)
np.savez_compressed(f'{FD}/results/post_samples_{tag}.npz', **{k_: np.asarray(v) for k_, v in D.items()})
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
keys = ['L3', 'k', 'b1', 'r1', 'J'] + (['M1', 'M2', 'fR', 'D_ratio'] if kind == 'ms' else ['rho1_sun_qMSk', 'esinw'])
n = len(keys)
fig, ax = plt.subplots(n, n, figsize=(1.9 * n, 1.9 * n))
for i in range(n):
    for j in range(n):
        a_ = ax[i, j]
        if j > i: a_.axis('off'); continue
        if i == j: a_.hist(D[keys[i]], 50, color='k', histtype='step')
        else: a_.plot(D[keys[j]], D[keys[i]], ',', alpha=0.3, color='k')
        if i == n - 1: a_.set_xlabel(keys[j], fontsize=8)
        if j == 0 and i > 0: a_.set_ylabel(keys[i], fontsize=8)
        a_.tick_params(labelsize=5)
plt.tight_layout(); plt.savefig(f'{FD}/plots/corner_{tag}.png', dpi=80)
