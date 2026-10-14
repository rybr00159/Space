"""Posterior summaries from an emcee chain (geo or ms/msell), derived quantities and the bound-vs-background photometric test.
Derived (distance independent): r1=R1/a, r2, b1 = cos i / r1, i, rho1 = 3 pi /(G P^2) / r1^3 / (1+q) [q from MS model or
a prior], f_EB = 1 - L3, undiluted depths, T2/T1 (MS model), A_ell prediction.
Bound test (ms kinds): photometric distance of the EB: m_T,EB = Tmag_tot - 2.5 log10 f_EB ; M_T,EB from Pecaut & Mamajek
MS relations (+0.2 mag scatter) ; A_T = 0.03 +- 0.02 ; D_EB = 10^((m-M-A+5)/5) ; compare with Gaia distance of the G star."""
import numpy as np, sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/scripts')
import fit3, fit2, lclib as L, ms_table as MS
from ebmodel import eb_flux
tag = sys.argv[1]; burn_frac = float(sys.argv[2]) if len(sys.argv) > 2 else 0.4; thin = int(sys.argv[3]) if len(sys.argv) > 3 else 10
z = np.load(L.HERE + f'/results/chain_{tag}.npz', allow_pickle=True)
ch = z['chain']; lnp = z['lnp']; names = list(z['names']); aell = z['aell']
ns = ch.shape[0]; b = int(burn_frac * ns)
X = ch[b::thin].reshape(-1, ch.shape[2]); LP = lnp[b::thin].reshape(-1); AE = aell[b::thin].reshape(-1)
ok = np.isfinite(LP); X, LP, AE = X[ok], LP[ok], AE[ok]
S = {n: X[:, i] for i, n in enumerate(names)}
kind = 'geo' if 'rsum' in names else 'ms'
G = 6.674e-8; P = L.P0 * 86400; RHO_SUN = 1.408
if kind == 'ms':
    M1 = S['M1']; M2 = S['M1'] * S['q']
    R1 = np.array([MS.radius(m) for m in M1]) * S['fR']; R2 = np.array([MS.radius(m) for m in M2]) * S['fR']
    a = fit2.a_over_R(M1 + M2, L.P0)
    rsum = (R1 + R2) / a; k = R2 / R1
    J = np.array([MS.SB_T(m2) / MS.SB_T(m1) for m1, m2 in zip(M1, M2)]) * np.exp(S['dlnJ'])
    q = S['q']
else:
    rsum, k, J = S['rsum'], S['k'], S['J']
    rng = np.random.default_rng(1); q = rng.uniform(0.2, 1.0, rsum.size)   # mass-ratio prior for the density
r1 = rsum / (1 + k); r2 = rsum * k / (1 + k); cosi = S['bn'] * rsum; inc = np.degrees(np.arccos(cosi))
b1 = cosi / r1
rho1 = 3 * np.pi / (G * P ** 2) / r1 ** 3 / (1 + q)
fEB = 1 - S['L3']
L2 = k ** 2 * J / (1 + k ** 2 * J)
D = dict(r1=r1, r2=r2, k=k, J=J, rsum=rsum, b1=b1, inc=inc, rho1_cgs=rho1, rho1_sun=rho1 / RHO_SUN, fEB=fEB, L3=S['L3'],
         sec_depth_undiluted=L2, esinw=S['esinw'], ecosw=S['ecosw'], e=np.hypot(S['esinw'], S['ecosw']), s95=S['s95'], s105=S['s105'], s106=S['s106'])
if kind == 'ms':
    T1 = np.array([MS.teff(m) for m in M1]); T2 = np.array([MS.teff(m) for m in M2])
    D.update(M1=M1, M2=M2, R1=R1, R2=R2, T1=T1, T2=T2, T2T1=T2 / T1, fR=S['fR'], dlnJ=S['dlnJ'], a_Rsun=a, Aell_ppm=AE * 1e6)
    # photometric distance of the EB
    rng = np.random.default_rng(2)
    MT1 = np.array([MS.MT(m) for m in M1]); MT2 = np.array([MS.MT(m) for m in M2])
    MTEB = -2.5 * np.log10(10 ** (-0.4 * MT1) + 10 ** (-0.4 * MT2)) + rng.normal(0, 0.2, M1.size)
    Tmag = 11.358 + rng.normal(0, 0.006, M1.size)
    mTEB = Tmag - 2.5 * np.log10(fEB)
    AT = np.clip(rng.normal(0.03, 0.02, M1.size), 0, None)
    DEB = 10 ** ((mTEB - MTEB - AT + 5) / 5)
    plx = rng.normal(3.3622 + 0.02, np.hypot(0.017, 0.015), M1.size)   # Lindegren+21 zero point ~ -0.02 mas (+-0.015 assumed)
    DG = 1000 / plx
    D.update(MT_EB=MTEB, mT_EB=mTEB, D_EB_pc=DEB, D_ratio=DEB / DG, dT_contrast=mTEB - (Tmag + 2.5 * np.log10(1 / (1 - fEB))) )
    # variant: luminosity scaled by fR^2 (radius inflation at the MS Teff used for J)
    DEB2 = 10 ** ((mTEB - (MTEB - 5 * np.log10(S['fR'])) - AT + 5) / 5)
    D.update(D_ratio_fR2=DEB2 / DG)
    # Gaia G of the pair if bound, and contrast with the G star (G_tot = 11.775)
    MG1 = np.array([MS.interp('MG', m) for m in M1]); MG2 = np.array([MS.interp('MG', m) for m in M2])
    GEB = -2.5 * np.log10(10 ** (-0.4 * MG1) + 10 ** (-0.4 * MG2)) + 5 * np.log10(DG / 10) + 1.0 * AT / 0.6 * 0.8
    GG = -2.5 * np.log10(np.clip(10 ** (-0.4 * 11.775) - 10 ** (-0.4 * GEB), 1e-30, None))
    D.update(G_EB_ifbound=GEB, dG_ifbound=GEB - GG)
    MK1 = np.array([MS.MKs(m) for m in M1]); MK2 = np.array([MS.MKs(m) for m in M2])
    MKEB = -2.5 * np.log10(10 ** (-0.4 * MK1) + 10 ** (-0.4 * MK2))
    # if bound (at the G star distance): predicted Ks of the pair and Ks contrast against the G star (K_tot = 10.341)
    KEB = MKEB + 5 * np.log10(DG / 10)
    Ktot = 10.341; KG = -2.5 * np.log10(np.clip(10 ** (-0.4 * Ktot) - 10 ** (-0.4 * KEB), 1e-30, None))
    D.update(Ks_EB_ifbound=KEB, dKs_ifbound=KEB - KG)
    # bound prediction of f_EB in TESS band
    MTtot = Tmag - 5 * np.log10(DG / 10) - AT
    D.update(fEB_pred_bound=10 ** (-0.4 * (MTEB - MTtot)))
# model eclipse depths (diluted, s=1 SPOC scale; undiluted) and T14, for 400 random posterior samples
rng = np.random.default_rng(9); idx = rng.choice(X.shape[0], min(400, X.shape[0]), replace=False)
dp, ds, t14 = [], [], []
for i in idx:
    u1a, u2a = fit3.q2u(S['q1a'][i], S['q2a'][i]); u1b, u2b = fit3.q2u(S['q1b'][i], S['q2b'][i])
    tt = np.linspace(-0.06, 0.06, 1201)
    args = (rsum[i], k[i], cosi[i], J[i], S['L3'][i], S['ecosw'][i], S['esinw'][i], u1a, u2a, u1b, u2b, 1e-6, 1, fit3.XG, fit3.WG)
    fp = eb_flux(tt, 0.0, L.P0, *args)
    fs = eb_flux(tt + L.P0 / 2 + 2 * L.P0 / np.pi * S['ecosw'][i], 0.0, L.P0, *args)
    dp.append(1 - fp.min()); ds.append(1 - fs.min()); w = tt[fp < 1 - 1e-9]; t14.append((w.max() - w.min()) * 24)
dp = np.array(dp); ds = np.array(ds); fe = fEB[idx]
D2 = dict(depth_prim_diluted_ppm=dp * 1e6, depth_sec_diluted_ppm=ds * 1e6, depth_prim_undiluted=dp / fe, depth_sec_undiluted=ds / fe,
          depth_ratio=ds / dp, T14_prim_h=np.array(t14))
def summ(x):
    q = np.percentile(x, [2.5, 16, 50, 84, 97.5]); return [float(v) for v in q]
out = {k: summ(v) for k, v in D.items()}
out.update({k: summ(v) for k, v in D2.items()})
out['n_samples'] = int(X.shape[0]); out['lnp_max'] = float(LP.max())
for k, v in out.items():
    if isinstance(v, list) and len(v) == 5: print(f'{k:22s} median {v[2]:.5g}  68%: [{v[1]:.5g}, {v[3]:.5g}]  95%: [{v[0]:.5g}, {v[4]:.5g}]')
if kind == 'ms':
    out['P_DEB_gt_DG'] = float(np.mean(D['D_ratio'] > 1)); out['P_DEB_gt_1.5DG'] = float(np.mean(D['D_ratio'] > 1.5))
    print('P(D_EB > D_G) = %.3f ; P(D_EB > 1.5 D_G) = %.3f' % (out['P_DEB_gt_DG'], out['P_DEB_gt_1.5DG']))
    # correlation summaries
    print('corr(L3, M1) = %.2f ; corr(L3,k) = %.2f ; corr(L3, b1) = %.2f' % (np.corrcoef(D['L3'], D['M1'])[0, 1], np.corrcoef(D['L3'], D['k'])[0, 1], np.corrcoef(D['L3'], D['b1'])[0, 1]))
else:
    print('corr(L3,k) = %.2f ; corr(L3, b1) = %.2f ; corr(L3, rsum) = %.2f; corr(k,b1) = %.2f' % (np.corrcoef(D['L3'], D['k'])[0, 1], np.corrcoef(D['L3'], D['b1'])[0, 1], np.corrcoef(D['L3'], D['rsum'])[0, 1], np.corrcoef(D['k'], D['b1'])[0, 1]))
json.dump(out, open(L.HERE + f'/results/post_{tag}.json', 'w'), indent=1)
np.savez_compressed(L.HERE + f'/results/post_samples_{tag}.npz', **{k: np.asarray(v) for k, v in D.items()})
# corner-style plot of key parameters
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
keys = ['L3', 'k', 'b1', 'r1', 'J'] + (['M1', 'M2', 'D_ratio'] if kind == 'ms' else ['rho1_sun'])
n = len(keys)
fig, ax = plt.subplots(n, n, figsize=(2 * n, 2 * n))
for i in range(n):
    for j in range(n):
        a_ = ax[i, j]
        if j > i: a_.axis('off'); continue
        if i == j: a_.hist(D[keys[i]], 50, color='k', histtype='step')
        else: a_.plot(D[keys[j]], D[keys[i]], ',', alpha=0.3, color='k')
        if i == n - 1: a_.set_xlabel(keys[j])
        if j == 0 and i > 0: a_.set_ylabel(keys[i])
        a_.tick_params(labelsize=6)
plt.tight_layout(); plt.savefig(L.HERE + f'/plots/corner_{tag}.png', dpi=80)
