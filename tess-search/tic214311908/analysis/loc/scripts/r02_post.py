"""Posterior / likelihood analysis of the attribution from the simulation outputs (out/sim_{tag}.npz).

(1) Model-based posterior, equal priors over the 93 catalogued stars (incl. the Mira at its measured position):
    P_j  ~  exp(-chi2_j / (2 kappa)),  chi2_j of the 'a*PSF_j + b' fit (per sector; combined = sum over sectors).
    kappa (temperature) is calibrated on the injections: maximum log-score of the true source over all injections
    (each of the 93 stars injected 1000x with the observed signal).  Reliability is checked by binning.
(2) Empirical likelihood of the per-sector free-fit position: for each star j and sector s the recovered free-fit
    positions of its injections define a Gaussian core (3-sigma-clipped mean/covariance) plus an outlier fraction
    (uniform over the grid); L_j = prod_s L_js(real free-fit position).  Built on odd injections, tested on even ones.
(3) Localisation: combined chi2 map (grid co-moving with the M3V, source colour BP-RP 2.76); the 68.3/95.4/99.73%
    thresholds are the quantiles of the combined delta-chi2 at the true position for M3V injections (and ring
    injections when present).
(4) Uncatalogued source: U = min_j chi2_j(combined) - min_grid; real value against injections on catalogued stars.
Usage: r02_post.py TAG [TAG_RING]
Output: out/post_{TAG}.json, figs/loc_map_{TAG}.png
"""
import sys, json
sys.path.insert(0, 'tic214311908/analysis/loc/scripts')
import numpy as np
from l0lib import BASE, to_json, SECTORS
tag = sys.argv[1]
Z = np.load(f'{BASE}/out/sim_{tag}.npz', allow_pickle=True)
names = Z['hyp_name']; nh = len(names); secs = list(Z['sectors'])
NS = len(secs)
inj = [k[4:-6] for k in Z.files if k.startswith('inj_') and k.endswith('__chiH')]
MAINS = ['M3V', 'pairA', 'pairB', 'Mira', '5925373272827260544', '5925373272827262208', '5925373272827256576']
idx = {n: i for i, n in enumerate(names)}
real = Z['real__chiH']                     # (NS, nh)


def softmax_post(chi, kappa=1.0):
    x = -0.5 * (chi - chi.min(-1, keepdims=True)) / kappa
    e = np.exp(x); return e / e.sum(-1, keepdims=True)


def group(P):
    """posterior summed into groups"""
    g = {'M3V': P[..., idx['M3V']], 'pair(A+B)': P[..., idx['pairA']] + P[..., idx['pairB']], 'Mira': P[..., idx['Mira']]}
    g['faint 3.2" from M3V (..260544)'] = P[..., idx['5925373272827260544']]
    g['faint 5.7" from M3V (..262208)'] = P[..., idx['5925373272827262208']]
    others = np.ones(nh, bool)
    for n in ('M3V', 'pairA', 'pairB', 'Mira', '5925373272827260544', '5925373272827262208'): others[idx[n]] = False
    g['all other catalogued'] = P[..., others].sum(-1)
    return g


out = dict(tag=tag, sectors=[int(s) for s in secs], n_inj_targets=len(inj))
# ---------- (1) model-based posterior, kappa calibration
chi_inj = {t: Z[f'inj_{t}__chiH'].astype(float) for t in inj}            # (ni, NS, nh)
truth = {t: idx[t] for t in inj}
def logscore(kappa, which):
    tot = 0.0; n = 0
    for t in inj:
        c = chi_inj[t][:, which, :].sum(1) if isinstance(which, list) else chi_inj[t][:, which, :]
        P = softmax_post(c, kappa)
        tot += np.log(np.clip(P[:, truth[t]], 1e-300, None)).sum(); n += len(P)
    return tot / n
kgrid = np.array([0.5, 0.7, 0.85, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0])
kap = {}
for k, s in enumerate(secs):
    ls = [logscore(kv, k) for kv in kgrid]; kap[str(s)] = float(kgrid[int(np.argmax(ls))])
lsC = [logscore(kv, list(range(NS))) for kv in kgrid]; kap['combined'] = float(kgrid[int(np.argmax(lsC))])
out['kappa'] = kap
print('calibrated kappa per sector / combined:', kap)
# reliability (combined, kappa=1 and calibrated), near-source subset
def reliability(kappa, near_only=True):
    bins = np.array([0, 0.01, 0.05, 0.2, 0.5, 0.8, 0.95, 0.99, 1.0001]); hit = np.zeros(len(bins) - 1); cnt = np.zeros(len(bins) - 1); psum = np.zeros(len(bins) - 1)
    near = [t for t in inj if np.hypot(Z['hyp_dE'][idx[t]] - Z['hyp_dE'][idx['M3V']], Z['hyp_dN'][idx[t]] - Z['hyp_dN'][idx['M3V']]) < 25 or
            np.hypot(Z['hyp_dE'][idx[t]], Z['hyp_dN'][idx[t]]) < 25] if near_only else inj
    for t in near:
        P = softmax_post(chi_inj[t].sum(1), kappa)
        cols = [idx[n] for n in near]
        for j in cols:
            b = np.digitize(P[:, j], bins) - 1
            for bb in range(len(bins) - 1):
                m = b == bb; cnt[bb] += m.sum(); psum[bb] += P[m, j].sum(); hit[bb] += (m.sum() if j == truth[t] else 0)
    return [dict(lo=float(bins[i]), hi=float(bins[i + 1]), n=int(cnt[i]), mean_P=float(psum[i] / max(cnt[i], 1)), frac_true=float(hit[i] / max(cnt[i], 1))) for i in range(len(cnt))], len(near)
rel1, nnear = reliability(1.0); relk, _ = reliability(kap['combined'])
out['reliability_kappa1_near'] = rel1; out['reliability_kappaC_near'] = relk; out['n_near_sources'] = nnear
print(f'reliability of combined posterior among {nnear} sources within 25" of the M3V or the Mira (kappa=1):')
for r in rel1: print(f"   P in [{r['lo']:.2f},{r['hi']:.2f}): n={r['n']:6d} mean P {r['mean_P']:.3f}  true fraction {r['frac_true']:.3f}")
# real posteriors
res_real = {}
for k, s in enumerate(secs):
    for kk, kv in (('k1', 1.0), ('kcal', kap[str(s)])):
        P = softmax_post(real[k], kv); res_real[f'S{s}_{kk}'] = {g: float(v) for g, v in group(P).items()}
for kk, kv in (('k1', 1.0), ('kcal', kap['combined'])):
    P = softmax_post(real.sum(0), kv); res_real[f'combined_{kk}'] = {g: float(v) for g, v in group(P).items()}
    order = np.argsort(-P)[:6]
    res_real[f'combined_{kk}_top'] = [(str(names[i]), float(P[i])) for i in order]
out['posterior_real'] = res_real
dchi = {}
for k, s in enumerate(secs):
    c = real[k]; dchi[str(s)] = {n: float(c[idx[n]] - c.min()) for n in MAINS}
    dchi[str(s)]['best'] = str(names[int(np.argmin(c))]); dchi[str(s)]['grid_min_vs_best'] = float(c.min() - Z['real__chiG'][k].min())
cc = real.sum(0); dchi['combined'] = {n: float(cc[idx[n]] - cc.min()) for n in MAINS}; dchi['combined']['best'] = str(names[int(np.argmin(cc))])
dchi['combined']['grid_min_vs_best'] = float(cc.min() - Z['real__gsum'].min())
out['dchi2_real'] = dchi
print('\nreal data: chi2 - min over catalogued stars (per sector, combined):')
for k, v in dchi.items():
    print(f'  {k:>8s}: ' + ', '.join(f'{n[-8:]} {v[n]:6.1f}' for n in MAINS) + f"; best {v['best']}; best - free-fit {v['grid_min_vs_best']:.1f}")
print('\nreal posteriors (equal priors):')
for k, v in res_real.items():
    if k.endswith('_top'): print(f'  {k}: {v}'); continue
    print(f'  {k:16s}: ' + ', '.join(f'{g} {p:.3g}' for g, p in v.items()))
# injection confusion matrix for the main sources (combined, kappa=1): mean posterior groups
conf = {}
for t in MAINS:
    if t not in chi_inj: continue
    P = softmax_post(chi_inj[t].sum(1), 1.0); g = group(P)
    conf[t] = {gn: dict(mean=float(v.mean()), frac_gt_0p5=float((v > 0.5).mean())) for gn, v in g.items()}
    conf[t]['P_true_median'] = float(np.median(P[:, idx[t]]))
    conf[t]['frac_true_is_MAP'] = float((np.argmax(P, 1) == idx[t]).mean())
    # per-sector + combined distributions of chi2(pair) - chi2(M3V) and chi2(Mira) - chi2(M3V)
    c = chi_inj[t]
    conf[t]['dchi_pair_minus_M3V_comb_q'] = np.percentile(c[:, :, idx['pairA']].sum(1) - c[:, :, idx['M3V']].sum(1), [2.5, 16, 50, 84, 97.5]).tolist()
    conf[t]['dchi_Mira_minus_M3V_comb_q'] = np.percentile(c[:, :, idx['Mira']].sum(1) - c[:, :, idx['M3V']].sum(1), [2.5, 16, 50, 84, 97.5]).tolist()
    conf[t]['dchi_260544_minus_M3V_comb_q'] = np.percentile(c[:, :, idx['5925373272827260544']].sum(1) - c[:, :, idx['M3V']].sum(1), [2.5, 16, 50, 84, 97.5]).tolist()
out['injection_confusion_combined_k1'] = conf
rc = real.sum(0)
out['real_dchi_pair_minus_M3V'] = float(rc[idx['pairA']] - rc[idx['M3V']]); out['real_dchi_Mira_minus_M3V'] = float(rc[idx['Mira']] - rc[idx['M3V']])
out['real_dchi_260544_minus_M3V'] = float(rc[idx['5925373272827260544']] - rc[idx['M3V']])
print('\ninjection test (combined): distribution of chi2(pairA)-chi2(M3V), chi2(Mira)-chi2(M3V), chi2(..260544)-chi2(M3V) [2.5,16,50,84,97.5%]:')
for t, v in conf.items():
    print(f"  true {t[-12:]:>12s}: pair-M3V {np.round(v['dchi_pair_minus_M3V_comb_q'], 1).tolist()}  Mira-M3V {np.round(v['dchi_Mira_minus_M3V_comb_q'], 1).tolist()}  "
          f"260544-M3V {np.round(v['dchi_260544_minus_M3V_comb_q'], 1).tolist()}  P(true) median {v['P_true_median']:.3f}, MAP correct {v['frac_true_is_MAP']:.3f}")
print(f"  REAL:              pair-M3V {out['real_dchi_pair_minus_M3V']:.1f}  Mira-M3V {out['real_dchi_Mira_minus_M3V']:.1f}  260544-M3V {out['real_dchi_260544_minus_M3V']:.1f}")
# empirical p-values: how often does a true-pair (true-Mira) injection give pair-M3V >= real value?
pv = {}
for t in ('pairA', 'pairB', 'Mira', '5925373272827260544', '5925373272827262208'):
    if t not in chi_inj: continue
    c = chi_inj[t]
    d_pm = c[:, :, idx[t]].sum(1) - c[:, :, idx['M3V']].sum(1)   # chi2(true) - chi2(M3V)
    real_d = rc[idx[t]] - rc[idx['M3V']]
    pv[t] = dict(real=float(real_d), frac_inj_ge_real=float((d_pm >= real_d).mean()), n=len(d_pm),
                 q=np.percentile(d_pm, [50, 84, 97.5, 99.85]).tolist())
    cm = chi_inj['M3V']; d_m = cm[:, :, idx[t]].sum(1) - cm[:, :, idx['M3V']].sum(1)
    pv[t]['frac_M3Vinj_le_real'] = float((d_m <= real_d).mean()); pv[t]['M3Vinj_q'] = np.percentile(d_m, [0.15, 2.5, 16, 50]).tolist()
out['pvalues_chi2_star_minus_M3V'] = pv
print('\nempirical tail probabilities of chi2(star)-chi2(M3V) [combined]:')
for t, v in pv.items():
    print(f"  {t:>20s}: real {v['real']:.1f}; if {t} were the source, P(>= real) = {v['frac_inj_ge_real']:.4f} (n={v['n']}), quantiles {np.round(v['q'],1).tolist()}; "
          f"if M3V were the source, P(<= real) = {v['frac_M3Vinj_le_real']:.3f}, quantiles {np.round(v['M3Vinj_q'],1).tolist()}")
# ---------- (2) empirical position likelihood
def robust_gauss(X):
    m = np.median(X, 0); C = np.cov(X.T) + 1e-6 * np.eye(2)
    for _ in range(10):
        d = np.einsum('ij,jk,ik->i', X - m, np.linalg.inv(C), X - m)
        keep = d < 11.8
        m = X[keep].mean(0); C = np.cov(X[keep].T) + 1e-6 * np.eye(2)
    eps = 1 - keep.mean()
    return m, C, eps
area = (Z['grid_axis'][-1] - Z['grid_axis'][0]) ** 2
def logN(x, m, C, eps):
    d = x - m; Ci = np.linalg.inv(C)
    q = np.einsum('...j,jk,...k->...', d, Ci, d)
    return np.log((1 - eps) * np.exp(-0.5 * q) / (2 * np.pi * np.sqrt(np.linalg.det(C))) + eps / area + 1e-300)
models = {}
for t in inj:
    F = Z[f'inj_{t}__free'].astype(float)        # (ni, NS, 2)
    models[t] = [robust_gauss(F[0::2, k]) for k in range(NS)]
realfree = []
for k in range(NS):
    # real free-fit position from the real grid map
    cg = Z['real__chiG'][k]; ax = Z['grid_axis']; n = len(ax); i = int(np.argmin(cg)); iy, ix = divmod(i, n)
    realfree.append([ax[ix], ax[iy]])
realfree = np.array(realfree)
out['real_free_fit_offset_from_M3V_arcsec'] = {str(s): realfree[k].tolist() for k, s in enumerate(secs)}
LL = np.zeros((NS, len(inj)))
for j, t in enumerate(inj):
    for k in range(NS):
        m, C, eps = models[t][k]; LL[k, j] = logN(realfree[k], m, C, eps)
def post_from_LL(L):
    x = L - L.max(); e = np.exp(x); return e / e.sum()
emp = {}
for k, s in enumerate(secs):
    P = post_from_LL(LL[k]); Pf = np.zeros(nh); Pf[[idx[t] for t in inj]] = P
    emp[f'S{s}'] = {g: float(v) for g, v in group(Pf).items()}
P = post_from_LL(LL.sum(0)); Pf = np.zeros(nh); Pf[[idx[t] for t in inj]] = P
emp['combined'] = {g: float(v) for g, v in group(Pf).items()}
emp['combined_top'] = [(str(names[idx[inj[i]]]), float(P[i])) for i in np.argsort(-P)[:6]]
out['posterior_empirical_position'] = emp
# calibration test of the empirical likelihood on even-numbered injections of the main sources
caltest = {}
for t in MAINS:
    if t not in inj: continue
    F = Z[f'inj_{t}__free'].astype(float)[1::2]
    Lt = np.zeros((len(F), len(inj)))
    for j, u in enumerate(inj):
        for k in range(NS):
            m, C, eps = models[u][k]; Lt[:, j] += logN(F[:, k], m, C, eps)
    Pt = np.exp(Lt - Lt.max(1, keepdims=True)); Pt /= Pt.sum(1, keepdims=True)
    jt = inj.index(t)
    caltest[t] = dict(mean_P_true=float(Pt[:, jt].mean()), frac_MAP_correct=float((np.argmax(Pt, 1) == jt).mean()),
                      mean_P_M3V=float(Pt[:, inj.index('M3V')].mean()))
out['empirical_position_calibration'] = caltest
print('\nempirical free-position likelihood posterior (equal priors):')
for k, v in emp.items(): print(f'  {k}: {v}')
print('  test on held-out injections:', caltest)
# ---------- (3) localisation
gs = Z['real__gsum']; ax = Z['grid_axis']; n = len(ax)
dmap = gs - gs.min()
thr = {}
dm3 = Z['inj_M3V__dchi_m3v'].astype(float)
thr['M3V_injections'] = np.percentile(dm3, [68.27, 95.45, 99.73]).tolist()
ring_keys = [k for k in Z.files if k.startswith('ring_') and k.endswith('__dchi_true')]
if len(sys.argv) > 2:
    R = np.load(f'{BASE}/out/sim_{sys.argv[2]}.npz', allow_pickle=True)
    ring_keys = [k for k in R.files if k.startswith('ring_') and k.endswith('__dchi_true')]
    allr = np.concatenate([R[k] for k in ring_keys]); thr['ring_injections'] = np.percentile(allr, [68.27, 95.45, 99.73]).tolist()
out['loc_thresholds'] = thr
T = thr['M3V_injections']
out['loc_thresholds_used'] = T
gE, gN = np.meshgrid(ax, ax); gE = gE.ravel(); gN = gN.ravel()
i0 = int(np.argmin(gs))
reg = {}
for lab, tv in zip(('68', '95', '99.7'), T):
    m = dmap <= tv
    reg[lab] = dict(threshold=float(tv), area_arcsec2=float(m.sum() * (ax[1] - ax[0]) ** 2),
                    max_dist_from_M3V=float(np.hypot(gE[m], gN[m]).max()), E_range=[float(gE[m].min()), float(gE[m].max())],
                    N_range=[float(gN[m].min()), float(gN[m].max())])
    # catalogued stars inside (by their own combined chi2 relative to the free-fit minimum)
    rr = rc - gs.min()
    reg[lab]['catalogued_inside'] = [(str(names[j]), round(float(rr[j]), 2)) for j in np.argsort(rr) if rr[j] <= tv]
out['localisation'] = dict(free_fit_offset_from_M3V=[float(gE[i0]), float(gN[i0])], regions=reg,
                           M3V_dchi=float(dmap[np.argmin(np.hypot(gE, gN))]))
print('\nlocalisation (combined, grid co-moving with M3V, colour 2.76): thresholds', np.round(T, 2).tolist(), 'ring:', np.round(thr.get('ring_injections', []), 2).tolist())
print(f'  free-fit minimum at ({gE[i0]:+.0f},{gN[i0]:+.0f})" from the M3V; M3V delta-chi2 {dmap[np.argmin(np.hypot(gE, gN))]:.2f}')
for lab, v in reg.items():
    print(f"  {lab}%: area {v['area_arcsec2']:.0f} arcsec2, max dist from M3V {v['max_dist_from_M3V']:.1f}\", E {v['E_range']}, N {v['N_range']}; stars inside: {v['catalogued_inside']}")
# ---------- (4) uncatalogued source statistic
U_real = float(rc.min() - gs.min())
Uinj = np.concatenate([Z[f'inj_{t}__chiH'].astype(float).sum(1).min(1) - Z[f'inj_{t}__gcomb_min'] for t in ('M3V',)])
Uall = np.concatenate([Z[f'inj_{t}__chiH'].astype(float).sum(1).min(1) - Z[f'inj_{t}__gcomb_min'] for t in MAINS if t in inj])
out['uncatalogued'] = dict(U_real=U_real, p_M3V_inj=float((Uinj >= U_real).mean()), p_main_inj=float((Uall >= U_real).mean()),
                           U_M3V_inj_q=np.percentile(Uinj, [50, 95, 99.7]).tolist())
print(f"\nuncatalogued-source statistic U = min_cat chi2 - free-fit chi2 = {U_real:.2f}; P(U >= real | catalogued M3V source) = {out['uncatalogued']['p_M3V_inj']:.3f}; quantiles {np.round(out['uncatalogued']['U_M3V_inj_q'],2).tolist()}")
to_json(out, f'{BASE}/out/post_{tag}.json')
# figure
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax_ = plt.subplots(1, 1, figsize=(7, 6.3))
M = dmap.reshape(n, n)
ax_.imshow(np.clip(M, 0, 60), origin='lower', extent=[ax[0] - .5, ax[-1] + .5, ax[0] - .5, ax[-1] + .5], cmap='viridis_r')
ax_.contour(ax, ax, M, levels=T, colors=['w', 'orange', 'r'], linewidths=1.2)
mira = -Z['m3v_off'].mean(0)
for j in range(nh):
    xE = Z['hyp_dE'][j] - Z['hyp_dE'][idx['M3V']]; xN = Z['hyp_dN'][j] - Z['hyp_dN'][idx['M3V']]
    ax_.plot(xE, xN, 'o', mfc='none', mec='k', ms=max(3, 14 - 0.8 * (Z['hyp_T'][j] - 11)))
for nm, c in (('M3V', 'c'), ('pairA', 'm'), ('Mira', 'r'), ('5925373272827260544', 'y'), ('5925373272827262208', 'y')):
    j = idx[nm]; ax_.plot(Z['hyp_dE'][j] - Z['hyp_dE'][idx['M3V']], Z['hyp_dN'][j] - Z['hyp_dN'][idx['M3V']], '+', color=c, ms=12, mew=2)
ax_.set_xlim(-40, 40); ax_.set_ylim(-35, 35); ax_.invert_xaxis()
ax_.set_xlabel('dE from M3V (arcsec, east left)'); ax_.set_ylabel('dN (arcsec)')
ax_.set_title(f'{tag}: combined delta-chi2 (contours 68/95/99.7%)\n(Gaia positions at 2016; grid co-moving with M3V)', fontsize=9)
fig.tight_layout(); fig.savefig(f'{BASE}/figs/loc_map_{tag}.png', dpi=90)
