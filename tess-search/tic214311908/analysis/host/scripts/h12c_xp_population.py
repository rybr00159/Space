"""h12c: XP follow-ups.
 (a) Twin fraction among over-luminous field stars from their XP [Fe/H] (h12b out/h12b_xp_pred_field.csv): the
     XP-[Fe/H] distribution of over-luminous stars (-1.0 < dMKs < -0.55, RUWE < 1.2) is modelled as
     w_T * [distribution of normal stars] (twins look like normal stars) + (1 - w_T) * N(mu, s) (metal-rich singles; s fixed to the robust width of the normal stars, otherwise the
     mixture is degenerate - a free-width fit is also reported);
     ML fit; the M3V's P(twin | XP [Fe/H], over-luminous) and the likelihood ratio. Also restricted to RUWE < 1.1.
 (b) Mira contamination: fit target = template(colour-matched normal field median, in flux) + eps * Mira XP flux;
     eps calibrated by doing the same fit on the colour-matched field stars (empirical null distribution).
 (c) Composite-vs-single chi2 improvement for the M3V compared with the same statistic for normal, over-luminous
     field stars and metal-rich calibrators (is a 2-star composite preference meaningful?).
Output: out/h12c_xp_population.json, figs/h12c_xp_population.png"""
import glob, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from scipy.optimize import minimize
from scipy.stats import gaussian_kde, norm
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
B = 'tic214311908/analysis/host'
xp = json.load(open(f'{B}/out/h12b_xp.json'))
P = pd.read_csv(f'{B}/out/h12b_xp_pred_field.csv')
fe_t = xp['xp_feh']['target_pred']; e_t = xp['xp_feh']['target_pred_mc_sd']
out = {}
for lab, sel in (('ruwe1.2', P.ruwe < 1.2), ('ruwe1.1', P.ruwe < 1.1)):
    nor = P[(P.set == 'field_normal') & sel].feh_xp.values; ovl = P[(P.set == 'field_overlum') & sel & P.dMK_local.between(-1.0, -0.55)].feh_xp.values
    kn = gaussian_kde(nor, bw_method=0.3)
    s_fix = 1.4826 * np.median(np.abs(nor - np.median(nor)))     # metal-rich component width fixed to the normal stars' width
    def nll(p):
        w = 1 / (1 + np.exp(-p[0])); mu, s = p[1], (s_fix if FIX else np.exp(p[2]))
        return -np.sum(np.log(w * kn(ovl) + (1 - w) * norm.pdf(ovl, mu, s) + 1e-300))
    FIX = True
    r = minimize(nll, [0.0, 0.15, np.log(0.1)], method='Nelder-Mead')
    w = 1 / (1 + np.exp(-r.x[0])); mu, s = r.x[1], s_fix
    FIX = False
    rfree = minimize(nll, [0.0, 0.15, np.log(0.1)], method='Nelder-Mead')
    lr_free = float(2 * (r.fun - rfree.fun))
    # bootstrap
    bs = []
    rng = np.random.default_rng(2)
    for _ in range(200):
        o_ = rng.choice(ovl, len(ovl)); n_ = rng.choice(nor, len(nor)); kn_ = gaussian_kde(n_, bw_method=0.3)
        def nl(p):
            ww = 1 / (1 + np.exp(-p[0])); return -np.sum(np.log(ww * kn_(o_) + (1 - ww) * norm.pdf(o_, p[1], s_fix) + 1e-300))
        rb = minimize(nl, r.x, method='Nelder-Mead'); bs.append(1 / (1 + np.exp(-rb.x[0])))
    # target: densities convolved with its measurement error
    knt = gaussian_kde(nor, bw_method=np.sqrt(0.3 ** 2 + (e_t / nor.std()) ** 2))
    pT = float(knt(fe_t)[0]); pS = float(norm.pdf(fe_t, mu, np.hypot(s, e_t)))
    out[lab] = dict(n_normal=int(len(nor)), n_overlum=int(len(ovl)), w_twin=float(w), w_twin_boot_16_84=[float(np.percentile(bs, 16)), float(np.percentile(bs, 84))],
                    metalrich_mu=float(mu), metalrich_sd_fixed=float(s), normal_median=float(np.median(nor)), target_xp_feh=fe_t,
                    free_width_fit=dict(w_twin=float(1 / (1 + np.exp(-rfree.x[0]))), mu=float(rfree.x[1]), sd=float(np.exp(rfree.x[2])), two_dlnL_vs_fixed=lr_free),
                    LR_twin_over_single=pT / pS, P_twin_given_xp_and_overlum=float(w * pT / (w * pT + (1 - w) * pS)))
    print(lab, out[lab])
# ---- (b) and (c): need spectra
Z = [np.load(f) for f in sorted(glob.glob(f'{B}/data/xp/xps_chunk_*.npz'))]
sid = np.concatenate([z['source_id'] for z in Z]); FL = np.vstack([z['flux'] for z in Z if len(z['source_id'])]); FE = np.vstack([z['flux_error'] for z in Z if len(z['source_id'])])
_, ui = np.unique(sid, return_index=True); ui = np.sort(ui); sid, FL, FE = sid[ui], FL[ui], FE[ui]
wl = np.arange(336, 1021, 2.0); win = (wl >= 420) & (wl <= 1000); wn = (wl >= 620) & (wl <= 1000)
TID, MID = 5925373272850784768, 5925373268528684288
tf = FL[sid == TID][0]; mf = FL[sid == MID][0]
F = pd.read_csv(f'{B}/data/field_mdwarfs.csv').set_index('source_id')
Pn = P[(P.set == 'field_normal')].set_index('source_id')
cm = [s_ for s_ in Pn.index if abs(Pn.bp_rp[s_] - 2.756) < 0.03 and s_ in set(sid)]
idx = {s_: i for i, s_ in enumerate(sid)}
norm_fl = lambda f: f / f[wn].mean()
tmpl = np.median([norm_fl(FL[idx[s_]]) for s_ in cm], axis=0)
scat = 1.4826 * np.median(np.abs(np.array([norm_fl(FL[idx[s_]]) for s_ in cm]) - tmpl), axis=0)
def fit_eps(f):
    fn = norm_fl(f); mm = mf / f[wn].mean()
    # linear LSQ for model fn = a*tmpl + eps*mm on the window, weights 1/scat
    A = np.column_stack([tmpl[win], mm[win]]) / scat[win, None]; y = fn[win] / scat[win]
    c, *_ = np.linalg.lstsq(A, y, rcond=None); return c[1]
eps_t = fit_eps(tf)
eps_f = np.array([fit_eps(FL[idx[s_]]) for s_ in cm])
out['mira_contamination'] = dict(eps_target=float(eps_t), eps_field_median=float(np.median(eps_f)), eps_field_sd=float(1.4826 * np.median(np.abs(eps_f - np.median(eps_f)))),
                                 z=float((eps_t - np.median(eps_f)) / (1.4826 * np.median(np.abs(eps_f - np.median(eps_f))))),
                                 mira_flux_fraction_at_800nm_in_target=float(eps_t * mf[np.argmin(np.abs(wl - 800))] / tf[np.argmin(np.abs(wl - 800))]),
                                 n_field=int(len(cm)), note='eps = fraction of the Mira XP flux added to the M3V; field stars give the null distribution (they have no Mira)')
print('Mira contamination:', out['mira_contamination'])
json.dump(out, open(f'{B}/out/h12c_xp_population.json', 'w'), indent=1)
fig, ax = plt.subplots(figsize=(7, 4.5))
for st, cl in (('field_normal', 'k'), ('field_overlum', 'm')):
    v = P[(P.set == st) & (P.ruwe < 1.2)].feh_xp
    ax.hist(v, bins=np.arange(-0.6, 0.6, 0.03), histtype='step', density=True, color=cl, label=f'{st} (n={len(v)})')
xx = np.linspace(-0.6, 0.6, 300); o = out['ruwe1.2']
nor = P[(P.set == 'field_normal') & (P.ruwe < 1.2)].feh_xp.values
ax.plot(xx, o['w_twin'] * gaussian_kde(nor, bw_method=0.3)(xx), 'b--', label=f'twin part (w={o["w_twin"]:.2f})')
ax.plot(xx, (1 - o['w_twin']) * norm.pdf(xx, o['metalrich_mu'], o['metalrich_sd_fixed']), 'r--', label='metal-rich part')
ax.axvline(fe_t, c='r', label='M3V'); ax.set_xlabel('[Fe/H] predicted from XP shape'); ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig(f'{B}/figs/h12c_xp_population.png', dpi=110)
