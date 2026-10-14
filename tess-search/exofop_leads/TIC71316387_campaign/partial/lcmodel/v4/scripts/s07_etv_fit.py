"""ETV analysis of etv_<tag>.json:
 * empirical per-event errors = formal x z-scale(sector, eclipse type) from the injections (floor 1); optional bias
   correction with the injection median offset;
 * 5-sigma outlier rejection about the sector/type median;
 * linear ephemeris (T0, P, secondary offset), quadratic, cubic, and circular-LTTE (sinusoid grid 0.5-60 yr) fits;
 * empirical null for Delta-chi2 (bootstrap of normalised residuals within sector and eclipse type), 95/99.7% thresholds;
 * timing-level positive control (inject c2 * tyr^2 into the timings with bootstrap noise, recover);
 * per-sector mean O-C (primary, secondary, combined).
Usage: s07_etv_fit.py <tag> [biascorr 0/1]"""
import numpy as np, json, sys
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
FD = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v4'
tag = sys.argv[1]; BIAS = int(sys.argv[2]) if len(sys.argv) > 2 else 0
D = json.load(open(f'{FD}/results/etv_{tag}.json'))
P, T0 = D['P'], D['T0']
rows = [r for r in D['rows'] if not r['edge']]; inj = [r for r in D['inj'] if not r['edge']]
cal = {}
for sec in sorted(set(r['sector'] for r in rows)):
    for ph in (0.0, 0.5):
        ii = [r for r in inj if r['sector'] == sec and r['kind'] == ph]
        x = np.array([r['dt_min'] for r in ii]); e = np.array([r['err_min'] for r in ii]); z = x / e
        rsd = 1.4826 * np.median(np.abs(z - np.median(z)))
        cal[(sec, ph)] = dict(n=len(ii), bias_min=float(np.median(x)), bias_err=float(1.2533 * np.std(x) / np.sqrt(len(x))),
                              sd_min=float(1.4826 * np.median(np.abs(x - np.median(x)))), zscale=float(max(1.0, rsd)), zscale_raw=float(rsd),
                              med_formal=float(np.median(e)))
E = np.array([r['n'] + r['ph'] for r in rows]); oc = np.array([r['dt_min'] for r in rows]); ef = np.array([r['err_min'] for r in rows])
sec = np.array([r['sector'] for r in rows]); ph = np.array([r['ph'] for r in rows])
es = np.array([ef[i] * cal[(sec[i], ph[i])]['zscale'] for i in range(len(rows))])
if BIAS:
    oc = oc - np.array([cal[(sec[i], ph[i])]['bias_min'] for i in range(len(rows))])
t = T0 + E * P
keep = np.ones(len(rows), bool)
for s in np.unique(sec):
    for f in (0.0, 0.5):
        m = (sec == s) & (ph == f)
        if m.sum() < 3: continue
        med = np.median(oc[m]); keep[m] &= np.abs(oc[m] - med) < 5 * es[m]
print(f'events {len(rows)}, kept {keep.sum()} after 5-sigma clipping')
E, oc, ef, es, sec, ph, t = E[keep], oc[keep], ef[keep], es[keep], sec[keep], ph[keep], t[keep]
TREF = 2900.0
tyr = (t - TREF) / 365.25
def design(kind, Pout=None):
    cols = [np.ones_like(E), E / 1000.0, (ph == 0.5).astype(float)]
    if kind in ('quad', 'cubic'): cols.append(tyr ** 2)
    if kind == 'cubic': cols.append(tyr ** 3)
    if kind == 'ltte': cols += [np.sin(2 * np.pi * tyr / Pout), np.cos(2 * np.pi * tyr / Pout)]
    return np.array(cols).T
def wls(A, y, e):
    Aw = A / e[:, None]; c, *_ = np.linalg.lstsq(Aw, y / e, rcond=None)
    r = (y - A @ c) / e; cov = np.linalg.inv(Aw.T @ Aw)
    return c, float(r @ r), cov
n = len(oc)
cl, chil, covl = wls(design('lin'), oc, es); cq, chiq, covq = wls(design('quad'), oc, es); cc, chic, covc = wls(design('cubic'), oc, es)
out = dict(cal={f'{k[0]}_{k[1]}': v for k, v in cal.items()}, n_events=int(n), bias_corrected=BIAS)
Pn = P + cl[1] / 1000 / 1440; Pe = np.sqrt(covl[1, 1]) / 1000 / 1440
T0n = T0 + cl[0] / 1440; T0e = np.sqrt(covl[0, 0]) / 1440
out['linear'] = dict(chi2=chil, dof=n - 3, T0=T0n, T0_err_d=T0e, P=Pn, P_err_d=Pe, sec_offset_min=float(cl[2]), sec_offset_err_min=float(np.sqrt(covl[2, 2])))
out['quad'] = dict(chi2=chiq, c2_min_per_yr2=float(cq[3]), c2_err=float(np.sqrt(covq[3, 3])), dchi2=chil - chiq)
out['cubic'] = dict(chi2=chic, dchi2_vs_lin=chil - chic)
print(f'linear: chi2 {chil:.1f} / {n-3} dof; T0 = {T0n:.6f} +- {T0e*1440:.3f} min (BTJD); P = {Pn:.8f} +- {Pe:.8f} d; secondary offset {cl[2]:+.2f} +- {np.sqrt(covl[2,2]):.2f} min')
print(f'quadratic: c2 = {cq[3]:+.4f} +- {np.sqrt(covq[3,3]):.4f} min/yr^2 (tref BTJD {TREF}); dchi2 = {chil-chiq:.2f};  cubic dchi2 vs lin = {chil-chic:.2f}')
res = oc - design('lin') @ cl
ps = []
for s in np.unique(sec):
    for f in (0.0, 0.5, None):
        m = (sec == s) if f is None else (sec == s) & (ph == f)
        if m.sum() == 0: continue
        w = 1 / es[m] ** 2; mu = np.sum(res[m] * w) / w.sum(); er = 1 / np.sqrt(w.sum())
        chi_in = float(np.sum((res[m] - mu) ** 2 * w))
        ps.append(dict(sector=int(s), ph=('all' if f is None else f), n=int(m.sum()), t=float(np.mean(t[m])), mean_res_min=float(mu), err_min=float(er),
                       chi2_within=chi_in, dof=int(m.sum() - 1)))
        print(f"S{s:3d} {('both' if f is None else ('prim' if f == 0 else 'sec ')):4s}: N={m.sum():2d} mean O-C {mu:+.3f} +- {er:.3f} min; within-sector chi2 {chi_in:.1f}/{m.sum()-1}")
out['per_sector'] = ps
psc = np.array([(p_['t'], p_['mean_res_min'], p_['err_min']) for p_ in ps if p_['ph'] == 'all'])
out['sector_means_chi2'] = dict(chi2=float(np.sum((psc[:, 1] / psc[:, 2]) ** 2)), n=len(psc))
print('chi2 of the 6 combined sector means about the linear ephemeris: %.2f (N=%d)' % (out['sector_means_chi2']['chi2'], len(psc)))
rng = np.random.default_rng(3)
Al = design('lin'); Aq = design('quad')
model_l = Al @ cl; zres = (oc - model_l) / es
groups = [np.where((sec == s) & (ph == f))[0] for s in np.unique(sec) for f in (0.0, 0.5)]
groups = [g for g in groups if g.size]
# null pool: normalised residuals with each (sector, type) group's own mean removed (so that no sector-level offset -- which
# is what an ETV signal would look like -- survives into the null), rescaled by sqrt(n/(n-1)); drawn independently per event
zpool = np.concatenate([(zres[g] - zres[g].mean()) * np.sqrt(g.size / (g.size - 1)) for g in groups if g.size > 1])
def boot(m=1):
    return model_l[:, None] + rng.choice(zpool, (len(oc), m), replace=True) * es[:, None]
NB = 4000; Y0 = boot(NB)
dq = np.empty(NB); c2n = np.empty(NB)
for b in range(NB):
    _, a, _ = wls(Al, Y0[:, b], es); cqb, bq, _ = wls(Aq, Y0[:, b], es); dq[b] = a - bq; c2n[b] = cqb[3]
out['quad_null'] = dict(p_dchi2=float(np.mean(dq >= out['quad']['dchi2'])), c2_null_sd=float(np.std(c2n)), dchi2_95=float(np.percentile(dq, 95)),
                        dchi2_997=float(np.percentile(dq, 99.7)), c2_95_twosided=float(np.percentile(np.abs(c2n), 95)))
print('quadratic null (bootstrap): p(dchi2 >= obs) = %.3f; sd(c2) = %.4f min/yr^2; dchi2 thresholds 95%% %.2f, 99.7%% %.2f' %
      (out['quad_null']['p_dchi2'], np.std(c2n), np.percentile(dq, 95), np.percentile(dq, 99.7)))
c2_lim = abs(cq[3]) + 2 * max(np.std(c2n), np.sqrt(covq[3, 3]))
out['quad']['c2_abs_95_limit'] = float(c2_lim)
print('|c2| < %.3f min/yr^2 (|c2_obs| + 2 sigma_emp)' % c2_lim)
Pg = np.exp(np.linspace(np.log(0.5), np.log(60), 300)); dl = []
for Po in Pg:
    c, ch, cov = wls(design('ltte', Po), oc, es); dl.append((chil - ch, np.hypot(c[3], c[4])))
dl = np.array(dl); j = np.argmax(dl[:, 0])
mx = []
for b in range(400):
    y = Y0[:, b]; _, a, _ = wls(Al, y, es)
    mx.append(max(a - wls(design('ltte', Po), y, es)[1] for Po in Pg[::3]))
out['ltte'] = dict(best_Pout_yr=float(Pg[j]), dchi2=float(dl[j, 0]), amp_min=float(dl[j, 1]), p_null=float(np.mean(np.array(mx) >= dl[j, 0])),
                   dchi2_997=float(np.percentile(mx, 99.7)))
print('LTTE sinusoid search 0.5-60 yr: best P %.2f yr, dchi2 %.2f, amplitude %.2f min; null p = %.3f (99.7%% threshold %.2f)' %
      (Pg[j], dl[j, 0], dl[j, 1], out['ltte']['p_null'], out['ltte']['dchi2_997']))
pc = []
for c2inj in [0.025, 0.05, 0.1, 0.2, 0.4]:
    Yi = boot(400) + c2inj * tyr[:, None] ** 2
    rec = []
    for b in range(400):
        cqb, bq, _ = wls(Aq, Yi[:, b], es); _, a, _ = wls(Al, Yi[:, b], es); rec.append((cqb[3], a - bq))
    rec = np.array(rec)
    span = tyr.max() - tyr.min()
    pc.append(dict(c2_inj=c2inj, c2_rec_med=float(np.median(rec[:, 0])), c2_rec_sd=float(np.std(rec[:, 0])),
                   frac_detected_997=float(np.mean(rec[:, 1] > out['quad_null']['dchi2_997'])),
                   sagitta_min=float(c2inj * span ** 2 / 8)))
    print('positive control (timing level): c2 = %.3f min/yr^2 (sagitta over %.1f yr = %.2f min) -> recovered %.4f +- %.4f; detected at 99.7%%: %.2f' %
          (c2inj, span, c2inj * span ** 2 / 8, np.median(rec[:, 0]), np.std(rec[:, 0]), pc[-1]['frac_detected_997']))
out['positive_control_timing'] = pc
json.dump(dict(out=out, E=E.tolist(), oc=oc.tolist(), es=es.tolist(), ef=ef.tolist(), sec=sec.tolist(), ph=ph.tolist(), t=t.tolist(),
               lin=cl.tolist(), quad=cq.tolist(), tref=TREF), open(f'{FD}/results/etvfit_{tag}_b{BIAS}.json', 'w'), indent=1, default=float)
fig, ax = plt.subplots(2, 1, figsize=(9, 7.5))
for f, col, lab in [(0.0, 'C0', 'primary (per event)'), (0.5, 'C3', 'secondary (per event)')]:
    m = ph == f
    ax[0].errorbar(t[m] + (0 if f == 0 else 8), res[m], es[m], fmt='.', color=col, alpha=0.35, label=lab)
for p_ in ps:
    if p_['ph'] == 'all':
        ax[0].errorbar(p_['t'], p_['mean_res_min'], p_['err_min'], fmt='s', color='k', ms=6, zorder=5)
tt = np.linspace(t.min(), t.max(), 300); ty = (tt - TREF) / 365.25
ax[0].plot(tt, (cq[0] - cl[0]) + (cq[1] - cl[1]) * ((tt - T0) / P) / 1000 + cq[3] * ty ** 2, 'g-', label='best quadratic')
ax[0].axhline(0, color='0.5'); ax[0].set_ylabel('O-C [min] (linear ephemeris removed)'); ax[0].set_xlabel('BTJD')
ax[0].legend(fontsize=8); ax[0].set_ylim(-12, 12); ax[0].set_title('TIC 71316387 EB: eclipse timings 2018-2026 (squares: sector means)')
ax[1].semilogx(Pg, dl[:, 0], 'k-'); ax[1].axhline(out['ltte']['dchi2_997'], color='r', ls=':', label='99.7% bootstrap null')
ax[1].set_xlabel('outer period [yr]'); ax[1].set_ylabel('Delta chi2 (sinusoid vs linear)'); ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig(f'{FD}/plots/etv_{tag}_b{BIAS}.png', dpi=110)
