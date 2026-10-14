"""ETV analysis from etv_<tag>.json: empirical per-event errors (injection-calibrated), outlier rejection,
linear / quadratic / LTTE fits, empirical null for the quadratic term, positive control (timing-level injection),
and outer-companion constraints (LTTE from a third body of mass M3 on a circular orbit)."""
import numpy as np, json, sys
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
HERE = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel'
tag = sys.argv[1] if len(sys.argv) > 1 else 'default'
D = json.load(open(f'{HERE}/results/etv_{tag}.json'))
P, T0 = D['P'], D['T0']
rows = [r for r in D['rows'] if not r['edge']]
inj = [r for r in D['inj'] if not r['edge']]
# ---- 1. empirical error scaling per sector & eclipse type
cal = {}
for sec in sorted(set(r['sector'] for r in rows)):
    for ph in (0.0, 0.5):
        ii = [r for r in inj if r['sector'] == sec and r['kind'] == ph]
        x = np.array([r['dt_min'] for r in ii]); e = np.array([r['err_min'] for r in ii])
        z = x / e
        rsd = 1.4826 * np.median(np.abs(z - np.median(z)))
        cal[(sec, ph)] = dict(n=len(ii), bias_min=float(np.median(x)), sd_min=float(1.4826 * np.median(np.abs(x - np.median(x)))),
                              zscale=float(max(1.0, rsd)), med_formal=float(np.median(e)))
E = np.array([r['n'] + r['ph'] for r in rows]); oc = np.array([r['dt_min'] for r in rows]); ef = np.array([r['err_min'] for r in rows])
sec = np.array([r['sector'] for r in rows]); ph = np.array([r['ph'] for r in rows])
es = np.array([ef[i] * cal[(sec[i], ph[i])]['zscale'] for i in range(len(rows))])
t = T0 + E * P
# outliers: > 5 sigma (empirical) from sector-type median
keep = np.ones(len(rows), bool)
for s in np.unique(sec):
    for f in (0.0, 0.5):
        m = (sec == s) & (ph == f)
        if m.sum() < 3: continue
        med = np.median(oc[m]); keep[m] &= np.abs(oc[m] - med) < 5 * es[m]
print('events', len(rows), 'kept', keep.sum())
E, oc, ef, es, sec, ph, t = E[keep], oc[keep], ef[keep], es[keep], sec[keep], ph[keep], t[keep]
tyr = (t - 2900.0) / 365.25    # years from ~2022.7 (centre of baseline)
def design(kind, Pout=None, extra=None):
    cols = [np.ones_like(E), E / 1000.0, (ph == 0.5).astype(float)]
    if kind in ('quad', 'cubic'): cols.append(tyr ** 2)
    if kind == 'cubic': cols.append(tyr ** 3)
    if kind == 'ltte':
        cols += [np.sin(2 * np.pi * tyr / Pout), np.cos(2 * np.pi * tyr / Pout)]
    return np.array(cols).T
def wls(A, y, e):
    Aw = A / e[:, None]; c, *_ = np.linalg.lstsq(Aw, y / e, rcond=None)
    r = (y - A @ c) / e; cov = np.linalg.inv(Aw.T @ Aw)
    return c, float(r @ r), cov
out = dict(cal={f'{k[0]}_{k[1]}': v for k, v in cal.items()}, n_events=int(keep.sum()))
cl, chil, covl = wls(design('lin'), oc, es)
cq, chiq, covq = wls(design('quad'), oc, es)
cc, chic, covc = wls(design('cubic'), oc, es)
n = len(oc)
out['linear'] = dict(chi2=chil, dof=n - 3, T0_shift_min=cl[0], dP_s=cl[1] / 1000 * 60, sec_offset_min=cl[2], sec_offset_err_min=float(np.sqrt(covl[2, 2])),
                     P_new=P + cl[1] / 1000 / 1440, P_err_d=float(np.sqrt(covl[1, 1]) / 1000 / 1440))
out['quad'] = dict(chi2=chiq, c2_min_per_yr2=float(cq[3]), c2_err=float(np.sqrt(covq[3, 3])), dchi2=chil - chiq, dBIC=chil - chiq - np.log(n))
out['cubic'] = dict(chi2=chic, dchi2_vs_lin=chil - chic)
print('linear chi2 %.1f dof %d; secondary offset %.2f +- %.2f min; P = %.8f +- %.8f' % (chil, n - 3, cl[2], np.sqrt(covl[2, 2]), out['linear']['P_new'], out['linear']['P_err_d']))
print('quadratic c2 = %.3f +- %.3f min/yr^2, dchi2 = %.2f' % (cq[3], np.sqrt(covq[3, 3]), chil - chiq))
# per-sector means (primary & secondary separately), residual about linear fit
res = oc - design('lin') @ cl
ps = []
for s in np.unique(sec):
    for f in (0.0, 0.5):
        m = (sec == s) & (ph == f)
        if m.sum() == 0: continue
        w = 1 / es[m] ** 2; mu = np.sum(res[m] * w) / w.sum(); er = 1 / np.sqrt(w.sum())
        chi_in = float(np.sum((res[m] - mu) ** 2 * w))
        ps.append(dict(sector=int(s), ph=f, n=int(m.sum()), t=float(np.mean(t[m])), mean_res_min=float(mu), err_min=float(er), chi2_within=chi_in, dof=int(m.sum() - 1)))
        print(f"S{s} ph{f}: N={m.sum():2d} mean resid {mu:+.3f} +- {er:.3f} min   within-sector chi2/dof {chi_in:.1f}/{m.sum()-1}")
out['per_sector'] = ps
# combined per-sector (prim+sec) means
psc = []
for s in np.unique(sec):
    m = sec == s; w = 1 / es[m] ** 2
    psc.append((float(np.mean(t[m])), float(np.sum(res[m] * w) / w.sum()), float(1 / np.sqrt(w.sum()))))
psc = np.array(psc)
chi_sec = float(np.sum((psc[:, 1] / psc[:, 2]) ** 2))
out['sector_means_chi2'] = dict(chi2=chi_sec, n=len(psc))
print('chi2 of 6 sector means about linear ephemeris: %.2f' % chi_sec)
# ---- 2. empirical null for quadratic term: parametric bootstrap with residuals resampled within sector
rng = np.random.default_rng(3)
Al = design('lin'); Aq = design('quad')
nb = 4000; dchi = np.empty(nb); c2n = np.empty(nb)
model_l = Al @ cl
resid_by = {s: (oc - model_l)[sec == s] / es[sec == s] for s in np.unique(sec)}
for b in range(nb):
    y = model_l.copy()
    for s in np.unique(sec):
        m = sec == s
        y[m] += rng.choice(resid_by[s], m.sum(), replace=True) * es[m]
    _, a, _ = wls(Al, y, es); cqb, bq, _ = wls(Aq, y, es); dchi[b] = a - bq; c2n[b] = cqb[3]
out['quad_null'] = dict(p_dchi2=float(np.mean(dchi >= out['quad']['dchi2'])), c2_null_sd=float(np.std(c2n)), dchi2_95=float(np.percentile(dchi, 95)),
                        dchi2_997=float(np.percentile(dchi, 99.7)))
print('null: p(dchi2 >= obs) = %.3f ; sd(c2) = %.3f min/yr^2 ; 95%% / 99.7%% dchi2 = %.2f / %.2f' % (out['quad_null']['p_dchi2'], np.std(c2n), np.percentile(dchi, 95), np.percentile(dchi, 99.7)))
# ---- 3. LTTE periodogram (circular), Pout 1-60 yr
Pg = np.exp(np.linspace(np.log(0.5), np.log(60), 400)); dl = []
for Po in Pg:
    c, ch, cov = wls(design('ltte', Po), oc, es); dl.append((chil - ch, np.hypot(c[3], c[4])))
dl = np.array(dl); j = np.argmax(dl[:, 0])
out['ltte'] = dict(best_Pout_yr=float(Pg[j]), dchi2=float(dl[j, 0]), amp_min=float(dl[j, 1]))
# null for max LTTE dchi2
mx = []
for b in range(300):
    y = model_l.copy()
    for s in np.unique(sec):
        m = sec == s; y[m] += rng.choice(resid_by[s], m.sum(), replace=True) * es[m]
    _, a, _ = wls(Al, y, es)
    mx.append(max(a - wls(design('ltte', Po), y, es)[1] for Po in Pg[::4]))
out['ltte']['p_null'] = float(np.mean(np.array(mx) >= dl[j, 0]))
print('LTTE best Pout %.2f yr dchi2 %.2f amp %.2f min, null p=%.3f' % (Pg[j], dl[j, 0], dl[j, 1], out['ltte']['p_null']))
# ---- 4. positive control: inject quadratic O-C into the timings, recover
pc = []
for c2inj in [0.05, 0.1, 0.2, 0.4, 0.8]:
    rec = []
    for b in range(300):
        y = model_l + c2inj * tyr ** 2
        for s in np.unique(sec):
            m = sec == s; y[m] += rng.choice(resid_by[s], m.sum(), replace=True) * es[m]
        cqb, bq, covb = wls(Aq, y, es); _, a, _ = wls(Al, y, es)
        rec.append((cqb[3], a - bq))
    rec = np.array(rec)
    pc.append(dict(c2_inj=c2inj, c2_rec_med=float(np.median(rec[:, 0])), c2_rec_sd=float(np.std(rec[:, 0])),
                   frac_detected_997=float(np.mean(rec[:, 1] > out['quad_null']['dchi2_997'])),
                   curvature_over_baseline_min=float(c2inj * (tyr.max() - tyr.min()) ** 2 / 4)))
    print('inject c2=%.2f min/yr^2 (peak-to-peak curvature after linear removal ~%.2f min): recovered %.3f +- %.3f, detected(99.7%%) %.2f' %
          (c2inj, c2inj * (tyr.max() - tyr.min()) ** 2 / 8, np.median(rec[:, 0]), np.std(rec[:, 0]), pc[-1]['frac_detected_997']))
out['positive_control_timing'] = pc
json.dump(dict(out=out, E=E.tolist(), oc=oc.tolist(), es=es.tolist(), ef=ef.tolist(), sec=sec.tolist(), ph=ph.tolist(), t=t.tolist(),
               lin=cl.tolist(), quad=cq.tolist()), open(f'{HERE}/results/etvfit_{tag}.json', 'w'), indent=1, default=float)
# ---- plot
fig, ax = plt.subplots(2, 1, figsize=(9, 7), sharex=False)
for f, col, lab in [(0.0, 'C0', 'primary'), (0.5, 'C3', 'secondary')]:
    m = ph == f
    ax[0].errorbar(t[m], res[m], es[m], fmt='.', color=col, alpha=0.4, label=lab)
for p_ in ps:
    ax[0].errorbar(p_['t'], p_['mean_res_min'], p_['err_min'], fmt='s', color='k' if p_['ph'] == 0 else 'm', ms=6, zorder=5)
tt = np.linspace(t.min(), t.max(), 300); ty = (tt - 2900) / 365.25
ax[0].plot(tt, (cq[0] - cl[0]) + (cq[1] - cl[1]) * ((tt - T0) / P) / 1000 + cq[3] * ty ** 2, 'g-', label='quadratic fit')
ax[0].axhline(0, color='0.5'); ax[0].set_ylabel('O-C (min) after linear ephemeris'); ax[0].set_xlabel('BTJD'); ax[0].legend(fontsize=8)
ax[0].set_ylim(-8, 8)
ax[1].semilogx(Pg, dl[:, 0]); ax[1].set_xlabel('outer period (yr)'); ax[1].set_ylabel('LTTE delta chi2 vs linear')
plt.tight_layout(); plt.savefig(f'{HERE}/plots/etv_{tag}.png', dpi=110)
