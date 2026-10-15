"""Eclipse depth scale per sector/product and eclipse type, relative to the best-fit template (shape fixed).
For each product: y = baseline_e(dt) * (1 - s * dep(t)), with s common to all events of that product and eclipse
type, per-event quadratic baselines profiled; chi2(s) on a grid -> s and formal error.  Empirical error: the same
measurement on fake-phase (0.25/0.75) injections of the template at s = 1 into the real light curves; the robust sd
of recovered s over 200 bootstrap resamplings of the injected events (and the scatter between halves) gives the
error.  Also QLP SAP (S95, S105) and own TESScut S1/S28/S68 for product comparison.
Usage: c06_depths.py <best json>"""
import numpy as np, sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/v5/lib')
import data as D, etv
V = D.V
best = json.load(open(sys.argv[1])); p = best['p']
tm = etv.Template(p, {})   # unscaled template
T0, P = p['T0'], p['P']
summ = json.load(open(f'{V}/results/data_summary.json'))
sig_of = {}
for ds in summ.values():
    for k, v in ds['info'].items(): sig_of[k] = v['sigma']
rng = np.random.default_rng(5)
SG = np.linspace(0.6, 1.4, 401)

def event_chi2_curves(lc, key, phases, inject=False):
    """returns list of (ph, chi2(s) array) per event"""
    sig = sig_of.get(key, lc['sig_pt'])
    out = []
    for w in D.windows(lc, P=P, T0=T0, phases=phases):
        i = w['idx']
        if inject:
            tcf = w['tc']; ph = 0.0 if w['ph'] == 0.25 else 0.5
            dep = tm.depth(lc['t'][i] - tcf, lc['texp'], ph, key)
            y = lc['f'][i] * (1 - dep)
        else:
            ph = w['ph']; tcf = T0 + w['n'] * P + (tm.dsec if ph == 0.5 else 0.0)
            dep = tm.depth(lc['t'][i] - tcf, lc['texp'], ph, key); y = lc['f'][i]
        dt = lc['t'][i] - tcf; wgt = np.full(i.size, 1 / sig ** 2)
        ch = np.array([etv._chi2_lin(y, 1 - s * dep, dt, wgt) for s in SG])
        out.append((ph, ch))
    return out

def best_s(curves):
    tot = np.sum(curves, axis=0)
    j = np.argmin(tot)
    if 2 <= j <= SG.size - 3:
        c = np.polyfit(SG[j - 2:j + 3], tot[j - 2:j + 3], 2); s0 = -c[1] / (2 * c[0]); e = 1 / np.sqrt(c[0])
    else:
        s0, e = SG[j], np.nan
    return float(s0), float(e)

prods = [('SPOC1', D.load_spoc(1)), ('SPOC28', D.load_spoc(28)), ('SPOC68', D.load_spoc(68)), ('TC95', D.load_tc(95)), ('TC105', D.load_tc(105)),
         ('TC106', D.load_tc(106)), ('TC1', D.load_tc(1)), ('TC28', D.load_tc(28)), ('TC68', D.load_tc(68)), ('QLP95', D.load_qlp(95)), ('QLP105', D.load_qlp(105)),
         ('QLP1', D.load_qlp(1)), ('QLP28', D.load_qlp(28)), ('QLP68', D.load_qlp(68))]
res = {}
for key, lc in prods:
    lc = D.clip(lc)
    real = event_chi2_curves(lc, key, (0.0, 0.5))
    fake = event_chi2_curves(lc, key, (0.25, 0.75), inject=True)
    r = {}
    for lab, sel in (('prim', 0.0), ('sec', 0.5), ('both', None)):
        cr = [c for ph, c in real if sel is None or ph == sel]
        cf = [c for ph, c in fake if sel is None or ph == sel]
        s, e = best_s(cr)
        sf, ef = best_s(cf)
        boots = []
        for b in range(200):
            idx = rng.integers(0, len(cf), len(cf)); boots.append(best_s([cf[i] for i in idx])[0])
        rb = []
        for b in range(200):
            idx = rng.integers(0, len(cr), len(cr)); rb.append(best_s([cr[i] for i in idx])[0])
        r[lab] = dict(s=s, err_formal=e, n=len(cr), inj_s=sf, inj_err_formal=ef, inj_boot_sd=float(np.std(boots)), real_boot_sd=float(np.std(rb)),
                      err_adopted=float(max(e, np.std(boots), np.std(rb))))
    res[key] = r
    print(f"{key:7s} prim s={r['prim']['s']:.3f}+-{r['prim']['err_adopted']:.3f} (formal {r['prim']['err_formal']:.3f}, inj-boot {r['prim']['inj_boot_sd']:.3f}, inj s {r['prim']['inj_s']:.3f}) | "
          f"sec s={r['sec']['s']:.3f}+-{r['sec']['err_adopted']:.3f} (inj s {r['sec']['inj_s']:.3f}) | both {r['both']['s']:.3f}+-{r['both']['err_adopted']:.3f}", flush=True)
json.dump(dict(template=best.get('source'), p=p, depths=res), open(f'{V}/results/depths_by_product.json', 'w'), indent=1)
