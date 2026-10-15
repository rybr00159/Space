"""Combine: measured depth ratios depth(band)/depth(TESS) vs hypotheses.

Measured ratio R = A_band x (1 + D_band) / A_TESS, with
  stat  : empirical null sd (max of fake-phase sd and fake-period beta-scaled error), deep-shallow correlation
          from the fake-phase null;
  sys   : sd of reasonable detrending/cleaning variants (survey-specific list below), dilution error,
          TESS reference error, template-shape (3%) and physical-shape (4%) terms, and an epoch term for
          possible depth evolution between the survey epoch and 2018-2026 (8% for 2000-2009 data, 3% for
          NEOWISE which overlaps TESS).
Hypotheses
  (a) bound K/M pair at 297 pc (predict3 pairs); 'TESS-consistent' = predicted TESS shallow depth within
      20% of the observed one;
  (b) background EB with a G/F eclipsed star: Teff1 in [5300, 7000] K, grazing/mid/deep geometry,
      secondary Teff from the TESS shallow/deep ratio (surface-brightness ratio 0.454);
  (c) generic eclipsed star of any Teff (small light fraction): Teff inference.
chi2 per dataset uses the 2x2 (deep, shallow) covariance; the background is given its most favourable
Teff/geometry (minimum chi2), as is the bound pair.
"""
import sys, json
import numpy as np
from scipy.stats import chi2 as CHI2, norm
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/v3/scripts')
import lib3 as L
import sed3 as S

R = lambda f: json.load(open(L.V3 + '/results/' + f))
T3, P3, W3, WD, X3, SP, DL, N3 = (R('tess3.json'), R('predict3.json'), R('wasp3.json'), R('wasp3_diag.json'), R('wise3.json'),
                                  R('sparse3.json'), R('dilution3.json'), R('null3.json'))
WE = R('wasp3_e2e.json')
AT = np.array([T3['adopted_deep']['mean'], T3['adopted_shallow']['mean']])
ATe = np.array([T3['adopted_deep']['err'], T3['adopted_shallow']['err']])
out = dict(TESS=dict(deep=AT[0], shallow=AT[1], deep_err=ATe[0], shallow_err=ATe[1]))

# ------------------------------------------------------------------ measurements
meas = {}


PASSBAND = {'WASP': 0.05, 'V': 0.02, 'W1': 0.03, 'W2': 0.03}   # prediction-side passband/LD-proxy uncertainty


def add(name, band, A, stat, corr, sys_abs, D, D_err, epoch_frac, note):
    A = np.asarray(A, float); stat = np.asarray(stat, float); sys_abs = np.asarray(sys_abs, float)
    fac = (1 + D) / AT
    r = A * fac
    # covariance: stat (correlated deep/shallow) + detrending sys (assume same correlation) + multiplicative terms
    C = np.outer(stat * fac, stat * fac) * np.array([[1, corr], [corr, 1]])
    C += np.outer(sys_abs * fac, sys_abs * fac) * np.array([[1, corr], [corr, 1]])
    mult = np.array([ATe[0] / AT[0], ATe[1] / AT[1]])
    frac2 = (D_err / (1 + D)) ** 2 + 0.03 ** 2 + 0.04 ** 2 + epoch_frac ** 2 + PASSBAND[band] ** 2
    # multiplicative terms: fully correlated between deep and shallow except the TESS reference errors
    C += np.outer(r, r) * frac2 + np.diag((r * mult) ** 2)
    meas[name] = dict(band=band, A=A.tolist(), stat=stat.tolist(), sys=sys_abs.tolist(), corr=corr, D=D, ratio=r.tolist(),
                      ratio_err=np.sqrt(np.diag(C)).tolist(), cov=C.tolist(), note=note)
    print('%-12s %-5s A_d %6.0f +- %5.0f (stat) +- %5.0f (sys) | A_s %6.0f +- %5.0f +- %5.0f | corr %.2f D %.3f | R_d %.3f +- %.3f  R_s %.3f +- %.3f' % (
        name, band, A[0] * 1e6, stat[0] * 1e6, sys_abs[0] * 1e6, A[1] * 1e6, stat[1] * 1e6, sys_abs[1] * 1e6, corr, D, r[0], np.sqrt(C[0, 0]), r[1], np.sqrt(C[1, 1])))


# SuperWASP
d = W3['default']
var = W3['variants']
vd = [var[k]['deep'] for k in ('TAM_o0', 'TAM_o2', 'TAM_clip3', 'TAM_clip6', 'TAM_minpts40', 'RAW_o1', 'RAW_o2', 'TAM_trim_night_ends')] + \
     [WD['airmass_o1']['deep'], WD['airmass_o2']['deep'], d['deep']]
vs = [var[k]['shallow'] for k in ('TAM_o0', 'TAM_o2', 'TAM_clip3', 'TAM_clip6', 'TAM_minpts40', 'RAW_o1', 'RAW_o2', 'TAM_trim_night_ends')] + \
     [WD['airmass_o1']['shallow'], WD['airmass_o2']['shallow'], d['shallow']]
stat_w = [max(d['deep_err'], N3['WASP']['err_scaled'][0]), max(d['shallow_err'], N3['WASP']['err_scaled'][1])]
out['WASP_variants'] = dict(deep=vd, shallow=vs, deep_sd=float(np.std(vd, ddof=1)), shallow_sd=float(np.std(vs, ddof=1)))
add('WASP', 'WASP', [d['deep'], d['shallow']], stat_w, d['null_corr'], [np.std(vd, ddof=1), np.std(vs, ddof=1)],
    DL['wasp_D']['adopted'], DL['wasp_D']['err'], 0.08, 'SuperWASP DR1 2006-2007, TAMFLUX2, linear night baselines')
# NEOWISE W1/W2 (strict cuts + AllWISE)
for b in ('W1', 'W2'):
    r = X3[f'{b}_strict_allwise']
    keys = [f'{b}_{c}_{a}' for c in ('strict', 'loose', 'none') for a in ('allwise', 'neowiseonly')]
    vd = [X3[k]['deep'] for k in keys]
    vs = [X3[k]['shallow'] for k in keys if 'strict' in k]
    stat = [max(r['deep_err'], N3[b]['err_scaled'][0]), max(r['shallow_err'], N3[b]['err_scaled'][1])]
    add('NEOWISE_' + b, b, [r['deep'], r['shallow']], stat, r['null_corr'], [0.5 * (max(vd) - min(vd)), 0.5 * (max(vs) - min(vs))],
        0.0, 0.01, 0.03, 'NEOWISE-R 2014-2024 + AllWISE 2010 single exposures, strict quality cuts')
# ASAS-SN V and ASAS-3 V (precision-limited)
r = SP['ASASSN_V_all']
add('ASASSN_V', 'V', [r['deep'], r['shallow']], [r['deep_err'], r['shallow_err']], r['null_corr'],
    [abs(r['deep'] - SP['ASASSN_V_all_o0']['deep']), abs(r['shallow'] - SP['ASASSN_V_all_o0']['shallow'])],
    DL['asassn_D']['adopted'], DL['asassn_D']['err'], 0.03, 'ASAS-SN photometry database V 2014-2018 (220 epochs)')
r = SP['ASAS3_V']
add('ASAS3_V', 'V', [r['deep'], r['shallow']], [r['deep_err'], r['shallow_err']], r['null_corr'], [0, 0],
    0.0065, 0.006, 0.08, 'ASAS-3 V 2001-2009, aperture MAG_2')
out['measurements'] = meas

# ------------------------------------------------------------------ hypotheses
bands_needed = sorted(set(m['band'] for m in meas.values()))
pairs = {}
for k, p in P3['pairs'].items():
    if not p.get('feasible'):
        continue
    cons = abs(p['TESS_shallow_pred'] / AT[1] - 1) < 0.20
    pairs[k] = dict(cons=bool(cons), sT=p['TESS_shallow_pred'], pred={b: (p['ratio_deep'][b], p['ratio_shallow'][b]) for b in bands_needed},
                    fpair=p['pair_frac_TESS'], dKs=p['dKs'], dI=p['dI'], dV=p['dV'])
# generic ratio(Teff) for each band and geometry (log-interpolated)
bg = P3['background']
TE = np.array(sorted(int(k) for k in bg))
GEO = ('grazing', 'mid', 'deep')


def generic(teff, band, geo):
    y = np.array([bg[str(t)][geo][band] for t in TE])
    return float(np.exp(np.interp(teff, TE, np.log(y))))


# secondary Teff for a background EB from the TESS surface-brightness ratio s/d
tab = S.mamajek()
_sb = {}


def sbT(teff):
    k = int(round(teff, -1))
    if k not in _sb:
        _sb[k] = S.band_flux(S.star_model(teff=k, tab=tab)['spec'], 'TESS')
    return _sb[k]


sod = AT[1] / AT[0]


def teff2_for(teff1):
    g = np.arange(2900, teff1, 20.0)
    v = np.array([sbT(x) / sbT(teff1) for x in g])
    return float(np.interp(sod, v, g))


bg_models = {}
for t1 in (5300, 5600, 5900, 6200, 6600, 7000):
    t2 = teff2_for(t1)
    for geo in GEO:
        bg_models[f'{t1}_{geo}'] = dict(Teff1=t1, Teff2=t2, pred={b: (generic(t1, b, geo), generic(t2, b, geo)) for b in bands_needed})
print('background EB secondary Teff from TESS s/d=%.3f: %s' % (sod, {t: round(teff2_for(t)) for t in (5300, 5900, 7000)}))


def chi2_of(m, pred):
    r = np.array(m['ratio']); C = np.array(m['cov'])
    dlt = r - np.array(pred)
    return float(dlt @ np.linalg.solve(C, dlt))


GF_STRICT = [k for k in bg_models if bg_models[k]['Teff1'] >= 5600]
res = {}
for name, m in meas.items():
    b = m['band']
    cb = {k: chi2_of(m, p['pred'][b]) for k, p in pairs.items()}
    cbc = {k: v for k, v in cb.items() if pairs[k]['cons']}
    cg = {k: chi2_of(m, v['pred'][b]) for k, v in bg_models.items()}
    kb = min(cbc, key=cbc.get); kg = min(cg, key=cg.get)
    # deep-only (1 dof) version
    r = np.array(m['ratio']); C = np.array(m['cov'])
    cbd = {k: (r[0] - pairs[k]['pred'][b][0]) ** 2 / C[0, 0] for k in cbc}
    cgd = {k: (r[0] - bg_models[k]['pred'][b][0]) ** 2 / C[0, 0] for k in cg}
    kbd = min(cbd, key=cbd.get); kgd = min(cgd, key=cgd.get)
    zg_d = (r[0] - bg_models[kgd]['pred'][b][0]) / np.sqrt(C[0, 0])
    kg2 = min(GF_STRICT, key=cg.get)
    kgd2 = min(GF_STRICT, key=cgd.get)
    zg_d2 = (r[0] - bg_models[kgd2]['pred'][b][0]) / np.sqrt(C[0, 0])
    res[name] = dict(best_bound=kb, chi2_bound=cbc[kb], chi2_bound_all=cbc, best_bg=kg, chi2_bg=cg[kg], dchi2=cg[kg] - cbc[kb],
                     best_bg_ge5600=kg2, chi2_bg_ge5600=cg[kg2], dchi2_ge5600=cg[kg2] - cbc[kb], z_bg_deep_ge5600=float(zg_d2),
                     p_bound=float(CHI2.sf(cbc[kb], 2)), p_bg=float(CHI2.sf(cg[kg], 2)),
                     deep_only=dict(best_bound=kbd, chi2_bound=cbd[kbd], best_bg=kgd, chi2_bg=cgd[kgd], z_bg=float(zg_d),
                                    z_bound=float((r[0] - pairs[kbd]['pred'][b][0]) / np.sqrt(C[0, 0]))))
    print('%-12s bound best %-8s chi2 %6.2f (p %.3g) | G/F bg best %-12s chi2 %6.2f (p %.3g) | dchi2 %.1f | deep-only: z_bound %+.2f (%s) z_bg %+.2f (%s)' % (
        name, kb, cbc[kb], res[name]['p_bound'], kg, cg[kg], res[name]['p_bg'], cg[kg] - cbc[kb],
        res[name]['deep_only']['z_bound'], kbd, zg_d, kgd))
    print('             (G/F restricted to Teff>=5600 K: best %s chi2 %.2f, dchi2 %.1f, deep-only z %+.2f)' % (kg2, cg[kg2], cg[kg2] - cbc[kb], zg_d2))
out['tests'] = res

# joint over the independent datasets (WASP, NEOWISE W1+W2 treated with their deep-deep correlation, ASAS-SN, ASAS-3)
wcorr = float(np.corrcoef(X3['W1_strict_allwise_null_deep'], X3['W2_strict_allwise_null_deep'])[0, 1])
out['W1_W2_null_corr_deep'] = wcorr
names = list(meas)
vec = np.concatenate([meas[n]['ratio'] for n in names])
Cbig = np.zeros((len(vec), len(vec)))
for i, n in enumerate(names):
    Cbig[2 * i:2 * i + 2, 2 * i:2 * i + 2] = meas[n]['cov']
iw1, iw2 = names.index('NEOWISE_W1'), names.index('NEOWISE_W2')
for j in (0, 1):
    a, bb = 2 * iw1 + j, 2 * iw2 + j
    Cbig[a, bb] = Cbig[bb, a] = wcorr * np.sqrt(Cbig[a, a] * Cbig[bb, bb])


def joint(predfun):
    pv = np.concatenate([predfun(meas[n]['band']) for n in names])
    dl = vec - pv
    return float(dl @ np.linalg.solve(Cbig, dl))


jb = {k: joint(lambda b, k=k: pairs[k]['pred'][b]) for k in pairs}
jg = {k: joint(lambda b, k=k: bg_models[k]['pred'][b]) for k in bg_models}
kb = min({k: v for k, v in jb.items() if pairs[k]['cons']}, key=jb.get)
kg = min(jg, key=jg.get)
kg2 = min(GF_STRICT, key=jg.get)
ndof = len(vec)
out['joint'] = dict(datasets=names, ndata=ndof, best_bound=kb, chi2_bound=jb[kb], best_bg=kg, chi2_bg=jg[kg], dchi2=jg[kg] - jb[kb],
                    chi2_bound_all=jb, chi2_bg_all=jg, p_bound=float(CHI2.sf(jb[kb], ndof)), p_bg=float(CHI2.sf(jg[kg], ndof)),
                    best_bg_ge5600=kg2, chi2_bg_ge5600=jg[kg2], dchi2_ge5600=jg[kg2] - jb[kb], p_bg_ge5600=float(CHI2.sf(jg[kg2], ndof)))
print('JOINT (%d numbers): bound best %s chi2 %.2f (p %.3g); G/F background best %s chi2 %.2f (p %.3g); dchi2 = %.1f (likelihood ratio %.2g)' % (
    ndof, kb, jb[kb], out['joint']['p_bound'], kg, jg[kg], out['joint']['p_bg'], jg[kg] - jb[kb], np.exp(0.5 * (jg[kg] - jb[kb]))))
print('  G/F restricted to >=5600 K: best %s chi2 %.2f (p %.3g); dchi2 %.1f' % (kg2, jg[kg2], out['joint']['p_bg_ge5600'], jg[kg2] - jb[kb]))
print('  bound pairs joint chi2:', {k: round(v, 1) for k, v in sorted(jb.items(), key=lambda x: x[1])})

emp = {}
for name, key, wd in (('WASP', 'WASP', W3['default']), ('NEOWISE_W1', 'W1', X3['W1_strict_allwise'])):
    smp = np.array(N3[key]['samples'])
    m = meas[name]; b = m['band']
    fac = (1 + m['D']) / AT[0]
    zs = smp[:, 0] / np.std(smp[:, 0])         # standardised fake-period null
    for lab, keys in (('ge5300', list(bg_models)), ('ge5600', GF_STRICT)):
        # most favourable background prediction for this deep amplitude
        preds = np.array([bg_models[k]['pred'][b][0] for k in keys]) / fac
        pbest = preds[np.argmin(np.abs(preds - m['A'][0]))]
        tot = np.sqrt(np.array(m['cov'])[0, 0]) / fac
        zobs = (m['A'][0] - pbest) / tot
        tail = float((np.sum(zs <= zobs) + 1) / (len(zs) + 1)) if zobs < 0 else float((np.sum(zs >= zobs) + 1) / (len(zs) + 1))
        emp[f'{name}_{lab}'] = dict(pred_amp=float(pbest), z=float(zobs), empirical_tail=tail, gauss_tail=float(norm.sf(abs(zobs))))
        print('%-11s vs most favourable G/F (%s): predicted amplitude %.0f ppm, z = %+.2f (total error), empirical fake-period tail %.2g (Gaussian %.2g)' % (
            name, lab, pbest * 1e6, zobs, tail, norm.sf(abs(zobs))))
out['empirical_tails'] = emp

# ------------------------------------------------------------------ Teff of the eclipsed (deep-eclipse) star, generic model
tg = np.arange(3300, 7001, 25.0)
post = {}
for name in ('WASP', 'NEOWISE_W1', 'NEOWISE_W2'):
    m = meas[name]; b = m['band']; r = m['ratio'][0]; e = m['ratio_err'][0]
    like = np.zeros_like(tg)
    for geo in GEO:
        pr = np.array([generic(t, b, geo) for t in tg])
        like += np.exp(-0.5 * ((r - pr) / e) ** 2) / 3
    post[name] = like
comb = post['WASP'] * post['NEOWISE_W1'] * post['NEOWISE_W2']
tinf = {}
for name, lk in list(post.items()) + [('combined', comb)]:
    c = np.cumsum(lk); c /= c[-1]
    q = np.interp([0.025, 0.16, 0.5, 0.84, 0.975], c, tg)
    p_gf = float(lk[tg >= 5300].sum() / lk.sum())
    tinf[name] = dict(q=q.tolist(), P_Teff_ge_5300=p_gf)
    print('Teff(eclipsed star) from %-11s: median %.0f K, 68%% [%.0f, %.0f], 95%% [%.0f, %.0f]; posterior mass at Teff>=5300 K (G/F): %.2g' % (
        name, q[2], q[1], q[3], q[0], q[4], p_gf))
out['teff_inference'] = tinf
out['pairs'] = pairs
out['bg_models'] = bg_models
json.dump(out, open(L.V3 + '/results/combine3.json', 'w'), indent=1)
