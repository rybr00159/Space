"""SuperWASP end-to-end injection-recovery (re-cleaned, paired) and template-shape systematics.

1. Inject (a) template-shaped and (b) physically shaped (K7V+M2V, WASP-band LD) eclipse pairs into the
   normalised TAMFLUX2 series at fake phases; re-run cleaning (with the injected phase modelled),
   weights and fit; paired gain = (A_inj - A_uninj)/A_expected.
2. Template-shape systematic: alternative TESS template (generalised V with p=2, q=1, T fitted) applied
   identically to TESS (all six sectors) and to WASP -> change in the WASP/TESS ratio.
"""
import sys, json
import numpy as np
from scipy.optimize import minimize
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/v3/scripts')
import lib3 as L
import wasp3 as W
import tess3 as TS
import predict3 as PR
import sed3 as S

T = json.load(open(L.V3 + '/results/tess3.json'))
E3 = json.load(open(L.V3 + '/results/ephem3.json'))
shp = L.survey_shapes(T['shapes'], E3)
sod = T['s_over_d']
EX = W.EXPT
out = {}
t, f, seg, cam, hjd = W.load(W.TARGET)
rng = np.random.default_rng(31)

# physical profile
P3 = json.load(open(L.V3 + '/results/predict3.json'))
pr = P3['pairs']['K7V+M2V']
tab = S.mamajek()
tg = np.linspace(-0.12, 0.12, 2401)
ld1 = S.ld_quadratic(tab['K7V']['Teff'], 4.65, 'WASP'); ld2 = S.ld_quadratic(tab['M2V']['Teff'], 4.8, 'WASP')
k1, k2 = PR.lightcurves(pr['R1'], pr['R2'], pr['a_Rsun'], pr['inc'], ld1, ld2, tg)
prof = (tg, k1 / k1.max(), k2 / k2.max())

for nm, pf, D in (('template_3000', None, 0.0030), ('physical_3000', prof, 0.0030), ('physical_6000', prof, 0.0060)):
    rec = L.e2e_injection(t, f, seg, shp, EX, 1, L.fake_offsets(80, rng), D, 0.35, minpts=10, floor=2e-4, profile=pf)
    ok = np.isfinite(rec[:, 2])
    out[nm] = dict(depth=D, n=int(ok.sum()), gain_deep=float(np.nanmean(rec[:, 2])), gain_deep_sem=float(np.nanstd(rec[:, 2]) / np.sqrt(ok.sum())),
                   gain_shallow=float(np.nanmean(rec[:, 3])), gain_shallow_sem=float(np.nanstd(rec[:, 3]) / np.sqrt(ok.sum())),
                   rec_deep_mean=float(np.nanmean(rec[:, 0])), rec_deep_sd=float(np.nanstd(rec[:, 0])))
    print('%-14s gain deep %.3f +- %.3f  shallow %.3f +- %.3f  (recovered deep %.0f +- %.0f ppm)' % (
        nm, out[nm]['gain_deep'], out[nm]['gain_deep_sem'], out[nm]['gain_shallow'], out[nm]['gain_shallow_sem'],
        out[nm]['rec_deep_mean'] * 1e6, out[nm]['rec_deep_sd'] * 1e6), flush=True)

# ---- alternative template: p=2, q=1 (parabolic V), T fitted on 200-s sectors, both eclipses
data = TS.load_all()
for sec, d in data.items():
    d['w'] = np.ones(len(d['t']))


def shp_alt(Td, Ts, dd, ds):
    return {'d': dict(T=Td, p=2.0, q=1.0, phc=dd), 's': dict(T=Ts, p=2.0, q=1.0, phc=0.5 + ds)}


def chi2(par):
    sh = shp_alt(*par)
    tot = 0
    for s in (68, 95, 105, 106):
        a, e, r = TS.fit_sector(data[s], sh)
        tot += np.sum(r ** 2) / L.mad_std(r) ** 2
    return tot


r = minimize(chi2, [0.07, 0.07, 0.0, 0.0], method='Nelder-Mead', options=dict(xatol=1e-6, fatol=1e-3))
alt = shp_alt(*r.x)
print('alt template T_d %.4f T_s %.4f d (%.2f, %.2f h); chi2 %.1f vs default-template chi2 %.1f' % (
    r.x[0], r.x[1], r.x[0] * 24, r.x[1] * 24, r.fun, chi2_def := sum(np.sum(TS.fit_sector(data[s], T['shapes'])[2] ** 2) / L.mad_std(TS.fit_sector(data[s], T['shapes'])[2]) ** 2 for s in (68, 95, 105, 106))))
amps = {}
for s, d in data.items():
    a, e, rr = TS.fit_sector(d, alt)
    a0, e0, r0 = TS.fit_sector(d, T['shapes'])
    amps[s] = (a[0], a[1], a0[0], a0[1])
amps_arr = np.array(list(amps.values()))
tess_alt = amps_arr[:, 0].mean(); tess_alt_s = amps_arr[:, 1].mean()
tess_def = amps_arr[:, 2].mean(); tess_def_s = amps_arr[:, 3].mean()
print('TESS mean (unweighted over 6 sectors): alt template deep %.0f shallow %.0f; default template deep %.0f shallow %.0f' % (
    tess_alt * 1e6, tess_alt_s * 1e6, tess_def * 1e6, tess_def_s * 1e6))
shp_alt_s = L.survey_shapes(alt, E3)
ra, _ = L.run_survey(t, f, seg, shp_alt_s, EX, 1, nfake=500, ninj=50, s_over_d=sod)
rd, _ = L.run_survey(t, f, seg, shp, EX, 1, nfake=500, ninj=50, s_over_d=sod)
out['alt_template'] = dict(shape=alt, tess_deep=float(tess_alt), tess_shallow=float(tess_alt_s), tess_deep_default_unw=float(tess_def),
                           wasp_deep=ra['deep'], wasp_deep_err=ra['deep_err'], wasp_shallow=ra['shallow'], wasp_shallow_err=ra['shallow_err'],
                           ratio_deep_alt=float(ra['deep'] / tess_alt), ratio_deep_default=float(rd['deep'] / tess_def),
                           ratio_shallow_alt=float(ra['shallow'] / tess_alt_s), ratio_shallow_default=float(rd['shallow'] / tess_def_s))
print('WASP alt template deep %.0f +- %.0f shallow %.0f +- %.0f ; ratio deep alt %.3f vs default %.3f ; shallow alt %.3f vs default %.3f' % (
    ra['deep'] * 1e6, ra['deep_err'] * 1e6, ra['shallow'] * 1e6, ra['shallow_err'] * 1e6, out['alt_template']['ratio_deep_alt'],
    out['alt_template']['ratio_deep_default'], out['alt_template']['ratio_shallow_alt'], out['alt_template']['ratio_shallow_default']))
json.dump(out, open(L.V3 + '/results/wasp3_e2e.json', 'w'), indent=1)
print('saved')
