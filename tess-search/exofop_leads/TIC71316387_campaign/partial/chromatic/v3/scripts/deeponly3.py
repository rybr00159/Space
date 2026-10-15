"""Deep-eclipse-only version of the MS-pair shape test (shallow eclipse modelled by the empirical template,
so an eccentric orbit / different shallow duration cannot drive the result)."""
import sys, json
import numpy as np
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/v3/scripts')
import lib3 as L
import tess3 as TS
import predict3 as PR
import sed3 as S
import importlib
src = open(L.V3 + '/scripts/pairs3.py').read()
src = src[:src.index("prims = [")]            # reuse the setup part of pairs3 (data, projectors, helpers)
exec(compile(src, 'pairs3_setup', 'exec'))
ps_t = L.template(tg, 0.0, shp['s']['T'], shp['s']['p'], shp['s']['q'], t0=0.0, per=L.P)
pd_t = L.template(tg, 0.0, shp['d']['T'], shp['d']['p'], shp['d']['q'], t0=0.0, per=L.P)
c_ref = chi2_free(pd_t, ps_t)[0]
out = {}
for p1, p2 in (('K7V', 'M2V'), ('M0V', 'M3V'), ('M0.5V', 'M3V'), ('M1V', 'M3.5V'), ('M1V', 'M4V')):
    r1, r2 = tab[p1], tab[p2]
    c1, c2 = comp(p1), comp(p2)
    F1V = S.ZP['V'] * 10 ** (-0.4 * r1['Mv']); F2V = S.ZP['V'] * 10 ** (-0.4 * r2['Mv'])
    FGV = FV_sys - F1V - F2V
    Ft = (FGV * rel(gm, 'TESS'), F1V * rel(c1, 'TESS'), F2V * rel(c2, 'TESS'))
    f1MS = Ft[1] / sum(Ft)
    a = (PR.GSUN * (r1['Msun'] + r2['Msun']) * PR.MSUN * (L.P * 86400) ** 2 / (4 * np.pi ** 2)) ** (1 / 3) / PR.RSUN
    R1, R2 = r1['R_Rsun'], r2['R_Rsun']
    ld1 = S.ld_quadratic(r1['Teff'], 4.7, 'TESS'); ld2 = S.ld_quadratic(r2['Teff'], 4.8, 'TESS')
    imin = np.degrees(np.arccos(min((R1 + R2) / a, 1))) + 0.1
    rows = []
    for inc in np.arange(imin, 90.01, 0.25):
        k1, k2 = PR.lightcurves(R1, R2, a, inc, ld1, ld2, tg)
        if k1.max() <= 0:
            continue
        c, Ad, As = chi2_free(k1, ps_t)
        rows.append((inc, c, Ad, k1.max()))
    rows = np.array(rows)
    i = np.argmin(rows[:, 1])
    # depth-matched (bound MS at 297 pc): inclination where f1MS * template-amp(k1) ~ observed deep amplitude
    amp_obs = chi2_free(pd_t, ps_t)[1]
    pred = []
    for inc in rows[:, 0]:
        k1, k2 = PR.lightcurves(R1, R2, a, inc, ld1, ld2, tg)
        pred.append(f1MS * PR.template_amp(tg, k1, shp, 'd'))
    pred = np.array(pred)
    j = np.argmin(np.abs(pred - amp_obs))
    d1 = 297.4 * np.sqrt(f1MS / rows[i, 2])
    out[p1 + '+' + p2] = dict(best_inc=float(rows[i, 0]), dchi2_best_vs_template=float(rows[i, 1] - c_ref), f1_fit=float(rows[i, 2]), f1_MS=float(f1MS),
                              d1_pc=float(d1), inc_bound=float(rows[j, 0]), dchi2_bound_vs_best=float(rows[j, 1] - rows[i, 1]))
    print('%-11s deep-only: best inc %.2f (dchi2 vs template %+.1f) f1 %.4f (MS@297 %.4f) -> d %.0f pc; bound depth-matched inc %.2f is worse by dchi2 %.1f' % (
        p1 + '+' + p2, rows[i, 0], rows[i, 1] - c_ref, rows[i, 2], f1MS, d1, rows[j, 0], rows[j, 1] - rows[i, 1]), flush=True)
json.dump(out, open(L.V3 + '/results/deeponly3.json', 'w'), indent=1)
