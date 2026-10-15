"""Chromatic predictions and chi2 for the MS pairs that the TESS eclipse SHAPE prefers (pairs3: small,
cool secondary, near-central eclipses), at the inclination fitted to TESS, with the full multi-band light
curves (LD per band, MS colours), fitted with the same TESS template as the data.  Joint chi2 over all
archival datasets uses the combine3 covariances (W1-W2 deep/shallow correlation included)."""
import sys, json
import numpy as np
from scipy.stats import chi2 as CHI2
sys.path.insert(0, '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/chromatic/v3/scripts')
import lib3 as L
import predict3 as PR
import sed3 as S

C3 = json.load(open(L.V3 + '/results/combine3.json'))
P3s = json.load(open(L.V3 + '/results/pairs3.json'))['pairs']
T3 = json.load(open(L.V3 + '/results/tess3.json'))
shp = T3['shapes']
tab = S.mamajek()
gm = S.star_model(teff=5883.0, tab=tab)
rel = lambda m, b: S.band_flux(m['spec'], b) / S.band_flux(m['spec'], 'V')
V_SYS, DIST, AV = 11.857, 297.4, 0.02
FV_sys = S.ZP['V'] * 10 ** (-0.4 * (V_SYS - 5 * np.log10(DIST / 10) - AV))
BANDS = ['TESS', 'WASP', 'V', 'W1', 'W2', 'Ks', 'I']
tg = np.linspace(-0.12, 0.12, 241)
meas = C3['measurements']
names = list(meas)
vec = np.concatenate([meas[n]['ratio'] for n in names])
Cbig = np.zeros((len(vec), len(vec)))
for i, n in enumerate(names):
    Cbig[2 * i:2 * i + 2, 2 * i:2 * i + 2] = meas[n]['cov']
iw1, iw2 = names.index('NEOWISE_W1'), names.index('NEOWISE_W2')
for j in (0, 1):
    a_, b_ = 2 * iw1 + j, 2 * iw2 + j
    Cbig[a_, b_] = Cbig[b_, a_] = C3['W1_W2_null_corr_deep'] * np.sqrt(Cbig[a_, a_] * Cbig[b_, b_])
out = {}
sel = sorted([k for k, v in P3s.items() if v['dchi2_shape'] < 15], key=lambda k: P3s[k]['dchi2_shape'])
sel += ['K7V+M2V', 'M0V+M3V']
print('%-12s dchi2_shape d1(pc) | ratios deep: WASP V W1 W2 | shallow: WASP V W1 W2 | dI dKs pairT | joint chi2 (10)' )
for key in sel:
    p = P3s[key]
    s1, s2 = key.split('+')
    r1, r2 = tab[s1], tab[s2]
    c1, c2 = S.star_model(s1, tab=tab), S.star_model(s2, tab=tab)
    F1V = S.ZP['V'] * 10 ** (-0.4 * r1['Mv']); F2V = S.ZP['V'] * 10 ** (-0.4 * r2['Mv'])
    FGV = FV_sys - F1V - F2V
    F = {b: (FGV * rel(gm, b), F1V * rel(c1, b), F2V * rel(c2, b)) for b in BANDS}
    frac = {b: (F[b][1] / sum(F[b]), F[b][2] / sum(F[b])) for b in BANDS}
    a = (PR.GSUN * (r1['Msun'] + r2['Msun']) * PR.MSUN * (L.P * 86400) ** 2 / (4 * np.pi ** 2)) ** (1 / 3) / PR.RSUN
    amps = {}
    for b in BANDS:
        ld1 = S.ld_quadratic(r1['Teff'], 4.7, b); ld2 = S.ld_quadratic(r2['Teff'], 4.8, b)
        k1, k2 = PR.lightcurves(r1['R_Rsun'], r2['R_Rsun'], a, p['inc'], ld1, ld2, tg)
        amps[b] = (PR.template_amp(tg, k1 * frac[b][0], shp, 'd'), PR.template_amp(tg, k2 * frac[b][1], shp, 's'))
    rd = {b: amps[b][0] / amps['TESS'][0] for b in BANDS}
    rs = {b: amps[b][1] / amps['TESS'][1] for b in BANDS}
    pv = np.concatenate([[rd[meas[n]['band']], rs[meas[n]['band']]] for n in names])
    dl = vec - pv
    jc = float(dl @ np.linalg.solve(Cbig, dl))
    dI = float(-2.5 * np.log10((F['I'][1] + F['I'][2]) / F['I'][0])); dK = float(-2.5 * np.log10((F['Ks'][1] + F['Ks'][2]) / F['Ks'][0]))
    pairT = float(frac['TESS'][0] + frac['TESS'][1])
    out[key] = dict(dchi2_shape=p['dchi2_shape'], d1=p['d1'], d2=p['d2'], inc=p['inc'], ratio_deep=rd, ratio_shallow=rs, joint_chi2=jc,
                    p_joint=float(CHI2.sf(jc, len(vec))), dI_MS297=dI, dKs_MS297=dK, pair_frac_TESS_MS297=pairT,
                    pair_frac_TESS_fit=float(p['f1_fit'] + p['f2_fit']))
    print('%-12s %+6.1f %5.0f | %.3f %.3f %.2f %.2f | %.3f %.3f %.2f %.2f | %.2f %.2f %.3f | %.2f' % (
        key, p['dchi2_shape'], p['d1'], rd['WASP'], rd['V'], rd['W1'], rd['W2'], rs['WASP'], rs['V'], rs['W1'], rs['W2'], dI, dK, pairT, jc), flush=True)
out['G/F_best_joint_chi2'] = C3['joint']['chi2_bg_ge5600']
print('best G/F (>=5600 K) joint chi2: %.2f' % C3['joint']['chi2_bg_ge5600'])
json.dump(out, open(L.V3 + '/results/ext_pairs3.json', 'w'), indent=1)
