"""Referee: planet quantities for depths 1.4-1.9% on the M3V (single) and on one twin (depth x2), with batman.
k solved so that the batman model depth (minimum flux) equals the given depth for quadratic LD (u1,u2) = (0.20, 0.38)
and also (0.15, 0.45) as an LD check; T14 for b = 0-0.9; a, S, Teq; K (Chen & Kipping not re-derived: use M_p grid);
TSM (Kempton+18, scale 1.28 for 4 < Rp < 10 Re; 1.26 for 2.75-4) with J = 11.363.
Output: out/r11_planet.json"""
import json, numpy as np, batman
from scipy.optimize import brentq
R = 'tic214311908/analysis/referee_star'
P = 5.6914968; Rsun_au = 0.00465047; Re_per_Rsun = 109.076
cases = dict(single=dict(R=0.3813, M=0.3702, L=0.0165, Teff=3350, dil=1.0), twin=dict(R=0.2888, M=0.2639, L=0.00825, Teff=3237, dil=2.0))
def depth_of(k, aR, b, u):
    pm = batman.TransitParams(); pm.t0 = 0; pm.per = P; pm.rp = k; pm.a = aR; pm.inc = np.degrees(np.arccos(b / aR)); pm.ecc = 0; pm.w = 90; pm.u = list(u); pm.limb_dark = 'quadratic'
    t = np.linspace(-0.002, 0.002, 21); m = batman.TransitModel(pm, t); return 1 - m.light_curve(pm).min()
def T14(k, aR, b): return P * 24 / np.pi * np.arcsin(np.sqrt((1 + k) ** 2 - b ** 2) / aR / np.sin(np.arccos(b / aR)))
out = {}
for lab, c in cases.items():
    a = (c['M'] * (P / 365.25) ** 2) ** (1 / 3); aR = a / (c['R'] * Rsun_au)
    S = c['L'] / a ** 2 if lab == 'single' else c['L'] / a ** 2
    Teq0 = c['Teff'] * np.sqrt(1 / aR) * 0.25 ** 0.25
    o = dict(a_au=a, a_over_R=aR, S=S, Teq_A0=Teq0, Teq_A03=Teq0 * 0.7 ** 0.25, grid={})
    for dep in (0.014, 0.0164, 0.0175, 0.019):
        for u in ((0.20, 0.38), (0.15, 0.45)):
            for b in (0.0, 0.3, 0.6, 0.8):
                d_ = dep * c['dil']
                k = brentq(lambda kk: depth_of(kk, aR, b, u) - d_, 0.03, 0.4)
                o['grid'][f'd{dep}_u{u[0]}_b{b}'] = dict(k=k, Rp=k * c['R'] * Re_per_Rsun, T14_h=T14(k, aR, b))
    out[lab] = o
    print(f"{lab}: a = {a:.4f} au, a/R = {aR:.2f}, S = {S:.2f} S_e, Teq(A=0) = {Teq0:.0f} K")
    for dep in (0.014, 0.0164, 0.019):
        print('  depth', dep, ' '.join(f"b={b}: Rp={o['grid'][f'd{dep}_u0.2_b{b}']['Rp']:.2f} T14={o['grid'][f'd{dep}_u0.2_b{b}']['T14_h']:.2f}h" for b in (0.0, 0.3, 0.6, 0.8)),
              '| LD (0.15,0.45) b=0.6: Rp', round(o['grid'][f'd{dep}_u0.15_b0.6']['Rp'], 2))
# K and TSM (single)
c = cases['single']; Rp = out['single']['grid']['d0.0164_u0.2_b0.6']['Rp']
for Mp in (10, 20, 23, 30):
    K = 28.4329 / np.sqrt(1) * (Mp / 317.83) * (c['M'] + Mp / 332946) ** (-2 / 3) * (P / 365.25) ** (-1 / 3)
    tsm = (1.28 if Rp > 4 else 1.26) * Rp ** 3 * out['single']['Teq_A0'] / (Mp * c['R'] ** 2) * 10 ** (-11.363 / 5)
    out['single'][f'K_Mp{Mp}'] = K; out['single'][f'TSM_Mp{Mp}'] = tsm
    print(f'  Mp {Mp}: K = {K:.2f} m/s, TSM = {tsm:.0f}')
# stellar density implied by T14 = 1.63 h with b = 0.6 (circular), for the single-star k
json.dump(out, open(f'{R}/out/r11_planet.json', 'w'), indent=1, default=float)
