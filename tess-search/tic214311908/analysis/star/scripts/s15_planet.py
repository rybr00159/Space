"""Planet parameters implied by a 1.4-1.9% transit on TIC 214311908 (P = 5.6914968 d).
Uses the stellar Monte Carlo of s11 (out/star_params_mc.npz; single-star case) and the equal-mass-binary case.
- Radius ratio k from the depth with quadratic limb darkening (TESS band, Teff 3350 K, logg 4.85: u1 = 0.20,
  u2 = 0.38, approximate Claret 2017 values) using batman, for b = 0-0.9.
- a, a/R*, insolation S (S_earth), Teq (A = 0 and 0.3, full redistribution), T14 and T23 for b = 0-0.9 (circular),
  the impact parameter matching T14 = 1.61-1.65 h (QLP) or 1.83 +- 0.12 h (deblended, phot/), and the stellar density
  implied by T14 for each b.
- RV semi-amplitude for M_p = 10, 20, 30, 50, 100 M_earth and for the Chen & Kipping (2017) mass-radius relation
  (Neptunian branch, with its 0.15 dex scatter).
- TSM and ESM (Kempton+18; TSM scale factor 1.28 for 4 < Rp < 10 R_earth; J = 11.363, Ks = 10.482).
Output: out/planet_params.json"""
import json, numpy as np, batman
B = 'tic214311908/analysis/star'
S = np.load(f'{B}/out/star_params_mc.npz'); sp = json.load(open(f'{B}/star_params.json'))
P = 5.6914968; Pd = P * 86400; G = 6.6743e-8; Msun, Rsun, Rearth, Mearth, AU = 1.98847e33, 6.957e10, 6.371e8, 5.972e27, 1.495978707e13
u1, u2 = 0.20, 0.38
def depth_of_k(k, b, aR=27.0):
    pm = batman.TransitParams(); pm.t0 = 0; pm.per = P; pm.rp = k; pm.a = aR; pm.inc = np.degrees(np.arccos(min(b / aR, 1))); pm.ecc = 0; pm.w = 90
    pm.u = [u1, u2]; pm.limb_dark = 'quadratic'
    t = np.array([0.0]); m = batman.TransitModel(pm, t); return 1 - m.light_curve(pm)[0]
kgrid = np.linspace(0.08, 0.25, 171)
def k_of_depth(dep, b):
    dd = np.array([depth_of_k(k, b) for k in kgrid]); return float(np.interp(dep, dd, kgrid))
cases = {'single': dict(R=S['R'], M=S['M'], L=S['L'], T=S['Teff']),
         'twin': dict(R=np.full(20000, sp['binary_twin_scenario']['radius_each_rsun']['median']) * (1 + np.random.default_rng(1).normal(0, 0.031, 20000)),
                      M=np.full(20000, sp['binary_twin_scenario']['mass_each_msun']['median']) * (1 + np.random.default_rng(2).normal(0, 0.031, 20000)),
                      L=S['L'] / 2, T=np.full(20000, sp['binary_twin_scenario']['teff_each_K']['median']))}
out = {}
def stat(x): return [float(np.median(x)), float(np.percentile(x, 15.87)), float(np.percentile(x, 84.13))]
rng = np.random.default_rng(9)
for cname, c in cases.items():
    R, M, L, T = c['R'], c['M'], c['L'], c['T']
    a = (G * M * Msun * Pd ** 2 / (4 * np.pi ** 2)) ** (1 / 3); aR = a / (R * Rsun)
    Sins = L / (a / AU) ** 2
    Teq0 = T * np.sqrt(R * Rsun / (2 * a)); Teq3 = Teq0 * 0.7 ** 0.25
    res = dict(a_au=stat(a / AU), a_over_Rstar=stat(aR), insolation_Searth=stat(Sins), Teq_A0=stat(Teq0), Teq_A03=stat(Teq3))
    dil = 2.0 if cname == 'twin' else 1.0
    tab = []
    for dep in (0.014, 0.0164, 0.0175, 0.019):
        for b in (0.0, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9):
            if dep * dil > 0.06: continue
            k = k_of_depth(dep * dil, b)
            Rp = k * R * Rsun / Rearth
            inc = np.arccos(b / aR)
            T14 = P / np.pi * np.arcsin(np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) / (aR * np.sin(inc))) * 24
            T23 = P / np.pi * np.arcsin(np.sqrt(np.clip((1 - k) ** 2 - b ** 2, 0, None)) / (aR * np.sin(inc))) * 24
            tab.append(dict(depth=dep, b=b, k=k, Rp_earth=stat(Rp), T14_h=stat(T14), T23_over_T14=stat(T23 / T14)))
    res['grid'] = tab
    # impact parameter and density from the observed T14
    for lab, (t14, et14) in {'QLP_1.63h': (1.63, 0.03), 'deblended_1.83h': (1.83, 0.12)}.items():
        bb = np.linspace(0, 0.95, 96); k = k_of_depth(0.0164 * dil, 0.6)
        t14s = rng.normal(t14, et14, len(aR))
        # for each MC star find b that gives T14
        bsol = []
        for i in range(0, len(aR), 20):
            inc = np.arccos(bb / aR[i]); tt = P / np.pi * np.arcsin(np.sqrt(np.clip((1 + k) ** 2 - bb ** 2, 0, None)) / (aR[i] * np.sin(inc))) * 24
            bsol.append(np.interp(t14s[i], tt[::-1], bb[::-1]) if (t14s[i] <= tt[0] and t14s[i] >= tt[-1]) else (0.0 if t14s[i] > tt[0] else np.nan))
        bsol = np.array(bsol)
        frac_long = float(np.mean(t14s[::20] > P / np.pi * np.arcsin((1 + k) / aR[::20]) * 24))
        res[f'b_from_T14_{lab}'] = dict(b=stat(bsol[np.isfinite(bsol)]), frac_T14_longer_than_central=frac_long)
    # implied stellar density vs b for T14 = 1.63 and 1.83 h (k from depth 1.64%)
    rho_b = {}
    for t14 in (1.63, 1.83):
        for b in (0.0, 0.3, 0.5, 0.6, 0.7, 0.8):
            k = k_of_depth(0.0164 * dil, b)
            # solve aR from T14 (circular)
            s = np.sin(np.pi * t14 / 24 / P)
            aR_ = np.sqrt(((1 + k) ** 2 - b ** 2) / s ** 2 + b ** 2)
            rho = 3 * np.pi * aR_ ** 3 / (G * Pd ** 2)
            rho_b[f'T14={t14}h,b={b}'] = float(rho)
    res['implied_density_gcc'] = rho_b
    res['star_density_gcc'] = stat(M * Msun / (4 / 3 * np.pi * (R * Rsun) ** 3))
    # RV and TSM/ESM
    kk = k_of_depth(0.0164 * dil, 0.6); Rp = kk * R * Rsun / Rearth
    rv = {}
    for mp in (10, 20, 30, 50, 100):
        K = 28.4329 / np.sqrt(1) * (mp * Mearth / 1.89813e30) * (M) ** (-2 / 3) * (P / 365.25) ** (-1 / 3)
        rv[f'{mp}Me'] = stat(K)
    # Chen & Kipping 2017 Neptunian branch: log10 R = C_N + 0.589 log10 M, continuous with terran branch at 2.04 Me (C_T = 0.00346, S_T = 0.2790)
    CN = 0.00346 + 0.2790 * np.log10(2.04) - 0.589 * np.log10(2.04)
    lgM = (np.log10(Rp) - CN) / 0.589 + rng.normal(0, 0.15 / 0.589, len(Rp))
    Mp = 10 ** lgM
    K_ck = 28.4329 * (Mp * Mearth / 1.89813e30) * M ** (-2 / 3) * (P / 365.25) ** (-1 / 3)
    rv['ChenKipping17'] = dict(Mp_earth=stat(Mp), K_ms=stat(K_ck))
    Teq = Teq0
    TSM = 1.28 * Rp ** 3 * Teq / (Mp * R ** 2) * 10 ** (-11.363 / 5)
    def Bl(T, lam=7.5e-4):
        return 1 / (np.exp(1.4387769 / (lam * T)) - 1)
    ESM = 4.29e6 * Bl(1.10 * Teq) / Bl(T) * (kk) ** 2 * 10 ** (-10.482 / 5)
    res.update(Rp_b06_depth1p64=stat(Rp), rv=rv, TSM=stat(TSM), ESM=stat(ESM))
    out[cname] = res
    print(f"== {cname}: a = {stat(a/AU)[0]:.4f} au, a/R* = {stat(aR)[0]:.1f}, S = {stat(Sins)[0]:.1f} S_E, Teq(A=0) = {stat(Teq0)[0]:.0f} K, Rp(1.64%, b=0.6) = {stat(Rp)}")
    print('   b from T14:', {k_: v for k_, v in res.items() if k_.startswith('b_from')})
    print('   rho* =', res['star_density_gcc'], ' implied by T14:', {k_: round(v, 2) for k_, v in rho_b.items()})
    print('   K:', {k_: (round(v[0], 1) if isinstance(v, list) else {kk_: [round(x, 1) for x in vv] for kk_, vv in v.items()}) for k_, v in rv.items()})
    print('   TSM', stat(TSM), 'ESM', stat(ESM))
    for t in tab:
        if t['b'] in (0.0, 0.6, 0.8): print(f"   depth {t['depth']:.4f} b {t['b']}: k {t['k']:.3f} Rp {t['Rp_earth'][0]:.2f} T14 {t['T14_h'][0]:.2f} h T23/T14 {t['T23_over_T14'][0]:.2f}")
json.dump(out, open(f'{B}/out/planet_params.json', 'w'), indent=1)
