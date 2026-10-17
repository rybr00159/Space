"""Referee: independent stellar parameters of TIC 214311908 (single-star and equal-twin cases).
Inputs: Gaia DR3 row from data/gaia_cone45.csv (own query); Lindegren+21 zero point via the zero_point package;
2MASS Ks = 10.482 +- 0.026 (2MASS PSC; note cc_flg = 'ddd' -- see r07); Mann+15 R(MKs) (eq. 4, 2.89% scatter),
Mann+15 R(MKs,[Fe/H]); Mann+19 M(MKs) (n=5 polynomial, 2.2% scatter); density; also an
independent L from a bolometric correction (Mann+15 BC_J vs (r-J)?? not used: we use the Mamajek EEM table's
BC_V and V-Ks to get BC_Ks as a function of BP-RP) to check the SED Fbol.
Parallax error inflated x1.3 (El-Badry+21 for G~14) plus 0.015 mas ZP systematic, as the investigator, and a
variant without inflation.
Output: out/r06_star_params.json"""
import json, numpy as np, pandas as pd
from zero_point import zpt
zpt.load_tables()
R = 'tic214311908/analysis/referee_star'
g = pd.read_csv(f'{R}/data/gaia_cone45.csv'); r = g[g.source_id == 5925373272850784768].iloc[0]
zp = zpt.get_zpt(r.phot_g_mean_mag, r.nu_eff_used_in_astrometry, r.pseudocolour, r.ecl_lat, r.astrometric_params_solved)
print('Gaia: plx', r.parallax, '+-', r.parallax_error, 'ZP', zp, 'G', r.phot_g_mean_mag, 'BP-RP', r.bp_rp, 'nu_eff', r.nu_eff_used_in_astrometry)
rng = np.random.default_rng(123); N = 200000
plx0 = r.parallax - zp
out = dict(parallax=float(r.parallax), parallax_error=float(r.parallax_error), zp=float(zp), plx_corr=float(plx0))
def R15(MK): return 1.9515 - 0.3520 * MK + 0.01680 * MK ** 2
def R15f(MK, fe): return (1.9305 - 0.3466 * MK + 0.01647 * MK ** 2) * (1 + 0.04458 * fe)
def M19(MK):
    x = MK - 7.5; a = [-0.642, -0.208, -8.43e-4, 7.87e-3, 1.42e-4, -2.13e-4]
    return 10 ** sum(ai * x ** i for i, ai in enumerate(a))
def st(x): return dict(median=float(np.median(x)), lo=float(np.percentile(x, 15.87)), hi=float(np.percentile(x, 84.13)), sd=float(np.std(x)))
for lab, eplx, Ks, eKs in (('inflated', np.hypot(1.3 * r.parallax_error, 0.015), 10.482, 0.026), ('formal', r.parallax_error, 10.482, 0.026),
                           ('Ks_DENIS_check', np.hypot(1.3 * r.parallax_error, 0.015), 10.823, 0.09)):
    plx = rng.normal(plx0, eplx, N); d = 1000 / plx
    MK = rng.normal(Ks, eKs, N) + 5 * np.log10(plx / 100)
    Rs = R15(MK) * (1 + rng.normal(0, 0.0289, N))
    Ms = M19(MK) * (1 + rng.normal(0, 0.022, N)) * (1 + rng.normal(0, 0.02, N))   # 2.2% scatter + 2% systematic
    rho = Ms / Rs ** 3 * 1.40987
    MKt = MK + 2.5 * np.log10(2); Rt = R15(MKt) * (1 + rng.normal(0, 0.0289, N)); Mt = M19(MKt) * (1 + rng.normal(0, 0.022, N)) * (1 + rng.normal(0, 0.02, N))
    fe = rng.normal(0.2, 0.15, N); Rf = R15f(MK, fe) * (1 + rng.normal(0, 0.027, N)); rhof = Ms / Rf ** 3 * 1.40987
    out[lab] = dict(d=st(d), MKs=st(MK), R=st(Rs), M=st(Ms), rho_gcc=st(rho), R_feh02=st(Rf), rho_feh02=st(rhof),
                    twin_R=st(Rt), twin_M=st(Mt), twin_rho=st(Mt / Rt ** 3 * 1.40987))
    o = out[lab]
    print(f"{lab}: d {o['d']['median']:.2f}+-{o['d']['sd']:.2f} MKs {o['MKs']['median']:.3f}+-{o['MKs']['sd']:.3f} R {o['R']['median']:.4f}+-{o['R']['sd']:.4f} "
          f"M {o['M']['median']:.4f}+-{o['M']['sd']:.4f} rho {o['rho_gcc']['median']:.2f} ({o['rho_gcc']['lo']:.2f}-{o['rho_gcc']['hi']:.2f}); "
          f"[Fe/H]=+0.2: R {o['R_feh02']['median']:.4f} rho {o['rho_feh02']['median']:.2f}; twin R {o['twin_R']['median']:.4f} M {o['twin_M']['median']:.4f} rho {o['twin_rho']['median']:.2f} ({o['twin_rho']['lo']:.2f}-{o['twin_rho']['hi']:.2f})")
# luminosity check via Mamajek EEM BC (BC_Ks = BC_V + V - Ks) interpolated at BP-RP, and L from SED Fbol (investigator)
try:
    rows = []
    for line in open(f'{R}/data/EEM_dwarf.txt'):
        p = line.split()
        if len(p) > 25 and p[0][:1] in 'KM' and p[0].endswith('V'):
            rows.append(p)
    # header columns located from the header line
    hdr = [l for l in open(f'{R}/data/EEM_dwarf.txt') if l.startswith('#SpT')][0].split()
    ci = {k: hdr.index(k) for k in ('#SpT', 'Teff', 'BCv', 'Bp-Rp', 'V-Ks', 'M_Ks', 'R_Rsun', 'Msun', 'M_G')}
    tab = []
    for p in rows:
        try: tab.append([float(p[ci[k]]) if p[ci[k]] not in ('...', '....', '.....') else np.nan for k in ('Teff', 'BCv', 'Bp-Rp', 'V-Ks', 'M_Ks', 'R_Rsun', 'Msun', 'M_G')])
        except Exception: pass
    tab = np.array(tab); k = np.isfinite(tab).all(1); tab = tab[k]; o_ = np.argsort(tab[:, 2]); tab = tab[o_]
    bprp = r.bp_rp
    BCK = np.interp(bprp, tab[:, 2], tab[:, 1] + tab[:, 3])
    MKs_t = 10.482 + 5 * np.log10(plx0 / 100)
    Mbol = MKs_t + BCK; L = 10 ** (-0.4 * (Mbol - 4.74))
    MG_ms = np.interp(bprp, tab[:, 2], tab[:, 7]); MK_ms = np.interp(bprp, tab[:, 2], tab[:, 4])
    out['mamajek'] = dict(BC_Ks=float(BCK), L_from_BCKs=float(L), MG_ms_at_bprp=float(MG_ms), MKs_ms_at_bprp=float(MK_ms),
                          dMG=float(r.phot_g_mean_mag + 5 * np.log10(plx0 / 100) - MG_ms), dMKs=float(MKs_t - MK_ms))
    print('Mamajek: BC_Ks', round(BCK, 3), 'L', round(L, 5), 'Lsun; MS at this BP-RP: MG', round(MG_ms, 2), 'MKs', round(MK_ms, 2),
          '-> offsets dMG', round(out['mamajek']['dMG'], 2), 'dMKs', round(out['mamajek']['dMKs'], 2))
except Exception as e:
    print('Mamajek check failed', e)
json.dump(out, open(f'{R}/out/r06_star_params.json', 'w'), indent=1)
