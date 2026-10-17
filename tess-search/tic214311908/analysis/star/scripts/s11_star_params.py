"""Adopted stellar parameters of TIC 214311908 (Gaia DR3 5925373272850784768), Monte Carlo (100k draws).
Inputs: Gaia DR3 parallax 15.4973 +- 0.0251 mas + Lindegren+21 zero point (-0.0521 mas; zero_point package),
error inflated x1.3 (El-Badry+21, G ~ 14) and 0.015 mas ZP systematic; 2MASS Ks = 10.482 +- 0.026;
Fbol from the SED fit (s06; +-2% systematic); [Fe/H]_phot = +0.19 +- 0.16 (s10, valid only if single).
Relations: R(MKs) Mann+15 eq. 4 (+-2.89%) and its [Fe/H] form (eq. 5, +-2.70%); M(MKs) Mann+19 n=5 polynomial
(+-2.0% scatter, +2% systematic); Teff from Stefan-Boltzmann with R_Mann and L_SED; Teff from Mann+15 colour
relations (V-J with V from Gaia G, BP-RP; Riello+21 G-V) and BT-Settl SED as checks.
Also: the 'equal-mass unresolved binary' scenario (each star 0.753 mag fainter in Ks, half of Fbol).
Output: star_params.json (top-level folder), out/star_params_mc.npz"""
import json, numpy as np
B = 'tic214311908/analysis/star'
rng = np.random.default_rng(42); N = 100000
sed = json.load(open(f'{B}/out/sed.json')); feh = json.load(open(f'{B}/out/metallicity.json'))
plx0 = 15.497306 + 0.052059; eplx = np.hypot(0.025122 * 1.3, 0.015)
plx = rng.normal(plx0, eplx, N); d = 1000 / plx
Ks = rng.normal(10.482, 0.026, N); MK = Ks + 5 * np.log10(plx / 100)
def R_mann15(MK): return 1.9515 - 0.3520 * MK + 0.01680 * MK ** 2
def R_mann15_feh(MK, fe): return (1.9305 - 0.3466 * MK + 0.01647 * MK ** 2) * (1 + 0.04458 * fe)
def M_mann19(MK):
    x = MK - 7.5; a = [-0.642, -0.208, -8.43e-4, 7.87e-3, 1.42e-4, -2.13e-4]
    return 10 ** sum(ai * x ** i for i, ai in enumerate(a))
R = R_mann15(MK) * (1 + rng.normal(0, 0.0289, N))
fe = rng.normal(feh['feh_phot_if_single'][0], feh['feh_phot_if_single'][1], N)
Rfe = R_mann15_feh(MK, fe) * (1 + rng.normal(0, 0.027, N))
M = M_mann19(MK) * (1 + rng.normal(0, 0.020, N)) * (1 + rng.normal(0, 0.02, N))
Fbol = sed['mc']['fbol_cgs'][0] * (1 + rng.normal(0, 0.02, N))
L = 4 * np.pi * (d * 3.0857e18) ** 2 * Fbol / 3.828e33
Teff_SB = 5772.0 * (L / R ** 2) ** 0.25
rho = M / R ** 3 * 1.40987   # g/cm3 (rho_sun = 1.40987 for nominal solar values)
logg = np.log10(6.6743e-8 * M * 1.98847e33 / (R * 6.957e10) ** 2)
# colour Teff (checks)
V = 14.281043 + 1.3269        # Riello+21 G-V at BP-RP = 2.756 (sigma 0.03; extrapolation edge for M dwarfs)
X = V - 11.363; Teff_VJ = 3500 * (2.840 - 1.3453 * X + 0.3906 * X ** 2 - 0.0546 * X ** 3 + 0.002913 * X ** 4)
X = 2.756149; Teff_BPRP = 3500 * (3.245 - 2.4309 * X + 1.043 * X ** 2 - 0.2127 * X ** 3 + 0.01649 * X ** 4)
# equal-mass binary scenario
MKb = MK + 2.5 * np.log10(2)
Rb = R_mann15(MKb) * (1 + rng.normal(0, 0.0289, N)); Mb = M_mann19(MKb) * (1 + rng.normal(0, 0.02, N)) * (1 + rng.normal(0, 0.02, N))
Lb = L / 2; Teffb = 5772.0 * (Lb / Rb ** 2) ** 0.25; rhob = Mb / Rb ** 3 * 1.40987
def st(x): return dict(median=float(np.median(x)), err_lo=float(np.median(x) - np.percentile(x, 15.87)), err_hi=float(np.percentile(x, 84.13) - np.median(x)), mean=float(np.mean(x)), sd=float(np.std(x)))
P = dict(
    star='TIC 214311908 = Gaia DR3 5925373272850784768 = 2MASS J17235292-5119271',
    notes=['AI-run analysis 2026-10-16; adopted values assume a SINGLE main-sequence star. The star lies 0.72 mag above the empirical MKs-(BP-RP) ridge '
           '(5.6% of comparison stars as bright); explained by [Fe/H] ~ +0.2-0.4 if single, or by a near-equal-mass companion within ~5 au (not excluded). '
           'See binary_twin_scenario.',
           'Teff adopted from Stefan-Boltzmann (L from SED Fbol, R from Mann+15 MKs relation).'],
    gaia_dr3=dict(ra_deg=260.969674, dec_deg=-51.323821, epoch=2016.0, parallax_mas=15.497306, parallax_err_mas=0.025122, parallax_zp_L21_mas=-0.052059,
                  pmra=-113.671152, pmra_err=0.031124, pmdec=85.54355, pmdec_err=0.0224, ruwe=1.076204, G=14.281043, BP=15.830741, RP=13.074592, bp_rp=2.756149,
                  rv_kms=16.283682, rv_err_kms=2.217367, rv_nb_transits=21, ipd_frac_multi_peak=0, ipd_gof_harmonic_amplitude=0.013223,
                  phot_bp_rp_excess_factor=1.452786, Cstar_sigma=0.71, non_single_star=0, phot_variable_flag='NOT_AVAILABLE'),
    photometry=dict(J_2MASS=[11.363, 0.028], H_2MASS=[10.740, 0.031], Ks_2MASS=[10.482, 0.026], SkyMapper_DR4_AB=dict(g=[15.817, 0.030], r=[14.908, 0.024], i=[13.387, 0.010], z=[12.768, 0.010], v=[18.192, 0.059]),
                    DECaPS2_AB=dict(g=16.177, r=14.737, i=13.355, z=12.761, Y=12.563), unWISE_W1=10.33, unWISE_W2=10.23,
                    AllWISE=dict(W1=[10.112, 0.022], W2=[9.884, 0.020], W3=[9.380, 0.037], W4=[7.152, 0.101], note='ccf=hhh0: halo of the Mira (W1 6.5) 17" away; W1/W2 0.28/0.33 mag too bright vs model, W3/W4 dominated by the Mira; not used'),
                    TESS_mag=dict(TIC=13.0135, SED_synthetic=sed['check_bands']['TESS']['model'], Stassun19_Gaia=13.134)),
    parallax_corrected_mas=st(plx), distance_pc=st(d), MKs=st(MK), MG=float(14.281043 + 5 * np.log10(plx0 / 100)),
    radius_rsun=st(R), radius_rsun_with_feh=st(Rfe), mass_msun=st(M), luminosity_lsun=st(L), fbol_cgs=st(Fbol),
    teff_K=st(Teff_SB), logg_cgs=st(logg), density_gcc=st(rho), density_solar=st(rho / 1.40987),
    teff_checks=dict(SED_BTSettl=[sed['mc']['teff'][0], 100.0], Mann15_VminusJ=[float(Teff_VJ), 60.0], Mann15_BPRP=[float(Teff_BPRP), 60.0], TIC_v8=[3301, 157], Gaia_msc1=3298),
    radius_checks=dict(SED_theta_d=[sed['mc']['R_Rsun'][0], 0.02], TIC_v8=[0.3807, 0.0117]),
    feh_phot_if_single=feh['feh_phot_if_single'], spectral_type='M3V (BP-RP 2.76 between Mamajek M3V 2.50 and M3.5V 2.78)',
    binary_twin_scenario=dict(description='if the light is shared equally by two unresolved stars (not excluded inside ~0.08\"/5 au)',
                              MKs_each=st(MKb), radius_each_rsun=st(Rb), mass_each_msun=st(Mb), teff_each_K=st(Teffb), density_each_gcc=st(rhob),
                              transit_dilution_factor=2.0),
)
json.dump(P, open(f'{B}/star_params.json', 'w'), indent=1)
np.savez_compressed(f'{B}/out/star_params_mc.npz', R=R[:20000], M=M[:20000], L=L[:20000], Teff=Teff_SB[:20000], rho=rho[:20000], d=d[:20000])
for k in ['distance_pc', 'MKs', 'radius_rsun', 'radius_rsun_with_feh', 'mass_msun', 'luminosity_lsun', 'teff_K', 'logg_cgs', 'density_gcc', 'density_solar']:
    v = P[k]; print(f"{k:22s} {v['median']:.4f} -{v['err_lo']:.4f} +{v['err_hi']:.4f}")
print('Teff checks', P['teff_checks'])
for k in ['radius_each_rsun', 'mass_each_msun', 'teff_each_K', 'density_each_gcc']:
    v = P['binary_twin_scenario'][k]; print(f"twin {k:18s} {v['median']:.4f} +- {v['sd']:.4f}")
