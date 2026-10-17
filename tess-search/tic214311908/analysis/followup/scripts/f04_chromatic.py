"""f04: chromatic predictions for the ground-based depth-ratio test, by synthetic photometry.
Spectra: BT-Settl CIFIST (star/data/models, [M/H]=0) for M dwarfs and Teff <= 4100 K; Castelli-Kurucz 2004 (shape/data/ck04)
for the background K dwarf; blackbody for a white dwarf. Filter curves: SVO (data/filters: Johnson V, Cousins I, SDSS g r i z,
PS1 z, 2MASS J Ks, TESS) and star/data/models/filters_all.npz (DECam g r i z); top-hat proxies (with a deep-depletion /
generic QE) for SPECULOOS I+z (750-1050 nm), ASTEP+ blue (410-700) / red (712-988), ExTrA (850-1550). Photon-weighted.
Hypotheses (depth in band X relative to the TESS depth 15.2 ppt):
 H1 planet on the M3V: k achromatic; raw mid-transit depth differs only by limb darkening (Claret & Bloemen 2011 PHOENIX
    quadratic LD, Teff 3300-3400 K; Claret 2017 for TESS; b = 0.5).
 H2 PTP: planet on an unresolved bound M-dwarf companion: depth_X ∝ F_c/(F_M3V + F_c) in X; companion Teff found from dT
    (ratio curves vs dT tabulated for use with the referee2_fpp posterior draws in f07).
 H3 WD totally occulted by the M3V: depth = F_WD/(F_M3V + F_WD); WD radius 0.0125 Rsun, Teff fitted to the TESS depth.
 H4 unresolved background K-dwarf EB (Teff 4500 K, E(B-V) 0.2-0.4, CCM R_V 3.1) diluted by the M3V.
 H5 the 3.1" star (G 17.48, BP-RP 1.32) is the eclipsing source: its own eclipse (51% in T) seen resolved, and the blended
    depth in apertures that include it.
 H6 systematic for a real planet: unocculted cool spots (Tspot 2950/3150 K, filling 5-20%) make the transit deeper in the
    blue relative to TESS (opposite sign to H2).
Outputs: out/chromatic.json, out/chromatic_ratio_vs_dT.csv, out/m3v_band_mags.json"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flib import *
import numpy as np, pandas as pd, glob, json
from astropy.io import fits
MOD = f'{B}/star/data/models'
lam = np.load(f'{MOD}/btsettl_3300_5.0.npy')[0].astype(float)
def bts(teff, logg=5.0):
    ts = sorted(set(int(os.path.basename(f).split('_')[1]) for f in glob.glob(f'{MOD}/btsettl_*_{logg:.1f}.npy')))
    lo = max(t for t in ts if t <= teff); hi = min(t for t in ts if t >= teff)
    a = np.load(f'{MOD}/btsettl_{lo}_{logg:.1f}.npy')[1].astype(float)
    if hi == lo: return a
    b = np.load(f'{MOD}/btsettl_{hi}_{logg:.1f}.npy')[1].astype(float)
    x = (teff - lo)/(hi - lo)
    return np.exp((1 - x)*np.log(np.clip(a, 1e-30, None)) + x*np.log(np.clip(b, 1e-30, None)))
def ck(teff, logg='g45'):
    h = fits.open(f'{B}/shape/data/ck04/ckp00_{teff}.fits')[1].data
    return np.interp(lam, h['WAVELENGTH'], h[logg])
def bb(T):
    h, c, k = 6.626e-27, 2.998e10, 1.381e-16
    l = lam*1e-8
    return 2*h*c**2/l**5/np.expm1(h*c/(l*k*T))*1e-8      # per Angstrom (arbitrary units are fine for ratios at fixed band)
def ccm(lam_A, Rv=3.1):
    x = 1e4/lam_A; a = np.zeros_like(x); b = np.zeros_like(x)
    ir = (x >= 0.3) & (x < 1.1); a[ir] = 0.574*x[ir]**1.61; b[ir] = -0.527*x[ir]**1.61
    op = (x >= 1.1) & (x <= 3.3); y = x[op] - 1.82
    a[op] = 1 + 0.17699*y - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 + 0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7
    b[op] = 1.41338*y + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4 - 0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7
    return a + b/Rv            # A_lambda / A_V
# ---- filters
FD = f'{F}/data/filters'
def svo(name):
    d = np.loadtxt(f'{FD}/{name}.dat'); return np.interp(lam, d[:, 0], d[:, 1], left=0, right=0)
fa = np.load(f'{MOD}/filters_all.npz')
def npzf(k): d = fa[k]; return np.interp(lam, d[0], d[1], left=0, right=0)
def tophat(l0, l1, qe=None):
    t = ((lam >= l0) & (lam <= l1)).astype(float)
    if qe == 'dd':   # generic deep-depletion CCD QE falling beyond 900 nm (approximate)
        t *= np.clip(np.interp(lam, [7000, 9000, 10000, 10500, 11000], [0.95, 0.9, 0.5, 0.25, 0.05]), 0, 1)
    return t
FILT = {'g': svo('SLOAN_SDSS.g'), 'r': svo('SLOAN_SDSS.r'), 'i': svo('SLOAN_SDSS.i'), 'z_s': svo('PAN-STARRS_PS1.z'),
        'z': svo('SLOAN_SDSS.z'), 'V': svo('Generic_Johnson.V'), 'Ic': svo('Generic_Cousins.I'), 'J': svo('2MASS_2MASS.J'),
        'Ks': svo('2MASS_2MASS.Ks'), 'TESS': svo('TESS_TESS.Red'), 'I+z': tophat(7500, 11000, 'dd'),
        'ASTEP_blue': tophat(4100, 7000), 'ASTEP_red': tophat(7120, 9880, 'dd'), 'ExTrA': tophat(8500, 15500),
        'G': npzf('GAIA__GAIA3.G'), 'BP': npzf('GAIA__GAIA3.Gbp'), 'RP': npzf('GAIA__GAIA3.Grp'), 'DECam_g': npzf('CTIO__DECam.g')}
def band(f, T):  # photon-weighted band flux
    return float(np.trapezoid(f*T*lam, lam))
# ---- M3V model normalised to R = 0.383; Mamajek radius vs Teff for companions
EEM = [(3430, 0.361), (3270, 0.300), (3210, 0.274), (3110, 0.217), (3060, 0.196), (2980, 0.156), (2940, 0.137), (2860, 0.126)]
ee = np.array(EEM)[::-1]
def R_of_T(T): return float(np.interp(T, ee[:, 0], ee[:, 1]))
fM = bts(3350)*RSTAR**2
FM = {b: band(fM, T) for b, T in FILT.items()}
# AB magnitudes of the M3V in every band: scale the model to the observed i (SDSS-like, gspc 13.57) -> check vs others
ab0 = {b: band(3.631e-20*2.998e18/lam**2, T) for b, T in FILT.items()}   # f_nu = 3631 Jy -> f_lambda [erg/s/cm2/A]
scale = 10**(-0.4*13.57)/(FM['i']/ab0['i'])
mAB = {b: float(-2.5*np.log10(scale*FM[b]/ab0[b])) for b in FILT}
obs_check = dict(g=16.36, r=14.98, i=13.57, z=12.81, V_AB=15.59 + 0.02, Ic_AB=12.86 + 0.45, J_AB=11.36 + 0.91, Ks_AB=10.48 + 1.85)
print('M3V model AB mags (scaled to i):', {k: round(v, 2) for k, v in mAB.items()}); print('observed:', obs_check)
# ---- H1 planet: limb darkening
ld = pd.read_csv(f'{F}/data/claret2011_ld_3300_3400.csv'); ld = ld[(ld.Met == 'L') & (ld.logg == 5.0)]
def ldq(filt):
    s = ld[ld.Filt == filt]; return float(s.a.mean()), float(s.b.mean())
cl17 = pd.read_csv(f'{B}/shape/data/claret17_tableab.csv'); cl17 = cl17[(cl17.Teff >= 3300) & (cl17.Teff <= 3400) & (cl17.logg == 5.0) & (cl17.Z == 0)]
LDc = {'TESS': (float(cl17.aLSM.mean()), float(cl17.bLSM.mean())), 'g': ldq('g*'), 'r': ldq('r*'), 'i': ldq('i*'), 'z': ldq('z*'), 'z_s': ldq('z*'),
       'V': ldq('V'), 'Ic': ldq('I'), 'J': ldq('J'), 'Ks': ldq('K'), 'I+z': ldq('z*'), 'ASTEP_blue': ldq('V'), 'ASTEP_red': ldq('I'), 'ExTrA': ldq('J')}
def mid_depth(k, b, u1, u2):
    # small-planet approximation of the mid-transit depth: k^2 I(mu)/<I>
    mu = np.sqrt(1 - b*b); I = 1 - u1*(1 - mu) - u2*(1 - mu)**2; Im = 1 - u1/3 - u2/6
    return k*k*I/Im
k = np.sqrt(DEPTH/ (mid_depth(1, 0.5, *LDc['TESS'])))
H1 = {b: mid_depth(k, 0.5, *LDc[b])/DEPTH for b in LDc}
# ---- H2 PTP: ratio vs dT
rows = []
for Tc in np.arange(2700, 3351, 10):
    Rc = R_of_T(Tc) if Tc < 3350 else RSTAR
    fc = bts(Tc)*Rc**2
    Fc = {b: band(fc, T) for b, T in FILT.items()}
    dT = -2.5*np.log10(Fc['TESS']/FM['TESS'])
    dil = {b: Fc[b]/(FM[b] + Fc[b]) for b in FILT}
    row = dict(Teff_c=Tc, R_c=Rc, dT=dT)
    for b in FILT:
        row[f'ratio_{b}'] = dil[b]/dil['TESS']; row[f'd_{b}'] = -2.5*np.log10(Fc[b]/FM[b])
    rows.append(row)
R2 = pd.DataFrame(rows); R2.to_csv(f'{OUT}/chromatic_ratio_vs_dT.csv', index=False)
def at_dT(dT, col): s = R2.sort_values('dT'); return float(np.interp(dT, s.dT, s[col]))
H2 = {q: {b: at_dT(dT, f'ratio_{b}') for b in FILT} for q, dT in (('dT0.88', 0.88), ('dT1.16', 1.16), ('dT1.45', 1.45))}
H2_VoverI = {q: H2[q]['V']/H2[q]['Ic'] for q in H2}
# ---- H3 WD occultation
def wd_depths(Twd, Rwd=0.0125):
    fw = bb(Twd)*1.0
    # normalise blackbody to the same units as BT-Settl surface flux: pi*B_lambda; BT-Settl files are surface flux [erg/s/cm2/A]
    fw = np.pi*fw*Rwd**2
    return {b: band(fw, T)/(FM[b] + band(fw, T)) for b, T in FILT.items()}
from scipy.optimize import brentq
units_check = band(np.pi*bb(3350)*RSTAR**2, FILT['TESS'])/FM['TESS']
Twd = brentq(lambda T: wd_depths(T)['TESS'] - DEPTH, 3000, 60000)
H3 = {b: v/DEPTH for b, v in wd_depths(Twd).items()}
H3_alt = {}
for Rw in (0.010, 0.015):
    Tw = brentq(lambda T: wd_depths(T, Rw)['TESS'] - DEPTH, 3000, 80000)
    H3_alt[f'R_wd={Rw}'] = dict(Teff=Tw, g_over_TESS=wd_depths(Tw, Rw)['g']/DEPTH, V_over_TESS=wd_depths(Tw, Rw)['V']/DEPTH)
# ---- H4 background K dwarf EB (small dilution)
H4 = {}
for ebv in (0.2, 0.4):
    fk = ck(4500)*10**(-0.4*3.1*ebv*ccm(lam))
    Fk = {b: band(fk, T) for b, T in FILT.items()}
    H4[f'E(B-V)={ebv}'] = {b: (Fk[b]/FM[b])/(Fk['TESS']/FM['TESS']) for b in FILT}
# ---- H5 the 3.1" star (observed magnitudes, AB/Vega as tabulated) as the source
m31 = dict(g=18.16, r=17.35, i=17.10, z=16.90, V=17.83, Ic=16.64, TESS=16.82)
mM = dict(g=16.36, r=14.98, i=13.57, z=12.81, V=15.59, Ic=12.86, TESS=13.01)
ecl31 = DEPTH*10**(-0.4*(mM['TESS'] - m31['TESS']))*(1 + 10**(-0.4*(m31['TESS'] - mM['TESS'])))   # eclipse depth in the 3.1" star's own light
H5 = dict(own_eclipse_depth_T=ecl31,
          blended_depth_if_both_in_aperture={b: ecl31*10**(-0.4*(m31[b] - mM[b]))/(1 + 10**(-0.4*(m31[b] - mM[b]))) for b in m31},
          note='own eclipse assumed grey (EB of similar stars); blended depth = depth diluted by the M3V when both are fully in the aperture')
# ---- H6 (systematic, not a false positive): unocculted cool spots on the active M3V deepen the transit more in the blue
H6 = {}
fphot = bts(3350)
for Tsp in (2950, 3150):
    fs = bts(Tsp)
    for ff in (0.05, 0.10, 0.20):
        dep = {b: 1/(1 - ff*(1 - band(fs, T)/band(fphot, T))) for b, T in FILT.items()}
        H6[f'Tspot={Tsp}_f={ff}'] = {b: dep[b]/dep['TESS'] for b in FILT}
res = dict(H6_unocculted_spots_depth_over_TESS=H6, M3V_model_AB_mags=mAB, M3V_observed=obs_check, LD_quadratic=LDc,
           H1_planet_raw_depth_over_TESS=H1, H2_PTP_depth_over_TESS=H2, H2_V_over_I=H2_VoverI,
           H3_WD=dict(Teff_wd=Twd, R_wd=0.0125, depth_over_TESS=H3, variants=H3_alt, units_check_bb_vs_btsettl_3350=units_check),
           H4_bg_Kdwarf_EB=H4, H5_star3p1=H5,
           notes=['ratios are depth_X / depth_TESS for the same event; for H1 they are LD effects only (k achromatic)',
                  'H2 uses BT-Settl [M/H]=0 with Mamajek radii; the M3V is metal-rich, so absolute colours are approximate but ratios are differential',
                  'z_s = PS1 z (MuSCAT zs proxy); I+z, ASTEP and ExTrA bands are top-hat proxies'])
jdump(res, f'{OUT}/chromatic.json')
jdump(dict(AB=mAB, note='BT-Settl 3350 K scaled to the XP-synthetic i = 13.57; use for count rates in custom bands'), f'{OUT}/m3v_band_mags.json')
show = ['g', 'r', 'i', 'z_s', 'V', 'Ic', 'I+z', 'J', 'ASTEP_blue', 'ASTEP_red', 'ExTrA']
print('H1 planet raw depth/TESS', {b: round(H1[b], 3) for b in show if b in H1})
for q in H2: print('H2', q, {b: round(H2[q][b], 3) for b in show}, 'V/I', round(H2_VoverI[q], 3))
print('H3 WD Teff', round(Twd), {b: round(H3[b], 2) for b in show}, H3_alt, 'units check', units_check)
for e in H4: print('H4', e, {b: round(H4[e][b], 2) for b in show})
print('H5', H5)
for k_ in H6: print('H6', k_, {b: round(H6[k_][b], 3) for b in show})
