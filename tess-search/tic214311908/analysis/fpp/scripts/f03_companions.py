"""f03: Monte-Carlo population of possible unresolved bound companions of TIC 214311908 and the likelihood of each
companion under the non-transit data. Used by f04 for the TP (diluted), TP-twin, PTP (planet on companion) and
companion-EB scenarios.

Prior population (Winters et al. 2019, M dwarfs within 25 pc): a companion exists with probability f_mult (f04);
log10(a/au) ~ N(1.30, 1.16) (log-normal peak 20 au); e ~ U(0, 0.8); isotropic orientation, uniform mean anomaly ->
projected separation at the Gaia/DECaPS2 epochs; q = m2/m1 drawn uniform in [0.075/m1, 1] here and re-weighted in f04
for other q distributions (power law, twin excess). Primary mass m1 is solved so that the system's total M_Ks equals
the observed 6.4406 (Mann+19 mass-M_Ks relation is metallicity-insensitive). Companion photometry from the Mamajek
sequence (T - Ks, M_V - (V - Ic) for I), radii from Mann+15 (M_Ks).

Likelihood of each configuration under non-transit data (stored separately so f04 can switch them):
* L_ipd: Gaia IPD/RUWE calibration cells (star/out/ipd_limits.json), LR per (separation, Delta I) cell for 0.08-2.5";
  cells with n < 5 or no data -> 1. Variant L_ipd84 uses the 84% upper bound of each LR (conservative).
* L_ast: RUWE + DR3-DR2 proper-motion simulation (star/out/astrometric_limits.json, L_combined on an M2 x a grid,
  sigma_AL as in the investigator's simulation); companions above the grid's 0.34 Msun use the 0.34 row; a > 31.6 au -> 1.
  The referee (referee_star r10) found the grid too strict for 0.3 Msun and for sigma_CCD = 0.35 mas: variant L_ast^0.4
  (errors x1.6) in f04.
* L_dec: DECaPS2 2018 residual images, 5-sigma contrast as corrected by the referee (r04: i 4.24/5.07/5.80/5.58,
  z 3.96/4.53/4.96/5.24 at 1/1.25/1.5/2"; investigator values beyond 2.5"; 2.0 mag at 0.75" (referee: ~10% of
  injected flux recovered there); none inside 0.75"):
  companions brighter than the limit get L = 0.02.
* L_cmd: the system's M_Ks offset from the main-sequence ridge at its BP-RP is -0.753 mag (referee_star r09). A companion
  brightens the system by dM_bin(q) = -2.5 log10(1 + 10^(-0.4 dKs)); the rest must be the primary's own offset eps
  (metallicity, age, errors). The CMD-offset distribution of ~1000 nearby M dwarfs with RUWE < 1.1 at the same colour
  (referee_star/data/mdw_cmd.csv; ridge = KDE mode) was examined first, with the idea of mirroring its faint side.
  NOTE: that faint side has a long tail (metal-poor stars, bad photometry: 97.5th percentile +1.4 mag) while the bright
  side is compact (16th percentile -0.27 mag, 1st percentile -0.86 mag including binaries), so the mirrored faint side
  would grossly overstate metal-rich singles. Adopted: Gaussian eps ~ N(0, sigma_s), sigma_s = 0.30 mag (bright-side
  width; f04 varies 0.25-0.40). The empirical percentiles are recorded in the summary.
Outputs: out/companions_mc.npz (draws + likelihood terms), out/companions_summary.json"""
import numpy as np, json, sys, pandas as pd
sys.path.insert(0, 'tic214311908/analysis/fpp/scripts')
sys.path.insert(0, 'tic214311908/analysis/shape/scripts')
from fpplib import *
import eem
rng = np.random.default_rng(20261017)
N = 400000
MKS_OBS = 6.4406; DM_OBS = -0.753

# ---------- companion draws
loga = rng.normal(1.30, 1.16, N); a = 10**loga
e = rng.uniform(0, 0.8, N)
cosi = rng.uniform(-1, 1, N); Om = rng.uniform(0, 2*np.pi, N); om = rng.uniform(0, 2*np.pi, N); Mano = rng.uniform(0, 2*np.pi, N)
# solve Kepler for projected separation
E = Mano.copy()
for _ in range(30): E = E - (E - e*np.sin(E) - Mano)/(1 - e*np.cos(E))
nu = 2*np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2))
r = a*(1 - e*np.cos(E))
sinI = np.sqrt(1 - cosi**2)
X = r*(np.cos(Om)*np.cos(om+nu) - np.sin(Om)*np.sin(om+nu)*cosi)
Y = r*(np.sin(Om)*np.cos(om+nu) + np.cos(Om)*np.sin(om+nu)*cosi)
s_au = np.hypot(X, Y); sep = s_au/DIST_PC
u = rng.uniform(0, 1, N)   # q drawn uniform in [qmin, 1] (qmin depends on m1: iterate)

def mks_of_m(m):  # Mann+19 inverse, vectorised
    mk = np.linspace(4.0, 11.0, 4000); mm = eem.mann_mass(mk)
    return np.interp(m, mm[::-1], mk[::-1])
m1 = np.full(N, 0.37)
for it in range(40):
    qmin = 0.075/m1
    q = qmin + u*(1 - qmin)
    m2 = q*m1
    mk1 = mks_of_m(m1); mk2 = mks_of_m(m2)
    tot = -2.5*np.log10(10**(-0.4*mk1) + 10**(-0.4*mk2))
    # adjust m1 so that tot -> MKS_OBS (primary fainter by the companion's share)
    mk1_new = mk1 + (MKS_OBS - tot)
    m1 = eem.mann_mass(mk1_new)
mk1 = mks_of_m(m1); mk2 = mks_of_m(m2)
dKs = mk2 - mk1
R1 = eem.mann_radius(mk1); R2 = eem.mann_radius(mk2)
rho1 = m1*MSUN/(4/3*np.pi*(R1*RSUN)**3); rho2 = m2*MSUN/(4/3*np.pi*(R2*RSUN)**3)
TK = lambda m: eem.ms(m, 'T_minus_Ks')
MI = lambda m: eem.ms(m, 'Mv') - eem.ms(m, 'V-Ic')
MKS_tab = lambda m: eem.ms(m, 'M_Ks')
T1 = mk1 + TK(m1); T2 = mk2 + TK(m2); dT = T2 - T1
dI = (MI(m2) - MKS_tab(m2) + mk2) - (MI(m1) - MKS_tab(m1) + mk1)
Dc = 1/(1 + 10**(0.4*dT))     # companion's share of the system TESS flux
Teff1 = eem.ms_Teff(m1); Teff2 = eem.ms_Teff(m2)
dM_bin = -2.5*np.log10(1 + 10**(-0.4*dKs))

# ---------- L_ipd
ipd = jload('star/out/ipd_limits.json')['cells']
def ipd_lr(sep, dm, key='LR'):
    out = np.ones(len(sep))
    for c in ipd:
        if c['n'] < 5 or not (c[key] == c[key]): continue
        m = (sep >= c['sep'][0]) & (sep < c['sep'][1]) & (dm >= c['dm'][0]) & (dm < c['dm'][1])
        out[m] = min(1.0, c[key])
    return out
L_ipd = ipd_lr(sep, dI, 'LR'); L_ipd84 = ipd_lr(sep, dI, 'LR_hi84')
# ---------- L_ast: the grid (star/s12c) assumed a 0.37 Msun primary and companion light from the Mamajek G sequence.
# The RUWE / PM signal scales with the photocentre wobble a*|B - beta| (B = mass fraction, beta = G-band light fraction
# of the companion), which vanishes for twins. Each draw is mapped to the grid M2 (low-mass branch) with the same
# |B - beta|, at the draw's a.
ast = jload('star/out/astrometric_limits.json')
AM = np.array(ast['M2']); AA = np.log10(np.array(ast['a_au'])); AL = np.array(ast['L_combined'])
FG = lambda m: 10**(-0.4*eem.ms(m, 'M_G'))
dg = AM/(0.37 + AM) - FG(AM)/(FG(0.37) + FG(AM))
ipk = int(np.argmax(dg))
from scipy.interpolate import RegularGridInterpolator
fi = RegularGridInterpolator((AM, AA), AL, bounds_error=False, fill_value=None)
Bm = m2/(m1 + m2); bG = FG(m2)/(FG(m1) + FG(m2)); dlt = np.abs(Bm - bG)
M2eq = np.interp(dlt, dg[:ipk+1], AM[:ipk+1])           # equivalent grid mass on the rising branch
M2eq = np.where(dlt > dg[ipk], AM[ipk], M2eq)
L_ast = fi(np.column_stack([np.clip(M2eq, AM[0], AM[-1]), np.clip(loga, AA[0], AA[-1])]))
L_ast[loga > AA[-1]] = 1.0
L_ast = np.clip(L_ast, 1e-4, None)
L_ast_raw = fi(np.column_stack([np.clip(m2, AM[0], AM[-1]), np.clip(loga, AA[0], AA[-1])])); L_ast_raw[loga > AA[-1]] = 1.0
# ---------- L_dec
dec_s = np.array([0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0])
dec_v = np.array([2.0, 4.24, 5.07, 5.80, 5.58, 6.17, 6.89, 7.53, 7.88, 8.59])
lim = np.where(sep >= 0.75, np.interp(sep, dec_s, dec_v), -99)
L_dec = np.where(dT < lim, 0.02, 1.0)
# ---------- L_cmd: empirical single-star offset distribution
d = pd.read_csv(f'{B}/referee_star/data/mdw_cmd.csv')
d = d[(d.ruwe < 1.1) & (np.abs(d.bp_rp - 2.756) < 0.06) & np.isfinite(d.ks_m) & (d.ph_qual.astype(str).str[2] == 'A')]
MKs_s = d.ks_m + 5*np.log10(d.parallax/100.)
# ridge = mode via KDE
from scipy.stats import gaussian_kde
xs = np.linspace(MKs_s.min(), MKs_s.max(), 2000)
mode = xs[np.argmax(gaussian_kde(MKs_s, 0.15)(xs))]
# colour-term: local slope from a linear fit of MKs vs bp_rp around the target colour, applied before taking offsets
sl = np.polyfit(d.bp_rp, MKs_s, 1)[0]
off = MKs_s - (mode + sl*(d.bp_rp - 2.756))
mode2 = xs[np.argmax(gaussian_kde(off, 0.12)(xs - xs.mean()))] if False else 0
faint = off[off > 0].values
frac_bright = float(np.mean(off <= DM_OBS))
pct = np.percentile(off, [1, 2.5, 16, 50, 84, 97.5, 99]).tolist()
SIG_S = 0.30   # adopted single-star scatter on the bright (metal-rich) side; f04 varies 0.25-0.40
def p_eps_gauss(x, s=SIG_S): return np.exp(-0.5*(x/s)**2)/np.sqrt(2*np.pi)/s
eps_needed = DM_OBS - dM_bin
L_cmd_g = p_eps_gauss(eps_needed)
L_cmd_single_g = float(p_eps_gauss(DM_OBS))
print(f'CMD sample n={len(off)}, slope {sl:.2f}, frac as bright as {DM_OBS}: {frac_bright:.4f}; percentiles 1/2.5/16/50/84/97.5/99: {np.round(pct,2)}')
np.savez(f'{F}/out/companions_mc.npz', a=a, loga=loga, e=e, sep=sep, s_au=s_au, q_u=u, q=q, m1=m1, m2=m2, R1=R1, R2=R2,
         rho1=rho1, rho2=rho2, dKs=dKs, dT=dT, dI=dI, Dc=Dc, Teff1=Teff1, Teff2=Teff2, dM_bin=dM_bin,
         L_ipd=L_ipd, L_ipd84=L_ipd84, L_ast=L_ast, L_ast_raw=L_ast_raw, M2eq=M2eq, L_dec=L_dec, eps_needed=eps_needed,
         L_cmd_g=L_cmd_g, L_cmd_single_g=L_cmd_single_g, sig_s=SIG_S)
# ---------- summary (uniform-q prior, all constraints, empirical CMD)
W = L_ipd*L_ast*L_dec
summ = dict(N=N, cmd=dict(n=int(len(off)), slope=float(sl), frac_as_bright=frac_bright, percentiles_1_2p5_16_50_84_97p5_99=pct,
                          sig_s_adopted=SIG_S, p_single_over_p_twin=float(p_eps_gauss(DM_OBS)/p_eps_gauss(0.0))),
            frac_prior_mass_surviving={})
for name, w in [('ipd', L_ipd), ('ast', L_ast), ('dec', L_dec), ('ipd*ast*dec', W), ('all_with_cmd', W*L_cmd_g/L_cmd_single_g)]:
    summ['frac_prior_mass_surviving'][name] = float(np.mean(w))
bins = [(0, 0.3), (0.3, 0.75), (0.75, 1.5), (1.5, 2.5), (2.5, 3.5), (3.5, 4.5), (4.5, 99)]
summ['by_dT'] = {}
for lo, hi in bins:
    m = (dT >= lo) & (dT < hi)
    summ['by_dT'][f'{lo}-{hi}'] = dict(prior_frac=float(m.mean()), surv_noCMD=float(np.mean(W*m)), surv_withCMD=float(np.mean(W*m*L_cmd_g/L_cmd_single_g)))
summ['by_a'] = {}
for lo, hi in [(0, 0.15), (0.15, 0.5), (0.5, 4), (4, 6), (6, 65), (65, 200), (200, 1e5)]:
    m = (a >= lo) & (a < hi)
    summ['by_a'][f'{lo}-{hi}'] = dict(prior_frac=float(m.mean()), surv_noCMD=float(np.mean(W*m)), surv_withCMD=float(np.mean(W*m*L_cmd_g/L_cmd_single_g)))
jdump(summ, 'out/companions_summary.json')
print(json.dumps(summ, indent=1))
