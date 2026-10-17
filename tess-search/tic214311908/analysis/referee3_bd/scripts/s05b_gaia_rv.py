"""s05b: own Gaia DR3 RVS-error constraint on a BD orbiting the M3V (rv_method_used = 2: combined CCF of 21 transits;
rv_amplitude_robust / rv_chisq_pvalue / rv_renormalised_gof are not provided for such stars - checked in data/gaia_target.csv).
Model: each transit CCF = C0 exp(-v^2/2w^2) shifted by the orbital RV at that transit; combined CCF = mean of the 21;
error from the Zucker (2003) curvature formula at the combined peak; ratio r = err(K)/err(0). Epoch sampling:
(i) 21 independent random phases; (ii) Gaia-like pairs (two FOV transits 106.5 min apart: 10 pairs + 1 single).
Intrinsic error density p0: KDE (log err) of my own comparison sample (s05a: rv_method_used = 2, |grvs - 12.55| < 0.3,
rv_template_teff <= 3700, parallax > 5, 2.4 < BP-RP < 3.1, 15-30 transits), additionally restricted to
|rv_expected_sig_to_noise/14.70 - 1| < 0.2. LR(K) = < p0(err_obs/r)/r > / p0(err_obs).
BD on the M3V: K = 28.4329 m/s (M/MJ)(P/yr)^-1/3 ((M*+M)/Msun)^-2/3 (1-e^2)^-1/2, sin i = 1; e = 0 and e = 0.55 with
w = 270 deg (the configuration a grazing BD needs: transit near apoastron). Mass priors: uniform 13-80 MJ; the
empirical masses of M-dwarf transiting BDs with P < 16 d (bd/out/b01_mdwarf_bds.csv values, checked against abstracts
for TOI-5610 b / TOI-5389Ab / TOI-6508 b in lit/).
Also: photocentre semi-axis alpha for the RUWE argument.
Output: out/s05_gaia_rv.json"""
import sys; sys.path.insert(0, 'tic214311908/analysis/referee3_bd/scripts')
from rlib3 import *
import pandas as pd
from scipy.stats import gaussian_kde
rng = np.random.default_rng(21)
tg = pd.read_csv(f'{R3}/data/gaia_target.csv'); tg = tg[tg.source_id == 5925373272850784768].iloc[0]
ERR, SN = float(tg.radial_velocity_error), float(tg.rv_expected_sig_to_noise)
cmp_ = pd.read_csv(f'{R3}/data/gaia_rvcomp.csv')
cmp_ = cmp_[cmp_.source_id != 5925373272850784768]
sel = cmp_[np.abs(cmp_.rv_expected_sig_to_noise/SN - 1) < 0.2]
out = dict(target=dict(err=ERR, sn=SN, grvs=float(tg.grvs_mag), nb_transits=int(tg.rv_nb_transits), template_teff=float(tg.rv_template_teff),
                       rv=float(tg.radial_velocity), rv_amplitude_robust=str(tg.rv_amplitude_robust), rv_chisq_pvalue=str(tg.rv_chisq_pvalue)),
           n_comp_all=int(len(cmp_)), n_comp_sel=int(len(sel)),
           pct_all={p: float(np.percentile(cmp_.radial_velocity_error, p)) for p in (1, 5, 16, 50, 84)},
           pct_sel={p: float(np.percentile(sel.radial_velocity_error, p)) for p in (1, 5, 16, 50, 84)},
           target_percentile_all=float(np.mean(cmp_.radial_velocity_error < ERR)*100), target_percentile_sel=float(np.mean(sel.radial_velocity_error < ERR)*100))
print(out, flush=True)
kde = gaussian_kde(np.log(sel.radial_velocity_error.values))
p0 = lambda e_: kde(np.log(e_))/e_
v = np.linspace(-200, 200, 2001); dv = v[1] - v[0]
def zerr(C):
    i = int(np.argmax(C)); c = C[i]; c2 = (C[i + 1] - 2*c + C[i - 1])/dv**2
    return np.sqrt(-1/(c2/c*c**2/(1 - c**2)))
def rv_curve(K, e, w, M):
    E = M.copy()
    for _ in range(40): E = E - (E - e*np.sin(E) - M)/(1 - e*np.cos(E))
    f = 2*np.arctan2(np.sqrt(1 + e)*np.sin(E/2), np.sqrt(1 - e)*np.cos(E/2))
    return K*(np.cos(f + w) + e*np.cos(w))
def phases(mode):
    if mode == 'indep': return rng.uniform(0, 2*np.pi, 21)
    m = rng.uniform(0, 2*np.pi, 11); d = 2*np.pi*106.5/1440/PC
    return np.r_[m, m[:10] + d]
Ks = np.array([2., 4., 6., 8., 10., 12., 14., 16., 20.])
LR = {}
ND = 400
for w_ in (13., 16., 20.):
    for C0 in (0.3, 0.6):
        base = zerr(C0*np.exp(-0.5*(v/w_)**2))
        for e, wdeg in ((0.0, 90.), (0.55, 270.)):
            for mode in ('indep', 'pairs'):
                lr, rmed = [], []
                for K in Ks:
                    rr = np.array([zerr(C0*np.mean(np.exp(-0.5*((v[:, None] - rv_curve(K, e, np.radians(wdeg), phases(mode))[None, :])/w_)**2), axis=1))/base for _ in range(ND)])
                    lr.append(float(np.mean(p0(ERR/rr)/rr)/p0(ERR)[0])); rmed.append(float(np.median(rr)))
                LR[f'w{int(w_)}_C{C0}_e{e}_{mode}'] = dict(K=Ks, LR=lr, r_median=rmed)
                print(w_, C0, e, mode, 'LR', np.round(lr, 3), 'r_med', np.round(rmed, 2), flush=True)
out['LR'] = LR
# mass-averaged factor for a BD on the M3V
def K_of(M_MJ, e=0.):
    return 28.4329*M_MJ*(PC/365.25)**(-1/3.)*((MSTAR + M_MJ*MJUP/MSUN))**(-2/3.)/np.sqrt(1 - e**2)/1e3
emp = np.array([33.7, 68.1, 54.6, 75.5, 61.6, 64.4, 77.0, 69.0, 40.4, 68.0, 62.1, 18.5])   # M-dwarf transiting BDs, P < 16 d
masses = dict(uniform=np.linspace(13, 80, 68), empirical_Pl16=emp)
fac = {}
for key, val in LR.items():
    e = 0.55 if '_e0.55_' in key else 0.0
    fk = {}
    for mn, Ms in masses.items():
        Kv = K_of(Ms, e)
        lrm = np.interp(Kv, np.r_[0., Ks], np.r_[1., val['LR']])
        fk[mn] = float(np.mean(lrm))
    fac[key] = fk
out['mass_averaged_LR'] = fac
out['K_examples_kms'] = {f'{m}MJ_e0': float(K_of(m)) for m in (13, 20, 30, 40, 50, 60, 70, 80)}
out['K_examples_kms'].update({f'{m}MJ_e0.55': float(K_of(m, 0.55)) for m in (13, 40, 80)})
summ = {mn: [fac[k][mn] for k in fac] for mn in masses}
out['mass_averaged_LR_range'] = {mn: [float(np.min(v_)), float(np.median(v_)), float(np.max(v_))] for mn, v_ in summ.items()}
print('K', out['K_examples_kms'], flush=True)
print('mass-averaged LR [min, median, max] over CCF width/peak/e/sampling', out['mass_averaged_LR_range'], flush=True)
# photocentre semi-axis
plx = 15.497
al = {}
for m in (13, 40, 60, 80):
    q = m*MJUP/MSUN
    a = (G*(MSTAR + q)*MSUN*(PC*DAY)**2/(4*np.pi**2))**(1/3.)/AU
    al[f'{m}MJ'] = dict(a_au=a, alpha_mas=a*q/(MSTAR + q)*plx)
out['photocentre'] = al
print(al)
jdump(dict(method=__doc__, **out), f'{R3}/out/s05_gaia_rv.json')
