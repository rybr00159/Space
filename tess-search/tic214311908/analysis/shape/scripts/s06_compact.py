"""s06: 'hot compact star totally occulted by a cooler, larger star' scenarios on the M3V system (bound).
The eclipsed source A is a white dwarf / compact star (radius R_A free in [0.005, 0.10] Rsun, i.e. WDs, ELM and
pre-ELM WDs; mass M_A free in [0.15, 1.35] Msun) that contributes a fraction f_A of the aperture light; it is
occulted once per orbit by the occulter B, with physical radius and mass set by the scenario:
  single : B = the M3V itself, R_B ~ N(0.381, 0.012) (allowed up to 0.40 for [Fe/H] = +0.4), M_B ~ N(0.370, 0.011)
  twin   : B = one component of an equal twin, R_B ~ N(0.289, 0.009), M_B ~ N(0.264, 0.008)
  lateM  : B = a bound late-M dwarf with dT 3.5-4.5 (M7-M9): R_B in [0.10, 0.13], M_B in [0.08, 0.10]
  any    : B = any M dwarf no more luminous than the observed system, R_B in [0.08, 0.40], M_B = R_B (M-dwarf M-R)
Kepler fixes a/R_A from M_A + M_B and P. Circular orbit by default ('_ecc': e, w free, e < 0.8).
The secondary eclipse (A transiting B: batman occultation of B's limb-darkened light by A, minus the self-lensing
excess 2 (R_E/R_B)^2, diluted by 1 - f_A) is part of the fit: circular orbits are compared with the measured
phase-0.5 depth (INV -0.16 +- 0.79 ppt empirical, s10; REF +0.29 +- 0.77 ppt, referee r18); eccentric orbits with
the phase scan at the predicted secondary phase (INV out/secondary_scan.csv with the 0.79 ppt empirical error; REF
../referee_phot/out/phase_scan_clean.csv). Fit by differential evolution (min chi2; Gaussian priors as penalties).
WD radius: by default tied to the mass, R_A = R_Nauenberg(M) x f, f in [1, 1.3] for M >= 0.45 Msun (finite
temperature, H envelope) and f in [1, 2.5] below (He-core / ELM white dwarfs, which are larger); scenario suffix
'_Rfree' leaves R_A free in [0.005, 0.10] Rsun (unphysical small radii allowed; most generous to the false positive).
usage: s06_compact.py DATASET SCENARIO   ->  out/compact_{DATASET}_{SCENARIO}.json"""
import sys, json, numpy as np
sys.path.insert(0, 'tic214311908/analysis/shape/scripts')
from shapelib import *
from scipy.optimize import differential_evolution
ds, sc = sys.argv[1], sys.argv[2]
rfree = '_Rfree' in sc; sc0 = sc.replace('_Rfree', ''); ecc = sc0.endswith('_ecc'); base = sc0.replace('_ecc', '')
lc = LC({'inv': f'{SD}/data/inv_win.npz', 'ref': f'{SD}/data/ref_win.npz'}[ds])
fit = json.load(open(f'{SD}/out/fit_{ds}_free.json'))
P = fit['post']['P']['med']; T0c = fit['post']['T0']['med']
Rb = dict(single=(0.37, 0.40), twin=(0.27, 0.31), lateM=(0.10, 0.13), any=(0.08, 0.40))[base]
Mb = dict(single=(0.34, 0.40), twin=(0.24, 0.29), lateM=(0.08, 0.10), any=(0.08, 0.40))[base]
pri = dict(single=((0.381, 0.012), (0.370, 0.011)), twin=((0.289, 0.009), (0.264, 0.008)), lateM=(None, None), any=(None, None))[base]
def r_nauenberg(M):
    mu = M/1.44
    return 0.0112*np.sqrt(mu**(-2/3.) - mu**(2/3.))
def unpack(x):
    T0, MA, lRA, bB, fA, RB, MB = x[:7]
    if not rfree:
        fmax = 2.5 if MA < 0.45 else 1.3
        lRA = np.log10(r_nauenberg(MA)*(1 + lRA*(fmax - 1)))
    if base == 'any': MB = RB
    e, w = (0., 90.) if not ecc else (x[7], x[8])
    return T0, MA, 10**lRA, bB, fA, RB, MB, e, w
def geom(MA, RA, RB, MB):
    a = (G*(MA+MB)*MSUN*(P*DAY)**2/(4*np.pi**2))**(1/3.)/RSUN
    return a, RB/RA, a/RA
def model(x):
    T0, MA, RA, bB, fA, RB, MB, e, w = unpack(x)
    a, k, aRA = geom(MA, RA, RB, MB)
    bA = bB*k
    m = lc.transit(T0, P, k, aRA, bA, 0.1, 0.2, e, w)
    return 1 + fA*(m - 1)
import pandas as pd, batman
SEC05 = {'inv': (-0.16e-3, 0.79e-3), 'ref': (0.29e-3, 0.77e-3)}[ds]
if ds == 'inv':
    _sc = pd.read_csv(f'{SD}/out/secondary_scan.csv'); SCAN = (_sc.phase.values, _sc.depth.values, np.full(len(_sc), 0.79e-3))
else:
    _sc = pd.read_csv('tic214311908/analysis/referee_phot/out/phase_scan_clean.csv')
    SCAN = (_sc.phase.values, _sc.depth.values*1e-3, _sc.err.values*1e-3)
def secondary(x):
    T0, MA, RA, bB, fA, RB, MB, e, w = unpack(x)
    a, k, aRA = geom(MA, RA, RB, MB)
    fac1 = (1 + e*np.sin(np.radians(w)))/(1 - e**2); cosi = bB*k/aRA*fac1
    fac2 = (1 + e*np.sin(np.radians(w+180)))/(1 - e**2)
    aRB = a/RB; bsec = cosi*aRB/fac2
    p = batman.TransitParams(); p.t0 = 0; p.per = P; p.rp = RA/RB; p.a = aRB; p.inc = np.degrees(np.arccos(min(1, cosi))); p.ecc = e; p.w = (w+180) % 360
    p.limb_dark = 'quadratic'; p.u = [0.2, 0.4]
    occ = 1 - batman.TransitModel(p, np.array([0.0])).light_curve(p)[0]
    RE = np.sqrt(4*G*MA*MSUN*a*RSUN/2.998e10**2)/RSUN
    lens = 2*(RE/RB)**2 if bsec < 1 else 0.
    # share of the total light that comes from the occulter B (the star the WD transits at the secondary)
    shareB = dict(single=1 - fA, twin=(1 - fA)/2, lateM=0.03, any=(1 - fA)*min(1., (RB/0.381)**3))[base]
    dep = shareB*(occ - lens)
    phs = (sec_offset(P, e, w)/P) % 1
    return dep, phs
def c2(x):
    T0, MA, RA, bB, fA, RB, MB, e, w = unpack(x)
    v = lc.chi2(model(x))
    dep, phs = secondary(x)
    if not ecc:
        v += ((dep - SEC05[0])/SEC05[1])**2
    else:
        j = np.argmin(np.abs(SCAN[0] - phs))
        if np.abs(SCAN[0][j] - phs) < 0.01:
            v += ((dep - SCAN[1][j])/SCAN[2][j])**2
        else:  # predicted secondary near the transit (not scanned): use the circular value as a proxy
            v += ((dep - SEC05[0])/SEC05[1])**2
    if pri[0] is not None:
        v += ((RB - pri[0][0])/pri[0][1])**2 + ((MB - pri[1][0])/pri[1][1])**2
    return v
bounds = [(T0c-0.008, T0c+0.008), (0.15, 1.35), (np.log10(0.005), np.log10(0.10)) if rfree else (0, 1), (0, 1.05), (0.005, 0.06), Rb, Mb]
if ecc: bounds += [(0, 0.8), (0, 360)]
r = differential_evolution(c2, bounds, seed=3, popsize=(15 if (ecc or base in ('lateM', 'any')) else 25), maxiter=(300 if ecc else 600), tol=1e-8, polish=True, init='sobol')
T0, MA, RA, bB, fA, RB, MB, e, w = unpack(r.x)
a, k, aRA = geom(MA, RA, RB, MB)
T14, T23 = t14_t23(P, k, aRA, bB*k, e, w)
RE = np.sqrt(4*G*MA*MSUN*a*RSUN/2.998e10**2)/RSUN
sec, sec_phase = secondary(r.x)
# implied TESS-band brightness temperature of A relative to B (blackbody, TESS response)
lam = np.loadtxt(f'{SD}/data/filt_TESS.Red.dat'); wl = lam[:, 0]*1e-8; tr = lam[:, 1]
def SB(T): return np.trapezoid(tr*wl/(wl**5*(np.exp(1.4388/(wl*T)) - 1)), wl)
TB = dict(single=3350, twin=3237, lateM=2600, any=3300)[base]
ratio_needed = fA/(1 - fA)/(RA/RB)**2
Tg = np.linspace(2000, 60000, 6000); s = np.array([SB(t) for t in Tg])/SB(TB)
TA = float(np.interp(ratio_needed, s, Tg)) if ratio_needed < s.max() else float('nan')
out = dict(dataset=ds, scenario=sc, chi2=float(r.fun), chi2_lc=float(lc.chi2(model(r.x))), chi2_planet_free=fit['chi2_best'],
           dchi2_vs_planet_free=float(r.fun - fit['chi2_best']), T0=T0, M_A=MA, R_A=RA, b_occ=bB, f_A=fA, R_B=RB, M_B=MB, e=e, w=w,
           a_Rsun=a, k=k, aR_A=aRA, T14_h=float(T14*24), T23_h=float(T23*24), ingress_min=float((T14-T23)/2*1440), secondary_ppt=float(sec*1e3), secondary_phase=float(sec_phase),
           einstein_radius_Rsun=float(RE), TA_needed_K_blackbody=TA, TB_K=TB)
json.dump(out, open(f'{SD}/out/compact_{ds}_{sc}.json', 'w'), indent=1, default=float)
print(json.dumps(out, default=float))
