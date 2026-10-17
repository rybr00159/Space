"""s12: physical eclipsing-binary fits for the two catalogued background stars, joint with their photometry.
Source: bg3 = Gaia DR3 5925373272827260544 (3.1", F = 30.3 e-/s in PRF units -> f = 0.031 of the M3V's 979 e-/s),
bg5 = Gaia DR3 5925373272827262208 (5.5", F = 22.1 e-/s -> f = 0.0225); the flux fraction gets a 10% Gaussian
nuisance (PRF/response uncertainty; a response < 1 would require even deeper eclipses).
Data: the M3V light curve (INV or REF windows) + the star's DECaPS2 mean AB magnitudes g, r, i, z (sigma 0.05 mag
floor) + Gaia DR3 parallax (+0.03 mas zero-point, sigma as published) + E(B-V) in [0, 0.30] (SFD 0.20 total).
Binary classes (P fixed at the planet value or twice it):
  msms5 / msms11 : two main-sequence stars (Mamajek sequence, 0.08-1.4 Msun), P or 2P (equal eclipses)
  cms            : hot compact star A (M 0.15-1.4, R 0.008-0.5 Rsun, Teff 4000-50000 K; ATLAS9 logg 5 colours)
                   totally/partially occulted by a main-sequence star B (WD, ELM/pre-ELM, sdB-like)
  cms_ecc        : same, eccentric (e < 0.6, any w) -- physically implausible for post-interaction binaries
  (cms has free M, R, Teff and finds unphysical objects, e.g. 1.4 Msun with R = 0.2 Rsun; the physical sub-classes
  below restrict A to real stellar loci:)
  cwd            : white dwarf, R tied to M (Nauenberg x 1-1.3, x 1-2.5 below 0.45 Msun), M 0.17-1.3, Teff 4-50 kK
  csdb           : hot subdwarf, M 0.40-0.55, R 0.10-0.30, Teff 20-40 kK
  celm           : (pre-)ELM white dwarf / EL CVn-type stripped star, M 0.15-0.30, R 0.03-0.50, Teff 6-25 kK
  msbd           : main-sequence star A transited by a dark body (BD / giant planet, R 0.06-0.25 Rsun)
Secondary eclipse at phase 0.5 (P = 5.69 d, circular) is constrained by the measured phase-0.5 depth; for 2P both
eclipses fall in the windows. Fit: differential evolution on the total chi2 = LC + secondary + photometry +
parallax + flux nuisance. 'lconly' variants drop the photometry/parallax terms (shape-only test).
usage: s12_bg_eb.py DATASET SOURCE CLASS [lconly]   -> out/bg_eb_{DATASET}_{SOURCE}_{CLASS}[_lconly].json"""
import sys, json, numpy as np
sys.path.insert(0, 'tic214311908/analysis/shape/scripts')
from shapelib import *
from eem import ms
from specs import S, absmag, EXT
from scipy.optimize import differential_evolution
import batman
ds, src, cls = sys.argv[1], sys.argv[2], sys.argv[3]
lconly = len(sys.argv) > 4 and sys.argv[4] == 'lconly'
lc = LC({'inv': f'{SD}/data/inv_win.npz', 'ref': f'{SD}/data/ref_win.npz'}[ds])
fit = json.load(open(f'{SD}/out/fit_{ds}_free.json'))
P0 = fit['post']['P']['med']; T0c = fit['post']['T0']['med']
SEC = {'inv': (-0.16e-3, 0.79e-3), 'ref': (0.29e-3, 0.77e-3)}[ds]  # s10 empirical (INV), referee r18 (REF)
OBS = dict(bg3=dict(f=30.3/979, mag=dict(g=18.164, r=17.352, i=17.100, z=16.900), plx=(0.655+0.03, 0.110)),
           bg5=dict(f=22.06/979, mag=dict(g=18.090, r=17.589, i=17.407, z=17.304), plx=(0.283+0.03, 0.115)))[src]
MSIG = 0.05
def ldc(T):
    return (0.13, 0.53) if T < 3900 else ((0.35, 0.25) if T < 5500 else ((0.25, 0.25) if T < 8000 else (0.12, 0.2)))
# parameter vectors per class
if cls in ('msms5', 'msms11'):
    names = ['mA', 'mB', 'bf', 'T0', 'logd', 'ebv', 'fs']
    bounds = [(0.08, 1.4), (0.08, 1.4), (0, 1), (T0c-0.006, T0c+0.006), (2.7, 4.3), (0, 0.3), (0.7, 1.3)]
elif cls in ('cms', 'cms_ecc'):
    names = ['MA', 'logRA', 'logTA', 'mB', 'bf', 'T0', 'logd', 'ebv', 'fs'] + (['e', 'w'] if cls == 'cms_ecc' else [])
    bounds = [(0.15, 1.4), (np.log10(0.008), np.log10(0.5)), (np.log10(4000), np.log10(50000)), (0.08, 1.4), (0, 1), (T0c-0.006, T0c+0.006), (2.7, 4.3), (0, 0.3), (0.7, 1.3)] + ([(0, 0.6), (0, 360)] if cls == 'cms_ecc' else [])
elif cls in ('cwd', 'csdb', 'celm'):
    names = ['MA', 'xRA', 'logTA', 'mB', 'bf', 'T0', 'logd', 'ebv', 'fs']
    mb, tb = dict(cwd=((0.17, 1.3), (4000, 50000)), csdb=((0.40, 0.55), (20000, 40000)), celm=((0.15, 0.30), (6000, 25000)))[cls]
    bounds = [mb, (0, 1), (np.log10(tb[0]), np.log10(tb[1])), (0.08, 1.4), (0, 1), (T0c-0.006, T0c+0.006), (2.7, 4.3), (0, 0.3), (0.7, 1.3)]
elif cls == 'msbd':
    names = ['mA', 'RB', 'MB', 'bf', 'T0', 'logd', 'ebv', 'fs']
    bounds = [(0.08, 1.4), (0.06, 0.25), (0.01, 0.08), (0, 1), (T0c-0.006, T0c+0.006), (2.7, 4.3), (0, 0.3), (0.7, 1.3)]
P = 2*P0 if cls == 'msms11' else P0
def comps(x):
    p = dict(zip(names, x))
    if cls in ('msms5', 'msms11'):
        A = dict(M=p['mA'], R=ms(p['mA'], 'R_Rsun'), T=ms(p['mA'], 'Teff'), compact=False)
        B = dict(M=p['mB'], R=ms(p['mB'], 'R_Rsun'), T=ms(p['mB'], 'Teff'), compact=False, dark=False)
    elif cls in ('cms', 'cms_ecc'):
        A = dict(M=p['MA'], R=10**p['logRA'], T=10**p['logTA'], compact=True)
        B = dict(M=p['mB'], R=ms(p['mB'], 'R_Rsun'), T=ms(p['mB'], 'Teff'), compact=False, dark=False)
    elif cls in ('cwd', 'csdb', 'celm'):
        if cls == 'cwd':
            mu = p['MA']/1.44; rn = 0.0112*np.sqrt(mu**(-2/3.) - mu**(2/3.)); R = rn*(1 + p['xRA']*((2.5 if p['MA'] < 0.45 else 1.3) - 1))
        elif cls == 'csdb':
            R = 0.10 + p['xRA']*0.20
        else:
            R = 10**(np.log10(0.03) + p['xRA']*(np.log10(0.5) - np.log10(0.03)))
        A = dict(M=p['MA'], R=R, T=10**p['logTA'], compact=True)
        B = dict(M=p['mB'], R=ms(p['mB'], 'R_Rsun'), T=ms(p['mB'], 'Teff'), compact=False, dark=False)
    else:
        A = dict(M=p['mA'], R=ms(p['mA'], 'R_Rsun'), T=ms(p['mA'], 'Teff'), compact=False)
        B = dict(M=p['MB'], R=p['RB'], T=1000., compact=False, dark=True)
    return p, A, B
def evaluate(x, detail=False):
    p, A, B = comps(x)
    e, w = (p.get('e', 0.), p.get('w', 90.))
    a = (G*(A['M']+B['M'])*MSUN*(P*DAY)**2/(4*np.pi**2))**(1/3.)/RSUN
    if a*(1-e) < 1.05*(A['R'] + B['R']):
        return 1e9 if not detail else None
    k = B['R']/A['R']; aR = a/A['R']; b = p['bf']*(1 + k)
    LA = S(A['T'], 'T', A['compact'])*A['R']**2
    LB = 0. if B.get('dark') else S(B['T'], 'T', False)*B['R']**2
    ftot = OBS['f']*p['fs']
    m = eb_flux(lc, p['T0'], P, k, aR, b, ldc(A['T']), ldc(B['T']), LB/LA, ftot, e, w)
    c_lc = lc.chi2(m)
    c_sec = 0.; sec = 0.
    if cls != 'msms11' and LB > 0:
        # predicted secondary depth (B behind A) at its conjunction
        fac1 = (1 + e*np.sin(np.radians(w)))/(1 - e**2); cosi = b/aR*fac1
        fac2 = (1 + e*np.sin(np.radians(w+180)))/(1 - e**2)
        bB = cosi*(aR/k)/fac2
        pp = batman.TransitParams(); pp.t0 = 0; pp.per = P; pp.rp = 1/k; pp.a = aR/k; pp.inc = np.degrees(np.arccos(min(1, cosi)))
        pp.ecc = e; pp.w = (w+180) % 360; pp.limb_dark = 'quadratic'; pp.u = list(ldc(B['T']))
        fB = batman.TransitModel(pp, np.array([0.0])).light_curve(pp)[0]
        sec = ftot*LB/(LA+LB)*(1 - fB)
        c_sec = ((sec - SEC[0])/SEC[1])**2
    c_ph = 0.; c_plx = 0.; mags = {}
    d = 10**p['logd']
    for bnd, mo in OBS['mag'].items():
        fA = 10**(-0.4*absmag(A['T'], A['R'], bnd, A['compact']))
        fB_ = 0. if B.get('dark') else 10**(-0.4*absmag(B['T'], B['R'], bnd, False))
        mags[bnd] = -2.5*np.log10(fA + fB_) + 5*np.log10(d/10) + EXT[bnd]*p['ebv']
        c_ph += ((mags[bnd] - mo)/MSIG)**2
    c_plx = ((1000/d - OBS['plx'][0])/OBS['plx'][1])**2
    c_fs = ((p['fs'] - 1)/0.1)**2
    tot = c_lc + c_sec + c_fs + (0 if lconly else c_ph + c_plx)
    if not detail: return tot
    T14, T23 = t14_t23(P, k, aR, b, e, w)
    return dict(total=float(tot), chi2_lc=float(c_lc), chi2_sec=float(c_sec), chi2_phot=float(c_ph), chi2_plx=float(c_plx), chi2_fs=float(c_fs),
                params={n: float(v) for n, v in zip(names, x)}, A={kk: (float(vv) if not isinstance(vv, bool) else vv) for kk, vv in A.items()},
                B={kk: (float(vv) if not isinstance(vv, bool) else vv) for kk, vv in B.items()}, a_Rsun=float(a), k=float(k), aR=float(aR), b=float(b),
                LB_LA_T=float(LB/LA), f_total=float(ftot), eclipsed_fraction_needed=float(0.0158/(ftot*LA/(LA+LB))) if LA > 0 else None,
                sec_ppt=float(sec*1e3), mags_model={kk: float(vv) for kk, vv in mags.items()}, mags_obs=OBS['mag'], d_pc=float(d),
                T14_h=float(T14*24), T23_h=float(T23*24))
def c2(x): return evaluate(x)
r = differential_evolution(c2, bounds, seed=5, popsize=14, maxiter=350, tol=1e-8, polish=True, init='sobol', mutation=(0.5, 1.0), recombination=0.7)
det = evaluate(r.x, detail=True)
det.update(dataset=ds, source=src, cls=cls, lconly=lconly, planet_free_chi2_lc=fit['chi2_best'], dchi2_lc_vs_planet=det['chi2_lc'] - fit['chi2_best'], P=P, nfev=int(r.nfev))
json.dump(det, open(f'{SD}/out/bg_eb_{ds}_{src}_{cls}{"_lconly" if lconly else ""}.json', 'w'), indent=1, default=float)
print(json.dumps(det, default=float))
