"""s02: FPP of the referee2_fpp framework (r05_fpp.py, imported READ-ONLY, nothing written to referee2_fpp/) with
(0) planet e priors: Kipping (2013) Beta(0.867,3.03), Van Eylen+19 singles half-normal 0.32, circular, and Sagear & Ballard
    (2023, 2023PNAS..12017398S, arXiv:2305.17157) M-dwarf single-transit Rayleigh sigma 0.19;
(1) the photo-eccentric sign corrected by my own patch: r05 evaluates the grid at rho_c = rho_host*ge^-3, so I make its
    ecc() return ge = 1/g (g = (1 + e sin w)/sqrt(1-e^2)), giving rho_c = rho_host*g^+3 (s01 verifies the sign with an
    independent Kepler-orbit code); draws increased from 2,500 to N_E (MC convergence check);
(2) a BD-specific eccentricity prior (Ferreira+26 short-period Beta(0.426+-0.092, 2.354+-0.73), long-period
    Beta(1.879+-0.92, 2.470+-1.25)) for the BD class only, planets keep Kipping / Van Eylen / circular;
(3) f_BD per star at 3-10 d: 1e-3 (referee2_fpp), 3.5e-4 (bd 95% upper), 1.3e-4 (bd central), 5e-5;
(4) an optional Gaia-RVS factor on the BD terms of the M3V hosts (BD on M3V single and BD on primary) from s05.
Output: out/s02_fpp_fixed.json, out/s02_fpp_fixed.csv"""
import sys, os
sys.path.insert(0, 'tic214311908/analysis/referee3_bd/scripts')
from rlib3 import *
import pandas as pd
sys.path.insert(0, f'{B}/referee2_fpp/scripts')
cwd = os.getcwd(); os.chdir(f'{B}/referee2_fpp')
import r05_fpp as R
os.chdir(cwd)
N_E = int(os.environ.get('S02_NE', 10000))
RVF = float(os.environ.get('S02_RVF', 1.0))
def edraw(kind, n, rng):
    if kind == 'kipping': return rng.beta(0.867, 3.03, n)
    if kind == 'vaneylen': return np.abs(rng.normal(0, 0.32, n))
    if kind == 'circ': return np.zeros(n)
    if kind == 'rayleigh19': return rng.rayleigh(0.19, n)
    if kind == 'bdsp':
        a = np.clip(rng.normal(0.426, 0.092, n), 0.05, None); b = np.clip(rng.normal(2.354, 0.73, n), 0.1, None); return rng.beta(a, b)
    if kind == 'bdlp':
        a = np.clip(rng.normal(1.879, 0.92, n), 0.05, None); b = np.clip(rng.normal(2.470, 1.25, n), 0.1, None); return rng.beta(a, b)
    raise ValueError(kind)
SIGN = {'mode': 'none'}
def ecc_patched(kind, n=2500):
    rng = np.random.default_rng(12345)
    nn = n if SIGN['mode'] == 'orig_n' else N_E
    e = np.clip(edraw(kind, nn, rng), 0, 0.9); w = rng.uniform(0, 2*np.pi, nn)
    g = (1 + e*np.sin(w))/np.sqrt(1 - e**2); geo = (1 + e*np.sin(w))/(1 - e**2)
    if SIGN['mode'] in ('wrong', 'orig_n'): return g, geo
    return 1/g, geo
orig_build = R.build
def build_patched(ds, kind, alpha):
    if '+' in kind:
        kp, kb = kind.split('+')
        tp = orig_build(ds, kp, alpha); tb = orig_build(ds, kb, alpha)
        out = dict(tp); out['BD'] = tb['BD']
        return out
    return orig_build(ds, kind, alpha)
R.ecc = ecc_patched; R.build = build_patched
rows = []
def run(tag, mode='fixed', rvf=1.0, **kw):
    if SIGN['mode'] != mode:
        R.ZT.clear(); R.ZCACHE.clear()
    SIGN['mode'] = mode
    R.rng = np.random.default_rng(11)
    p = dict(R.BASE); p.update(kw)
    r = R.evaluate(p, detail=True)
    po = dict(r['post'])
    if rvf != 1.0:     # Gaia RV factor on BDs orbiting the M3V (single or primary of a binary): rescale and renormalise
        po['BD_on_M3V_single'] *= rvf; po['BD_on_primary_of_binary'] *= rvf
        s = sum(po.values()); po = {k: v/s for k, v in po.items()}
    planet = po['TP_single'] + po['TP_diluted'] + po['TP_on_twin']
    bd = po['BD_on_M3V_single'] + po['BD_on_primary_of_binary'] + po['BDHEB_on_companion']
    row = dict(tag=tag, sign=mode, rvf=rvf, **{k: str(v) for k, v in kw.items()}, FPP=1 - planet, FPP_system=1 - planet - po['PTP_companion'],
               PTP=po['PTP_companion'], BD_on_M3V=po['BD_on_M3V_single'], BD_on_primary=po['BD_on_primary_of_binary'],
               BD_on_companion=po['BDHEB_on_companion'], BD_total=bd, NTP_pairA=po['NTP_pairA'], WD=po['EB_WD_on_M3V'],
               HEB=po['HEB_MSMS'], BG=po['BG_all'], TP_single=po['TP_single'], TP_diluted=po['TP_diluted'], TP_twin=po['TP_on_twin'],
               P_planet_any_star=planet + po['PTP_companion'])
    rows.append(row)
    print(f"{tag:58s} FPP {row['FPP']:.5f} sys {row['FPP_system']:.5f} PTP {row['PTP']:.2e} BDM3V {row['BD_on_M3V']:.2e} BDcomp {row['BD_on_companion']:.2e} BDtot {bd:.2e} NTP {row['NTP_pairA']:.1e}", flush=True)
    return row
# (a) reproduction of r05 as published (wrong sign, 2,500 draws) and the bd/b10 fix at 2,500 draws
for ds in ('ref_beta', 'inv_beta'):
    run(f'{ds} kipping WRONG sign n=2500 (r05 as published)', mode='orig_n', ds=ds, ecc='kipping')
    run(f'{ds} kipping WRONG sign n={N_E}', mode='wrong', ds=ds, ecc='kipping')
    run(f'{ds} kipping FIXED n={N_E}', ds=ds, ecc='kipping')
# (b) full grid
for ds in ('ref_beta', 'inv_beta', 'ref_white', 'inv_white'):
    for pe in ('kipping', 'vaneylen', 'circ', 'rayleigh19'):
        for be in ('same', 'bdsp', 'bdlp'):
            kind = pe if be == 'same' else f'{pe}+{be}'
            for fbd in (1e-3, 3.5e-4, 1.3e-4, 5e-5):
                run(f'{ds} planet-e {pe} BD-e {be} f_BD {fbd:g}', ds=ds, ecc=kind, f_BD=fbd)
pd.DataFrame(rows).to_csv(f'{R3}/out/s02_fpp_fixed.csv', index=False)
jdump(dict(method=__doc__, N_E=N_E, rows=rows), f'{R3}/out/s02_fpp_fixed.json')
